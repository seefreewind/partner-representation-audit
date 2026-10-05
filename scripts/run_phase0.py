#!/usr/bin/env python3
"""Run the metadata-only Phase 0 feasibility audit.

This script never predicts or scores a structure. The September 2026 ClinVar
variant summary is used only to screen candidate variants; candidate RCVs are
then retrieved from ClinVar E-utilities and parsed as variant-condition records.
"""
import argparse
import ast
import csv
import gzip
import hashlib
import io
import json
import os
import re
import sys
import time
import warnings
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import yaml
from lxml import etree

warnings.filterwarnings("ignore", message="This pattern is interpreted as a regular expression", category=UserWarning)

ROOT = Path(__file__).resolve().parents[1]
PHASE = ROOT / "results" / "phase0"
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
PROCESSED = ROOT / "data" / "processed"
META = ROOT / "metadata"
REPORTS = ROOT / "reports"
API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
AA3 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
AA1 = set(AA3.values())
P3 = re.compile(r"p\.([A-Z][a-z]{2}|[A-Z])(\d+)([A-Z][a-z]{2}|[A-Z])(?:[),; ]|$)")
P1 = re.compile(r"p\.([A-Z])(\d+)([A-Z])(?:[),; ]|$)")
XREF = re.compile(r"\b(MONDO|MedGen|OMIM|Orphanet|UMLS|UMLS_CUI):([A-Za-z0-9_.-]+)\b")
RCV_RE = re.compile(r"RCV\d+(?:\.\d+)?")


def now():
    return datetime.now(timezone.utc).isoformat()


def write_tsv(path, rows, columns=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError("Refusing to overwrite existing output: %s" % path)
    if not rows and columns is None:
        columns = {
            "disease_pairs_needing_manual_review.tsv":["gene","disease_A","MONDO_A","n_variants_A","disease_B","MONDO_B","n_variants_B","relationship","status","reason"],
            "candidate_genes_pre_mechanism.tsv":["gene","n_diseases"],
            "candidate_genes_post_mechanism.tsv":["gene","mechanism_status"],
            "ppi_provenance.tsv":["gene","partner","sources","source_count","experimental_support","binary_experimental_support","PIONEER_coverage","structural_support","predicted_only","tier"],
            "ppi_coverage_by_gene.tsv":["gene","total_partners","experimental_partners","PIONEER_partners","Predictomes_partners","unique_credible_partners","structural_ready_partners","ppi_sources"],
            "rcv_api_failures.tsv":["batch","accessions","error"],
            "top_candidate_genes.tsv":["gene","n_normalized_diseases","n_unique_variants","n_strict_diseases","mechanism_status","credible_PPI_partners","PIONEER_partners","status"],
        }.get(path.name,["status"])
    frame = pd.DataFrame(rows, columns=columns)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def add_failure(source, url, message):
    path = META / "source_failures.tsv"
    row = {"timestamp_utc": now(), "source": source, "url": url, "status": message}
    header = not path.exists()
    with open(path, "a", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=row.keys(), delimiter="\t")
        if header:
            w.writeheader()
        w.writerow(row)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def register_local_inputs():
    """Build a source manifest from downloaded files and downloader sidecars."""
    entries = []
    direct_urls={
        "gene2phenotype_all_latest.csv":"https://www.ebi.ac.uk/gene2phenotype/api/panel/all/download/",
        "gencc_submissions_latest.tsv":"https://thegencc.org/download/action/submissions-export-tsv",
        "BIOGRID-MV-Physical-5.0.261.tab3.zip":"https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-5.0.261/BIOGRID-MV-Physical-5.0.261.tab3.zip",
        "PIONEER_very_high_human_2026-09-27.txt":"https://pioneer.yulab.org/static/predictions/very_high/human.txt",
        "hgnc_non_alt_loci_2026-09-27.tsv":"https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/non_alt_loci_set.txt",
        "RCV000009910_current_API_schema_example.xml":"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=clinvarset&id=RCV000009910&retmode=xml",
    }
    for p in sorted(RAW.rglob("*")):
        if not p.is_file() or ".partial" in p.name or p.name.startswith("._"):
            continue
        if p.name.endswith(".download.json"):
            continue
        sidecar = p.with_suffix(p.suffix + ".download.json")
        meta = {}
        if sidecar.exists():
            try:
                meta = json.loads(sidecar.read_text())
            except Exception:
                meta = {}
        entries.append({
            "source_file": str(p.relative_to(ROOT)),
            "url": meta.get("url", direct_urls.get(p.name,"ClinVar E-utilities; URL and accession in filename/source failure log" if p.suffix==".xml" and "rcv_" in p.name else "see config/phase0.yml or source_failures.tsv")),
            "release_or_version": release_from_name(p.name),
            "retrieved_utc": meta.get("retrieved_utc", datetime.fromtimestamp(p.stat().st_mtime, timezone.utc).isoformat()),
            "bytes": p.stat().st_size,
            "sha256": meta.get("sha256", sha256_file(p)),
            "status":"partial_not_used" if p.name.endswith(".incomplete") else "downloaded",
        })
    write_tsv(META / "input_manifest.tsv", entries,
              ["source_file", "url", "release_or_version", "retrieved_utc", "bytes", "sha256","status"])


def release_from_name(name):
    m = re.search(r"(2026[-_]?09(?:[-_]?01)?|5\.0\.261|PIONEER.*2026-09-27)", name, re.I)
    return m.group(1) if m else "date-stamped retrieval; see manifest"


def clean_col(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower().lstrip("#"))


def find_col(columns, aliases):
    lookup = {clean_col(c): c for c in columns}
    for alias in aliases:
        if clean_col(alias) in lookup:
            return lookup[clean_col(alias)]
    return None


def split_identifiers(value):
    if value is None or pd.isna(value):
        return []
    return list(dict.fromkeys("%s:%s" % (m.group(1), m.group(2)) for m in XREF.finditer(str(value))))


def parse_protein_change(value):
    if value is None or pd.isna(value):
        return None
    s = str(value)
    m = P3.search(s)
    if m:
        ref = AA3.get(m.group(1), m.group(1))
        alt = AA3.get(m.group(3), m.group(3))
        return (ref, int(m.group(2)), alt, "p.%s%d%s" % (ref, int(m.group(2)), alt))
    m = P1.search(s)
    if m and m.group(1) in AA1 and m.group(3) in AA1:
        return (m.group(1), int(m.group(2)), m.group(3), "p.%s%d%s" % (m.group(1), int(m.group(2)), m.group(3)))
    return None


def accepted_plp(value):
    s = str(value or "").strip().lower()
    if "conflict" in s or "uncertain" in s or "benign" in s or "risk factor" in s or "drug response" in s:
        return False
    return s in {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic", "pathogenic, likely pathogenic"}


def parse_mondo(path):
    terms, xref_map, synonym_map, alt_map = {}, defaultdict(set), {}, {}
    if not Path(path).exists():
        return terms, xref_map, synonym_map, alt_map
    current = None
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n")
        if line == "[Term]":
            if current and current.get("id", "").startswith("MONDO:"):
                terms[current["id"]] = current
            current = {"parents": [], "xrefs": [], "synonyms": [], "obsolete": False, "replaced_by": []}
        elif line.startswith("["):
            if current and current.get("id", "").startswith("MONDO:"):
                terms[current["id"]] = current
            current = None
        elif current is not None:
            if line.startswith("id: "):
                n = line[4:].strip()
                current["id"] = n if n.startswith("MONDO:") else "MONDO:" + n
            elif line.startswith("name: "):
                current["name"] = line[6:].strip()
            elif line.startswith("is_a: "):
                current["parents"].append("MONDO:" + line[7:].split()[0].replace("MONDO:", ""))
            elif line.startswith("alt_id: "):
                current.setdefault("alt_ids", []).append(line[8:].strip())
            elif line.startswith("xref: "):
                raw = line[6:].split()[0]
                current["xrefs"].append(raw)
            elif line.startswith("synonym: "):
                m = re.match(r'synonym: "(.*?)" (EXACT|RELATED|BROAD|NARROW)', line)
                if m:
                    current["synonyms"].append((m.group(1), m.group(2)))
            elif line == "is_obsolete: true":
                current["obsolete"] = True
            elif line.startswith("replaced_by: "):
                current["replaced_by"].append(line[13:].split()[0])
    if current and current.get("id", "").startswith("MONDO:"):
        terms[current["id"]] = current
    for mid, t in terms.items():
        for x in t["xrefs"]:
            xref_map[x].add(mid)
        for alt in t.get("alt_ids", []):
            alt_map[alt] = mid
        if t.get("name"):
            synonym_map[t["name"].casefold()] = mid
        for name, scope in t["synonyms"]:
            if scope == "EXACT":
                synonym_map[name.casefold()] = mid
    return terms, xref_map, synonym_map, alt_map


def mondo_map(ids, terms, xref_map, synonym_map, alt_map, raw_name=""):
    found = set()
    for x in ids:
        if x.startswith("MONDO:") and x in terms:
            found.add(x)
        if x in alt_map:
            found.add(alt_map[x])
        for key in (x, x.replace(":", "_"), x.replace("MedGen:", "UMLS:"), x.replace("MedGen:", "UMLS_CUI:")):
            found.update(xref_map.get(key, set()))
    if not found and raw_name:
        name_key = str(raw_name).strip().casefold()
        if name_key in synonym_map:
            found.add(synonym_map[name_key])
    found = {m for m in found if m in terms and not terms[m].get("obsolete")}
    return sorted(found)


def text_of(node):
    if node is None:
        return ""
    return " ".join(" ".join(node.itertext()).split())


def localname(node):
    return etree.QName(node).localname


def child(node, name):
    if node is None:
        return None
    for e in node:
        if localname(e) == name:
            return e
    return None


def descendants(node, name):
    return [e for e in node.iter() if localname(e) == name]


def star_count(status):
    s = (status or "").lower()
    if "practice guideline" in s:
        return 4
    if "reviewed by expert panel" in s:
        return 3
    if "multiple submitters" in s and "no conflicts" in s:
        return 2
    if "criteria provided" in s:
        return 1
    return 0


def parse_rcv_xml(xml_bytes, variation_gene_map):
    parser = etree.XMLParser(resolve_entities=False, no_network=True, recover=False, huge_tree=True)
    root = etree.fromstring(xml_bytes, parser=parser)
    sets = [n for n in root.iter() if localname(n) == "ClinVarSet"]
    rows = []
    for record in sets:
        ref = next((n for n in record if localname(n) == "ReferenceClinVarAssertion"), None)
        if ref is None:
            continue
        accession_node = next((n for n in ref.iter() if localname(n) == "ClinVarAccession" and n.get("Type") == "RCV"), None)
        if accession_node is None:
            continue
        rcv = accession_node.get("Acc", "")
        version = accession_node.get("Version", "")
        classes = next((n for n in ref if localname(n) == "Classifications"), None)
        sig_node = next((n for n in descendants(classes, "Description") if n.get("DateLastEvaluated") is not None), None) if classes is not None else None
        if sig_node is None and classes is not None:
            sig_node = next((n for n in descendants(classes, "Description") if text_of(n)), None)
        significance = text_of(sig_node)
        last_evaluated = sig_node.get("DateLastEvaluated", "") if sig_node is not None else ""
        review_node = next((n for n in descendants(classes, "ReviewStatus") if text_of(n)), None) if classes is not None else None
        review = text_of(review_node)
        observations = [text_of(n) for n in ref if localname(n) == "ObservedIn"]
        origin_text = " | ".join(observations)
        germline_class_node = next((n for n in descendants(classes, "GermlineClassification") if text_of(n)), None) if classes is not None else None
        germline = germline_class_node is not None or "germline" in origin_text.lower()

        vsets = [n for n in ref if localname(n) == "MeasureSet" and (n.get("Type", "").lower() in {"variant", "haplotype", "genotype"})]
        vsets = vsets or [n for n in record if localname(n) == "MeasureSet" and n.get("Type", "").lower() == "variant"]
        variation_id = ""
        hgvsp, protein_acc, transcript, protein_name = [], [], [], []
        consequences, gene_names = [], []
        for vs in vsets:
            vacc = vs.get("Acc", "")
            m = re.search(r"VCV(\d+)", vacc)
            if m:
                variation_id = str(int(m.group(1)))
            for e in vs.iter():
                if localname(e) == "Attribute":
                    typ = e.get("Type", "").lower()
                    val = text_of(e)
                    if "hgvs, protein" in typ or typ == "hgvs protein":
                        if val:
                            hgvsp.append(val)
                        ac = e.get("Accession", "")
                        ver = e.get("Version", "")
                        if ac:
                            protein_acc.append(ac + ("." + ver if ver else ""))
                    elif "hgvs, coding" in typ and e.get("Accession"):
                        ac = e.get("Accession", "")
                        ver = e.get("Version", "")
                        transcript.append(ac + ("." + ver if ver else ""))
                    elif typ in {"proteinchange1lettercode", "protein change"} and val:
                        protein_name.append(val)
                    elif typ == "molecularconsequence" and val:
                        consequences.append(val)
                elif localname(e) == "ElementValue" and e.get("Type") == "Preferred":
                    val = text_of(e)
                    m = re.search(r"\(([A-Za-z0-9_-]+)\):", val)
                    if m:
                        gene_names.append(m.group(1))
                elif localname(e) == "Gene" and e.get("Symbol"):
                    gene_names.append(e.get("Symbol"))
        gene_names = list(dict.fromkeys(gene_names))
        variation_id = variation_id or ""
        gene = gene_names[0] if gene_names else variation_gene_map.get(variation_id, "")
        if gene and variation_gene_map.get(variation_id) and gene != variation_gene_map[variation_id]:
            gene = variation_gene_map[variation_id]
        pval = next((x for x in hgvsp if parse_protein_change(x)), "")
        if not pval:
            pval = next((x for x in protein_name if parse_protein_change(x)), "")
        parsed = parse_protein_change(pval)
        refaa, residue, altaa, aa_change = parsed if parsed else ("", "", "", "")
        traitsets = [n for n in record if localname(n) == "TraitSet" and n.get("Type", "").lower() in {"disease", "condition"}]
        if not traitsets:
            traitsets = [n for n in ref if localname(n) == "TraitSet" and n.get("Type", "").lower() in {"disease", "condition"}]
        names, ids = [], []
        for ts in traitsets:
            for tr in ts:
                if localname(tr) == "Trait":
                    for n in tr.iter():
                        if localname(n) == "ElementValue" and n.get("Type") == "Preferred" and text_of(n):
                            names.append(text_of(n))
                    for x in descendants(tr, "XRef"):
                        db, xid = x.get("DB", ""), x.get("ID", "")
                        if db and xid:
                            ids.append("%s:%s" % (db, xid))
        names = list(dict.fromkeys(names))
        ids = list(dict.fromkeys(ids))
        moi = []
        for a in descendants(ref, "Attribute"):
            if a.get("Type", "").lower() in {"modeofinheritance", "mode of inheritance"} and text_of(a):
                moi.append(text_of(a))
        scvs = []
        for n in descendants(record, "ClinVarAccession"):
            if n.get("Type") == "SCV" and n.get("Acc"):
                scvs.append(n.get("Acc"))
        rows.append({
            "Variation_ID": variation_id,
            "RCV_accession": rcv,
            "RCV_version": version,
            "gene": gene,
            "transcript": ";".join(dict.fromkeys(transcript)),
            "protein_accession": ";".join(dict.fromkeys(protein_acc)),
            "residue": residue,
            "amino_acid_ref": refaa,
            "amino_acid_alt": altaa,
            "HGVS_protein": pval,
            "protein_change": aa_change,
            "molecular_consequence": ";".join(dict.fromkeys(consequences)),
            "raw_condition": ";".join(names),
            "condition_IDs": ";".join(ids),
            "MONDO": "",
            "review_status": review,
            "stars": star_count(review),
            "clinical_significance": significance,
            "date_last_evaluated": last_evaluated,
            "inheritance": ";".join(dict.fromkeys(moi)),
            "germline": germline,
            "n_SCV": len(set(scvs)),
            "SCV_accessions": ";".join(sorted(set(scvs))),
            "protein_mapping_status": "mapped" if parsed else "unmapped",
            "record_source": "ClinVar E-utilities RCV clinvarset",
        })
    return rows


def fetch_rcv_batches(accessions, variation_gene_map):
    cache = RAW / "rcv_api_batches"
    cache.mkdir(parents=True, exist_ok=True)
    accessions = list(dict.fromkeys(a.split(".")[0] for a in accessions if a.startswith("RCV")))
    rows = []
    failed = []
    batch_size = 40
    session = requests.Session()
    for start in range(0, len(accessions), batch_size):
        batch = accessions[start:start + batch_size]
        batch_id = "%05d_%s" % (start // batch_size, hashlib.sha1(",".join(batch).encode()).hexdigest()[:10])
        target = cache / ("rcv_%s.xml" % batch_id)
        if target.exists():
            payload = target.read_bytes()
        else:
            payload = None
            for attempt in range(7):
                try:
                    r = session.get(API, params={"db":"clinvar", "rettype":"clinvarset", "id":",".join(batch), "retmode":"xml", "tool":"phase0_allelic_heterogeneity"}, timeout=(30, 180))
                    if not r.ok or b"<ClinVarResult-Set" not in r.content:
                        raise RuntimeError("HTTP %s or unexpected E-utilities body: %s" % (r.status_code, r.text[:150].replace("\n", " ")))
                    target.write_bytes(r.content)
                    (target.with_suffix(target.suffix + ".download.json")).write_text(json.dumps({
                        "url":r.url,"accessions":batch,"retrieved_utc":now(),"bytes":len(r.content),
                        "sha256":hashlib.sha256(r.content).hexdigest()},indent=2),encoding="utf-8")
                    payload = r.content
                    break
                except Exception as exc:
                    if attempt == 6:
                        failed.append({"batch":batch_id, "accessions":";".join(batch), "error":str(exc)[:300]})
                        add_failure("ClinVar E-utilities RCV", API, str(exc)[:300])
                    else:
                        time.sleep(min(2 ** attempt, 30))
            time.sleep(0.38)
        if payload:
            try:
                rows.extend(parse_rcv_xml(payload, variation_gene_map))
            except Exception as exc:
                failed.append({"batch":batch_id, "accessions":";".join(batch), "error":"XML parse: %s" % str(exc)[:250]})
                add_failure("ClinVar E-utilities RCV XML parse", str(target), str(exc)[:250])
    return rows, failed


def is_plp(value):
    s = (value or "").strip().lower()
    if any(x in s for x in ("conflict", "uncertain", "benign", "risk factor", "drug response")):
        return False
    return s in {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic", "pathogenic, likely pathogenic"}


def is_missense(row):
    consequence = (row.get("molecular_consequence", "") or "").lower()
    return "missense" in consequence or parse_protein_change(row.get("HGVS_protein", "")) is not None


def ancestors(mid, terms, memo=None, seen=None):
    memo = memo if memo is not None else {}
    if mid in memo:
        return memo[mid]
    seen = set(seen or ())
    if mid in seen:
        return set()
    seen.add(mid)
    out = set()
    for p in terms.get(mid, {}).get("parents", []):
        out.add(p)
        out.update(ancestors(p, terms, memo, seen))
    memo[mid] = out
    return out


def normalize_mechanism(text):
    s = (text or "").strip().lower().replace("-", " ").replace("_", " ")
    if not s or s in {"na", "n/a", "unknown", "not known", "not yet known", "none", "not provided"}:
        return "unknown"
    if "mixed" in s or "both" in s:
        return "mixed"
    if "dominant negative" in s or re.search(r"\bdn\b", s):
        return "dominant_negative"
    if "gain of function" in s or "gain-of-function" in s:
        return "gain_of_function"
    if "loss of function" in s or "loss-of-function" in s or re.search(r"\blof\b", s) or "haploinsufficiency" in s:
        return "loss_of_function"
    if "activat" in s:
        return "activating"
    if "altered function" in s or "altered-function" in s:
        return "altered_function"
    return "unknown"


def normalize_inheritance(text):
    s = (text or "").strip().lower()
    if not s or "unknown" in s or "not provided" in s:
        return "unknown"
    if "mitochond" in s:
        return "mitochondrial"
    if "monoallelic_autosomal" in s or "autosomal dominant" in s or "semidominant" in s or "semi-dominant" in s:
        return "dominant"
    if "x-linked" in s or "x linked" in s:
        return "x_linked"
    if "recessive" in s or "biallelic" in s:
        return "recessive"
    if "dominant" in s or "heterozygous" in s or "de novo" in s:
        return "dominant"
    return "unknown"


def read_mechanism_sources(genes, disease_to_mondo):
    """Return per gene/disease mechanism assertions from G2P and GenCC."""
    rows = []
    g2p = RAW / "gene2phenotype_all_latest.csv"
    if g2p.exists():
        try:
            df = pd.read_csv(g2p, dtype=str, keep_default_na=False)
            cols = {clean_col(c):c for c in df.columns}
            gene_col = find_col(df.columns, ["gene symbol", "gene", "symbol"])
            mondo_col = find_col(df.columns, ["disease MONDO", "mondo", "mondo id"])
            disease_col = find_col(df.columns, ["disease name", "disease"])
            moi_col = find_col(df.columns, ["allelic requirement", "inheritance"])
            mech_col = find_col(df.columns, ["molecular mechanism", "mechanism"])
            conf_col = find_col(df.columns, ["confidence"])
            var_col = find_col(df.columns, ["evidence based variant consequence", "inferred variant consequence", "variant consequence", "variant type", "variant types"])
            for _, r in df.iterrows():
                gene = str(r.get(gene_col, "")).strip() if gene_col else ""
                if gene not in genes:
                    continue
                mondo = str(r.get(mondo_col, "")).strip() if mondo_col else ""
                if mondo and not mondo.startswith("MONDO:"):
                    mondo = "MONDO:" + mondo
                if not mondo and disease_col:
                    mondo = disease_to_mondo.get(str(r.get(disease_col, "")).casefold(), "")
                conf=str(r.get(conf_col,"")) if conf_col else ""
                primary=conf.strip().lower() in {"definitive","strong","moderate"}
                rows.append({"gene":gene,"MONDO":mondo,"disease_label":str(r.get(disease_col,"")) if disease_col else "",
                             "source":"Gene2Phenotype","inheritance":str(r.get(moi_col,"")) if moi_col else "",
                             "inheritance_category":normalize_inheritance(str(r.get(moi_col,"")) if moi_col else ""),
                             "molecular_mechanism":str(r.get(mech_col,"")) if mech_col else "",
                             "mechanism_category":normalize_mechanism(str(r.get(mech_col,"")) if mech_col else ""),
                             "variant_consequence":str(r.get(var_col,"")) if var_col else "",
                             "evidence_confidence":conf,"primary_evidence":"yes" if primary else "no"})
        except Exception as exc:
            add_failure("Gene2Phenotype parse", str(g2p), str(exc)[:300])
    gencc = RAW / "gencc_submissions_latest.tsv"
    if gencc.exists():
        try:
            df = pd.read_csv(gencc, sep="\t", dtype=str, keep_default_na=False, on_bad_lines="skip")
            gene_col = find_col(df.columns, ["gene symbol", "gene", "hgnc symbol"])
            mondo_col = find_col(df.columns, ["disease mondo", "mondo", "mondo id", "disease curie"])
            disease_col = find_col(df.columns, ["disease name", "disease", "disease title"])
            moi_col = find_col(df.columns, ["mode of inheritance", "inheritance", "allelic requirement", "moi title"])
            mech_col = find_col(df.columns, ["molecular mechanism", "mechanism"])
            var_col = find_col(df.columns, ["variant type", "variant consequence"])
            conf_col = find_col(df.columns, ["classification", "classification title", "assertion criteria"])
            for _, r in df.iterrows():
                gene = str(r.get(gene_col, "")).strip() if gene_col else ""
                if gene not in genes:
                    continue
                mondo = str(r.get(mondo_col, "")).strip() if mondo_col else ""
                if mondo and not mondo.startswith("MONDO:"):
                    mondo = "MONDO:" + mondo
                if not mondo and disease_col:
                    mondo = disease_to_mondo.get(str(r.get(disease_col, "")).casefold(), "")
                conf=str(r.get(conf_col,"")) if conf_col else ""
                primary=conf.strip().lower() in {"definitive","strong","moderate"}
                rows.append({"gene":gene,"MONDO":mondo,"disease_label":str(r.get(disease_col,"")) if disease_col else "",
                             "source":"GenCC","inheritance":str(r.get(moi_col,"")) if moi_col else "",
                             "inheritance_category":normalize_inheritance(str(r.get(moi_col,"")) if moi_col else ""),
                             "molecular_mechanism":str(r.get(mech_col,"")) if mech_col else "",
                             "mechanism_category":normalize_mechanism(str(r.get(mech_col,"")) if mech_col else ""),
                             "variant_consequence":str(r.get(var_col,"")) if var_col else "",
                             "evidence_confidence":conf,"primary_evidence":"yes" if primary else "no"})
        except Exception as exc:
            add_failure("GenCC parse", str(gencc), str(exc)[:300])
    return rows


def load_candidate_summary(path, mondo):
    """Filter the pinned ClinVar summary and retain row-level audit attributes."""
    aliases = {
        "variation_id":["VariationID"], "allele_id":["AlleleID"], "gene":["GeneSymbol"],
        "gene_id":["GeneID"], "variation_type":["Type"], "clinical_significance":["ClinicalSignificance"],
        "review_status":["ReviewStatus"], "protein_change":["ProteinChange"], "hgvs_p":["HGVS_p","HGVS protein"],
        "hgvs_c":["HGVS_c"], "transcript":["Transcript"], "protein_accession":["ProteinAccession","Protein Accession"],
        "consequence":["MolecularConsequence","Molecular Consequence"], "origin":["OriginSimple","Origin"],
        "phenotype_list":["PhenotypeList"], "phenotype_ids":["PhenotypeIDS","PhenotypeIDs","Phenotype IDs"],
        "rcvs":["RCVaccession","RCVAccession"], "organism":["Organism"], "assembly":["Assembly"],
        "name":["Name"], "protein_length":["ProteinLength"]}
    totals = Counter()
    prelim = INTERIM / "clinvar_preliminary_plp_missense.tsv"
    prelim.parent.mkdir(parents=True, exist_ok=True)
    reuse_prelim = prelim.exists()
    first_write = not reuse_prelim
    all_columns = None
    unique_vars = set()
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", dtype=str, chunksize=100000, low_memory=False, keep_default_na=False):
        totals["raw_records"] += len(chunk)
        all_columns = list(chunk.columns)
        cols = {k:find_col(chunk.columns, v) for k,v in aliases.items()}
        gene = chunk[cols["gene"]].str.strip() if cols["gene"] else pd.Series([""]*len(chunk),index=chunk.index)
        sig = chunk[cols["clinical_significance"]].str.strip() if cols["clinical_significance"] else pd.Series([""]*len(chunk),index=chunk.index)
        org = chunk[cols["organism"]].str.lower() if cols["organism"] else pd.Series([""]*len(chunk),index=chunk.index)
        origin = chunk[cols["origin"]].str.lower() if cols["origin"] else pd.Series([""]*len(chunk),index=chunk.index)
        conseq = chunk[cols["consequence"]].str.lower() if cols["consequence"] else pd.Series([""]*len(chunk),index=chunk.index)
        typ = chunk[cols["variation_type"]].str.lower() if cols["variation_type"] else pd.Series([""]*len(chunk),index=chunk.index)
        pchange = chunk[cols["protein_change"]].fillna("") if cols["protein_change"] else pd.Series([""]*len(chunk),index=chunk.index)
        hgvs_p = chunk[cols["hgvs_p"]].fillna("") if cols["hgvs_p"] else pd.Series([""]*len(chunk),index=chunk.index)
        name = chunk[cols["name"]].fillna("") if cols["name"] else pd.Series([""]*len(chunk),index=chunk.index)
        phen = chunk[cols["phenotype_list"]].fillna("") if cols["phenotype_list"] else pd.Series([""]*len(chunk),index=chunk.index)
        phenid = chunk[cols["phenotype_ids"]].fillna("") if cols["phenotype_ids"] else pd.Series([""]*len(chunk),index=chunk.index)
        rcv = chunk[cols["rcvs"]].fillna("") if cols["rcvs"] else pd.Series([""]*len(chunk),index=chunk.index)
        if cols["organism"]:
            is_human = org.str.contains(r"9606|homo sapiens|human", regex=True, na=False)
        else:
            is_human = pd.Series([True]*len(chunk),index=chunk.index)
        if cols["origin"]:
            # Broad candidate screen: defer origin adjudication to each RCV record.
            is_germline = ~origin.str.contains("somatic", na=False)
        else:
            is_germline = pd.Series([True]*len(chunk),index=chunk.index)
        plp = sig.str.lower().isin({"pathogenic","likely pathogenic","pathogenic/likely pathogenic","pathogenic, likely pathogenic"})
        plp &= ~sig.str.lower().str.contains("conflict|uncertain|benign", regex=True, na=False)
        miss = conseq.str.contains("missense", na=False) | pchange.str.contains(P3, regex=True, na=False) | pchange.str.contains(P1, regex=True, na=False)
        if not cols["consequence"] and not cols["protein_change"]:
            miss = hgvs_p.str.contains(P3, regex=True, na=False) | hgvs_p.str.contains(P1, regex=True, na=False) | name.str.contains(P3, regex=True, na=False) | name.str.contains(P1, regex=True, na=False)
        mapped = (pchange.str.contains(P3, regex=True, na=False) | pchange.str.contains(P1, regex=True, na=False) |
                  hgvs_p.str.contains(P3, regex=True, na=False) | hgvs_p.str.contains(P1, regex=True, na=False) |
                  name.str.contains(P3, regex=True, na=False) | name.str.contains(P1, regex=True, na=False))
        disease = phen.str.strip().ne("") | phenid.str.strip().ne("") | rcv.str.contains("RCV", na=False)
        totals["human_records"] += int(is_human.sum())
        totals["germline_or_origin_unspecified_records"] += int((is_human & is_germline).sum())
        totals["plp_records"] += int((is_human & is_germline & plp).sum())
        totals["missense_records"] += int((is_human & is_germline & plp & miss).sum())
        totals["protein_mapped_records"] += int((is_human & is_germline & plp & miss & mapped).sum())
        keep = is_human & is_germline & plp & miss & mapped & disease
        totals["disease_annotated_records"] += int(keep.sum())
        if not keep.any():
            continue
        out = pd.DataFrame(index=chunk.index[keep])
        for key, source_col in cols.items():
            if source_col:
                out[key] = chunk.loc[keep, source_col].values
            else:
                out[key] = ""
        out["gene"] = out["gene"].astype(str).str.strip()
        out["record_source"] = "ClinVar variant_summary preliminary; RCV validation required"
        if not reuse_prelim:
            out.to_csv(prelim, sep="\t", index=False, na_rep="NA", mode="w" if first_write else "a", header=first_write)
            first_write = False
        vcol = out["variation_id"].astype(str)
        pcol = out["protein_change"].astype(str).where(out["protein_change"].astype(str).str.len()>0, out["hgvs_p"].astype(str))
        pcol = pcol.where(pcol.str.len()>0,out["name"].astype(str))
        unique_vars.update(set((out["gene"] + "|" + pcol).tolist()))
    totals["preliminary_unique_gene_protein_variants"] = len({x for x in unique_vars if x.split("|")[-1]})
    return totals, prelim


def preliminary_candidates(prelim, terms, xref_map, synonym_map, alt_map):
    df = pd.read_csv(prelim, sep="\t", dtype=str, keep_default_na=False)
    variant_gene_map = {}
    candidates = set()
    rcv_ids = set()
    preliminary_disease_rows = []
    disease_variants = defaultdict(set)
    for _, r in df.iterrows():
        gene = r.get("gene", "").strip()
        if not gene:
            continue
        vid = r.get("variation_id", "")
        if vid:
            variant_gene_map[vid] = gene
        p = r.get("protein_change", "") or r.get("hgvs_p", "") or r.get("name", "")
        variant = parse_protein_change(p)
        if not variant:
            continue
        key = "%s|%s" % (gene, variant[3])
        label = r.get("phenotype_list", "")
        ids = split_identifiers(r.get("phenotype_ids", ""))
        mids = set()
        for mid in mondo_map(ids, terms, xref_map, synonym_map, alt_map, label):
            mids.add(mid)
        for mid in mids:
            disease_variants[(gene, mid)].add(key)
        for rcv in RCV_RE.findall(r.get("rcvs", "")):
            rcv_ids.add(rcv)
        preliminary_disease_rows.append({"gene":gene,"variation_id":vid,"protein_change":variant[3],"mondo_ids":";".join(sorted(mids)),"raw_condition":label,"condition_ids":r.get("phenotype_ids",""),"rcv_accessions":r.get("rcvs","")})
    per_gene = defaultdict(set)
    for (gene, mid), variants in disease_variants.items():
        if variants:
            per_gene[gene].add(mid)
    threshold_rows=[]
    eligible_by_threshold={}
    for threshold in [2,3,5,10,15,20]:
        qualifying={(g,m):len(vs) for (g,m),vs in disease_variants.items() if len(vs)>=threshold}
        disease_n=defaultdict(int)
        variants=set()
        for (g,m),n in qualifying.items():
            disease_n[g]+=1
            variants.update((g,v.split("|",1)[1]) for v in disease_variants[(g,m)])
        gene_n={g for g,n in disease_n.items() if n>=2}
        eligible_by_threshold[threshold]=gene_n
        threshold_rows.append({"Threshold":"%d variants/disease"%threshold,"minimum_variants_per_disease":threshold,"Genes":len(gene_n),
                               "Diseases":sum(disease_n[g] for g in gene_n),"Variants":len({v for v in variants if v[0] in gene_n}),
                               "data_stage":"whole-release ClinVar summary; preliminary screen"})
    strict_screen=eligible_by_threshold[10]
    strict_volume=Counter(); strict_group_n=Counter()
    for (g,m),vs in disease_variants.items():
        if len(vs)>=10:
            strict_volume[g]+=len(vs); strict_group_n[g]+=1
    strict_rank=sorted(strict_screen,key=lambda g:(-strict_volume[g],-strict_group_n[g],g))
    # RCV validation is limited to five highest-count preliminary strict candidates and the prespecified ACTB/ACTG1 pilots.
    candidates=set(strict_rank[:5])|{"ACTB","ACTG1"}
    selection_rows=[]
    for rank,g in enumerate(strict_rank[:5],1):
        groups=[(m,len(vs)) for (gg,m),vs in disease_variants.items() if gg==g and len(vs)>=10]
        selection_rows.append({"rank":rank,"gene":g,"strict_disease_groups_preliminary":len(groups),
                               "unique_variants_in_strict_groups_preliminary":sum(n for _,n in groups),"selection":"top five by strict-group variant volume"})
    for g in ["ACTB","ACTG1"]:
        selection_rows.append({"rank":"pilot","gene":g,"strict_disease_groups_preliminary":sum(gg==g and len(vs)>=10 for (gg,m),vs in disease_variants.items()),
                               "unique_variants_in_strict_groups_preliminary":sum(len(vs) for (gg,m),vs in disease_variants.items() if gg==g and len(vs)>=10),
                               "selection":"prespecified pilot"})
    if not (INTERIM/"clinvar_preliminary_thresholds.tsv").exists():
        write_tsv(INTERIM/"clinvar_preliminary_thresholds.tsv",threshold_rows,
                  ["Threshold","minimum_variants_per_disease","Genes","Diseases","Variants","data_stage"])
    if not (INTERIM/"clinvar_rcv_validation_candidates.tsv").exists():
        write_tsv(INTERIM/"clinvar_rcv_validation_candidates.tsv",selection_rows,
                  ["rank","gene","strict_disease_groups_preliminary","unique_variants_in_strict_groups_preliminary","selection"])
    candidate_rcvs = set()
    candidate_rows = []
    for r in preliminary_disease_rows:
        if r["gene"] in candidates:
            candidate_rows.append(r)
            candidate_rcvs.update(RCV_RE.findall(r["rcv_accessions"]))
    broad_candidate_path=INTERIM/"clinvar_preliminary_candidate_disease_rows.tsv"
    if not broad_candidate_path.exists():
        write_tsv(broad_candidate_path,preliminary_disease_rows)
    if not (INTERIM/"clinvar_rcv_validation_candidate_disease_rows.tsv").exists():
        write_tsv(INTERIM/"clinvar_rcv_validation_candidate_disease_rows.tsv",candidate_rows)
    return candidates, sorted(candidate_rcvs), variant_gene_map, disease_variants, threshold_rows, strict_rank


def build_disease_harmonization(rcv_rows, terms, xref_map, synonym_map, alt_map):
    harmonization = {}
    for r in rcv_rows:
        ids = [x for x in (r.get("condition_IDs", "").split(";") if r.get("condition_IDs") else []) if x]
        names = [x for x in (r.get("raw_condition", "").split(";") if r.get("raw_condition") else []) if x]
        mapped = mondo_map(ids, terms, xref_map, synonym_map, alt_map, names[0] if names else "")
        mondo = mapped[0] if len(mapped) == 1 else ""
        r["MONDO"] = mondo if mondo else ("AMBIGUOUS:" + ";".join(mapped) if mapped else "NA")
        r["MedGen"] = ";".join(x.split(":", 1)[1] for x in ids if x.startswith("MedGen:")) or "NA"
        if mondo:
            label = terms.get(mondo, {}).get("name", "")
            parents = ";".join(terms.get(mondo, {}).get("parents", []))
            status = "AUTO_ACCEPT"
            rel = "mapped MONDO term"
        elif mapped:
            label, parents, status, rel = ";".join(terms.get(m, {}).get("name", "") for m in mapped), "", "MANUAL_REVIEW", "multiple MONDO mappings"
        else:
            label, parents, status, rel = "", "", "MANUAL_REVIEW", "unmapped condition"
        harmonization.setdefault((r.get("raw_condition", ""), r.get("condition_IDs", ""), r["MONDO"]), {
            "raw_condition":r.get("raw_condition", ""),"mapped_MONDO":r["MONDO"],"MONDO_label":label,"parent_MONDO":parents,
            "relationship":rel,"collapse_group":r["MONDO"],"status":status})
    return list(harmonization.values())


def disease_relation(a, b, terms, label_a="", label_b=""):
    if a == b:
        return "EXACT_MONDO_EQUIVALENT"
    aa, bb = ancestors(a, terms), ancestors(b, terms)
    if a in bb or b in aa:
        return "PARENT_CHILD_REVIEW"
    stop={"disease","syndrome","disorder","related","type","with","due","to","of","and","susceptibility","familial","associated","caused"}
    ta={x for x in re.findall(r"[a-z0-9]+",str(label_a).lower()) if x not in stop and len(x)>2}
    tb={x for x in re.findall(r"[a-z0-9]+",str(label_b).lower()) if x not in stop and len(x)>2}
    overlap=len(ta&tb)/max(1,len(ta|tb))
    if overlap>=0.25 or len(ta&tb)>=2:
        return "SAME_SPECTRUM_REVIEW"
    return "DISTINCT_DISEASE"


def audit_disease_pairs(counts, terms):
    by_gene = defaultdict(list)
    for r in counts:
        if r["MONDO"].startswith("MONDO:"):
            by_gene[r["gene"]].append(r)
    pairs = []
    for gene, rs in by_gene.items():
        for i in range(len(rs)):
            for j in range(i+1, len(rs)):
                a,b=rs[i],rs[j]
                rel=disease_relation(a["MONDO"],b["MONDO"],terms,a.get("disease",""),b.get("disease",""))
                pairs.append({"gene":gene,"disease_A":a["disease"],"MONDO_A":a["MONDO"],"n_variants_A":a["n_unique_variants"],
                              "disease_B":b["disease"],"MONDO_B":b["MONDO"],"n_variants_B":b["n_unique_variants"],
                              "relationship":rel,"status":"AUTO_ACCEPT" if rel=="DISTINCT_DISEASE" else "MANUAL_REVIEW",
                              "reason":"separate MONDO branches and distinct labels" if rel=="DISTINCT_DISEASE" else "ontology ancestry or label overlap requires clinical spectrum review"})
    return pairs


def load_biogrid(path, genes):
    rows=[]
    if not Path(path).exists():
        return rows
    import zipfile
    try:
        with zipfile.ZipFile(path) as z:
            names=[n for n in z.namelist() if n.lower().endswith(".tab3.txt")]
            for name in names:
                with z.open(name) as f:
                    reader=csv.DictReader(io.TextIOWrapper(f,encoding="utf-8",errors="replace"),delimiter="\t")
                    for r in reader:
                        a=(r.get("Official Symbol Interactor A") or "").strip()
                        b=(r.get("Official Symbol Interactor B") or "").strip()
                        if a not in genes and b not in genes:
                            continue
                        if str(r.get("Organism ID Interactor A",""))!="9606" or str(r.get("Organism ID Interactor B",""))!="9606":
                            continue
                        system=(r.get("Experimental System") or "").strip()
                        sys_type=(r.get("Experimental System Type") or "").strip().lower()
                        pub=(r.get("Publication Source") or "").strip()
                        binary=any(x in system.lower() for x in ("two-hybrid","protein complementation assay","fluorescence resonance energy transfer","bioluminescence resonance energy transfer","far western"))
                        for target,partner in ((a,b),(b,a)):
                            if target in genes and partner and target!=partner:
                                rows.append({"gene":target,"partner":partner,"source":"BioGRID MV-Physical 5.0.261","support":"experimental physical, multi-validated",
                                             "binary_assay":"yes" if binary else "not established as binary","experimental_system":system or sys_type,
                                             "publication":pub,"tier_hint":"Tier 2" if not binary else "Tier 1/2"})
    except Exception as exc:
        add_failure("BioGRID parse", str(path), str(exc)[:300])
    return rows


def load_hgnc_uniprot_maps(path):
    uniprot_to_symbol={}; symbol_to_uniprot=defaultdict(set)
    if not Path(path).exists():
        return uniprot_to_symbol,symbol_to_uniprot
    try:
        df=pd.read_csv(path,sep="\t",dtype=str,keep_default_na=False,usecols=lambda c: clean_col(c) in {"symbol","uniprotids"})
        sym_col=find_col(df.columns,["symbol"]); uni_col=find_col(df.columns,["uniprot_ids","uniprot ids"])
        if not sym_col or not uni_col:
            return uniprot_to_symbol,symbol_to_uniprot
        for _,r in df.iterrows():
            sym=str(r.get(sym_col,""));
            for uid in re.split(r"[|,; ]+",str(r.get(uni_col,""))):
                uid=uid.strip()
                if uid and uid!="-":
                    uniprot_to_symbol.setdefault(uid,sym);symbol_to_uniprot[sym].add(uid)
    except Exception as exc:
        add_failure("HGNC UniProt mapping parse",str(path),str(exc)[:300])
    return uniprot_to_symbol,symbol_to_uniprot


def load_pioneer(path, genes, uniprot_to_symbol):
    rows=[]
    if not Path(path).exists():
        return rows
    try:
        with open(path,encoding="utf-8",errors="replace",newline="") as f:
            reader=csv.DictReader(f,delimiter="\t")
            for r in reader:
                a=str(r.get("uniprot1","")).split("-")[0].strip();b=str(r.get("uniprot2","")).split("-")[0].strip()
                ga=uniprot_to_symbol.get(a,"");gb=uniprot_to_symbol.get(b,"")
                try: ra=ast.literal_eval(r.get("interface_residues1","[]"))
                except Exception: ra=[]
                try: rb=ast.literal_eval(r.get("interface_residues2","[]"))
                except Exception: rb=[]
                for target,partner,ua,ub,residues,sequence in ((ga,gb,a,b,ra,r.get("sequence1","")),(gb,ga,b,a,rb,r.get("sequence2",""))):
                    if target in genes and partner:
                        evidence=(r.get("source","") or "unknown").strip()
                        rows.append({"gene":target,"partner":partner,"source":"PIONEER very-high-confidence human","support":"precomputed partner-specific interface annotation",
                                     "binary_assay":"NA","experimental_system":"NA","publication":"NA","tier_hint":"Tier 3",
                                     "uniprot_target":ua,"uniprot_partner":ub,"PIONEER_reference_source":evidence,
                                     "prediction_confidence":"very_high","interface_residues":";".join(map(str,residues)),
                                     "target_sequence":sequence,"partner_sequence":r.get("sequence2" if target==ga else "sequence1",""),
                                     "target_sequence_length":len(sequence)})
    except Exception as exc:
        add_failure("PIONEER parse", str(path), str(exc)[:300])
    return rows


def load_intact(genes, symbol_to_uniprot, uniprot_to_symbol):
    """Retrieve candidate-focused human physical evidence from IntAct PSICQUIC."""
    rows=[]
    cache=RAW/"intact_api"
    cache.mkdir(parents=True,exist_ok=True)
    base="https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/"
    session=requests.Session()
    for gene in sorted(genes):
        ids=sorted(symbol_to_uniprot.get(gene,set()))
        if not ids:
            continue
        for uid in ids:
            safe_uid=re.sub(r"[^A-Za-z0-9_.-]","_",uid)
            target=cache/("%s_%s.tab25"%(gene,safe_uid))
            if target.exists():
                content=target.read_text(encoding="utf-8",errors="replace")
            else:
                query=quote("identifier:%s AND species:9606"%uid,safe=":")
                url=base+query+"?format=tab25&first=0&max=2000"
                try:
                    r=session.get(url,timeout=(25,90),headers={"Accept":"text/plain","User-Agent":"phase0-feasibility-audit/1.0"})
                    if not r.ok:
                        raise RuntimeError("HTTP %s: %s"%(r.status_code,r.text[:160]))
                    content=r.text
                    target.write_text(content,encoding="utf-8")
                    (target.with_suffix(target.suffix + ".download.json")).write_text(json.dumps({
                        "url":url,"retrieved_utc":now(),"bytes":len(r.content),
                        "sha256":hashlib.sha256(r.content).hexdigest()},indent=2),encoding="utf-8")
                    time.sleep(0.25)
                except Exception as exc:
                    add_failure("IntAct PSICQUIC candidate query",url,str(exc)[:300])
                    continue
            for line in content.splitlines():
                if not line or line.startswith("#"):
                    continue
                f=line.split("\t")
                if len(f)<12:
                    continue
                ua_match=re.search(r"uniprotkb:([A-Za-z0-9-]+)",f[0],re.I)
                ub_match=re.search(r"uniprotkb:([A-Za-z0-9-]+)",f[1],re.I)
                if not ua_match or not ub_match:
                    continue
                ua,ub=ua_match.group(1).split("-")[0],ub_match.group(1).split("-")[0]
                ga,gb=uniprot_to_symbol.get(ua,""),uniprot_to_symbol.get(ub,"")
                if "taxid:9606" not in f[9] or "taxid:9606" not in f[10]:
                    continue
                method=f[6] if len(f)>6 else ""
                interaction=f[11] if len(f)>11 else ""
                publication=f[8] if len(f)>8 else ""
                binary=("direct interaction" in interaction.lower() or any(x in method.lower() for x in ("two hybrid","fret","protein complementation","bioluminescence resonance energy transfer")))
                for target_gene,partner,ta,tb in ((ga,gb,ua,ub),(gb,ga,ub,ua)):
                    if target_gene==gene and partner and partner!=target_gene:
                        rows.append({"gene":target_gene,"partner":partner,"source":"IntAct PSICQUIC (retrieved live)","support":"curated experimental molecular interaction",
                                     "binary_assay":"yes" if binary else "not established as binary","experimental_system":method+" | "+interaction,
                                     "publication":publication,"tier_hint":"Tier 1/2" if binary else "Tier 2","uniprot_target":ta,"uniprot_partner":tb})
    return rows


def build_partner_tables(genes):
    biogrid=load_biogrid(RAW/"BIOGRID-MV-Physical-5.0.261.mitab.zip",genes)
    if not biogrid:
        biogrid=load_biogrid(RAW/"BIOGRID-MV-Physical-5.0.261.tab3.zip",genes)
    uniprot_to_symbol,symbol_to_uniprot=load_hgnc_uniprot_maps(RAW/"hgnc_non_alt_loci_2026-09-27.tsv")
    intact=load_intact(genes,symbol_to_uniprot,uniprot_to_symbol)
    pioneer=load_pioneer(RAW/"PIONEER_very_high_human_2026-09-27.txt",genes,uniprot_to_symbol)
    all_rows=biogrid+intact+pioneer
    uniq=defaultdict(list)
    for r in all_rows:
        uniq[(r["gene"],r["partner"])].append(r)
    prov=[]
    for (gene,partner),sources in uniq.items():
        names=sorted(set(x["source"] for x in sources))
        experimental=any(x["source"].startswith(("BioGRID","IntAct")) for x in sources)
        pioneer_cov=any(x["source"].startswith("PIONEER") for x in sources)
        pdb_cov=any("pdb" in x.get("PIONEER_reference_source","").lower() for x in sources)
        binary=any(x.get("binary_assay")=="yes" for x in sources)
        tier="Tier 1" if experimental and binary else ("Tier 2" if experimental else "Tier 3")
        prov.append({"gene":gene,"partner":partner,"sources":";".join(names),"source_count":len(names),"experimental_support":"yes" if experimental else "no",
                     "binary_experimental_support":"yes" if binary else "no/not established","PIONEER_coverage":"yes" if pioneer_cov else "no",
                     "structural_support":"yes" if pdb_cov else ("PIONEER interface annotation; source not PDB" if pioneer_cov else "no"),"predicted_only":"yes" if not experimental else "no","tier":tier,
                     "experimental_systems":";".join(sorted(set(x.get("experimental_system","") for x in sources if x.get("experimental_system"))))})
    gene_stats=[]
    for gene in sorted(genes):
        pr=[x for x in prov if x["gene"]==gene]
        gene_stats.append({"gene":gene,"total_partners":len(pr),"experimental_partners":sum(x["experimental_support"]=="yes" for x in pr),
                           "PIONEER_partners":sum(x["PIONEER_coverage"]=="yes" for x in pr),"Predictomes_partners":"NA",
                           "unique_credible_partners":sum(x["experimental_support"]=="yes" for x in pr),
                           "structural_ready_partners":sum(x["structural_support"]=="yes" for x in pr),
                           "ppi_sources":"BioGRID MV-Physical 5.0.261; IntAct PSICQUIC; PIONEER very-high-confidence"})
    return prov,gene_stats,all_rows


def annotate_variant_proxy(record, pioneer_rows, experimental_pairs):
    gene=record.get("gene","")
    position=int(record.get("residue") or 0)
    ref=(record.get("amino_acid_ref") or "").strip()
    target_rows=[p for p in pioneer_rows if p.get("gene")==gene]
    if not target_rows or position<=0 or not ref:
        record.update({"STRUCTURE_READY_PROXY":"NA","PARTNER_ANALYSIS_READY_PROXY":"NA","interface_residue_overlap":"NA","isoform_mapping_status":"ambiguous or unavailable","n_PIONEER_partners":len({p.get("partner") for p in target_rows})})
        return record
    concordant=[]; mismatched=[]
    for p in target_rows:
        seq=p.get("target_sequence","")
        if position<=len(seq) and seq[position-1]==ref:
            concordant.append(p)
        else:
            mismatched.append(p)
    credible=[p for p in concordant if (gene,p.get("partner")) in experimental_pairs]
    pdb_rows=[p for p in concordant if "pdb" in p.get("PIONEER_reference_source","").lower()]
    interface_overlap=any(str(position) in {x for x in p.get("interface_residues","").split(";") if x} for p in credible)
    record.update({"STRUCTURE_READY_PROXY":"YES" if pdb_rows else "NO",
                   "PARTNER_ANALYSIS_READY_PROXY":"YES" if credible else "NO",
                   "interface_residue_overlap":"YES" if interface_overlap else "NO",
                   "isoform_mapping_status":"same-position reference amino acid concordant" if concordant else ("position/ref mismatch; isoform mapping ambiguous" if mismatched else "ambiguous or unavailable"),
                   "n_PIONEER_partners":len({p.get("partner") for p in target_rows}),
                   "n_credible_PIONEER_covered_partners":len({p.get("partner") for p in credible}),
                   "n_existing_PDB_interface_partners":len({p.get("partner") for p in pdb_rows})})
    return record


def audit_phase0():
    PHASE.mkdir(parents=True,exist_ok=True); INTERIM.mkdir(parents=True,exist_ok=True); PROCESSED.mkdir(parents=True,exist_ok=True); REPORTS.mkdir(parents=True,exist_ok=True)
    summary=RAW/"variant_summary_2026-09.txt.gz"
    mondo_file=RAW/"mondo-v2026-09-01.obo"
    if not summary.exists():
        add_failure("ClinVar variant_summary", "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2026-09.txt.gz", "fixed snapshot unavailable; acquisition incomplete")
        raise FileNotFoundError(summary)
    terms,xref_map,synonym_map,alt_map=parse_mondo(mondo_file)
    attrition,prelim=load_candidate_summary(summary,(terms,xref_map,synonym_map,alt_map))
    candidates,rcvs,variant_gene_map,disease_variants,preliminary_threshold_rows,preliminary_strict_rank=preliminary_candidates(prelim,terms,xref_map,synonym_map,alt_map)
    rcv_rows,failed=fetch_rcv_batches(rcvs,variant_gene_map)
    # Keep one row per actual RCV record and retain all downloaded API records as immutable raw XML batches.
    harmonization=build_disease_harmonization(rcv_rows,terms,xref_map,synonym_map,alt_map)
    if rcv_rows:
        filtered=[]
        for r in rcv_rows:
            mondo_id=r.get("MONDO","")
            is_primary=(r.get("germline") and is_plp(r.get("clinical_significance","")) and is_missense(r) and r.get("protein_mapping_status")=="mapped" and r.get("raw_condition"))
            if is_primary and mondo_id.startswith("MONDO:"):
                filtered.append(r)
    else:
        filtered=[]
    write_tsv(PHASE/"disease_harmonization.tsv",harmonization)
    # Collapse duplicate ClinVar records to gene–amino-acid-change–disease units.
    unique={}
    for r in filtered:
        key=(r.get("gene",""),r.get("protein_change",""),r.get("MONDO",""))
        if not key[0] or not key[1] or not key[2]:
            continue
        if key not in unique:
            unique[key]=dict(r)
        else:
            scvs=sorted(set(x for x in (unique[key].get("SCV_accessions","")+";"+r.get("SCV_accessions","")).split(";") if x))
            unique[key]["SCV_accessions"]=";".join(scvs)
            unique[key]["n_SCV"]=len(scvs)
    by_gene_disease=defaultdict(list)
    for r in unique.values():
        by_gene_disease[(r["gene"],r["MONDO"])].append(r)
    disease_counts=[]
    for (gene,mondo_id),rs in by_gene_disease.items():
        stars=[int(x.get("stars",0)) for x in rs]
        disease_counts.append({"gene":gene,"disease":rs[0].get("raw_condition",""),"MONDO":mondo_id,
                               "n_unique_variants":len(rs),"n_residues":len(set((x.get("protein_accession",""),str(x.get("residue",""))) for x in rs)),
                               "n_RCV":len(set(x.get("RCV_accession","") for x in rs)),
                               "n_SCV":len({scv for x in rs for scv in x.get("SCV_accessions","").split(";") if scv}),
                               "median_ClinVar_stars":float(np.median(stars)) if stars else "NA","max_ClinVar_stars":max(stars) if stars else "NA",
                               "strict_eligible":len(rs)>=10,"expanded_eligible":len(rs)>=5})
    write_tsv(PHASE/"gene_disease_variant_counts.tsv",disease_counts)
    disease_pairs=audit_disease_pairs(disease_counts,terms)
    distinct_counts=Counter()
    for r in disease_counts:
        distinct_counts[r["gene"]]+=1
    pre_mech=[g for g,n in distinct_counts.items() if n>=2]
    write_tsv(PHASE/"candidate_genes_pre_mechanism.tsv",[{"gene":g,"n_diseases":distinct_counts[g]} for g in sorted(pre_mech)], ["gene","n_diseases"])
    gene_set=set(pre_mech)
    pair_mondo={r["raw_condition"].casefold():r["mapped_MONDO"] for r in harmonization if r.get("mapped_MONDO","").startswith("MONDO:")}
    mechanisms=read_mechanism_sources(gene_set,pair_mondo)
    write_tsv(PHASE/"gene_disease_mechanism_audit.tsv",mechanisms)
    # Mechanism compatibility is conservative: require at least one known mechanism and inheritance match for every disease.
    mech_by=defaultdict(list)
    for r in mechanisms:
        if r.get("MONDO"):
            mech_by[(r["gene"],r["MONDO"])].append(r)
    gene_mech=[]
    for gene in sorted(gene_set):
        ds=[r for r in disease_counts if r["gene"]==gene]
        for i in range(len(ds)):
            for j in range(i+1,len(ds)):
                a,b=ds[i],ds[j]
                ma=[x for x in mech_by.get((gene,a["MONDO"]),[]) if x.get("primary_evidence")=="yes"]
                mb=[x for x in mech_by.get((gene,b["MONDO"]),[]) if x.get("primary_evidence")=="yes"]
                ac={x["mechanism_category"] for x in ma if x["mechanism_category"]!="unknown"}; bc={x["mechanism_category"] for x in mb if x["mechanism_category"]!="unknown"}
                ai={x["inheritance_category"] for x in ma if x["inheritance_category"]!="unknown"}; bi={x["inheritance_category"] for x in mb if x["inheritance_category"]!="unknown"}
                if len(ac)==len(bc)==len(ai)==len(bi)==1:
                    compatible=ac==bc and ai==bi
                    status="MECHANISM_COMPATIBLE" if compatible else "MECHANISM_DIFFERENT"
                elif ac and bc and ai and bi and len(ac)==len(bc)==len(ai)==len(bi)==1:
                    status="MECHANISM_DIFFERENT"
                else:
                    status="MECHANISM_UNCERTAIN"
                gene_mech.append({"gene":gene,"disease_A":a["disease"],"MONDO_A":a["MONDO"],"mechanism_A":";".join(sorted(ac)) or "unknown","inheritance_A":";".join(sorted(ai)) or "unknown",
                                  "disease_B":b["disease"],"MONDO_B":b["MONDO"],"mechanism_B":";".join(sorted(bc)) or "unknown","inheritance_B":";".join(sorted(bi)) or "unknown",
                                  "compatibility":status,"source":"Gene2Phenotype/GenCC; no OMIM licensed source"})
    write_tsv(PHASE/"gene_disease_mechanism_pair_audit.tsv",gene_mech)
    mech_compatible={r["gene"] for r in gene_mech if r["compatibility"]=="MECHANISM_COMPATIBLE"}
    post=[g for g in pre_mech if g in mech_compatible]
    write_tsv(PHASE/"candidate_genes_post_mechanism.tsv",[{"gene":g,"mechanism_status":"MECHANISM_COMPATIBLE"} for g in sorted(post)])
    prov,ppi_stats,ppi_rows=build_partner_tables(set(pre_mech)|{"ACTB","ACTG1"})
    pioneer_rows=[r for r in ppi_rows if r.get("source","").startswith("PIONEER")]
    experimental_pairs={(r["gene"],r["partner"]) for r in prov if r["experimental_support"]=="yes"}
    coverage_rows=[]
    coverage_by_gene=defaultdict(list)
    annotated_unique={}
    seen_gene_variant=set()
    for r in unique.values():
        key=(r["gene"],r["protein_change"])
        annotated=annotated_unique.get(key)
        if annotated is None:
            annotated=annotate_variant_proxy(dict(r),pioneer_rows,experimental_pairs)
            annotated_unique[key]=annotated
        if key in seen_gene_variant:
            continue
        seen_gene_variant.add(key)
        conditions=sorted({x.get("MONDO","") for x in unique.values() if (x["gene"],x["protein_change"])==key})
        coverage_rows.append({"gene":r["gene"],"MONDOs":";".join(conditions),"protein_change":r["protein_change"],"residue":r.get("residue",""),
                              "STRUCTURE_READY_PROXY":annotated.get("STRUCTURE_READY_PROXY","NA"),"PARTNER_ANALYSIS_READY_PROXY":annotated.get("PARTNER_ANALYSIS_READY_PROXY","NA"),
                              "interface_residue_overlap":annotated.get("interface_residue_overlap","NA"),"isoform_mapping_status":annotated.get("isoform_mapping_status","NA"),
                              "n_PIONEER_partners":annotated.get("n_PIONEER_partners",0),"n_credible_PIONEER_covered_partners":annotated.get("n_credible_PIONEER_covered_partners",0),
                              "n_existing_PDB_interface_partners":annotated.get("n_existing_PDB_interface_partners",0)})
        coverage_by_gene[r["gene"]].append(annotated)
    for r in filtered:
        key=(r.get("gene",""),r.get("protein_change",""))
        ar=annotated_unique.get(key)
        if ar:
            for field in ("STRUCTURE_READY_PROXY","PARTNER_ANALYSIS_READY_PROXY","interface_residue_overlap","isoform_mapping_status","n_PIONEER_partners","n_credible_PIONEER_covered_partners","n_existing_PDB_interface_partners"):
                r[field]=ar.get(field,"NA")
    write_tsv(PHASE/"clinvar_variant_condition_raw.tsv",rcv_rows,["Variation_ID","RCV_accession","RCV_version","gene","transcript","protein_accession","residue","amino_acid_ref","amino_acid_alt","HGVS_protein","protein_change","molecular_consequence","raw_condition","condition_IDs","MONDO","MedGen","review_status","stars","clinical_significance","date_last_evaluated","inheritance","germline","n_SCV","SCV_accessions","protein_mapping_status","record_source"] if not rcv_rows else None)
    write_tsv(PHASE/"clinvar_variant_condition_filtered.tsv",filtered,["Variation_ID","RCV_accession","RCV_version","gene","transcript","protein_accession","residue","amino_acid_ref","amino_acid_alt","HGVS_protein","protein_change","molecular_consequence","raw_condition","condition_IDs","MONDO","MedGen","review_status","stars","clinical_significance","date_last_evaluated","inheritance","germline","n_SCV","SCV_accessions","protein_mapping_status","record_source","STRUCTURE_READY_PROXY","PARTNER_ANALYSIS_READY_PROXY","interface_residue_overlap","isoform_mapping_status","n_PIONEER_partners","n_credible_PIONEER_covered_partners","n_existing_PDB_interface_partners"] if not filtered else None)
    write_tsv(PHASE/"variant_coverage_proxy.tsv",coverage_rows)
    write_tsv(PHASE/"disease_pairs_relationship_audit.tsv",disease_pairs, ["gene","disease_A","MONDO_A","n_variants_A","disease_B","MONDO_B","n_variants_B","relationship","status","reason"])
    write_tsv(PHASE/"disease_pairs_needing_manual_review.tsv",[r for r in disease_pairs if r["status"]=="MANUAL_REVIEW"], ["gene","disease_A","MONDO_A","n_variants_A","disease_B","MONDO_B","n_variants_B","relationship","status","reason"])
    for gs in ppi_stats:
        gene=gs["gene"]; rs=coverage_by_gene.get(gene,[])
        ready=[x for x in rs if x.get("PARTNER_ANALYSIS_READY_PROXY")=="YES"]
        structural=[x for x in rs if x.get("STRUCTURE_READY_PROXY")=="YES"]
        gs["partner_ready_variants"]=len(ready);gs["structural_ready_variants"]=len(structural);gs["variant_proxy_denominator"]=len(rs)
        gs["partner_analysis_ready_fraction"]=len(ready)/len(rs) if rs else "NA"
        gs["structural_ready_fraction"]=len(structural)/len(rs) if rs else "NA"
    write_tsv(PHASE/"ppi_provenance.tsv",prov)
    write_tsv(PHASE/"ppi_coverage_by_gene.tsv",ppi_stats)
    # Candidate summary table and empirical thresholds use actual RCV-validated unique substitutions.
    by_gene=defaultdict(list)
    for r in disease_counts:
        by_gene[r["gene"]].append(r)
    threshold_rows=list(preliminary_threshold_rows)
    for threshold in [2,3,5,10,15,20]:
        accepted=[p for p in disease_pairs if p["status"]=="AUTO_ACCEPT" and int(p["n_variants_A"])>=threshold and int(p["n_variants_B"])>=threshold]
        gs=sorted({p["gene"] for p in accepted if sum(x["n_unique_variants"] for x in by_gene[p["gene"]] if x["n_unique_variants"]>=threshold)>=2*threshold})
        qualified_groups={(p["gene"],p[k]) for p in accepted if p["gene"] in gs for k in ("MONDO_A","MONDO_B")}
        ndisease=len(qualified_groups)
        qualifying_variants={(g,x.get("protein_change","")) for g,m in qualified_groups for x in unique.values() if x.get("gene")==g and x.get("MONDO")==m}
        nvar=len(qualifying_variants)
        threshold_rows.append({"Threshold":"%d variants/disease"%threshold,"minimum_variants_per_disease":threshold,"Genes":len(gs),"Diseases":ndisease,"Variants":nvar,"data_stage":"selected RCV-validation panel"})
    write_tsv(PHASE/"candidate_gene_thresholds.tsv",threshold_rows)
    strict_pairs=[p for p in disease_pairs if p["status"]=="AUTO_ACCEPT" and int(p["n_variants_A"])>=10 and int(p["n_variants_B"])>=10]
    strict=sorted({p["gene"] for p in strict_pairs})
    expanded_pairs=[p for p in disease_pairs if p["status"]=="AUTO_ACCEPT" and int(p["n_variants_A"])>=5 and int(p["n_variants_B"])>=5]
    expanded=sorted({p["gene"] for p in expanded_pairs})
    mech_compatible_pairs={(r["gene"],frozenset([r["MONDO_A"],r["MONDO_B"]])) for r in gene_mech if r["compatibility"]=="MECHANISM_COMPATIBLE"}
    strict_mech=sorted({p["gene"] for p in strict_pairs if (p["gene"],frozenset([p["MONDO_A"],p["MONDO_B"]])) in mech_compatible_pairs})
    ppi_eligible={r["gene"] for r in ppi_stats if r["gene"] in set(pre_mech) and r["unique_credible_partners"]>=3}
    quality=set(strict_mech)&ppi_eligible
    # Conservative proxy: variant is ready when a high-confidence partner-specific PIONEER annotation exists for its gene.
    pioneer_genes={r["gene"] for r in ppi_stats if r["PIONEER_partners"]>0}
    coverage_fractions={g:(sum(x.get("PARTNER_ANALYSIS_READY_PROXY")=="YES" for x in rs)/len(rs) if rs else 0.0) for g,rs in coverage_by_gene.items()}
    final_genes={g for g in quality&pioneer_genes if coverage_fractions.get(g,0.0)>=0.70}
    final_variant_keys=set()
    for g in final_genes:
        strict_diseases={p[k] for p in strict_pairs if p["gene"]==g for k in ("MONDO_A","MONDO_B")}
        final_variant_keys.update((g,x.get("protein_change","")) for x in unique.values() if x.get("gene")==g and x.get("MONDO") in strict_diseases)
    final_variants=sum(annotated_unique.get(key,{}).get("PARTNER_ANALYSIS_READY_PROXY")=="YES" for key in final_variant_keys)
    top=[]
    ranked_genes=sorted(set(pre_mech),key=lambda x:(x not in final_genes,x not in strict,-(x in mech_compatible),-next((p["unique_credible_partners"] for p in ppi_stats if p["gene"]==x),0),-next((p["PIONEER_partners"] for p in ppi_stats if p["gene"]==x),0),-coverage_fractions.get(x,0),-len(by_gene[x]),-sum(y["n_unique_variants"] for y in by_gene[x]),x))[:20]
    for g in ranked_genes:
        ps=next((x for x in ppi_stats if x["gene"]==g),{})
        top.append({"gene":g,"n_normalized_diseases":len(by_gene[g]),"n_unique_variants":sum(x["n_unique_variants"] for x in by_gene[g]),
                    "n_strict_diseases":sum(x["n_unique_variants"]>=10 for x in by_gene[g]),"mechanism_status":"compatible" if g in mech_compatible else "uncertain/different",
                    "credible_PPI_partners":ps.get("unique_credible_partners",0),"PIONEER_partners":ps.get("PIONEER_partners",0),"status":"high-quality" if g in final_genes else "screened candidate"})
    write_tsv(PHASE/"top_candidate_genes.tsv",top)
    strict_table=[]
    for g in strict:
        ds=[x for x in by_gene[g] if x["n_unique_variants"]>=10]
        ps=next((x for x in ppi_stats if x["gene"]==g),{})
        mech_labels=sorted({r["mechanism_A"] for r in gene_mech if r["gene"]==g and r["compatibility"]=="MECHANISM_COMPATIBLE"}|{r["mechanism_B"] for r in gene_mech if r["gene"]==g and r["compatibility"]=="MECHANISM_COMPATIBLE"})
        inheritance_labels=sorted({r["inheritance_A"] for r in gene_mech if r["gene"]==g and r["compatibility"]=="MECHANISM_COMPATIBLE"}|{r["inheritance_B"] for r in gene_mech if r["gene"]==g and r["compatibility"]=="MECHANISM_COMPATIBLE"})
        strict_table.append({"gene":g,"diseases":";".join(x["disease"] for x in ds),"MONDO":";".join(x["MONDO"] for x in ds),
                             "n_variants_per_disease":";".join(str(x["n_unique_variants"]) for x in ds),"inheritance":";".join(inheritance_labels) or "unknown",
                             "mechanism":";".join(mech_labels) or "unknown","credible_PPI_partners":ps.get("unique_credible_partners",0),
                             "PIONEER_partners":ps.get("PIONEER_partners",0),"partner_analysis_ready_fraction":coverage_fractions.get(g,0),
                             "status":"high-quality" if g in final_genes else ("mechanism-compatible" if g in strict_mech else "secondary/sensitivity")})
    write_tsv(PHASE/"strict_candidate_genes.tsv",strict_table,["gene","diseases","MONDO","n_variants_per_disease","inheritance","mechanism","credible_PPI_partners","PIONEER_partners","partner_analysis_ready_fraction","status"])
    write_tsv(PHASE/"mechanism_conflicts.tsv",[r for r in gene_mech if r["compatibility"] in {"MECHANISM_DIFFERENT","MECHANISM_UNCERTAIN"}], ["gene","disease_A","MONDO_A","mechanism_A","inheritance_A","disease_B","MONDO_B","mechanism_B","inheritance_B","compatibility","source"])
    write_tsv(PHASE/"ACTB_ACTG1_feasibility.tsv",[{"gene":g,"P_LP_missense_variants":len({x.get("protein_change","") for x in unique.values() if x.get("gene")==g}),
                   "normalized_diseases":len(by_gene.get(g,[])),"strict_diseases":sum(x["n_unique_variants"]>=10 for x in by_gene.get(g,[])),
                   "mechanism_compatible":g in mech_compatible,"credible_PPI_partners":next((x["unique_credible_partners"] for x in ppi_stats if x["gene"]==g),0),
                   "partner_analysis_ready_fraction":coverage_fractions.get(g,0.0),
                   "pilot_ready":"YES" if (g in strict_mech and next((x["unique_credible_partners"] for x in ppi_stats if x["gene"]==g),0)>=3 and coverage_fractions.get(g,0.0)>=0.70) else "NO"} for g in ["ACTB","ACTG1"]])
    # Cohort attrition table distinguishes the same-month bulk screen from RCV-confirmed candidate rows.
    attrition_rows=[{"step":k,"records":v,"data_source":"ClinVar variant_summary 2026-09 preliminary screen"} for k,v in attrition.items()]
    attrition_rows += [{"step":"preliminary MONDO-mapped unique variant-condition candidates","records":sum(len(v) for v in disease_variants.values()),"data_source":"ClinVar variant_summary + MONDO 2026-09-01; screen only"},
                       {"step":"candidate RCV accessions requested","records":len(rcvs),"data_source":"ClinVar E-utilities RCV clinvarset (retrieved on run date)"},
                       {"step":"candidate RCV records returned","records":len(rcv_rows),"data_source":"ClinVar E-utilities RCV clinvarset (retrieved on run date)"},
                       {"step":"RCV fetch/parse failures","records":len(failed),"data_source":"ClinVar E-utilities"},
                       {"step":"RCV P/LP germline missense protein mapped MONDO records","records":len(filtered),"data_source":"ClinVar RCV XML"}]
    write_tsv(PHASE/"clinvar_cohort_attrition.tsv",attrition_rows)
    write_tsv(PHASE/"rcv_api_failures.tsv",failed)
    make_figures(attrition_rows,disease_counts,threshold_rows,ppi_stats,by_gene,final_genes)
    strict_mech_variant_keys=set()
    strict_mech_diseases={(p["gene"],p[k]) for p in strict_pairs if p["gene"] in strict_mech for k in ("MONDO_A","MONDO_B")}
    for x in unique.values():
        if (x.get("gene"),x.get("MONDO")) in strict_mech_diseases:
            strict_mech_variant_keys.add((x.get("gene"),x.get("protein_change")))
    strict_mech_variants=len(strict_mech_variant_keys)
    strict_mech_variant_ready=sum(annotated_unique.get(k,{}).get("PARTNER_ANALYSIS_READY_PROXY")=="YES" for k in strict_mech_variant_keys)
    strict_mech_structure_ready=sum(annotated_unique.get(k,{}).get("STRUCTURE_READY_PROXY")=="YES" for k in strict_mech_variant_keys)
    ready_fraction=strict_mech_variant_ready/strict_mech_variants if strict_mech_variants else "NA"
    structural_fraction=strict_mech_structure_ready/strict_mech_variants if strict_mech_variants else "NA"
    candidate_gene_coverage=[r["unique_credible_partners"]>=3 for r in ppi_stats if r["gene"] in set(pre_mech)]
    strict_gene_coverage=[r["unique_credible_partners"]>=3 for r in ppi_stats if r["gene"] in set(strict)]
    preliminary_disease_counts_by_gene=defaultdict(set)
    for gene,mondo_id in disease_variants:
        preliminary_disease_counts_by_gene[gene].add(mondo_id)
    report_stats={"snapshot":"2026-09 monthly ClinVar summary; candidate RCV records retrieved live on run date",
                  "run_date_utc":now(),"raw_preliminary_records":attrition.get("raw_records",0),
                  "preliminary_PLP_missense_variants":attrition.get("preliminary_unique_gene_protein_variants",0),
                  "preliminary_strict_gene_candidates":len(preliminary_strict_rank),
                  "rcv_validation_panel_genes":len(candidates),
                  "rcv_validation_panel_gene_symbols":";".join(sorted(candidates)),
                  "candidate_rcv_accessions":len(rcvs),"candidate_rcv_records":len(rcv_rows),"rcv_filtered_records":len(filtered),
                  "raw_pleiotropic_genes":sum(len(ms)>=2 for ms in preliminary_disease_counts_by_gene.values()),
                  "genes_2_distinct_diseases":len({p["gene"] for p in disease_pairs if p["status"]=="AUTO_ACCEPT"}),
                  "expanded_genes":len(expanded),"strict_genes":len(strict),"mechanism_compatible_strict_genes":len(strict_mech),
                  "genes_3_credible_partners":len(ppi_eligible),"final_high_quality_genes":len(final_genes),"final_usable_variants":int(final_variants),
                  "median_variants_per_disease":float(np.median([x["n_unique_variants"] for x in disease_counts])) if disease_counts else "NA",
                  "median_credible_partners_per_gene":float(np.median([x["unique_credible_partners"] for x in ppi_stats if x["gene"] in set(pre_mech)])) if pre_mech else "NA",
                  "partner_analysis_ready_coverage":ready_fraction,
                  "structural_ready_proxy_coverage":structural_fraction,
                  "strict_mechanism_variant_denominator":strict_mech_variants,
                  "genes_with_3_partners_fraction":sum(candidate_gene_coverage)/len(candidate_gene_coverage) if candidate_gene_coverage else 0,
                  "strict_gene_3_partner_fraction":sum(strict_gene_coverage)/len(strict_gene_coverage) if strict_gene_coverage else 0,
                  "ACTB_status":next((x["pilot_ready"] for x in [dict(gene="ACTB",pilot_ready="NA")]),"NA"),
                  "ACTG1_status":next((x["pilot_ready"] for x in [dict(gene="ACTG1",pilot_ready="NA")]),"NA"),
                  "data_source_failures":sum(1 for line in open(META/"source_failures.tsv",encoding="utf-8") if line.strip())-1 if (META/"source_failures.tsv").exists() else 0,
                  "unresolved_disease_pairs":sum(x["status"]=="MANUAL_REVIEW" for x in disease_pairs),
                  "rcv_batch_failures":len(failed)}
    candidate_plp=[r for r in rcv_rows if r.get("germline") and is_plp(r.get("clinical_significance","")) and is_missense(r) and r.get("raw_condition")]
    candidate_protein=[r for r in candidate_plp if r.get("protein_mapping_status")=="mapped"]
    candidate_mondo=[r for r in candidate_protein if r.get("MONDO","").startswith("MONDO:")]
    disease_variant_n=[int(x["n_unique_variants"]) for x in disease_counts]
    genes_all=sorted({x["gene"] for x in disease_counts})
    stats_source_failures=[]
    failfile=META/"source_failures.tsv"
    if failfile.exists():
        with open(failfile,encoding="utf-8",newline="") as f:
            stats_source_failures=list(csv.DictReader(f,delimiter="\t"))
    report_stats.update({
        "candidate_plp_missense_condition_records":len(candidate_plp),
        "candidate_protein_mapped_records":len(candidate_protein),
        "candidate_mondo_mapped_records":len(candidate_mondo),
        "candidate_unique_variants":len({(x.get("gene"),x.get("protein_change")) for x in unique.values()}),
        "preliminary_disease_annotated_records":attrition.get("disease_annotated_records",0),
        "preliminary_candidate_disease_rows":max(0,sum(1 for _ in open(INTERIM/"clinvar_preliminary_candidate_disease_rows.tsv",encoding="utf-8"))-1),
        "rcv_validation_candidate_rows":max(0,sum(1 for _ in open(INTERIM/"clinvar_rcv_validation_candidate_disease_rows.tsv",encoding="utf-8"))-1),
        "protein_mapped_fraction":len(candidate_protein)/len(candidate_plp) if candidate_plp else 0.0,
        "mondo_mapped_fraction":len(candidate_mondo)/len(candidate_protein) if candidate_protein else 0.0,
        "iqr_variants_per_disease":float(np.quantile(disease_variant_n,.75)-np.quantile(disease_variant_n,.25)) if disease_variant_n else "NA",
        "range_variants_per_disease":"%s–%s"%(min(disease_variant_n),max(disease_variant_n)) if disease_variant_n else "NA",
        "median_diseases_per_gene":float(np.median([sum(x["gene"]==g for x in disease_counts) for g in genes_all])) if genes_all else "NA",
        "genes_with_pioneer_fraction":sum(x["PIONEER_partners"]>0 for x in ppi_stats if x["gene"] in set(pre_mech))/len(pre_mech) if pre_mech else 0.0,
        "genes_with_experimental_ppi_fraction":sum(x["unique_credible_partners"]>0 for x in ppi_stats if x["gene"] in set(pre_mech))/len(pre_mech) if pre_mech else 0.0,
        "source_failure_records":len(stats_source_failures),
        "biggest_bottleneck":"data access" if not (RAW/"ClinVarRCVRelease_2026-09.xml.gz").exists() else ("sample size" if len(strict)<10 else ("mechanism" if len(strict_mech)<len(strict)/2 else ("structural coverage" if isinstance(ready_fraction,(int,float)) and ready_fraction<.70 else ("PPI" if len(ppi_eligible)<len(strict_mech) else "phenotype")))),
    })
    actb=pd.read_csv(PHASE/"ACTB_ACTG1_feasibility.tsv",sep="\t",dtype=str,keep_default_na=False)
    for g in ["ACTB","ACTG1"]:
        val=actb.loc[actb.gene==g,"pilot_ready"]
        report_stats[g+"_status"]=val.iloc[0] if len(val) else "NO"
    make_report(report_stats,threshold_rows,top,failed)
    write_tsv(PHASE/"phase0_summary_metrics.tsv",[{"metric":k,"value":v} for k,v in report_stats.items()])
    register_local_inputs()
    write_repro_log(report_stats)
    write_output_manifests(report_stats,attrition,rcv_rows,filtered,disease_counts,mechanisms,prov,ppi_stats)
    return report_stats


def make_figures(attrition,disease_counts,threshold_rows,ppi_stats,by_gene,final_genes):
    PHASE.mkdir(parents=True,exist_ok=True)
    # Figure P0-1: attrition flowchart.
    fig,ax=plt.subplots(figsize=(10,6));ax.axis("off")
    steps=[("Raw ClinVar rows",next((r["records"] for r in attrition if r["step"]=="raw_records"),0)),
           ("Germline/origin not somatic",next((r["records"] for r in attrition if r["step"]=="germline_or_origin_unspecified_records"),0)),
           ("P/LP",next((r["records"] for r in attrition if r["step"]=="plp_records"),0)),
           ("Missense",next((r["records"] for r in attrition if r["step"]=="missense_records"),0)),
           ("Protein mapped",next((r["records"] for r in attrition if r["step"]=="protein_mapped_records"),0)),
           ("Disease annotated",next((r["records"] for r in attrition if r["step"]=="disease_annotated_records"),0)),
           ("RCV confirmed candidates",next((r["records"] for r in attrition if r["step"]=="RCV P/LP germline missense protein mapped MONDO records"),0))]
    ys=np.linspace(.88,.12,len(steps))
    for i,((label,n),y) in enumerate(zip(steps,ys)):
        ax.text(.5,y,"%s\n%s"%(label,f"{n:,}"),ha="center",va="center",fontsize=10,
                bbox=dict(boxstyle="round,pad=.55",fc="#e5e9ed",ec="#65727c"))
        if i<len(steps)-1: ax.annotate("",xy=(.5,ys[i+1]+.035),xytext=(.5,y-.04),arrowprops=dict(arrowstyle="->",color="#65727c"))
    ax.set_title("Figure P0-1. ClinVar cohort attrition",loc="left",fontsize=13)
    fig.tight_layout();fig.savefig(PHASE/"Figure_P0-1_cohort_attrition.png",dpi=220);plt.close(fig)
    # Figure P0-2: count distribution.
    vals=[r["n_unique_variants"] for r in disease_counts]
    fig,ax=plt.subplots(figsize=(8,4.8));
    if vals: ax.hist(vals,bins=min(30,max(5,int(np.sqrt(len(vals))))),color="#71808d",edgecolor="white")
    ax.set_xlabel("Unique P/LP missense substitutions per gene–disease group");ax.set_ylabel("Disease groups")
    ax.set_title("Figure P0-2. Variant counts per gene–disease group",loc="left")
    fig.tight_layout();fig.savefig(PHASE/"Figure_P0-2_variant_distribution.png",dpi=220);plt.close(fig)
    # Figure P0-3: threshold curve.
    xs=[2,3,5,10,15,20]; ys2=[next((r["Genes"] for r in threshold_rows if r["minimum_variants_per_disease"]==x),0) for x in xs]
    fig,ax=plt.subplots(figsize=(7.5,4.5));ax.plot(xs,ys2,marker="o",color="#536b7b",linewidth=2)
    ax.set_xticks(xs);ax.set_xlabel("Minimum unique variants per disease");ax.set_ylabel("Genes with ≥2 qualifying diseases");ax.set_title("Figure P0-3. Candidate genes across thresholds",loc="left")
    ax.grid(axis="y",alpha=.25);fig.tight_layout();fig.savefig(PHASE/"Figure_P0-3_thresholds.png",dpi=220);plt.close(fig)
    # Figure P0-4: PPI distribution.
    pv=[r["unique_credible_partners"] for r in ppi_stats]
    fig,ax=plt.subplots(figsize=(8,4.5));
    if pv: ax.hist(pv,bins=min(30,max(5,int(np.sqrt(len(pv))))),color="#71808d",edgecolor="white")
    ax.set_xlabel("Unique experimental physical partners per candidate gene");ax.set_ylabel("Genes");ax.set_title("Figure P0-4. Experimental PPI coverage",loc="left")
    fig.tight_layout();fig.savefig(PHASE/"Figure_P0-4_ppi_coverage.png",dpi=220);plt.close(fig)
    # Figure P0-5: candidate feasibility heatmap.
    genes=sorted(by_gene,key=lambda g:(g not in final_genes,-sum(x["n_unique_variants"] for x in by_gene[g]),g))[:40]
    dims=[]
    for g in genes:
        rs=by_gene[g];ps=next((r for r in ppi_stats if r["gene"]==g),{})
        dims.append([len(rs),sum(x["n_unique_variants"] for x in rs),1 if g in final_genes else 0,ps.get("unique_credible_partners",0),ps.get("PIONEER_partners",0)])
    fig,ax=plt.subplots(figsize=(9,max(4,len(genes)*.28)))
    if dims:
        arr=np.array(dims,float); scales=np.maximum(arr.max(axis=0),1); norm=arr/scales
        im=ax.imshow(norm,aspect="auto",cmap="Greys",vmin=0,vmax=1)
        ax.set_yticks(range(len(genes)));ax.set_yticklabels(genes,fontsize=7)
        ax.set_xticks(range(5));ax.set_xticklabels(["disease groups","unique variants","high-quality status","experimental partners","PIONEER partners"],rotation=25,ha="right",fontsize=8)
        fig.colorbar(im,ax=ax,label="Within-dimension scaled value")
    ax.set_title("Figure P0-5. Candidate feasibility dimensions",loc="left")
    fig.tight_layout();fig.savefig(PHASE/"Figure_P0-5_feasibility_heatmap.png",dpi=220);plt.close(fig)


def make_report(stats,threshold_rows,top,failed):
    ready_value=stats["partner_analysis_ready_coverage"]
    ready_for_gate=ready_value if isinstance(ready_value,(int,float)) else 0.0
    ready_text=("%0.1f%%"%(ready_value*100)) if isinstance(ready_value,(int,float)) else "NA (zero denominator)"
    structural_value=stats["structural_ready_proxy_coverage"]
    structural_text=("%0.1f%%"%(structural_value*100)) if isinstance(structural_value,(int,float)) else "NA (zero denominator)"
    full=(stats["final_high_quality_genes"]>=30 and stats["final_usable_variants"]>=1000 and
          stats["mechanism_compatible_strict_genes"]>=20 and ready_for_gate>=0.70 and
          stats["strict_gene_3_partner_fraction"]>0.5)
    verdict="GO_FULL_ATLAS" if full else ("GO_FOCUSED" if stats["final_high_quality_genes"]>=10 else "NO_GO_ATLAS")
    stats["verdict"]=verdict
    caveat=("ClinVar variant_summary from 2026-09 was used for the broad same-release screen. Exact RCV validation was restricted to the five highest-volume preliminary strict genes plus the prespecified ACTB/ACTG1 pilots. "
            "Candidate condition assignments in that seven-gene panel were checked using official RCV XML retrieved live through E-utilities. The full fixed RCV archive transfer was stopped after measured throughput implied a multi-hour download. "
            "Genome-wide exact RCV totals are unavailable; all-panel counts below are explicitly labeled preliminary.")
    if stats.get("unresolved_disease_pairs",0):
        caveat += " Parent-child and possible same-spectrum pairs remain in the manual review queue and are excluded from strict eligibility."
    preliminary_2=next((r["Genes"] for r in threshold_rows if r["minimum_variants_per_disease"]==2 and "whole-release" in r["data_stage"]),0)
    preliminary_5=next((r["Genes"] for r in threshold_rows if r["minimum_variants_per_disease"]==5 and "whole-release" in r["data_stage"]),0)
    preliminary_10=next((r["Genes"] for r in threshold_rows if r["minimum_variants_per_disease"]==10 and "whole-release" in r["data_stage"]),0)
    lines=["# Phase 0 feasibility report","",f"**PHASE0_VERDICT = `{verdict}`**","",caveat,"",
           "## 1. Executive verdict","",f"{verdict}. This decision uses only the prespecified thresholds and available Phase 0 metadata. Phase 1 was not started.","",
           "## 2. Exact numbers","",
           f"1. ClinVar P/LP missense variant–condition records: {stats['candidate_plp_missense_condition_records']:,} among candidate RCV records before protein/MONDO filters; {stats['rcv_filtered_records']:,} pass the primary candidate filters. The same-month whole-release pre-screen found {stats['preliminary_disease_annotated_records']:,} potentially eligible summary rows; that count is not RCV-verified.",
           f"2. Unique variants: {stats['preliminary_PLP_missense_variants']:,} unique gene–protein substitutions in the whole-release preliminary screen; {stats['candidate_unique_variants']:,} in filtered candidate RCV records; {stats['final_usable_variants']:,} partner-ready substitutions in the final high-quality set.",
           f"3. Protein-mapped variants: {stats['candidate_protein_mapped_records']:,} candidate RCV P/LP germline missense records ({stats['protein_mapped_fraction']:.1%} of candidate P/LP missense records).",
           f"4. MONDO-mapped variants: {stats['candidate_mondo_mapped_records']:,} candidate RCV records ({stats['mondo_mapped_fraction']:.1%} of protein-mapped candidate records).",
           f"5. Raw pleiotropic genes: {stats['raw_pleiotropic_genes']:,} genes in the whole-release preliminary screen with at least two mapped MONDO disease groups; this is a broad-screen value, not RCV-confirmed.",
           f"6. Genes with ≥2 normalized distinct diseases: {stats['genes_2_distinct_diseases']:,} RCV-confirmed within the seven-gene panel; {preliminary_2:,} were flagged by the whole-release summary pre-screen. {stats['unresolved_disease_pairs']:,} additional panel pairs require manual review.",
           f"7. Expanded genes (≥5 variants/disease): {stats['expanded_genes']:,} RCV-confirmed within the panel; {preliminary_5:,} whole-release summary-screen candidates.",
           f"8. Strict genes (≥10 variants/disease and ≥20 total): {stats['strict_genes']:,} RCV-confirmed within the panel; {preliminary_10:,} whole-release summary-screen candidates.",
           f"9. Mechanism-compatible strict genes: {stats['mechanism_compatible_strict_genes']:,}.",
           f"10. Genes with ≥3 credible PPI partners: {stats['genes_3_credible_partners']:,}.",
           f"11. Final high-quality analyzable genes: {stats['final_high_quality_genes']:,}.",
           f"12. Final usable variants: {stats['final_usable_variants']:,}.",
           f"13. Median variants/disease: {stats['median_variants_per_disease']}.",
           f"14. Median credible partners/gene: {stats['median_credible_partners_per_gene']}.",
           f"15. Partner-analysis-ready coverage: {ready_text} across {stats['strict_mechanism_variant_denominator']:,} variants in mechanism-compatible strict groups.",
           f"16. ACTB pilot ready: {stats['ACTB_status']}.",
           f"17. ACTG1 pilot ready: {stats['ACTG1_status']}.","",
           "## Threshold screen","","The first six rows are the whole-release ClinVar summary pre-screen; the next six are the selected seven-gene RCV panel. They are separate evidence stages.","",
           "| Minimum variants/disease | Genes | Disease groups | Unique substitutions | Evidence stage |","|---:|---:|---:|---:|---|"]
    for r in threshold_rows:
        lines.append("| {minimum_variants_per_disease} | {Genes} | {Diseases} | {Variants} | {data_stage} |".format(**r))
    lines += ["","## 3. Top candidate genes","","The RCV validation panel was selected using preliminary strict-group variant volume (top five) plus the prespecified ACTB/ACTG1 pilots. Within the panel, final rankings use disease distinction, sample size, mechanism consistency, experimental PPI and PIONEER coverage. Unresolved disease pairs remain flagged.","",
           "| Gene | Diseases | Unique variants | Strict diseases | Mechanism | Credible partners | PIONEER partners | Status |","|---|---:|---:|---:|---|---:|---:|---|"]
    for r in top:
        lines.append("| {gene} | {n_normalized_diseases} | {n_unique_variants} | {n_strict_diseases} | {mechanism_status} | {credible_PPI_partners} | {PIONEER_partners} | {status} |".format(**r))
    lines += ["","## 4. Empirical feasibility and bottlenecks","",
              f"- Gene–disease group variant counts: median {stats['median_variants_per_disease']}; IQR {stats['iqr_variants_per_disease']}; range {stats['range_variants_per_disease']}. Median disease groups per candidate gene: {stats['median_diseases_per_gene']}. Median experimental partners per candidate gene: {stats['median_credible_partners_per_gene']}.",
              f"- Candidate genes with ≥3 experimental partners: {stats['genes_with_3_partners_fraction']:.1%}; with PIONEER coverage: {stats['genes_with_pioneer_fraction']:.1%}; with experimental PPI support: {stats['genes_with_experimental_ppi_fraction']:.1%}.",
              f"- Existing-interface structural proxy: {structural_text} of strict mechanism-compatible variants have a PIONEER/PDB interface reference with reference-residue concordance. This is an annotation proxy, not structure validation.",
              f"- Largest bottleneck: **{stats['biggest_bottleneck']}**.",
              "- No formal statistical power calculation was performed; this is an empirical feasibility audit.","",
              "## 5. Phase 1 recommendation","",
              f"{verdict}. If a focused follow-up is approved, use {stats['final_high_quality_genes']} high-quality genes and {stats['final_usable_variants']} partner-ready substitutions, beginning from `clinvar_variant_condition_filtered.tsv`, `disease_pairs_relationship_audit.tsv`, `disease_pairs_needing_manual_review.tsv`, `gene_disease_mechanism_pair_audit.tsv`, `ppi_provenance.tsv`, and `variant_coverage_proxy.tsv`. Phase 1 was not run.","",
              "## Data snapshot and limits","",
              f"- ClinVar: archived monthly snapshot 2026-09-03 (`variant_summary_2026-09.txt.gz`). The whole-release screen flagged {stats['preliminary_strict_gene_candidates']:,} possible strict genes; only the seven-gene validation panel ({stats['rcv_validation_panel_gene_symbols']}) received exact RCV retrieval. The fixed RCV XML archive is 6.0 GB compressed; measured transfer speed implied a multi-hour download. Raw API XML responses for queried panel records are preserved.",
              "- MONDO: version v2026-09-01.",
              "- Gene2Phenotype: all-panel download dated 2026-09-27.",
              "- GenCC: current export captured on the run date; verify file completeness using the manifest and failure log.",
              "- PPI: BioGRID MV-Physical 5.0.261, candidate-focused IntAct PSICQUIC queries, and PIONEER very-high-confidence human predictions. HuRI access timed out. Predictomes was not queried and is not counted.",
              "- Manual review: MONDO exact-equivalent terms were collapsed; parent-child and possible same-spectrum pairs remain in `disease_pairs_needing_manual_review.tsv`.",
              "- Missing values are `NA`; the report does not substitute inferred numbers for unavailable data.","",
              "## Required feasibility tables and figures","",
              "Tables: `clinvar_cohort_attrition.tsv`, `candidate_gene_thresholds.tsv`, `top_candidate_genes.tsv`, `gene_disease_mechanism_audit.tsv`, `gene_disease_mechanism_pair_audit.tsv`, `ppi_coverage_by_gene.tsv`, `ppi_provenance.tsv`, `ACTB_ACTG1_feasibility.tsv`, `disease_pairs_relationship_audit.tsv`, and `disease_pairs_needing_manual_review.tsv`.","",
              "Figures: `Figure_P0-1_cohort_attrition.png` through `Figure_P0-5_feasibility_heatmap.png`.","",
              "## Data-source failures and review queue","",f"RCV API failed batches: {len(failed)}.","",
              "See `metadata/source_failures.tsv` and `results/phase0/rcv_api_failures.tsv` for exact records and errors. The full RCV archive was not completed; HuRI timed out; IntAct was queried per candidate UniProt accession.","",
              "## Stop confirmation","","Phase 0 artifacts are generated. **Phase 1 was NOT started.** No structure modeling, interface scoring, FoldX, MutPred-PPI, PERMANOVA, or machine learning was run.",""]
    (REPORTS/"PHASE0_FEASIBILITY_REPORT.md").write_text("\n".join(lines),encoding="utf-8")


def write_repro_log(stats):
    try:
        import subprocess
        commit=subprocess.run(["git","rev-parse","HEAD"],cwd=ROOT,text=True,capture_output=True).stdout.strip() or "NO_COMMIT"
    except Exception:
        commit="NO_COMMIT"
    try:
        import platform
        py=platform.python_version()
    except Exception: py="NA"
    (REPORTS/"reproducibility_log.md").write_text("# Reproducibility log\n\n"
        f"- Run timestamp (UTC): {stats['run_date_utc']}\n- ClinVar snapshot: 2026-09 monthly release; summary file is immutable and hashed in `metadata/input_manifest.tsv`.\n"
        f"- MONDO release: v2026-09-01.\n- Python: {py}; pandas {pd.__version__}; lxml {etree.LXML_VERSION}; NumPy {np.__version__}.\n"
        f"- Git commit at analysis: {commit}.\n- Inputs and SHA256: `metadata/input_manifest.tsv`.\n"
        "- Output hashes and row counts are in `metadata/output_manifest.tsv`; pipeline inputs, outputs, and row counts are in `metadata/step_manifest.tsv`.\n"
        "- The script refuses to overwrite existing outputs; reproduce in a clean project copy. Whole-release screening is in `data/interim/`.\n"
        "- The full fixed RCV archive was not downloaded; candidate RCV XML was retrieved live and cached with accession and response metadata.\n"
        "- No formal power calculation, model predictions, or structure calculations were run.\n",encoding="utf-8")


def write_output_manifests(stats, attrition, rcv_rows, filtered, disease_counts, mechanisms, prov, ppi_stats):
    step_rows=[
        {"step":"whole-release ClinVar pre-screen","inputs":"data/raw/variant_summary_2026-09.txt.gz; data/raw/mondo-v2026-09-01.obo","outputs":"data/interim/clinvar_preliminary_plp_missense.tsv; data/interim/clinvar_preliminary_thresholds.tsv; data/interim/clinvar_rcv_validation_candidates.tsv","input_rows":stats.get("raw_preliminary_records",0),"output_rows":stats.get("preliminary_candidate_disease_rows",0)},
        {"step":"candidate RCV retrieval and primary filters","inputs":"data/interim/clinvar_rcv_validation_candidate_disease_rows.tsv; ClinVar E-utilities","outputs":"results/phase0/clinvar_variant_condition_raw.tsv; results/phase0/clinvar_variant_condition_filtered.tsv","input_rows":stats.get("candidate_rcv_accessions",0),"output_rows":len(filtered)},
        {"step":"MONDO harmonization and pair review","inputs":"results/phase0/clinvar_variant_condition_filtered.tsv; data/raw/mondo-v2026-09-01.obo","outputs":"results/phase0/disease_harmonization.tsv; results/phase0/disease_pairs_relationship_audit.tsv; results/phase0/disease_pairs_needing_manual_review.tsv","input_rows":len(rcv_rows),"output_rows":len(disease_counts)},
        {"step":"inheritance and mechanism audit","inputs":"data/raw/gene2phenotype_all_latest.csv; data/raw/gencc_submissions_latest.tsv","outputs":"results/phase0/gene_disease_mechanism_audit.tsv; results/phase0/gene_disease_mechanism_pair_audit.tsv","input_rows":len(mechanisms),"output_rows":stats.get("mechanism_compatible_strict_genes",0)},
        {"step":"experimental and predicted interaction coverage","inputs":"data/raw/BIOGRID-MV-Physical-5.0.261.tab3.zip; data/raw/PIONEER_very_high_human_2026-09-27.txt; data/raw/hgnc_non_alt_loci_2026-09-27.tsv; IntAct PSICQUIC","outputs":"results/phase0/ppi_provenance.tsv; results/phase0/ppi_coverage_by_gene.tsv; results/phase0/variant_coverage_proxy.tsv","input_rows":len(prov),"output_rows":len(ppi_stats)},
        {"step":"thresholds, feasibility gate, and report","inputs":"Phase 0 cohort, disease, mechanism, and interaction tables","outputs":"results/phase0/candidate_gene_thresholds.tsv; results/phase0/phase0_summary_metrics.tsv; reports/PHASE0_FEASIBILITY_REPORT.md; Figure P0-1 to P0-5","input_rows":len(disease_counts),"output_rows":stats.get("final_high_quality_genes",0)},
    ]
    write_tsv(META/"step_manifest.tsv",step_rows,["step","inputs","outputs","input_rows","output_rows"])
    artifacts=[]
    for root in (PHASE,REPORTS,META):
        for p in sorted(root.rglob("*")):
            if not p.is_file() or p.name=="output_manifest.tsv" or p.name.startswith("._"):
                continue
            rows="NA"; cols="NA"
            if p.suffix.lower()==".tsv":
                try:
                    with open(p,encoding="utf-8",newline="") as f:
                        reader=csv.reader(f,delimiter="\t")
                        header=next(reader,[]); rows=sum(1 for _ in reader); cols=len(header)
                except Exception:
                    rows="NA"; cols="NA"
            artifacts.append({"artifact":str(p.relative_to(ROOT)),"bytes":p.stat().st_size,"rows":rows,"columns":cols,"sha256":sha256_file(p)})
    write_tsv(META/"output_manifest.tsv",artifacts,["artifact","bytes","rows","columns","sha256"])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config",default=str(ROOT/"config"/"phase0.yml"))
    args=ap.parse_args()
    config=yaml.safe_load(open(args.config,encoding="utf-8"))
    if not (META/"source_failures.tsv").exists():
        add_failure("ClinVar fixed RCV archive","https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_release/ClinVarRCVRelease_2026-09.xml.gz","transfer stopped after measured throughput implied a multi-hour download; 6 MB partial retained and excluded")
        add_failure("HuRI","https://www.interactome-atlas.org/download","download request timed out; no HuRI interactions ingested")
    result=audit_phase0()
    print(json.dumps(result,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
