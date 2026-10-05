#!/usr/bin/env python3
"""Build the frozen Phase 1A partner-resolved representation/QC package.

Consumes the completed Phase 0R tables and frozen source files. Network-derived
UniProt, RCSB and NCBI records are cached under data/interim/phase1A_*.
Disease-label association testing is deliberately out of scope.
"""
from __future__ import annotations

import ast
import csv
import gzip
import hashlib
import json
import math
import re
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import requests
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw"
INTERIM = ROOT / "data/interim"
PHASE0 = ROOT / "results/phase0R"
OUT = ROOT / "results/phase1A"
REPORT = ROOT / "reports/PHASE1A_PARTNER_REPRESENTATION_FEASIBILITY.md"
STRICT_DISEASES = PHASE0 / "gene_disease_variant_counts_exact_rcv.tsv"
MECH_GENES = PHASE0 / "mechanism_gene_states_exact_strict.tsv"
MECH_PAIRS = PHASE0 / "mechanism_pair_states_exact_strict.tsv"
PPI_COVERAGE = PHASE0 / "ppi_coverage_phase0R_candidates.tsv"
PPI_PROVENANCE = PHASE0 / "ppi_provenance_phase0R_candidates.tsv"
SEED = PHASE0 / "global_plp_missense_rcv_seed.tsv"
RECOVERED = PHASE0 / "global_rcv_recovered.tsv"
VARIANT_SUMMARY = RAW / "variant_summary_2026-09.txt.gz"
PIONEER = RAW / "PIONEER_very_high_human_2026-09-27.txt"
HGNC = RAW / "hgnc_non_alt_loci_2026-09-27.tsv"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_FILES = "https://files.rcsb.org/download/{}-assembly{}.cif.gz"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{}"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{}.fasta"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
DISTANCE_CUTOFF_A = 5.0
OVERLAP_HIGH = 0.50
OVERLAP_LOW = 0.10
DIFFERENTIAL_COVERAGE_PP = 20.0
LOW_COVERAGE_CUTOFF = 70.0
AA3 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H",
    "Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V",
}


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def split_ids(value: Any) -> list[str]:
    return [x.strip() for x in re.split(r"[|,;\s]+", str(value or "")) if x.strip() and x.strip() not in {"-", ".", "NA"}]


def base_accession(value: Any) -> str:
    return str(value or "").strip().split("-")[0]


def as_list(value: Any) -> list[Any]:
    if isinstance(value, list):
        return value
    if value is None:
        return []
    return [value]


def cell(d: dict[str, Any], key: str, n: int) -> Any:
    vals = d.get(key, [])
    return vals[n] if n < len(vals) else ""


def write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        columns = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in columns})


def stream_tsv(path: Path, columns: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in columns})


def load_inputs():
    strict = pd.read_csv(STRICT_DISEASES, sep="\t", dtype=str, keep_default_na=False)
    strict = strict[(strict["strict_disease_ge10"] == "YES") & strict["MONDO"].str.startswith("MONDO:")].copy()
    mg = pd.read_csv(MECH_GENES, sep="\t", dtype=str, keep_default_na=False)
    mg = mg[mg["analysis_set"] == "STRICT"].copy()
    strict = strict[strict["gene"].isin(set(mg["gene"]))].copy()
    mp = pd.read_csv(MECH_PAIRS, sep="\t", dtype=str, keep_default_na=False)
    mp = mp[mp["analysis_set"] == "STRICT"].copy()
    cov = pd.read_csv(PPI_COVERAGE, sep="\t", dtype=str, keep_default_na=False)
    prov = pd.read_csv(PPI_PROVENANCE, sep="\t", dtype=str, keep_default_na=False)
    seed = pd.read_csv(SEED, sep="\t", dtype=str, keep_default_na=False)
    rec = pd.read_csv(RECOVERED, sep="\t", dtype=str, keep_default_na=False)
    hgnc = pd.read_csv(HGNC, sep="\t", dtype=str, keep_default_na=False)
    return strict, mg, mp, cov, prov, seed, rec, hgnc


def make_cohort(strict, mg, mp, cov):
    strict_genes = set(strict["gene"])
    mech = mg.set_index("gene")["mechanism_status"].to_dict()
    ppi = cov.set_index("gene").to_dict("index")
    partner_n = {g: int(ppi.get(g, {}).get("unique_credible_partners", 0) or 0) for g in strict_genes}
    incompatible = {g for g in strict_genes if mech.get(g) == "CONFIRMED_INCOMPATIBLE"}
    primary = {g for g in strict_genes if g not in incompatible and partner_n[g] >= 3}
    strata = {}
    for gene in strict_genes:
        pair_states = set(mp.loc[mp["gene"] == gene, "mechanism_state"].tolist())
        if pair_states and pair_states <= {"CONFIRMED_COMPATIBLE"}:
            strata[gene] = "KNOWN_COMPATIBLE"
        elif "CONFIRMED_INCOMPATIBLE" in pair_states:
            strata[gene] = "KNOWN_INCOMPATIBLE"
        else:
            strata[gene] = "UNKNOWN_COMPATIBILITY"
    rows = []
    for gene in sorted(strict_genes):
        disease_n = int((strict["gene"] == gene).sum())
        partners = partner_n[gene]
        base = {
            "gene": gene, "gene_mechanism_stratum": strata[gene],
            "gene_level_phase0R_mechanism_status": mech.get(gene, "NA"),
            "strict_disease_group_n": disease_n, "credible_partner_n": partners,
            "partner_band": "<2" if partners < 2 else "2" if partners == 2 else "3-5" if partners <= 5 else "6-10" if partners <= 10 else ">10",
            "primary_cohort_eligible": "YES" if gene in primary else "NO",
            "primary_exclusion_reason": "" if gene in primary else "CONFIRMED_INCOMPATIBLE" if gene in incompatible else "FEWER_THAN_3_CREDIBLE_PARTNERS",
        }
        gene_pairs = mp[mp["gene"] == gene]
        if gene_pairs.empty:
            rows.append({**base, "MONDO_A":"", "disease_A":"", "MONDO_B":"", "disease_B":"", "mechanism_state":"NO_DISEASE_PAIR_ROW"})
        else:
            for r in gene_pairs.to_dict("records"):
                rows.append({**base, **{k:r.get(k,"") for k in ["MONDO_A","disease_A","mechanism_A","inheritance_A","evidence_rows_A","MONDO_B","disease_B","mechanism_B","inheritance_B","evidence_rows_B","mechanism_state"]}})
    meta = {
        "strict_genes": len(strict_genes), "confirmed_incompatible_genes": len(incompatible), "primary_cohort_genes": len(primary),
        "strict_ge3_credible_partners": sum(n >= 3 for n in partner_n.values()),
        "partner_bands_all_strict": {
            "exactly_2":sum(n==2 for n in partner_n.values()), "3_to_5":sum(3<=n<=5 for n in partner_n.values()),
            "6_to_10":sum(6<=n<=10 for n in partner_n.values()), "over_10":sum(n>10 for n in partner_n.values()), "under_2":sum(n<2 for n in partner_n.values())},
        "partner_bands_after_mechanism_exclusion": {
            "exactly_2":sum(n==2 for g,n in partner_n.items() if g not in incompatible), "3_to_5":sum(3<=n<=5 for g,n in partner_n.items() if g not in incompatible),
            "6_to_10":sum(6<=n<=10 for g,n in partner_n.items() if g not in incompatible), "over_10":sum(n>10 for g,n in partner_n.items() if g not in incompatible),
            "under_2":sum(n<2 for g,n in partner_n.items() if g not in incompatible)},
        "primary_genes": sorted(primary), "incompatible_genes": sorted(incompatible),
    }
    return rows, primary, incompatible, partner_n, strata, meta


def parse_protein_change(value: str) -> tuple[str, int, str] | None:
    s = str(value or "").strip().replace("p.", "").replace("(", "").replace(")", "")
    m = re.fullmatch(r"([A-Z][a-z]{2}|[A-Z])(\d+)([A-Z][a-z]{2}|[A-Z])", s)
    if not m:
        return None
    return AA3.get(m.group(1),m.group(1)), int(m.group(2)), AA3.get(m.group(3),m.group(3))


def parse_clinvar_transcripts(variant_keys: set[tuple[str,str]]) -> dict[tuple[str,str],set[str]]:
    transcripts: dict[tuple[str,str],set[str]] = defaultdict(set)
    names: dict[tuple[str,str],set[str]] = defaultdict(set)
    with gzip.open(VARIANT_SUMMARY,"rt",encoding="utf-8",errors="replace",newline="") as f:
        reader=csv.reader(f,delimiter="\t")
        header=next(reader); idx={x:i for i,x in enumerate(header)}
        i_name,i_gene,i_var=idx["Name"],idx["GeneSymbol"],idx["VariationID"]
        for row in reader:
            if len(row)<=max(i_name,i_gene,i_var): continue
            key=(row[i_gene],row[i_var])
            if key not in variant_keys: continue
            names[key].add(row[i_name])
            transcripts[key].update(re.findall(r"\b(?:NM|NR)_\d+\.\d+",row[i_name]))
    for key in variant_keys:
        if not transcripts.get(key):
            for name in names.get(key,set()):
                transcripts[key].update(re.findall(r"\b(?:NM|NR)_\d+(?:\.\d+)?",name))
    return transcripts


def fetch_refseq_proteins(transcripts: set[str]) -> dict[str,dict[str,str]]:
    cache_json=INTERIM/"phase1A_refseq_protein_map.json"
    cached={}
    if cache_json.exists():
        try: cached=json.loads(cache_json.read_text(encoding="utf-8"))
        except Exception: cached={}
    records={k:cached[k] for k in transcripts if k in cached and cached[k].get("sequence")}
    session=requests.Session(); to_fetch=sorted(transcripts-set(records))
    for start in range(0,len(to_fetch),40):
        batch=to_fetch[start:start+40]
        for attempt in range(4):
            try:
                resp=session.get(NCBI_EFETCH,params={"db":"nuccore","id":",".join(batch),"rettype":"gb","retmode":"text"},timeout=(15,120),headers={"User-Agent":"phase1A-allelic-interface/1.0"})
                resp.raise_for_status()
                for rec in SeqIO.parse(StringIO(resp.text),"genbank"):
                    cds=next((ft for ft in rec.features if ft.type=="CDS" and ft.qualifiers.get("translation")),None)
                    if not cds: continue
                    seq=re.sub(r"[^ACDEFGHIKLMNPQRSTVWYBXZJUO]","",cds.qualifiers["translation"][0].upper())
                    records[rec.id]={"protein_accession":cds.qualifiers.get("protein_id",[""])[0],"sequence":seq,"definition":rec.description}
                break
            except Exception as exc:
                if attempt==3:
                    for tx in batch: records.setdefault(tx,{"protein_accession":"","sequence":"","error":str(exc)[:180]})
                else: time.sleep((attempt+1)*2)
        cache_json.parent.mkdir(parents=True,exist_ok=True)
        cache_json.write_text(json.dumps(records,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        time.sleep(0.34)
    return records


def fasta_sequence(text: str) -> str:
    lines=[x.strip() for x in text.splitlines() if x and not x.startswith(">")]
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWYBXZJUO]","", "".join(lines).upper())


def fetch_uniprot_sequences(accessions: set[str]) -> dict[str,str]:
    seq_dir=INTERIM/"phase1A_uniprot_sequences"; seq_dir.mkdir(parents=True,exist_ok=True)
    wanted=sorted({base_accession(x) for x in accessions if x})
    sequences={}; session=requests.Session()
    missing=[]
    for acc in wanted:
        cache=seq_dir/f"{acc}.fasta"
        if cache.exists():
            seq=fasta_sequence(cache.read_text(encoding="utf-8",errors="replace"))
            if seq: sequences[acc]=seq; continue
        missing.append(acc)
    for start in range(0,len(missing),50):
        batch=missing[start:start+50]
        query="accession:("+" OR ".join(batch)+")"
        try:
            r=session.get("https://rest.uniprot.org/uniprotkb/stream",params={"format":"fasta","query":query},timeout=(10,90),headers={"User-Agent":"phase1A-allelic-interface/1.0"})
            r.raise_for_status()
            for record in SeqIO.parse(StringIO(r.text),"fasta"):
                parts=record.id.split("|"); acc=base_accession(parts[1] if len(parts)>1 else record.id)
                seq=fasta_sequence(str(record.seq))
                if not seq: continue
                cache=seq_dir/f"{acc}.fasta"
                cache.write_text(f">{record.description}\n"+"\n".join(seq[i:i+60] for i in range(0,len(seq),60))+"\n",encoding="utf-8")
                sequences[acc]=seq
        except requests.RequestException:
            # Only fall back to individual lookups for a failed bulk request.
            for acc in batch:
                try:
                    r=session.get(UNIPROT_FASTA.format(acc),timeout=(10,25),headers={"User-Agent":"phase1A-allelic-interface/1.0"})
                    if r.ok:
                        seq=fasta_sequence(r.text)
                        if seq:
                            (seq_dir/f"{acc}.fasta").write_text(r.text,encoding="utf-8"); sequences[acc]=seq
                except requests.RequestException:
                    continue
    return sequences


def sequence_map(source: str,target: str,residue_position: int|None=None) -> dict[str,Any]:
    """Map 1-based positions in source to target without guessing tied mappings."""
    source=re.sub(r"[^ACDEFGHIKLMNPQRSTVWYBXZJUO]","",str(source or "").upper())
    target=re.sub(r"[^ACDEFGHIKLMNPQRSTVWYBXZJUO]","",str(target or "").upper())
    if not source or not target:
        return {"status":"UNMAPPABLE","mapping":{},"identity":0.0,"coverage":0.0,"ambiguous_position":False}
    if source==target:
        return {"status":"EXACT_CANONICAL","mapping":{i:i for i in range(1,len(source)+1)},"identity":1.0,"coverage":1.0,"ambiguous_position":False}
    aligner=PairwiseAligner(); aligner.mode="global"; aligner.match_score=2.0; aligner.mismatch_score=-1.0; aligner.open_gap_score=-8.0; aligner.extend_gap_score=-0.5
    try: als=aligner.align(source,target); best=als[0]
    except Exception: return {"status":"UNMAPPABLE","mapping":{},"identity":0.0,"coverage":0.0,"ambiguous_position":False}
    mapping={}; aligned=identical=0
    for si,ti in zip(best.indices[0],best.indices[1]):
        if si>=0 and ti>=0:
            mapping[int(si)+1]=int(ti)+1; aligned+=1; identical+=source[si]==target[ti]
    identity=identical/aligned if aligned else 0.0; coverage=aligned/len(source) if source else 0.0
    ambiguous=False
    if residue_position and residue_position in mapping:
        try:
            top=best.score
            for idx in range(1, min(12, len(als))):
                other = als[idx]
                if abs(other.score-top)>1e-8: break
                mapped=None
                for si,ti in zip(other.indices[0],other.indices[1]):
                    if si==residue_position-1 and ti>=0: mapped=int(ti)+1; break
                if mapped!=mapping[residue_position]: ambiguous=True; break
        except Exception: ambiguous=True
    status="SAFE_ISOFORM_REMAP" if identity>=0.90 and coverage>=0.70 else "UNMAPPABLE"
    return {"status":status,"mapping":mapping,"identity":identity,"coverage":coverage,"ambiguous_position":ambiguous}


def map_variants(strict,mg,seed,rec,hgnc,transcripts_by_key,primary):
    strict_genes=set(mg["gene"]); strict_groups=set(map(tuple,strict[["gene","MONDO"]].values))
    rec=rec.copy(); rec["MONDO_PRIMARY"]=rec["MONDO"].apply(lambda x:next((v for v in str(x).split(";") if v.startswith("MONDO:")),""))
    seed_rec=seed.merge(rec[["RCVaccession","VariationID","MONDO_PRIMARY","aggregate_germline_classification","review_status","preferred_condition_name"]],on=["RCVaccession","VariationID"],how="inner",suffixes=("","_recovered"))
    seed_rec=seed_rec[(seed_rec["GeneSymbol"].isin(strict_genes))&seed_rec.apply(lambda r:(r["GeneSymbol"],r["MONDO_PRIMARY"]) in strict_groups,axis=1)]
    variant_rows={}
    for r in seed_rec.to_dict("records"):
        key=(r["GeneSymbol"],str(r["VariationID"]))
        v=variant_rows.setdefault(key,{"gene":key[0],"VariationID":key[1],"protein_change":r.get("ProteinChange",""),"MONDOs":set(),"RCVs":set(),"clinical_significances":set(),"review_statuses":set()})
        for col,field in [("MONDO_PRIMARY","MONDOs"),("RCVaccession","RCVs"),("ClinicalSignificance","clinical_significances"),("ReviewStatus","review_statuses")]:
            if r.get(col): v[field].add(r[col])
        if not v["protein_change"] and r.get("ProteinChange"): v["protein_change"]=r["ProteinChange"]
    keys=set(variant_rows); txset=set(t for k in keys for t in transcripts_by_key.get(k,set()))
    refseq=fetch_refseq_proteins(txset)
    hg_by_gene={}
    for r in hgnc.to_dict("records"):
        if r.get("symbol") in strict_genes: hg_by_gene[r["symbol"]]=split_ids(r.get("uniprot_ids",""))
    target_ids={base_accession(x) for ids in hg_by_gene.values() for x in ids if x}
    uni_seq=fetch_uniprot_sequences(target_ids)
    gene_accs={g:[base_accession(x) for x in ids if base_accession(x) in uni_seq] for g,ids in hg_by_gene.items()}
    incompatible=set(mg.loc[mg.mechanism_status=="CONFIRMED_INCOMPATIBLE","gene"])
    map_rows=[]; variant_map={}
    for key,v in sorted(variant_rows.items()):
        gene,var=key; parsed=parse_protein_change(v["protein_change"]); ref,pos,alt=parsed if parsed else ("",0,"")
        txs=sorted(transcripts_by_key.get(key,set())); candidates=[]
        for tx in txs:
            prot=refseq.get(tx,{}); pseq=prot.get("sequence","")
            if not pseq or not pos or pos>len(pseq) or (ref and pseq[pos-1]!=ref): continue
            for uacc in gene_accs.get(gene,[]):
                useq=uni_seq.get(uacc,""); m=sequence_map(pseq,useq,pos); cpos=m.get("mapping",{}).get(pos)
                if cpos and (not ref or useq[cpos-1]==ref) and not m.get("ambiguous_position"):
                    candidates.append({"transcript":tx,"protein_accession":prot.get("protein_accession",""),"uniprot_accession":uacc,"canonical_residue":cpos,"canonical_ref":useq[cpos-1],
                        "map_status":"EXACT_CANONICAL" if pseq==useq else "SAFE_ISOFORM_REMAP","identity":m.get("identity",0),"coverage":m.get("coverage",0)})
        unique={(x["uniprot_accession"],x["canonical_residue"],x["canonical_ref"]) for x in candidates}
        if len(unique)==1:
            chosen=sorted(candidates,key=lambda x:(x["map_status"]!="EXACT_CANONICAL",-x["coverage"],-x["identity"],x["transcript"]))[0]; status=chosen["map_status"]
        elif len(unique)>1: chosen={}; status="AMBIGUOUS_ISOFORM"
        else:
            chosen={}; status="AMBIGUOUS_ISOFORM" if any(refseq.get(t,{}).get("sequence") for t in txs) and pos else "UNMAPPABLE"
        row={"gene":gene,"VariationID":var,"protein_change":v["protein_change"],"reference_aa":ref,"source_residue":pos,"alternate_aa":alt,
            "ClinVar_transcripts":";".join(txs),"RefSeq_protein_accessions":";".join(sorted({refseq.get(t,{}).get("protein_accession","") for t in txs if refseq.get(t,{}).get("protein_accession")})),
            "canonical_UniProt_accession":chosen.get("uniprot_accession",""),"canonical_residue":chosen.get("canonical_residue",""),"canonical_reference_aa":chosen.get("canonical_ref",""),
            "mapping_status":status,"mapping_identity":f"{chosen.get('identity',0):.4f}" if chosen else "","mapping_source_coverage":f"{chosen.get('coverage',0):.4f}" if chosen else "",
            "MONDOs":";".join(sorted(v["MONDOs"])),"RCVs":";".join(sorted(v["RCVs"])),"clinical_significances":";".join(sorted(v["clinical_significances"])),
            "review_statuses":";".join(sorted(v["review_statuses"])),"primary_cohort_eligible":"YES" if gene in primary else "NO"}
        map_rows.append(row); variant_map[key]=row
    columns=["gene","VariationID","protein_change","reference_aa","source_residue","alternate_aa","ClinVar_transcripts","RefSeq_protein_accessions","canonical_UniProt_accession","canonical_residue","canonical_reference_aa","mapping_status","mapping_identity","mapping_source_coverage","MONDOs","RCVs","clinical_significances","review_statuses","primary_cohort_eligible"]
    write_tsv(OUT/"protein_mapping.tsv",map_rows,columns)
    return map_rows,variant_map,uni_seq,gene_accs,seed_rec,refseq


def build_pioneer(strict_genes:set[str],hgnc:pd.DataFrame):
    uniprot_to_symbol={}
    for r in hgnc.to_dict("records"):
        for uid in split_ids(r.get("uniprot_ids","")): uniprot_to_symbol[base_accession(uid)]=r.get("symbol","")
    by_edge=defaultdict(list)
    with PIONEER.open("r",encoding="utf-8",errors="replace",newline="") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            a,b=base_accession(r.get("uniprot1")),base_accession(r.get("uniprot2")); ga,gb=uniprot_to_symbol.get(a,""),uniprot_to_symbol.get(b,"")
            try: ra=ast.literal_eval(r.get("interface_residues1","[]"))
            except Exception: ra=[]
            try: rb=ast.literal_eval(r.get("interface_residues2","[]"))
            except Exception: rb=[]
            for gene,partner,ua,ub,rr,seq,pseq in ((ga,gb,a,b,ra,r.get("sequence1",""),r.get("sequence2","")),(gb,ga,b,a,rb,r.get("sequence2",""),r.get("sequence1",""))):
                if gene in strict_genes and partner:
                    by_edge[(gene,partner)].append({"gene":gene,"partner":partner,"target_uniprot":ua,"partner_uniprot":ub,"source":(r.get("source","") or "unknown").strip(),
                        "interface_residues":sorted({int(x) for x in rr if str(x).isdigit()}),"target_sequence":seq,"partner_sequence":pseq,"confidence":"very_high"})
    return by_edge,uniprot_to_symbol


def rcsb_query_one(pair:tuple[str,str,str,str],cache_dir:Path)->dict[str,Any]:
    gene,partner,ua,ub=pair; pair_key="__".join(sorted([base_accession(ua),base_accession(ub)])); cache=cache_dir/f"{pair_key}.json"
    if cache.exists():
        try:
            d=json.loads(cache.read_text(encoding="utf-8")); return {"pair":pair,"status":d.get("status","CACHED"),"ids":d.get("ids",[]),"error":d.get("error","")}
        except Exception: pass
    query={"query":{"type":"group","logical_operator":"and","nodes":[
        {"type":"terminal","service":"text","parameters":{"attribute":"rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession","operator":"exact_match","value":base_accession(ua)}},
        {"type":"terminal","service":"text","parameters":{"attribute":"rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession","operator":"exact_match","value":base_accession(ub)}}]},
        "return_type":"entry","request_options":{"results_content_type":["experimental"]}}
    result={"status":"QUERY_FAILED","ids":[],"error":""}
    for attempt in range(5):
        try:
            r=requests.post(RCSB_SEARCH,json=query,timeout=(10,50),headers={"User-Agent":"phase1A-allelic-interface/1.0","Content-Type":"application/json"})
            if r.status_code==204: result={"status":"NO_HIT","ids":[],"error":""}; break
            r.raise_for_status(); j=r.json()
            ids=sorted({x.get("identifier","").upper() for x in j.get("result_set",[]) if x.get("identifier")})
            result={"status":"HIT" if ids else "NO_HIT","ids":ids,"total_count":j.get("total_count",len(ids)),"error":""}; break
        except Exception as exc:
            result["error"]=str(exc)[:250]
            if attempt<4: time.sleep(2**attempt)
    cache_dir.mkdir(parents=True,exist_ok=True); cache.write_text(json.dumps({**result,"query":query,"retrieved_utc":utcnow()},indent=2)+"\n",encoding="utf-8")
    return {"pair":pair,**result}


def map_entity_sequence(entity_seq:str,reference_seq:str)->dict[int,int]:
    return sequence_map(entity_seq,reference_seq).get("mapping",{})


def rcsb_assembly_ids(pid:str,cache_dir:Path)->tuple[list[str],str]:
    cache=cache_dir/f"{pid}.json"
    if cache.exists():
        try:
            d=json.loads(cache.read_text(encoding="utf-8"))
            return d.get("assembly_ids",["1"]),d.get("status","CACHED")
        except Exception:
            pass
    try:
        r=requests.get(RCSB_ENTRY.format(pid),timeout=(10,40),headers={"User-Agent":"phase1A-allelic-interface/1.0"})
        r.raise_for_status(); j=r.json()
        ids=j.get("rcsb_entry_container_identifiers",{}).get("assembly_ids",[]) or []
        status="OK" if ids else "ASSEMBLY_LIST_EMPTY_FALLBACK_1"
        ids=[str(x) for x in ids] or ["1"]
    except Exception as exc:
        ids=["1"]; status="ASSEMBLY_METADATA_FAILED_FALLBACK_1:"+str(exc)[:120]
    cache_dir.mkdir(parents=True,exist_ok=True)
    cache.write_text(json.dumps({"PDB_ID":pid,"assembly_ids":ids,"status":status,"retrieved_utc":utcnow()},indent=2)+"\n",encoding="utf-8")
    return ids,status


def parse_cif_contacts(cif_path:Path,target_accession:str,partner_accession:str,target_reference:str,partner_reference:str):
    """Return mapped target-residue nearest distances for matching experimental chains."""
    try: d=MMCIF2Dict(str(cif_path))
    except Exception as exc: return [],f"CIF_PARSE_FAILED:{str(exc)[:160]}"
    refs=defaultdict(set)
    accs=as_list(d.get("_struct_ref_seq.pdbx_db_accession",[])); entities=as_list(d.get("_struct_ref_seq.entity_id",[]))
    for acc,eid in zip(accs,entities): refs[base_accession(acc)].add(str(eid))
    poly_ids=as_list(d.get("_entity_poly.entity_id",[])); poly_seq=as_list(d.get("_entity_poly.pdbx_seq_one_letter_code_can",[]))
    entity_seq={str(eid):re.sub(r"[^A-Za-z]","",str(seq)).upper() for eid,seq in zip(poly_ids,poly_seq)}
    t_entities=refs.get(base_accession(target_accession),set()); p_entities=refs.get(base_accession(partner_accession),set())
    entity_tmap={}; entity_pmap={}
    if not t_entities or not p_entities:
        # Biological-assembly mmCIFs can omit struct_ref_seq. In that case map
        # polymer entities by sequence to the exact queried UniProt references.
        # Ambiguous equal-quality mappings are left unmapped.
        for eid,seq in entity_seq.items():
            tm=sequence_map(seq,target_reference) if target_reference else {"mapping":{},"identity":0,"coverage":0}
            pm=sequence_map(seq,partner_reference) if partner_reference else {"mapping":{},"identity":0,"coverage":0}
            tq=float(tm.get("identity",0))*float(tm.get("coverage",0))
            pq=float(pm.get("identity",0))*float(pm.get("coverage",0))
            tgood=tm.get("identity",0)>=0.85 and tm.get("coverage",0)>=0.70
            pgood=pm.get("identity",0)>=0.85 and pm.get("coverage",0)>=0.70
            if base_accession(target_accession)==base_accession(partner_accession):
                if tgood: t_entities.add(eid); p_entities.add(eid); entity_tmap[eid]=tm["mapping"]; entity_pmap[eid]=pm["mapping"]
            elif tgood and (not pgood or tq>=pq+0.03):
                t_entities.add(eid); entity_tmap[eid]=tm["mapping"]
            elif pgood and (not tgood or pq>=tq+0.03):
                p_entities.add(eid); entity_pmap[eid]=pm["mapping"]
        if not t_entities or not p_entities: return [],"NO_MATCHING_UNIPROT_ENTITY_OR_SEQUENCE"
    else:
        for eid in t_entities:
            entity_tmap[eid]=map_entity_sequence(entity_seq.get(eid,""),target_reference)
        for eid in p_entities:
            entity_pmap[eid]=map_entity_sequence(entity_seq.get(eid,""),partner_reference)
    asym_ids=as_list(d.get("_struct_asym.id",[])); asym_entities=as_list(d.get("_struct_asym.entity_id",[]))
    asym_entity={str(a):str(e) for a,e in zip(asym_ids,asym_entities)}
    t_chains=[c for c,e in asym_entity.items() if e in t_entities]; p_chains=[c for c,e in asym_entity.items() if e in p_entities]
    t_maps={c:entity_tmap.get(asym_entity.get(c,""),{}) for c in t_chains}
    p_maps={c:entity_pmap.get(asym_entity.get(c,""),{}) for c in p_chains}
    atom_asym=as_list(d.get("_atom_site.label_asym_id",[]))
    atoms=defaultdict(lambda:defaultdict(list)); valid_chains=set(t_chains)|set(p_chains)
    for i,chain in enumerate(atom_asym):
        if chain not in valid_chains: continue
        if str(cell(d,"_atom_site.group_PDB",i)).upper() not in {"ATOM","HETATM"}: continue
        if str(cell(d,"_atom_site.type_symbol",i)).upper() in {"H","D"}: continue
        try: pos=int(cell(d,"_atom_site.label_seq_id",i)); xyz=tuple(float(cell(d,k,i)) for k in ["_atom_site.Cartn_x","_atom_site.Cartn_y","_atom_site.Cartn_z"])
        except (ValueError,TypeError): continue
        atoms[str(chain)][pos].append(xyz)
    rows=[]
    for tc in t_chains:
        tm=t_maps.get(tc,{})
        if not tm: continue
        tcoords={tm[pos]:np.asarray(coords,dtype=np.float32) for pos,coords in atoms.get(tc,{}).items() if pos in tm}
        if not tcoords: continue
        # Biological assemblies may contain many repeated copies of each
        # protein. Query all mapped partner-chain atoms together once per target
        # chain, then retain the nearest distance for each canonical target
        # residue. This gives the union interface across partner copies while
        # avoiding a quadratic chain-pair scan in large assemblies.
        partner_coords=[]; partner_chain_by_atom=[]; partner_residue_count={}
        for pc in p_chains:
            if tc==pc: continue
            pm=p_maps.get(pc,{})
            chain_atoms=atoms.get(pc,{})
            if not pm or not any(pos in pm for pos in chain_atoms): continue
            mapped_count=sum(pos in pm for pos in chain_atoms)
            partner_residue_count[pc]=len(pm)
            for pos,coords in chain_atoms.items():
                if pos not in pm: continue
                partner_coords.extend(coords)
                partner_chain_by_atom.extend([pc]*len(coords))
        if not partner_coords: continue
        tree=cKDTree(np.asarray(partner_coords,dtype=np.float32))
        target_points=[]; target_positions=[]
        for cpos,coords in tcoords.items():
            target_points.extend(coords)
            target_positions.extend([int(cpos)]*len(coords))
        distances,nearest=tree.query(np.asarray(target_points,dtype=np.float32),k=1)
        target_positions=np.asarray(target_positions,dtype=np.int32)
        unique_positions,inverse=np.unique(target_positions,return_inverse=True)
        min_distances=np.full(len(unique_positions),np.inf,dtype=np.float32)
        np.minimum.at(min_distances,inverse,np.asarray(distances,dtype=np.float32))
        min_by_pos={int(pos):float(dist) for pos,dist in zip(unique_positions,min_distances)}
        interface_positions=sorted(pos for pos,dst in min_by_pos.items() if dst<=DISTANCE_CUTOFF_A)
        contacting_partner_chains=sorted({partner_chain_by_atom[int(nearest[i])] for i,dst in enumerate(distances) if dst<=DISTANCE_CUTOFF_A})
        rows.append({"target_chain":tc,"partner_chain":";".join(sorted(partner_residue_count)),
            "contacting_partner_chains":";".join(contacting_partner_chains),
            "target_mapped_residues":len(tcoords),"partner_mapped_residues":sum(partner_residue_count.values()),
            "min_distance_A":min(min_by_pos.values()),"interface_distances":min_by_pos,"interface_positions":interface_positions})
    return rows,"OK" if rows else "NO_MAPPED_CHAIN_PAIR"


def structure_audit(primary,strict_genes,ppi_rows,pioneer_by_edge,uni_seq,gene_accs):
    edges={(r["gene"],r["partner"]) for r in ppi_rows if r["gene"] in strict_genes}
    candidates=[]
    for edge,records in pioneer_by_edge.items():
        if edge not in edges: continue
        for r in records:
            if r["source"].upper()=="PDB": candidates.append((*edge,r["target_uniprot"],r["partner_uniprot"]))
    candidates=sorted(set(candidates)); cache_dir=INTERIM/"phase1A_rcsb_pair_cache"; cache_dir.mkdir(parents=True,exist_ok=True)
    query_results=[]
    if candidates:
        with ThreadPoolExecutor(max_workers=4) as pool:
            futures=[pool.submit(rcsb_query_one,p,cache_dir) for p in candidates]
            for fut in as_completed(futures): query_results.append(fut.result())
    query_by_edge=defaultdict(list)
    for q in query_results:
        g,p,ua,ub=q["pair"]; query_by_edge[(g,p)].append(q)
    all_edge_to_ids={e:sorted({pid for q in qs for pid in q.get("ids",[])}) for e,qs in query_by_edge.items()}
    # Pair queries can return many structures containing the same UniProt-linked
    # chains. Review a fixed, data-blind three-entry cap per edge to keep the
    # coordinate audit bounded; retain the complete query hit list in coverage.
    edge_to_ids={e:ids[:3] for e,ids in all_edge_to_ids.items()}
    required_accs={base_accession(x) for p in candidates for x in p[2:4]}
    all_uni={**uni_seq,**fetch_uniprot_sequences(required_accs)}
    pdb_dir=RAW/"phase1A_rcsb_structures"; pdb_dir.mkdir(parents=True,exist_ok=True)
    parsed_dir=INTERIM/"phase1A_parsed_cif"; parsed_dir.mkdir(parents=True,exist_ok=True)
    entry_cache=INTERIM/"phase1A_rcsb_entry_cache"; entry_cache.mkdir(parents=True,exist_ok=True)
    all_ids=sorted({pid for ids in edge_to_ids.values() for pid in ids})
    assembly_meta={}
    with ThreadPoolExecutor(max_workers=6) as pool:
        futures={pool.submit(rcsb_assembly_ids,pid,entry_cache):pid for pid in all_ids}
        for fut in as_completed(futures): assembly_meta[futures[fut]]=fut.result()
    edge_structures=defaultdict(list); structure_rows=[]; errors=[]
    for pid in all_ids:
        accessions_by_edge={}
        for edge,qs in query_by_edge.items():
            if pid not in edge_to_ids.get(edge,[]): continue
            accessions_by_edge[edge]={tuple(q["pair"][2:4]) for q in qs if pid in q.get("ids",[])}
        for assembly_id in assembly_meta.get(pid,(["1"],"ASSEMBLY_METADATA_UNAVAILABLE"))[0]:
            cif=pdb_dir/f"{pid}-assembly{assembly_id}.cif.gz"
            if not cif.exists():
                try:
                    r=requests.get(RCSB_FILES.format(pid,assembly_id),timeout=(15,120),headers={"User-Agent":"phase1A-allelic-interface/1.0"}); r.raise_for_status(); cif.write_bytes(r.content)
                except Exception as exc: errors.append({"PDB_ID":pid,"status":"ASSEMBLY_DOWNLOAD_FAILED","error":str(exc)[:180]}); continue
            parsed_path=parsed_dir/f"{pid}-assembly{assembly_id}.cif"
            if not parsed_path.exists():
                try: parsed_path.write_bytes(gzip.decompress(cif.read_bytes()))
                except Exception as exc: errors.append({"PDB_ID":pid,"status":"ASSEMBLY_GZIP_FAILED","error":str(exc)[:180]}); continue
            for edge,access_pairs in accessions_by_edge.items():
                for ua,ub in sorted(access_pairs):
                    uase=all_uni.get(base_accession(ua),""); ubse=all_uni.get(base_accession(ub),"")
                    if not uase: uase=next((r["target_sequence"] for r in pioneer_by_edge.get(edge,[]) if r["target_uniprot"]==base_accession(ua) and r["target_sequence"]),"")
                    if not ubse: ubse=next((r["partner_sequence"] for r in pioneer_by_edge.get(edge,[]) if r["partner_uniprot"]==base_accession(ub) and r["partner_sequence"]),"")
                    contacts,status=parse_cif_contacts(parsed_path,ua,ub,uase,ubse)
                    if not contacts: errors.append({"PDB_ID":pid,"gene":edge[0],"partner":edge[1],"status":status,"error":f"assembly {assembly_id}"}); continue
                    if not any(c["interface_positions"] for c in contacts):
                        errors.append({"PDB_ID":pid,"gene":edge[0],"partner":edge[1],"status":"NO_HEAVY_ATOM_CONTACT_WITHIN_5A","error":f"assembly {assembly_id}"}); continue
                    for c in contacts:
                        row={"gene":edge[0],"partner":edge[1],"PDB_ID":pid,"assembly_ID":assembly_id,"target_uniprot":base_accession(ua),"partner_uniprot":base_accession(ub),
                            "target_chain":c["target_chain"],"partner_chain":c["partner_chain"],"target_mapped_residues":c["target_mapped_residues"],
                            "partner_mapped_residues":c["partner_mapped_residues"],"min_distance_A":c["min_distance_A"],"interface_positions":c["interface_positions"],"interface_distances":c["interface_distances"]}
                        structure_rows.append(row); edge_structures[edge].append(row)
    # Map the preferred PIONEER interface source to the ClinVar gene's canonical
    # reference; retain per-source residues separately for traceability.
    pioneer_mapped={}
    for edge,records in pioneer_by_edge.items():
        if edge not in edges: continue
        gene=edge[0]; accs=gene_accs.get(gene,[]); canon_acc=accs[0] if accs else ""; canon_seq=uni_seq.get(canon_acc,"")
        mapped=[]
        for r in sorted(records,key=lambda x:(0 if x["source"].upper()=="PDB" else 1 if x["source"].upper()=="PIONEER" else 2,x["target_uniprot"],x["partner_uniprot"])):
            sm=sequence_map(r["target_sequence"],canon_seq) if canon_seq else {"mapping":{},"status":"UNMAPPABLE","identity":0,"coverage":0}
            cpos={int(sm["mapping"][x]) for x in r["interface_residues"] if x in sm.get("mapping",{})}
            mapped.append({**r,"canonical_positions":cpos,"isoform_map_status":sm.get("status","UNMAPPABLE"),"sequence_identity":sm.get("identity",0),"sequence_coverage":sm.get("coverage",0)})
        pioneer_mapped[edge]=mapped
    prov_by_edge={(r["gene"],r["partner"]):r for r in ppi_rows}
    provenance=[]; coverage=[]; iface_rows=[]
    for edge,records in pioneer_mapped.items():
        gene,partner=edge; canon=gene_accs.get(gene,[""])[0] if gene_accs.get(gene) else ""
        preferred=next((r for r in records if r["canonical_positions"]),None)
        for r in records:
            for pos in sorted(r["canonical_positions"]):
                iface_rows.append({"gene":gene,"partner":partner,"interface_source":"PIONEER-derived","PIONEER_reference_source":r["source"],"structure_id":"","canonical_UniProt_accession":canon,
                    "canonical_interface_residue":pos,"residue_mapping_status":r["isoform_map_status"],"confidence":"very_high","min_distance_A":"","distance_cutoff_A":DISTANCE_CUTOFF_A})
    for edge,rows in edge_structures.items():
        gene,partner=edge; canon=gene_accs.get(gene,[""])[0] if gene_accs.get(gene) else ""
        for s in rows:
            for pos in s["interface_positions"]:
                iface_rows.append({"gene":gene,"partner":partner,"interface_source":"experimental PDB coordinates","PIONEER_reference_source":"","structure_id":f"{s['PDB_ID']}:assembly{s['assembly_ID']}","canonical_UniProt_accession":canon,
                    "canonical_interface_residue":pos,"residue_mapping_status":"sequence-aligned to canonical UniProt","confidence":"experimental coordinates",
                    "min_distance_A":f"{s['interface_distances'].get(pos,''):.3f}" if pos in s["interface_distances"] else "","distance_cutoff_A":DISTANCE_CUTOFF_A})
    for edge in sorted(edges):
        r=prov_by_edge[edge]; gene,partner=edge; records=pioneer_mapped.get(edge,[]); has_pioneer=bool(records)
        actual=edge_structures.get(edge,[])
        pdb_ids=sorted({s["PDB_ID"] for s in actual})
        observed_complex=bool(pdb_ids and actual)
        experimental=r.get("experimental_support")=="yes" or observed_complex
        binary=r.get("binary_experimental_support")=="yes" or observed_complex
        if experimental and binary and observed_complex: tier="Tier 1"
        elif experimental and (has_pioneer or observed_complex): tier="Tier 2"
        elif experimental: tier="Tier 3"
        else: tier="Tier 4"
        source_text=r.get("sources",""); pred="yes" if not experimental else "no"
        preferred=next((x for x in records if x["canonical_positions"]),None)
        iface_n=len(preferred["canonical_positions"]) if preferred else 0
        query_status="CONFIRMED_PDB_COMPLEX" if pdb_ids else "QUERIED_NO_PAIR_COMPLEX" if edge in query_by_edge else "NOT_QUERIED_OUTSIDE_PDB_TAGGED_PIONEER_SET"
        level="LEVEL_1_EXPERIMENTAL_COMPLEX" if pdb_ids and actual else "LEVEL_2_PIONEER_INTERFACE" if has_pioneer else "LEVEL_3_PREDICTION_ONLY" if "Predictomes" in source_text else "LEVEL_4_NONE"
        provenance.append({"gene":gene,"partner":partner,"sources":source_text,"source_count":r.get("source_count",""),"experimental_support":"yes" if experimental else "no",
            "binary_experimental_support":"yes (observed PDB complex)" if observed_complex else r.get("binary_experimental_support","no/not established"),"experimental_systems":r.get("experimental_systems",""),
            "PIONEER_coverage":"yes" if has_pioneer else "no","PIONEER_reference_sources":";".join(sorted({x["source"] for x in records})),
            "PIONEER_interface_residue_n":iface_n,"PIONEER_interface_confidence":"very_high" if has_pioneer else "","PDB_complex_search_status":query_status,
            "experimental_PDB_IDs":";".join(pdb_ids),"experimental_PDB_assemblies":";".join(sorted({f"{s['PDB_ID']}:assembly{s['assembly_ID']}" for s in actual})),
            "structure_level":level,"tier":tier,"predicted_only":pred,"included_in_primary_cohort":"YES" if gene in primary else "NO"})
        queried_ids=all_edge_to_ids.get(edge,[])
        reviewed_ids=edge_to_ids.get(edge,[])
        edge_errors=[x for x in errors if x.get("gene")==gene and x.get("partner")==partner]
        coverage.append({"gene":gene,"partner":partner,"tier":tier,"PPI_source":source_text,"structure_level":level,"PDB_pair_search_scope":"PIONEER_source_PDB_pair_subset" if edge in query_by_edge else "not queried",
            "PDB_pair_search_result":";".join(queried_ids),"coordinate_reviewed_PDB_IDs":";".join(reviewed_ids),
            "observed_experimental_complex_IDs":";".join(sorted({f"{s['PDB_ID']}:assembly{s['assembly_ID']}" for s in actual})),"PDB_target_chain_observation_n":len(actual),
            "target_residues_with_coordinates_n":max((s["target_mapped_residues"] for s in actual),default=0),
            "observed_interface_residues_lt_5A_n":len({p for s in actual for p in s["interface_positions"]}),"PIONEER_interface_residues_n":iface_n,
            "canonical_UniProt_accession":gene_accs.get(gene,[""])[0] if gene_accs.get(gene) else "",
            "annotation":"no PDB result is inferred only for the queried subset; unqueried pairs remain unknown; coordinate audit reviewed up to the first three lexically ordered query hits per edge; assembly IDs follow RCSB entry metadata; "+("; ".join(str(x.get("PDB_ID",""))+":"+str(x.get("status","")) for x in edge_errors) if edge_errors else "")})
    return provenance,coverage,iface_rows,structure_rows,errors,edge_structures,pioneer_mapped,query_results


def build_overlap(pioneer_mapped,edge_structures,primary,ppi_rows):
    t12={(r["gene"],r["partner"]) for r in ppi_rows if r.get("tier") in {"Tier 1","Tier 2"}}
    sets={}; sources={}
    for edge in set(pioneer_mapped)|set(edge_structures):
        if edge[0] not in primary: continue
        observed={p for s in edge_structures.get(edge,[]) for p in s["interface_positions"]}
        predicted=set()
        for r in pioneer_mapped.get(edge,[]): predicted.update(r.get("canonical_positions",set()))
        use=observed if observed else predicted
        if use:
            sets[edge]=use
            sources[edge]="experimental PDB coordinates" if observed else "PIONEER high-confidence interface"
    rows=[]; by_gene=defaultdict(list)
    for gene in sorted(primary):
        partners=sorted(p for g,p in sets if g==gene and (g,p) in t12)
        for i,a in enumerate(partners):
            for b in partners[i+1:]:
                A,B=sets[(gene,a)],sets[(gene,b)]
                inter=len(A&B); union=len(A|B); j=inter/union if union else float("nan")
                fa=inter/len(A) if A else float("nan"); fb=inter/len(B) if B else float("nan")
                if min(len(A),len(B))<5: cls="INSUFFICIENT_INTERFACE_RESIDUES"
                elif j>=OVERLAP_HIGH: cls="HIGHLY_OVERLAPPING"
                elif j<=OVERLAP_LOW: cls="SEPARABLE_INTERFACES"
                else: cls="PARTIALLY_OVERLAPPING"
                rows.append({"gene":gene,"partner_A":a,"partner_B":b,"interface_source_A":sources[(gene,a)],"interface_source_B":sources[(gene,b)],
                    "interface_residues_A_n":len(A),"interface_residues_B_n":len(B),"shared_residues_n":inter,"union_residues_n":union,
                    "Jaccard":f"{j:.4f}" if not math.isnan(j) else "","shared_fraction_A":f"{fa:.4f}" if not math.isnan(fa) else "",
                    "shared_fraction_B":f"{fb:.4f}" if not math.isnan(fb) else "","classification":cls,
                    "frozen_thresholds":"HIGH>=0.50; SEPARABLE<=0.10; each interface >=5 mapped residues"})
                by_gene[gene].append(cls)
    gene_class={}
    for gene in primary:
        vals=[x for x in by_gene.get(gene,[]) if x!="INSUFFICIENT_INTERFACE_RESIDUES"]
        if len(vals)<3: gene_class[gene]="INSUFFICIENT_INTERFACE_COMPARISONS"; continue
        n=len(vals)
        if sum(x=="SEPARABLE_INTERFACES" for x in vals)/n>=0.80: gene_class[gene]="SEPARABLE_INTERFACES"
        elif sum(x=="HIGHLY_OVERLAPPING" for x in vals)/n>=0.80: gene_class[gene]="HIGHLY_OVERLAPPING"
        else: gene_class[gene]="PARTIALLY_OVERLAPPING"
    return rows,gene_class,sets


def write_structural_matrix(variants,variant_map,primary,ppi_rows,edge_structures,interface_sets):
    ppi_by_gene=defaultdict(list)
    for r in ppi_rows:
        if r["gene"] in primary: ppi_by_gene[r["gene"]].append(r)
    cols=["gene","VariationID","protein_change","MONDO_groups","canonical_UniProt_accession","canonical_residue","mapping_status","partner","tier","interaction_sources","structure_source","interface_status","interface_min_distance_A","interface_confidence","experimental_PDB_IDs","distance_cutoff_A"]
    def rows_iter():
        for v in variants:
            g=v["gene"]; mp=variant_map.get((g,v["VariationID"]),{})
            if g not in primary or mp.get("mapping_status") not in {"EXACT_CANONICAL","SAFE_ISOFORM_REMAP"}: continue
            try: pos=int(mp.get("canonical_residue") or 0)
            except (ValueError,TypeError): pos=0
            for edge_row in ppi_by_gene.get(g,[]):
                partner=edge_row["partner"]; edge=(g,partner); observed=[]
                for s in edge_structures.get(edge,[]):
                    if pos in s["interface_distances"]: observed.append(float(s["interface_distances"][pos]))
                if observed:
                    dist=min(observed); iface="INTERFACE" if dist<=DISTANCE_CUTOFF_A else "NON_INTERFACE_IN_OBSERVED_COMPLEX"
                    source="experimental PDB coordinates"; confidence="experimental"
                else:
                    pred_positions=interface_sets.get(edge,set())
                    if pred_positions:
                        iface="INTERFACE" if pos in pred_positions else "NON_INTERFACE_IN_PIONEER_INTERFACE_MODEL"
                        dist=""; source="PIONEER interface annotation"; confidence="very_high"
                    else:
                        iface="NO_INTERFACE_MODEL"; dist=""; source=""; confidence=""
                yield {"gene":g,"VariationID":v["VariationID"],"protein_change":v["protein_change"],"MONDO_groups":mp.get("MONDOs",""),
                    "canonical_UniProt_accession":mp.get("canonical_UniProt_accession",""),"canonical_residue":pos,"mapping_status":mp.get("mapping_status",""),
                    "partner":partner,"tier":edge_row.get("tier","Tier 4"),"interaction_sources":edge_row.get("sources",""),
                    "structure_source":source,"interface_status":iface,"interface_min_distance_A":f"{dist:.3f}" if isinstance(dist,(int,float)) else dist,
                    "interface_confidence":confidence,"experimental_PDB_IDs":";".join(sorted({f"{s['PDB_ID']}:assembly{s['assembly_ID']}" for s in edge_structures.get(edge,[])})),
                    "distance_cutoff_A":DISTANCE_CUTOFF_A}
    stream_tsv(OUT/"partner_specific_matrix_structural.tsv",cols,rows_iter())


def per_variant_structural(gene,variant_ids,variant_map,ppi_rows,interface_sets,edge_structures):
    partner_rows=[r for r in ppi_rows if r["gene"]==gene]
    n_partner=len(partner_rows); n_t12=sum(r["tier"] in {"Tier 1","Tier 2"} for r in partner_rows); out={}
    for var in variant_ids:
        v=variant_map.get((gene,var),{})
        try: pos=int(v.get("canonical_residue") or 0)
        except (ValueError,TypeError): pos=0
        models=[]; t12=[]; observed=0
        for r in partner_rows:
            edge=(gene,r["partner"])
            if interface_sets.get(edge):
                models.append(edge)
                if r["tier"] in {"Tier 1","Tier 2"}: t12.append(edge)
            if any(pos in s["interface_distances"] for s in edge_structures.get(edge,[])): observed+=1
        mapped=v.get("mapping_status") in {"EXACT_CANONICAL","SAFE_ISOFORM_REMAP"}
        out[var]={"mapping_status":v.get("mapping_status","UNMAPPABLE"),"structural_partner_n":len(models),"tier1_2_interface_partner_n":len(t12),
            "observed_PDB_partner_n":observed,"usable_ge3_tier1_2":mapped and len(t12)>=3,
            "complete_partner_vector":n_partner>0 and len(models)==n_partner,"at_least_one_interface_partner":mapped and len(models)>0}
    return n_partner,n_t12,out


def build_coverage(strict,primary,incompatible,partner_n,variant_map,ppi_rows,interface_sets,edge_structures,gene_class):
    groups=defaultdict(list)
    for r in strict.to_dict("records"): groups[r["gene"]].append(r)
    variants=defaultdict(list)
    strict_mondo={g:{r["MONDO"] for r in rs} for g,rs in groups.items()}
    for (g,var),r in variant_map.items():
        if set(r.get("MONDOs","").split(";")) & strict_mondo.get(g,set()): variants[g].append(var)
    ppi_by_gene=defaultdict(list)
    for r in ppi_rows: ppi_by_gene[r["gene"]].append(r)
    cov_rows=[]; disease_rows=[]
    for gene in sorted(groups):
        vars_=sorted(set(variants.get(gene,[]))); n_partner,n_t12,perv=per_variant_structural(gene,vars_,variant_map,ppi_rows,interface_sets,edge_structures)
        nvar=len(vars_); mapped=sum(v["mapping_status"] in {"EXACT_CANONICAL","SAFE_ISOFORM_REMAP"} for v in perv.values())
        ge1=sum(v["at_least_one_interface_partner"] for v in perv.values()); ge3=sum(v["usable_ge3_tier1_2"] for v in perv.values())
        complete=sum(v["complete_partner_vector"] for v in perv.values()); pdb_addr=sum(v["observed_PDB_partner_n"]>0 for v in perv.values())
        mapping_pct=100*mapped/nvar if nvar else 0; usable_pct=100*ge3/nvar if nvar else 0
        pass_core=(len(groups[gene])>=2 and gene not in incompatible and partner_n.get(gene,0)>=3 and mapping_pct>=70 and usable_pct>=70)
        cov_rows.append({"gene":gene,"primary_cohort_eligible":"YES" if gene in primary else "NO",
            "mechanism_stratum":"KNOWN_INCOMPATIBLE" if gene in incompatible else "UNKNOWN_COMPATIBILITY","strict_disease_group_n":len(groups[gene]),
            "unique_strict_group_variants_n":nvar,"protein_mappable_variants_n":mapped,"protein_mapping_rate_pct":f"{mapping_pct:.2f}",
            "credible_partner_n":partner_n.get(gene,0),"Tier1_partner_n":sum(r["tier"]=="Tier 1" for r in ppi_by_gene.get(gene,[])),
            "Tier2_partner_n":sum(r["tier"]=="Tier 2" for r in ppi_by_gene.get(gene,[])),"Tier1_2_partner_n":n_t12,
            "PIONEER_interface_partner_n":sum(bool(interface_sets.get((gene,r["partner"]))) for r in ppi_by_gene.get(gene,[])),
            "variants_with_ge1_partner_interface_n":ge1,"variants_with_ge1_partner_interface_pct":f"{100*ge1/nvar:.2f}" if nvar else "0.00",
            "variants_with_ge3_Tier1_2_partner_interfaces_n":ge3,"variants_with_ge3_Tier1_2_partner_interfaces_pct":f"{usable_pct:.2f}",
            "variants_with_complete_partner_vector_n":complete,"variants_with_complete_partner_vector_pct":f"{100*complete/nvar:.2f}" if nvar else "0.00",
            "variants_with_observed_PDB_residue_coordinates_n":pdb_addr,"fraction_MutPred_PPI_structurally_addressable_pct":f"{100*pdb_addr/nvar:.2f}" if nvar else "0.00",
            "fraction_structurally_addressable_pct":f"{usable_pct:.2f}","interface_overlap_class":gene_class.get(gene,"INSUFFICIENT_INTERFACE_COMPARISONS"),
            "passes_ge70_representation_threshold":"YES" if pass_core else "NO",
            "failure_reasons":";".join(x for x,c in [("FEWER_THAN_2_DISEASE_GROUPS",len(groups[gene])<2),("KNOWN_INCOMPATIBLE",gene in incompatible),
                ("FEWER_THAN_3_CREDIBLE_PARTNERS",partner_n.get(gene,0)<3),("PROTEIN_MAPPING_LT70",mapping_pct<70),("TIER1_2_STRUCTURAL_REPRESENTATION_LT70",usable_pct<70)] if c)})
        for d in groups[gene]:
            mondo=d["MONDO"]; ids={var for (g,var),v in variant_map.items() if g==gene and mondo in v.get("MONDOs","").split(";")}
            n=len(ids); m=sum(perv.get(var,{}).get("mapping_status") in {"EXACT_CANONICAL","SAFE_ISOFORM_REMAP"} for var in ids)
            ge3d=sum(perv.get(var,{}).get("usable_ge3_tier1_2",False) for var in ids)
            disease_rows.append({"gene":gene,"MONDO":mondo,"MONDO_label":d.get("MONDO_label",""),"strict_group_variant_n":n,"mappable_variant_n":m,
                "mappable_rate_pct":f"{100*m/n:.2f}" if n else "0.00","variants_with_ge3_Tier1_2_interface_partners_n":ge3d,
                "ge3_partner_representation_rate_pct":f"{100*ge3d/n:.2f}" if n else "0.00"})
    return cov_rows,disease_rows


def build_differential_missingness(disease_rows,primary):
    by_gene=defaultdict(list)
    for r in disease_rows:
        by_gene[r["gene"]].append(r)
    out=[]
    for gene,rows in sorted(by_gene.items()):
        for metric,vals in [("protein_mappability",[float(x["mappable_rate_pct"]) for x in rows]),
                            ("ge3_Tier1_2_partner_representation",[float(x["ge3_partner_representation_rate_pct"]) for x in rows])]:
            spread=max(vals)-min(vals) if vals else 0.0
            flag=(spread>=DIFFERENTIAL_COVERAGE_PP or (min(vals)<LOW_COVERAGE_CUTOFF and max(vals)>=LOW_COVERAGE_CUTOFF)) if vals else False
            out.append({"gene":gene,"coverage_metric":metric,"strict_disease_group_n":len(rows),
                "minimum_disease_coverage_pct":f"{min(vals):.2f}" if vals else "0.00","maximum_disease_coverage_pct":f"{max(vals):.2f}" if vals else "0.00",
                "absolute_differential_percentage_points":f"{spread:.2f}","flag":"DIFFERENTIAL_COVERAGE" if flag else "NO_DIFFERENTIAL_COVERAGE",
                "primary_cohort_eligible":"YES" if gene in primary else "NO",
                "flag_rule":"range >=20 percentage points OR one disease <70% while another >=70%",
                "disease_group_rates":";".join(f"{x['MONDO']}={x['mappable_rate_pct'] if metric=='protein_mappability' else x['ge3_partner_representation_rate_pct']}%" for x in rows)})
    return out


def build_exemplars(strict,coverage,overlap_classes,primary):
    by_gene=defaultdict(list)
    for r in strict.to_dict("records"): by_gene[r["gene"]].append(r)
    cov={r["gene"]:r for r in coverage}; sep_rank={"SEPARABLE_INTERFACES":0,"PARTIALLY_OVERLAPPING":1,"HIGHLY_OVERLAPPING":2,"INSUFFICIENT_INTERFACE_COMPARISONS":3}
    rows=[]
    for gene in primary:
        c=cov.get(gene,{}); group_counts=[int(x["n_unique_variants"]) for x in by_gene.get(gene,[])]
        rows.append({"gene":gene,"strict_disease_group_n":len(group_counts),"unique_strict_group_variants_n":c.get("unique_strict_group_variants_n",0),
            "minimum_variants_per_strict_disease":min(group_counts) if group_counts else 0,"maximum_variants_per_strict_disease":max(group_counts) if group_counts else 0,
            "disease_balance_min_over_max":f"{min(group_counts)/max(group_counts):.4f}" if group_counts and max(group_counts)>0 else "0.0000",
            "Tier1_2_partner_n":c.get("Tier1_2_partner_n",0),"protein_mapping_rate_pct":c.get("protein_mapping_rate_pct","0"),
            "ge3_Tier1_2_interface_variant_rate_pct":c.get("variants_with_ge3_Tier1_2_partner_interfaces_pct","0"),
            "interface_overlap_class":overlap_classes.get(gene,"INSUFFICIENT_INTERFACE_COMPARISONS"),
            "selection_rule":"rank by >=3 Tier1-2 interface variant coverage, protein mapping, Tier1-2 partner count, interface separability, minimum strict-disease variant count; no disease-separation statistic"})
    rows.sort(key=lambda r:(-float(r["ge3_Tier1_2_interface_variant_rate_pct"]),-float(r["protein_mapping_rate_pct"]),-int(r["Tier1_2_partner_n"]),
        sep_rank.get(r["interface_overlap_class"],4),-int(r["minimum_variants_per_strict_disease"]),r["gene"]))
    for i,r in enumerate(rows[:10],1): r["rank"]=i
    columns=["rank","gene","strict_disease_group_n","unique_strict_group_variants_n","minimum_variants_per_strict_disease","maximum_variants_per_strict_disease","disease_balance_min_over_max",
        "Tier1_2_partner_n","protein_mapping_rate_pct","ge3_Tier1_2_interface_variant_rate_pct","interface_overlap_class","selection_rule"]
    return rows[:10],columns


def training_overlap_audit(variant_map,primary):
    root=INTERIM/"phase1A_mutpred_training"; files=sorted(p for p in root.rglob("*_rows.csv.gz") if not p.name.startswith("._")) if root.exists() else []
    train_variant=set(); train_pair=set(); train_protein=set(); train_pair_lookup=defaultdict(set); datasets=set(); row_n=0
    for path in files:
        try:
            for frame in pd.read_csv(path,compression="gzip",dtype=str,keep_default_na=False,chunksize=100000):
                row_n+=len(frame); datasets.add(path.name.replace("_rows.csv.gz",""))
                if not {"interactor","partner","mutation"}.issubset(frame.columns): continue
                for a,b,mut in frame[["interactor","partner","mutation"]].itertuples(index=False,name=None):
                    a=base_accession(a); b=base_accession(b); mut=str(mut).strip()
                    train_pair.add((a,b,mut)); train_pair_lookup[(a,mut)].add(b); train_variant.add((a,mut)); train_protein.update([a,b])
        except Exception: continue
    out=[]
    for (gene,var),v in sorted(variant_map.items()):
        if gene not in primary: continue
        acc=base_accession(v.get("canonical_UniProt_accession",""))
        try: mut=f"{v.get('reference_aa','')}{int(v.get('canonical_residue') or 0)}{v.get('alternate_aa','')}"
        except (ValueError,TypeError): mut=""
        variant_hit=bool(acc and mut and (acc,mut) in train_variant)
        training_partners=sorted(train_pair_lookup.get((acc,mut),set()))
        prot_hit=bool(acc and acc in train_protein); flags=[]
        if variant_hit: flags.append("TRAIN_VARIANT_OVERLAP")
        if training_partners: flags.append("TRAIN_PAIR_OVERLAP")
        if prot_hit: flags.append("TRAIN_PROTEIN_OVERLAP")
        if not flags: flags.append("NO_KNOWN_OVERLAP")
        out.append({"gene":gene,"VariationID":var,"canonical_UniProt_accession":acc,"canonical_variant":mut,"mapping_status":v.get("mapping_status",""),
            "training_overlap_flags":";".join(flags),"TRAIN_VARIANT_OVERLAP":"YES" if variant_hit else "NO","TRAIN_PAIR_OVERLAP":"YES" if training_partners else "NO",
            "TRAIN_PROTEIN_OVERLAP":"YES" if prot_hit else "NO","training_partner_accessions":";".join(training_partners),"training_partner_n":len(training_partners),
            "training_rows_available_n":row_n,"training_datasets":";".join(sorted(datasets)) or "NOT_AVAILABLE",
            "audit_scope":"public Zenodo training/evaluation row tables; unpublished VarChAMP rows are not auditable",
            "model_generation":"current repository checkpoint; score generation is a separate step"})
    cols=["gene","VariationID","canonical_UniProt_accession","canonical_variant","mapping_status","training_overlap_flags","TRAIN_VARIANT_OVERLAP","TRAIN_PAIR_OVERLAP","TRAIN_PROTEIN_OVERLAP",
        "training_partner_accessions","training_partner_n","training_rows_available_n","training_datasets","audit_scope","model_generation"]
    write_tsv(OUT/"mutpred_training_overlap.tsv",out,cols)
    return out,{"training_dataset_archive_present":bool(files),"training_datasets":sorted(datasets),"training_row_n":row_n,
        "training_files":[str(x.relative_to(ROOT)) for x in files],
        "training_file_sha256":{str(x.relative_to(ROOT)):sha256(x) for x in files},
        "unpublished_varchamp_auditable":False,"all_known_data_audited":bool(files)}


def build_disruption_matrix():
    # No score is inferred from structural occupancy. Rows will be added only
    # after a model-specific input can be prepared and inference completes.
    write_tsv(OUT/"partner_specific_matrix_disruption.tsv",[],
        ["gene","VariationID","partner","model_version_primary","score","score_status","complex_structure_ID","training_overlap_flags","note"])


def actb_actg1_summary():
    path=PHASE0/"ACTB_ACTG1_feasibility_exact_rcv.tsv"
    if not path.exists(): return []
    return pd.read_csv(path,sep="\t",dtype=str,keep_default_na=False).to_dict("records")


def write_report(meta,strict,coverage,disease_rows,differential,exemplars,provenance,structure_cov,structure_rows,training_meta,actb):
    OUT.mkdir(parents=True,exist_ok=True); REPORT.parent.mkdir(parents=True,exist_ok=True)
    primary=set(meta["primary_genes"]); incompatible=set(meta["incompatible_genes"])
    pass_genes=sorted(r["gene"] for r in coverage if r.get("passes_ge70_representation_threshold")=="YES" and r["gene"] in primary)
    map_total=sum(int(r.get("unique_strict_group_variants_n") or 0) for r in coverage if r["gene"] in primary)
    mapped=sum(int(r.get("protein_mappable_variants_n") or 0) for r in coverage if r["gene"] in primary)
    mapped_pct=100*mapped/map_total if map_total else 0
    ge3=sum(int(r.get("variants_with_ge3_Tier1_2_partner_interfaces_n") or 0) for r in coverage if r["gene"] in primary)
    structural_pct=100*ge3/map_total if map_total else 0
    mutpred_addressable=sum(int(r.get("variants_with_observed_PDB_residue_coordinates_n") or 0) for r in coverage if r["gene"] in primary)
    mutpred_pct=100*mutpred_addressable/map_total if map_total else 0
    t12={r["gene"]:int(r["Tier1_2_partner_n"]) for r in coverage if r["gene"] in primary}; ge3_t12=sum(x>=3 for x in t12.values())
    differential_genes=sorted({r["gene"] for r in differential if r["flag"]=="DIFFERENTIAL_COVERAGE" and r.get("primary_cohort_eligible")=="YES"})
    sep=sorted(r["gene"] for r in coverage if r.get("interface_overlap_class")=="SEPARABLE_INTERFACES" and r["gene"] in primary)
    if len(pass_genes)>=30: verdict="GO_PHASE1B_FULL"
    elif len(pass_genes)>=10: verdict="GO_PHASE1B_FOCUSED"
    else: verdict="TECHNICAL_HOLD"
    exemplar_text="\n".join(f"{r['rank']}. **{r['gene']}** — {r['unique_strict_group_variants_n']} variants across {r['strict_disease_group_n']} diseases; {r['Tier1_2_partner_n']} Tier 1–2 partners; mapping {r['protein_mapping_rate_pct']}%; ≥3 Tier 1–2 interface representation {r['ge3_Tier1_2_interface_variant_rate_pct']}%; {r['interface_overlap_class']}." for r in exemplars) or "No candidates."
    diff_summary=f"{len(differential_genes)} genes: "+(", ".join(differential_genes) if differential_genes else "none")
    t1=sum(r["tier"]=="Tier 1" and r["gene"] in primary for r in provenance); t2=sum(r["tier"]=="Tier 2" and r["gene"] in primary for r in provenance)
    t3=sum(r["tier"]=="Tier 3" and r["gene"] in primary for r in provenance); t4=sum(r["tier"]=="Tier 4" and r["gene"] in primary for r in provenance)
    pdb_q=sum(1 for r in structure_cov if r["gene"] in primary and r["PDB_pair_search_scope"]!="not queried")
    pdb_hits=len({r["PDB_ID"] for r in structure_rows if r["gene"] in primary})
    pdb_assemblies=len({(r["PDB_ID"],r["assembly_ID"]) for r in structure_rows if r["gene"] in primary})
    mutpred_pair_n=0
    matrix_path=OUT/"partner_specific_matrix_structural.tsv"
    if matrix_path.exists():
        with matrix_path.open(encoding="utf-8",newline="") as handle:
            mutpred_pair_n=sum(1 for r in csv.DictReader(handle,delimiter="\t") if r.get("structure_source")=="experimental PDB coordinates")
    score_path=OUT/"partner_specific_matrix_disruption.tsv"
    scored_pair_n=0
    if score_path.exists():
        with score_path.open(encoding="utf-8",newline="") as handle:
            scored_pair_n=max(0,sum(1 for _ in handle)-1)
    mutpred_pair_pct=100*scored_pair_n/mutpred_pair_n if mutpred_pair_n else 0.0
    mapping_fail=sorted(r["gene"] for r in coverage if r["gene"] in primary and float(r["protein_mapping_rate_pct"])<70)
    low_ppi=sorted(r["gene"] for r in coverage if r["gene"] in primary and int(r["credible_partner_n"])<3)
    low_struct=sorted(r["gene"] for r in coverage if r["gene"] in primary and float(r["variants_with_ge3_Tier1_2_partner_interfaces_pct"])<70)
    actb_text="Phase 0R exact-RCV pilot records: "+json.dumps(actb,ensure_ascii=False)
    body=f"""# Phase 1A — partner-resolved representation feasibility

Run date (UTC): {utcnow()}

## Verdict

**{verdict}**

The frozen Phase 0R baseline was reused without revision. The true Phase 1A intersection is **{len(primary)} genes**: {meta['strict_genes']} strict genes, excluding {len(incompatible)} explicitly CONFIRMED_INCOMPATIBLE genes, followed by the ≥3 credible-partner criterion. Phase 1A is a representation and missingness audit only. No disease-label association, classifier, disease-specific enrichment, UMAP interpretation, or phenotype prediction was run.

## Cohort and mechanism

- Strict starting genes: {meta['strict_genes']}; strict ClinVar disease groups: {len(strict)}.
- Explicitly incompatible genes excluded from primary: {len(incompatible)} ({", ".join(sorted(incompatible))}).
- Primary cohort: {len(primary)} genes; {meta['strict_ge3_credible_partners']} of all strict genes have ≥3 credible partners before mechanism exclusion.
- Partner-count bands among non-incompatible genes: {json.dumps(meta['partner_bands_after_mechanism_exclusion'],ensure_ascii=False)}. All-strict bands: {json.dumps(meta['partner_bands_all_strict'],ensure_ascii=False)}.
- Pairwise mechanism states remain in cohort_intersection.tsv; UNKNOWN_COMPATIBILITY has not been recoded as compatible.

## PPI and structure reuse

Primary-cohort credible edges by frozen hierarchy: Tier 1 **{t1}**, Tier 2 **{t2}**, Tier 3 **{t3}**, Tier 4 **{t4}**. Thus {ge3_t12} genes have ≥3 Tier 1–2 partners. Tier 4 remains explicitly prediction-only.

RCSB pairwise queries covered {pdb_q} edges from the subset already tagged source=PDB in the frozen very-high-confidence PIONEER table; the coordinate audit reviewed up to three lexically ordered search hits per edge and used RCSB-listed biological assemblies. {pdb_hits} unique PDB entry IDs and {pdb_assemblies} biological assemblies yielded mapped target-partner contacts. Within each assembly, mapped partner-chain copies were combined and each target residue was assigned its nearest heavy-atom distance to any mapped partner copy. `PDB_target_chain_observation_n` records target-chain observations; assembly IDs are retained per edge. The remaining pairs outside that query set are marked unqueried, not as having no structure. Interface membership uses a frozen **5.0 Å heavy-atom** cutoff for experimental coordinates; PIONEER interface residues are separately labeled as very-high-confidence predictions. New complex modelling and FoldX were not run.

## Protein mapping and representation

- Unique ClinVar variants across primary genes: {map_total}; EXACT_CANONICAL + SAFE_ISOFORM_REMAP: {mapped} ({mapped_pct:.2f}%).
- Variants with ≥3 Tier 1–2 partner-specific interface representations: {ge3} ({structural_pct:.2f}%).
- Variants with ≥1 observed experimental PDB residue mapping suitable for a future interaction-specific model input: {mutpred_addressable} ({mutpred_pct:.2f}%).
- MutPred-PPI scores generated: {scored_pair_n}/{mutpred_pair_n} structurally addressable variant-partner rows ({mutpred_pair_pct:.2f}%). The score matrix is header-only; structural addressability is not treated as a disruption score.
- Genes passing the ≥70% protein-mapping and ≥70% ≥3-partner-representation thresholds: {len(pass_genes)}.
- Genes below the protein-mapping threshold: {len(mapping_fail)}; genes below ≥3 credible partners within the primary cohort: {len(low_ppi)}; genes below the Tier 1–2 structural representation threshold: {len(low_struct)}.
- Differential coverage flagged by the frozen rule: {diff_summary}.

## MutPred-PPI

The official MutPred-PPI repository was pinned at commit 01a53eaf5661382a97a9b05a0b41e7f699cd0974; the current repository checkpoint is the primary generation. Its published v1.0 ensemble is a separate sensitivity generation and scores are not comparable. The training overlap table audits recoverable public rows; the repository documents an unpublished VarChAMP component, so an absent match means no known overlap in the public audit, not proof of no overlap in all training data. MutPred-PPI inference was not run in this pass: the {mutpred_pair_n} residue-addressable variant-partner rows remain candidates, and each requires a matching mmCIF-derived contact graph and ProtT5 embedding before the frozen current checkpoint can score it. The score matrix is header-only; no score is imputed from interface occupancy.

Training rows recovered: {training_meta.get('training_row_n',0)} across {", ".join(training_meta.get('training_datasets',[])) or "none"}. Public training audit files: {", ".join(training_meta.get('training_files',[])) or "not recovered"}.

## Interface overlap and future exemplars

Overlap thresholds were frozen before any phenotype analysis: Jaccard ≥{OVERLAP_HIGH:.2f} is highly overlapping; Jaccard ≤{OVERLAP_LOW:.2f} is separable if both interfaces have ≥5 mapped residues; intermediate values are partial. Gene-level labels require ≥3 comparisons and ≥80% of pairwise comparisons in the corresponding class. Separable genes: {len(sep)}.

Top 10 data-blind candidates use only variant count/balance, Tier 1–2 PPI quality, mapping coverage, ≥3-partner structural coverage, and interface overlap:

{exemplar_text}

## ACTB and ACTG1 secondary cases

{actb_text}

These remain historical/biological comparison cases and were not forced into the strict primary cohort. Their low variant counts within individual disease groups do not support quantitative primary disease-label testing without further data.

## Methods, thresholds, provenance, and limitations

1. Exact cohort source: results/phase0R/gene_disease_variant_counts_exact_rcv.tsv, restricted to strict_disease_ge10=YES; mechanism state from the frozen Phase 0R gene/pair tables; PPI counts/provenance from the Phase 0R frozen PPI outputs.
2. Variant-to-RCV-to-MONDO linkage uses the Phase 0R seed and recovered RCV table, requiring exact gene, VariationID, RCV accession, and one of the strict MONDO groups. Transcript HGVS was parsed from the frozen ClinVar variant_summary_2026-09.txt.gz. RefSeq protein sequence was obtained from NCBI E-utilities and aligned to HGNC-associated UniProt reference sequences. Equal-best alignments that map a residue differently are not forced.
3. PIONEER snapshot: PIONEER_very_high_human_2026-09-27.txt. Experimental evidence: BioGRID MV-Physical 5.0.261 and cached IntAct PSICQUIC evidence recorded in Phase 0R. HuRI remained unavailable as recorded in Phase 0R.
4. RCSB Search API exact UniProt accession pair query, entry assembly metadata, and biological assembly mmCIF files are cached under data/interim/phase1A_rcsb_pair_cache/, data/interim/phase1A_rcsb_entry_cache/, and data/raw/phase1A_rcsb_structures/. A focused pair-query scope is stated per edge in structure_coverage.tsv. The first three lexically ordered hits per edge were reviewed to bound coordinate processing while retaining the full search result list. Partner-chain copies were aggregated within an assembly; `PDB_target_chain_observation_n` records target-chain observations and the 5.0 Å distance is the minimum to any mapped partner-chain copy.
5. Coverage hard flag: within a gene, flag when disease-group coverage spans ≥{DIFFERENTIAL_COVERAGE_PP:.0f} percentage points, or at least one group is <{LOW_COVERAGE_CUTOFF:.0f}% while another is ≥{LOW_COVERAGE_CUTOFF:.0f}%. Reported separately for protein mappability and ≥3 Tier 1–2 interface partners.
6. Interface definitions do not estimate functional effect. MutPred-PPI predictions are interaction-specific functional scores; none were substituted by structural occupancy.
7. Public source endpoints: [ClinVar monthly variant summary](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2026-09.txt.gz), [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/), [UniProt REST API](https://www.uniprot.org/help/api), [RCSB Search API](https://search.rcsb.org/), [PIONEER frozen human predictions](https://pioneer.yulab.org/static/predictions/very_high/human.txt), [MutPred-PPI source and model documentation](https://github.com/rosstewart/mutpred-ppi), and [public MutPred-PPI training-data deposit](https://doi.org/10.5281/zenodo.22774603).
8. The main limitation is incomplete structure search outside the PIONEER source=PDB subset and unavailable unpublished VarChAMP training rows. The RCSB no-hit category is therefore limited to queried edges. Phase 1B admission should preserve this uncertainty and the differential-missingness flags.

## Output files

All twelve required result tables are in results/phase1A/. Re-run with python3 scripts/run_phase1a.py. The frozen inputs' SHA-256 values and analysis date are recorded in metadata/phase1A_manifest.json.
"""
    REPORT.write_text(body,encoding="utf-8")


def main():
    OUT.mkdir(parents=True,exist_ok=True); INTERIM.mkdir(parents=True,exist_ok=True)
    strict,mg,mp,cov,phase0_prov,seed,rec,hgnc=load_inputs()
    cohort_rows,primary,incompatible,partner_n,strata,meta=make_cohort(strict,mg,mp,cov)
    write_tsv(OUT/"cohort_intersection.tsv",cohort_rows)
    print(f"cohort: strict={meta['strict_genes']} incompatible={len(incompatible)} primary={len(primary)}",flush=True)
    strict_genes=set(strict["gene"]); strict_groups=set(map(tuple,strict[["gene","MONDO"]].values))
    rec_copy=rec.copy(); rec_copy["MONDO_PRIMARY"]=rec_copy["MONDO"].apply(lambda x:next((v for v in str(x).split(";") if v.startswith("MONDO:")),""))
    seed_join=seed.merge(rec_copy[["RCVaccession","VariationID","MONDO_PRIMARY"]],on=["RCVaccession","VariationID"],how="inner")
    seed_join=seed_join[(seed_join["GeneSymbol"].isin(strict_genes))&seed_join.apply(lambda r:(r["GeneSymbol"],r["MONDO_PRIMARY"]) in strict_groups,axis=1)]
    variant_keys=set((r.GeneSymbol,str(r.VariationID)) for r in seed_join.itertuples())
    print(f"strict exact variants: {len(variant_keys)}; scanning frozen ClinVar summary for transcripts",flush=True)
    transcripts_by_key=parse_clinvar_transcripts(variant_keys)
    map_rows,variant_map,uni_seq,gene_accs,seed_rec,refseq=map_variants(strict,mg,seed,rec,hgnc,transcripts_by_key,primary)
    print(f"protein mapping rows: {len(map_rows)}; fetched transcript proteins: {sum(bool(x.get('sequence')) for x in refseq.values())}",flush=True)
    pioneer_by_edge,uniprot_to_symbol=build_pioneer(strict_genes,hgnc)
    ppi_rows=[]
    ppi_df=phase0_prov[phase0_prov["gene"].isin(strict_genes)]
    for r in ppi_df.to_dict("records"):
        ppi_rows.append({"gene":r["gene"],"partner":r["partner"],"sources":r.get("sources",""),"source_count":r.get("source_count",""),
            "experimental_support":r.get("experimental_support","no"),"binary_experimental_support":r.get("binary_experimental_support","no/not established"),
            "experimental_systems":r.get("experimental_systems",""),"phase0R_tier_field":r.get("tier","")})
    ppi_by_edge={(r["gene"],r["partner"]):r for r in ppi_rows}
    for edge,r in ppi_by_edge.items(): r["PIONEER_coverage"]="yes" if edge in pioneer_by_edge else "no"
    print(f"PPI edges: {len(ppi_rows)}; PIONEER-annotated edges: {sum(r['PIONEER_coverage']=='yes' for r in ppi_rows)}; querying focused RCSB subset",flush=True)
    ppi_final,structure_cov,iface_rows,structure_rows,structure_errors,edge_structures,pioneer_mapped,query_results=structure_audit(primary,strict_genes,ppi_rows,pioneer_by_edge,uni_seq,gene_accs)
    prov_cols=["gene","partner","sources","source_count","experimental_support","binary_experimental_support","experimental_systems","PIONEER_coverage","PIONEER_reference_sources",
        "PIONEER_interface_residue_n","PIONEER_interface_confidence","PDB_complex_search_status","experimental_PDB_IDs","experimental_PDB_assemblies","structure_level","tier","predicted_only","included_in_primary_cohort"]
    write_tsv(OUT/"ppi_provenance.tsv",ppi_final,prov_cols)
    structure_cols=["gene","partner","tier","PPI_source","structure_level","PDB_pair_search_scope","PDB_pair_search_result","coordinate_reviewed_PDB_IDs","observed_experimental_complex_IDs","PDB_target_chain_observation_n",
        "target_residues_with_coordinates_n","observed_interface_residues_lt_5A_n","PIONEER_interface_residues_n","canonical_UniProt_accession","annotation"]
    write_tsv(OUT/"structure_coverage.tsv",structure_cov,structure_cols)
    interface_cols=["gene","partner","interface_source","PIONEER_reference_source","structure_id","canonical_UniProt_accession","canonical_interface_residue","residue_mapping_status","confidence","min_distance_A","distance_cutoff_A"]
    write_tsv(OUT/"interface_mapping.tsv",iface_rows,interface_cols)
    overlap_rows,gene_class,interface_sets=build_overlap(pioneer_mapped,edge_structures,primary,ppi_final)
    write_tsv(OUT/"interface_overlap.tsv",overlap_rows)
    variants=[{"gene":r["gene"],"VariationID":r["VariationID"],"protein_change":r["protein_change"]} for r in map_rows]
    t12edges={(r["gene"],r["partner"]) for r in ppi_final if r["tier"] in {"Tier 1","Tier 2"}}
    gate_interface_sets={e:s for e,s in interface_sets.items() if e in t12edges}
    write_structural_matrix(variants,variant_map,primary,ppi_final,edge_structures,interface_sets)
    cov_rows,disease_rows=build_coverage(strict,primary,incompatible,partner_n,variant_map,ppi_final,gate_interface_sets,edge_structures,gene_class)
    write_tsv(OUT/"coverage_by_gene.tsv",cov_rows)
    differential=build_differential_missingness(disease_rows,primary)
    write_tsv(OUT/"differential_missingness.tsv",differential)
    exemplar_rows,exemplar_cols=build_exemplars(strict,cov_rows,gene_class,primary)
    write_tsv(OUT/"exemplar_candidates.tsv",exemplar_rows,exemplar_cols)
    training_rows,training_meta=training_overlap_audit(variant_map,primary)
    build_disruption_matrix()
    actb=actb_actg1_summary()
    manifest={"phase":"PHASE 1A","created_utc":utcnow(),"verdict_baseline":"GO_PHASE1_CONDITIONAL","phase0R_frozen":True,"code":"scripts/run_phase1a.py","python":__import__("platform").python_version(),
        "packages":{"pandas":pd.__version__,"numpy":np.__version__,"biopython":__import__("Bio").__version__,"requests":requests.__version__},
        "thresholds":{"interface_heavy_atom_distance_A":DISTANCE_CUTOFF_A,"overlap_jaccard_high":OVERLAP_HIGH,"overlap_jaccard_low":OVERLAP_LOW,
            "differential_coverage_pp":DIFFERENTIAL_COVERAGE_PP,"low_coverage_pct":LOW_COVERAGE_CUTOFF,"gene_minimum_mappable_pct":70,"gene_minimum_ge3_interface_pct":70},
        "cohort":meta,"strict_disease_group_n":len(strict),"mapped_variant_row_n":len(map_rows),"clinvar_transcript_accession_n":len({x for vals in transcripts_by_key.values() for x in vals}),
        "rcsb_pair_query_n":len(query_results),"rcsb_confirmed_pdb_entry_n":len({r["PDB_ID"] for r in structure_rows}),"structure_parse_error_n":len(structure_errors),"training_overlap":training_meta,
        "frozen_inputs":{str(p.relative_to(ROOT)):{"bytes":p.stat().st_size,"sha256":sha256(p)} for p in [STRICT_DISEASES,MECH_GENES,MECH_PAIRS,PPI_COVERAGE,PPI_PROVENANCE,SEED,RECOVERED,VARIANT_SUMMARY,PIONEER,HGNC]}}
    (ROOT/"metadata/phase1A_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    write_report(meta,strict,cov_rows,disease_rows,differential,exemplar_rows,ppi_final,structure_cov,structure_rows,training_meta,actb)
    print(f"Phase1A tables written to {OUT}; report={REPORT}",flush=True)


if __name__=="__main__":
    main()
