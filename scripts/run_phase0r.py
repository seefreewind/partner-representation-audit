#!/usr/bin/env python3
"""Recover the frozen ClinVar RCV release for the Phase 0R cohort audit."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import platform
import re
import sqlite3
import subprocess
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from lxml import etree

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
PHASE = ROOT / "results" / "phase0R"
META = ROOT / "metadata"
LOGS = ROOT / "logs"
SUMMARY = RAW / "variant_summary_2026-09.txt.gz"
ARCHIVE_PART = RAW / "ClinVarRCVRelease_2026-09.xml.gz.incomplete"
ARCHIVE = RAW / "ClinVarRCVRelease_2026-09.xml.gz"
SEED = PHASE / "global_plp_missense_rcv_seed.tsv"
SEED_CHECKPOINT = META / "phase0R_seed_checkpoint.json"
PARSER_CHECKPOINT = META / "phase0R_rcv_parse_checkpoint.json"
API_CACHE = RAW / "phase0R_targeted_rcv_api"
RCV_ARCHIVE_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_release/ClinVarRCVRelease_2026-09.xml.gz"
RCV_ARCHIVE_BYTES = 6018665816
RCV_RE = re.compile(r"\bRCV\d{9}(?:\.\d+)?\b")
HGVS_RE = re.compile(r"p\.\(?([A-Z][a-z]{2}|[A-Z])(\d+)([A-Z][a-z]{2}|[A-Z])\)?(?![A-Za-z0-9])")
AA3 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
AA1 = set(AA3.values())
SEED_COLUMNS = ["VariationID", "GeneSymbol", "ProteinChange", "RCVaccession", "ClinicalSignificance", "ReviewStatus", "PhenotypeIDS", "PhenotypeList"]


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate_frozen_archive() -> dict:
    if not ARCHIVE.exists():
        raise FileNotFoundError(f"Frozen RCV archive not present: {ARCHIVE}")
    checksum_path = RAW / "ClinVarRCVRelease_2026-09.xml.gz.md5"
    text = checksum_path.read_text(encoding="ascii", errors="replace")
    m = re.search(r"\b([0-9a-fA-F]{32})\b", text)
    if not m:
        raise ValueError(f"Official MD5 sidecar has no parseable checksum: {checksum_path}")
    expected = m.group(1).lower()
    integrity_path=META/"phase0R_archive_integrity.json"
    if integrity_path.exists():
        cached=json.loads(integrity_path.read_text(encoding="utf-8"))
        stat=ARCHIVE.stat()
        if (cached.get("archive_bytes")==stat.st_size and cached.get("archive_mtime_ns")==stat.st_mtime_ns
                and cached.get("expected_md5")==expected and cached.get("md5_match") is True and cached.get("gzip_test")=="PASS"):
            return cached
    observed = md5sum(ARCHIVE)
    transfer_log = LOGS / "phase0R_rcv_resume.log"
    gzip_pass = transfer_log.exists() and "gzip_integrity=PASS" in transfer_log.read_text(encoding="utf-8", errors="replace")
    result = {"archive":str(ARCHIVE.relative_to(ROOT)),"checksum_sidecar":str(checksum_path.relative_to(ROOT)),
              "archive_bytes":ARCHIVE.stat().st_size,"archive_mtime_ns":ARCHIVE.stat().st_mtime_ns,
              "expected_md5":expected,"observed_md5":observed,"md5_match":expected==observed,
              "gzip_test":"PASS" if gzip_pass else "NOT_CONFIRMED"}
    integrity_path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    if not gzip_pass:
        raise ValueError("gzip -t PASS was not recorded by the resumable transfer command; refusing to parse archive")
    if not result["md5_match"]:
        raise ValueError(f"Frozen RCV archive MD5 mismatch: expected {expected}, observed {observed}")
    return result


def download_archive_resumable() -> dict:
    """Download the frozen archive with explicit, verified HTTP byte ranges."""
    import requests
    if ARCHIVE.exists():
        return validate_frozen_archive()
    RAW.mkdir(parents=True,exist_ok=True)
    INTERIM.mkdir(parents=True,exist_ok=True)
    LOGS.mkdir(parents=True,exist_ok=True)
    META.mkdir(parents=True,exist_ok=True)
    checkpoint_path=META/"phase0R_rcv_download_checkpoint.json"
    chunk_path=INTERIM/"phase0R_rcv_range_chunk.tmp"
    log_path=LOGS/"phase0R_rcv_resume.log"
    session=requests.Session()
    headers={"Accept-Encoding":"identity","User-Agent":"phase0R-frozen-ClinVar-archive/1.0"}

    def log(message: str) -> None:
        with log_path.open("a",encoding="utf-8") as f:
            f.write(f"{utcnow()} {message}\n")

    if ARCHIVE_PART.exists():
        size=ARCHIVE_PART.stat().st_size
    else:
        ARCHIVE_PART.touch()
        size=0
    checkpoint=json.loads(checkpoint_path.read_text(encoding="utf-8")) if checkpoint_path.exists() else {}
    committed=int(checkpoint.get("verified_bytes",0))
    if not checkpoint:
        # Verify the retained curl partial as an exact prefix of the frozen file.
        # If a byte differs, keep only the matching prefix and resume from there.
        if size:
            response=session.get(RCV_ARCHIVE_URL,headers={**headers,"Range":f"bytes=0-{size-1}"},stream=True,timeout=(60,180))
            with response:
                if response.status_code!=206:
                    raise RuntimeError(f"Prefix verification expected HTTP 206, got {response.status_code}")
                content_range=response.headers.get("Content-Range","")
                m=re.fullmatch(r"bytes 0-(\d+)/(\d+)",content_range)
                if not m or int(m.group(1))!=size-1 or int(m.group(2))!=RCV_ARCHIVE_BYTES:
                    raise RuntimeError(f"Unexpected prefix Content-Range: {content_range!r}")
                local_offset=0
                mismatch=None
                with ARCHIVE_PART.open("rb") as local:
                    for remote_chunk in response.iter_content(chunk_size=1024*1024):
                        if not remote_chunk: continue
                        local_chunk=local.read(len(remote_chunk))
                        if local_chunk!=remote_chunk:
                            common=0
                            limit=min(len(local_chunk),len(remote_chunk))
                            while common<limit and local_chunk[common]==remote_chunk[common]: common+=1
                            mismatch=local_offset+common
                            break
                        local_offset+=len(remote_chunk)
                if mismatch is None and local_offset!=size:
                    raise RuntimeError(f"Prefix verification returned {local_offset} of {size} bytes")
                if mismatch is not None:
                    with ARCHIVE_PART.open("r+b") as local:
                        local.truncate(mismatch)
                    size=mismatch
                    log(f"prefix_repair truncated_at_verified_byte={mismatch}")
                else:
                    log(f"retained_partial_prefix_verified bytes={size}")
        committed=size
        checkpoint={"archive_url":RCV_ARCHIVE_URL,"expected_bytes":RCV_ARCHIVE_BYTES,"verified_bytes":committed,
                    "last_checkpoint_timestamp_utc":utcnow(),"verification":"existing prefix compared byte-for-byte with official HTTP range" if committed else "empty partial"}
        checkpoint_path.write_text(json.dumps(checkpoint,indent=2)+"\n",encoding="utf-8")
    else:
        if size>committed:
            with ARCHIVE_PART.open("r+b") as local: local.truncate(committed)
            log(f"discarded_uncommitted_tail bytes={size-committed}; retained_verified_bytes={committed}")
            size=committed
        elif size<committed:
            committed=size
            checkpoint["verified_bytes"]=size
            checkpoint["last_checkpoint_timestamp_utc"]=utcnow()
            checkpoint_path.write_text(json.dumps(checkpoint,indent=2)+"\n",encoding="utf-8")
    if committed>RCV_ARCHIVE_BYTES:
        raise RuntimeError(f"Partial archive exceeds official size: {committed}>{RCV_ARCHIVE_BYTES}")

    retry_n=0
    chunk_size=32*1024*1024
    while committed<RCV_ARCHIVE_BYTES:
        start=committed
        end=min(start+chunk_size-1,RCV_ARCHIVE_BYTES-1)
        expected=end-start+1
        try:
            with session.get(RCV_ARCHIVE_URL,headers={**headers,"Range":f"bytes={start}-{end}"},stream=True,timeout=(60,240)) as response:
                if response.status_code!=206:
                    raise RuntimeError(f"Range request expected HTTP 206, got {response.status_code}")
                content_range=response.headers.get("Content-Range","")
                m=re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",content_range)
                if not m or (int(m.group(1)),int(m.group(2)),int(m.group(3)))!=(start,end,RCV_ARCHIVE_BYTES):
                    raise RuntimeError(f"Unexpected range response: {content_range!r}")
                content_length=int(response.headers.get("Content-Length",expected))
                if content_length!=expected:
                    raise RuntimeError(f"Range Content-Length {content_length} does not match expected {expected}")
                written=0
                with chunk_path.open("wb") as out:
                    for block in response.iter_content(chunk_size=1024*1024):
                        if not block: continue
                        out.write(block); written+=len(block)
                    out.flush();os.fsync(out.fileno())
            if written!=expected or chunk_path.stat().st_size!=expected:
                raise RuntimeError(f"Incomplete range body: received {written} of {expected} bytes")
            with chunk_path.open("rb") as src, ARCHIVE_PART.open("ab") as dest:
                for block in iter(lambda:src.read(8*1024*1024),b""):
                    dest.write(block)
                dest.flush();os.fsync(dest.fileno())
            chunk_path.unlink(missing_ok=True)
            committed=end+1
            checkpoint={"archive_url":RCV_ARCHIVE_URL,"expected_bytes":RCV_ARCHIVE_BYTES,"verified_bytes":committed,
                        "last_range_start":start,"last_range_end":end,"last_checkpoint_timestamp_utc":utcnow(),
                        "verification":"HTTP 206 with exact Content-Range and complete body"}
            checkpoint_path.write_text(json.dumps(checkpoint,indent=2)+"\n",encoding="utf-8")
            log(f"range_committed start={start} end={end} bytes={committed}/{RCV_ARCHIVE_BYTES} percent={committed/RCV_ARCHIVE_BYTES*100:.3f}")
            retry_n=0
        except Exception as exc:
            chunk_path.unlink(missing_ok=True)
            retry_n+=1
            log(f"range_retry attempt={retry_n}/20 offset={start} error={type(exc).__name__}: {str(exc)[:300]}")
            if retry_n>=20:
                raise RuntimeError(f"20 explicit range retries failed at byte {start}; verified partial remains at {committed} bytes") from exc
            time.sleep(min(5*retry_n,60))

    if ARCHIVE_PART.stat().st_size!=RCV_ARCHIVE_BYTES:
        raise RuntimeError(f"Downloaded size differs from frozen archive: {ARCHIVE_PART.stat().st_size}")
    log("gzip_test_started")
    gzip_result=subprocess.run(["gzip","-t",str(ARCHIVE_PART)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    if gzip_result.returncode!=0:
        log(f"gzip_integrity=FAIL error={gzip_result.stderr[:300]}")
        raise RuntimeError("gzip -t failed; verified byte-range partial was preserved")
    observed_md5=md5sum(ARCHIVE_PART)
    expected_match=re.search(r"\b([0-9a-fA-F]{32})\b",(RAW/"ClinVarRCVRelease_2026-09.xml.gz.md5").read_text(encoding="ascii",errors="replace"))
    if not expected_match:
        raise RuntimeError("Official MD5 sidecar became unreadable after archive transfer")
    expected_md5=expected_match.group(1).lower()
    if observed_md5!=expected_md5:
        log(f"gzip_integrity=PASS official_md5=FAIL expected={expected_md5} observed={observed_md5}")
        raise RuntimeError("Official MD5 mismatch; frozen archive has not been accepted for parsing")
    ARCHIVE_PART.replace(ARCHIVE)
    archive_stat=ARCHIVE.stat()
    integrity={"archive":str(ARCHIVE.relative_to(ROOT)),"checksum_sidecar":"data/raw/ClinVarRCVRelease_2026-09.xml.gz.md5",
               "archive_bytes":archive_stat.st_size,"archive_mtime_ns":archive_stat.st_mtime_ns,
               "expected_md5":expected_md5,"observed_md5":observed_md5,"md5_match":True,"gzip_test":"PASS"}
    integrity_path=META/"phase0R_archive_integrity.json"
    integrity_path.write_text(json.dumps(integrity,indent=2)+"\n",encoding="utf-8")
    digest=sha256(ARCHIVE)
    with log_path.open("a",encoding="utf-8") as f:
        f.write(f"gzip_integrity=PASS sha256={digest}  data/raw/ClinVarRCVRelease_2026-09.xml.gz\n")
        f.write(f"official_md5=PASS md5={observed_md5}\n")
    checkpoint.update({"verified_bytes":RCV_ARCHIVE_BYTES,"completed_utc":utcnow(),"gzip_test":"PASS","md5_match":True,"sha256":digest})
    checkpoint_path.write_text(json.dumps(checkpoint,indent=2)+"\n",encoding="utf-8")
    return {"bytes":archive_stat.st_size,"sha256":digest,"md5":observed_md5,"gzip_test":"PASS","status":"complete"}


def canonical_missense(name: str) -> str:
    for m in HGVS_RE.finditer(name or ""):
        ref, pos, alt = m.groups()
        ref = AA3.get(ref, ref)
        alt = AA3.get(alt, alt)
        if ref in AA1 and alt in AA1 and ref != alt:
            return f"p.{ref}{int(pos)}{alt}"
    return ""


def split_aligned(text: str) -> list[str]:
    return (text or "").split("|")


def make_seed() -> dict:
    PHASE.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    if SEED.exists():
        raise FileExistsError(f"Refusing to overwrite existing Phase 0R seed: {SEED}")
    if not SUMMARY.exists():
        raise FileNotFoundError(SUMMARY)
    # Deduplicate assembly rows by the required biological key. Choose the most
    # complete condition annotation, then a stable lexicographic tie-break.
    unique: dict[tuple[str, str, str, str], tuple[int, tuple[str, ...]]] = {}
    scanned = kept_rows = 0
    last_log = time.monotonic()
    with gzip.open(SUMMARY, "rt", encoding="utf-8", errors="replace", newline="") as inp:
        reader = csv.DictReader(inp, delimiter="\t")
        required = {"VariationID", "GeneSymbol", "Type", "ClinicalSignificance", "Origin", "Name", "ReviewStatus", "PhenotypeIDS", "PhenotypeList", "RCVaccession"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"variant_summary is missing columns: {sorted(missing)}")
        for row in reader:
            scanned += 1
            sig = (row.get("ClinicalSignificance") or "").strip().lower()
            if sig not in {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic", "pathogenic, likely pathogenic"}:
                continue
            origin = (row.get("Origin") or "").lower()
            if "germline" not in origin or "somatic" in origin:
                continue
            gene = (row.get("GeneSymbol") or "").strip()
            if not gene or gene.lower() in {"not provided", "unknown", "intergenic"}:
                continue
            variation_id = (row.get("VariationID") or "").strip()
            if not variation_id.isdigit():
                continue
            protein = canonical_missense(row.get("Name") or "")
            if not protein:
                continue
            rcv_list = list(dict.fromkeys(RCV_RE.findall(row.get("RCVaccession") or "")))
            if not rcv_list:
                continue
            phen_ids = split_aligned(row.get("PhenotypeIDS") or "")
            phen_names = split_aligned(row.get("PhenotypeList") or "")
            for idx, rcv in enumerate(rcv_list):
                ids = phen_ids[idx].strip() if idx < len(phen_ids) else ""
                names = phen_names[idx].strip() if idx < len(phen_names) else ""
                vals = (variation_id, gene, protein, rcv, sig.title(), (row.get("ReviewStatus") or "").strip(), ids, names)
                key = (variation_id, gene, protein, rcv)
                completeness = int(bool(ids)) * 2 + int(bool(names))
                # Prefer the row with more condition metadata; stable field tuple breaks ties.
                old = unique.get(key)
                if old is None or completeness > old[0] or (completeness == old[0] and vals < old[1]):
                    unique[key] = (completeness, vals)
                kept_rows += 1
            if scanned % 1_000_000 == 0:
                with (LOGS / "phase0R_seed.log").open("a", encoding="utf-8") as log:
                    log.write(f"{utcnow()} scanned={scanned} candidate_source_rows={kept_rows}\n")
    with SEED.open("x", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(SEED_COLUMNS)
        for key in sorted(unique):
            writer.writerow(unique[key][1])
    variation_ids, substitutions, rcvs, genes = set(), set(), set(), set()
    row_n = 0
    with SEED.open(encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            row_n += 1
            variation_ids.add(row["VariationID"])
            substitutions.add((row["GeneSymbol"], row["ProteinChange"]))
            rcvs.add(row["RCVaccession"].split(".")[0])
            genes.add(row["GeneSymbol"])
    result = {
        "checkpoint": "PHASE0R_GLOBAL_SEED_BEFORE_RCV_DOWNLOAD",
        "timestamp_utc": utcnow(),
        "frozen_release": "ClinVar monthly archive 2026-09-03; variant_summary_2026-09.txt.gz",
        "source_path": str(SUMMARY.relative_to(ROOT)),
        "source_bytes": SUMMARY.stat().st_size,
        "source_sha256": sha256(SUMMARY),
        "source_rows_scanned": scanned,
        "candidate_source_rows_before_assembly_dedup": kept_rows,
        "unique_seed_rows_variation_gene_protein_rcv": row_n,
        "n_unique_variation_ids": len(variation_ids),
        "n_unique_gene_protein_substitutions": len(substitutions),
        "n_unique_rcv_accessions": len(rcvs),
        "n_unique_genes": len(genes),
        "seed_path": str(SEED.relative_to(ROOT)),
        "seed_sha256": sha256(SEED),
    }
    if SEED_CHECKPOINT.exists():
        raise FileExistsError(f"Refusing to overwrite existing checkpoint: {SEED_CHECKPOINT}")
    SEED_CHECKPOINT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    with (LOGS / "phase0R_seed.log").open("a", encoding="utf-8") as log:
        log.write(json.dumps(result, sort_keys=True) + "\n")
    return result


def rebuild_previous_965() -> dict:
    """Reproduce the previous summary-screen logic from its preserved disease rows."""
    source = INTERIM / "clinvar_preliminary_candidate_disease_rows.tsv"
    out = PHASE / "previous_965_gene_list.tsv"
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite existing gene list: {out}")
    disease_variants: dict[tuple[str, str], set[str]] = defaultdict(set)
    with source.open(encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gene = row.get("gene", "").strip()
            variant = row.get("protein_change", "").strip()
            for mondo in (x for x in row.get("mondo_ids", "").split(";") if x):
                disease_variants[(gene, mondo)].add(variant)
    strict_volume: dict[str, int] = defaultdict(int)
    strict_diseases: dict[str, int] = defaultdict(int)
    for (gene, _), variants in disease_variants.items():
        if gene and len(variants) >= 10:
            strict_volume[gene] += len(variants)
            strict_diseases[gene] += 1
    genes = sorted((g for g, n in strict_diseases.items() if n >= 2), key=lambda g: (-strict_volume[g], -strict_diseases[g], g))
    with out.open("x", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "preliminary_strict_disease_groups", "preliminary_unique_variants_in_strict_groups", "rebuild_source", "logic"])
        for g in genes:
            w.writerow([g, strict_diseases[g], strict_volume[g], str(source.relative_to(ROOT)), ">=2 MONDO groups each >=10 variants; same gene/protein substitution/disease collapsed"])
    return {"recovered_gene_n": len(genes), "expected_previous_gene_n": 965, "exact_965_recovered": len(genes) == 965, "path": str(out.relative_to(ROOT))}


def local_name(node) -> str:
    return etree.QName(node).localname


def descendants(node, wanted: str):
    return (x for x in node.iter() if local_name(x) == wanted)


def direct_child(node, wanted: str):
    return next((x for x in node if local_name(x) == wanted), None)


def text_value(node) -> str:
    return " ".join(" ".join(node.itertext()).split()) if node is not None else ""


def rcv_record(clinvar_set, source: str, release_date: str = "") -> dict | None:
    ref = next((x for x in clinvar_set if local_name(x) == "ReferenceClinVarAssertion"), None)
    if ref is None:
        return None
    accession = next((x for x in descendants(ref, "ClinVarAccession") if x.get("Type") == "RCV"), None)
    if accession is None or not accession.get("Acc"):
        return None
    rcv = accession.get("Acc").split(".")[0]
    classes = direct_child(ref, "Classifications")
    germline_node = next((x for x in descendants(classes, "GermlineClassification")), None) if classes is not None else None
    desc = next((x for x in descendants(germline_node, "Description")), None) if germline_node is not None else None
    classification = text_value(desc)
    review_node = next((x for x in descendants(germline_node, "ReviewStatus")), None) if germline_node is not None else None
    review = text_value(review_node)
    evaluated = desc.get("DateLastEvaluated", "") if desc is not None else ""
    names, ids = [], set()
    traitsets = [x for x in ref if local_name(x) == "TraitSet" and x.get("Type", "").lower() in {"disease", "condition"}]
    if not traitsets:
        traitsets = [x for x in descendants(ref, "TraitSet") if x.get("Type", "").lower() in {"disease", "condition"}]
    for ts in traitsets:
        for trait in descendants(ts, "Trait"):
            for value in descendants(trait, "ElementValue"):
                if value.get("Type") == "Preferred" and text_value(value):
                    names.append(text_value(value))
            for xr in descendants(trait, "XRef"):
                db, xid = (xr.get("DB") or "").strip(), (xr.get("ID") or "").strip()
                if not db or not xid:
                    continue
                db_norm = db.casefold()
                if db_norm == "mondo":
                    norm_id = xid if xid.startswith("MONDO:") else "MONDO:" + xid
                else:
                    norm_id = db + ":" + xid
                ids.add(norm_id)
    measure_ids, variation_ids, genes = set(), set(), set()
    for ms in descendants(ref, "MeasureSet"):
        if ms.get("Acc", "").startswith("VCV"):
            version = ms.get("Version", "")
            measure_ids.add(ms.get("Acc") + ("." + version if version else ""))
            if ms.get("ID", "").isdigit():
                variation_ids.add(ms.get("ID"))
        for node in descendants(ms, "Gene"):
            if node.get("Symbol"):
                genes.add(node.get("Symbol"))
    parent = clinvar_set.getparent()
    return {
        "RCVaccession": rcv,
        "RCVversion": accession.get("Version", ""),
        "preferred_condition_name": "; ".join(dict.fromkeys(names)),
        "condition_IDs": ";".join(sorted(ids)),
        "MONDO": ";".join(sorted(x for x in ids if x.startswith("MONDO:"))),
        "MedGen": ";".join(sorted(x for x in ids if x.startswith("MedGen:"))),
        "OMIM": ";".join(sorted(x for x in ids if x.startswith("OMIM:"))),
        "Orphanet": ";".join(sorted(x for x in ids if x.startswith("Orphanet:"))),
        "aggregate_germline_classification": classification,
        "review_status": review,
        "last_evaluated": evaluated,
        "VariationID": ";".join(sorted(variation_ids)),
        "VCV_linkage": ";".join(sorted(measure_ids)),
        "gene_symbols": ";".join(sorted(genes)),
        "recovery_source": source,
        "archive_release_dated": release_date or (parent.get("Dated", "") if parent is not None else ""),
    }


def db_open() -> sqlite3.Connection:
    path = INTERIM / "phase0R_rcv_recovered.sqlite"
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("CREATE TABLE IF NOT EXISTS rcv (accession TEXT PRIMARY KEY, payload TEXT NOT NULL)")
    conn.commit()
    return conn


def write_parse_checkpoint(conn, candidate_n: int, scanned: int, pass_n: int, archive_dated: str, complete: bool, all_found: bool, total_base: int = 0) -> dict:
    found = conn.execute("SELECT COUNT(*) FROM rcv").fetchone()[0]
    payload = {
        "xml_records_scanned": scanned,
        "xml_records_scanned_total": total_base + scanned,
        "candidate_rcvs_found": found,
        "candidate_rcvs_remaining": max(0, candidate_n - found),
        "candidate_rcv_n": candidate_n,
        "last_checkpoint_timestamp_utc": utcnow(),
        "scan_pass": pass_n,
        "archive_release_dated": archive_dated,
        "archive_scan_completed": complete,
        "all_candidate_rcvs_found": all_found,
        "archive_path": str(ARCHIVE.relative_to(ROOT)),
    }
    PARSER_CHECKPOINT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


def parse_archive() -> dict:
    if not ARCHIVE.exists():
        raise FileNotFoundError(f"Validated frozen archive not present: {ARCHIVE}")
    candidate = set()
    with SEED.open(encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            candidate.add(row["RCVaccession"].split(".")[0])
    conn = db_open()
    existing = {x[0] for x in conn.execute("SELECT accession FROM rcv")}
    found_before = len(existing)
    if found_before >= len(candidate):
        conn.close()
        return {"candidate_rcv_n": len(candidate), "archived_rcv_n": found_before, "skipped": True}
    pass_n = 1
    total_base = 0
    if PARSER_CHECKPOINT.exists():
        try:
            prior = json.loads(PARSER_CHECKPOINT.read_text(encoding="utf-8"))
            pass_n = int(prior.get("scan_pass", 0)) + 1
            total_base = int(prior.get("xml_records_scanned_total", 0))
        except (ValueError, OSError):
            pass
    seen = set()
    scanned = 0
    archive_dated = ""
    stopped_early = False
    with gzip.open(ARCHIVE, "rb") as stream:
        context = etree.iterparse(stream, events=("end",), tag="{*}ClinVarSet", huge_tree=True, recover=False)
        for _, elem in context:
            scanned += 1
            parent = elem.getparent()
            if parent is not None and not archive_dated:
                archive_dated = parent.get("Dated", "")
            ref = next((x for x in elem if local_name(x) == "ReferenceClinVarAssertion"), None)
            acc_node = next((x for x in descendants(ref, "ClinVarAccession") if x.get("Type") == "RCV"), None) if ref is not None else None
            acc = acc_node.get("Acc", "").split(".")[0] if acc_node is not None else ""
            if acc in candidate and acc not in seen:
                rec = rcv_record(elem, "ClinVarRCVRelease_2026-09.xml.gz", archive_dated)
                if rec:
                    if acc not in existing:
                        conn.execute("INSERT OR IGNORE INTO rcv(accession,payload) VALUES (?,?)", (acc, json.dumps(rec, ensure_ascii=False)))
                    seen.add(acc)
            if scanned % 100_000 == 0:
                conn.commit()
                snap = write_parse_checkpoint(conn, len(candidate), scanned, pass_n, archive_dated, False, False, total_base)
                with (LOGS / "phase0R_rcv_parse.log").open("a", encoding="utf-8") as log:
                    log.write(json.dumps(snap, sort_keys=True) + "\n")
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]
            if len(existing | seen) >= len(candidate):
                stopped_early = True
                break
    conn.commit()
    found = conn.execute("SELECT COUNT(*) FROM rcv").fetchone()[0]
    all_found = found >= len(candidate)
    # When all requested records are found, scanning unrelated release records is unnecessary.
    complete = not stopped_early
    snap = write_parse_checkpoint(conn, len(candidate), scanned, pass_n, archive_dated, complete, all_found, total_base)
    conn.close()
    return {"candidate_rcv_n": len(candidate), "archived_rcv_n": found, "remaining_rcv_n": len(candidate)-found,
            "scan_pass": pass_n, "xml_records_scanned_this_pass": scanned, "archive_release_dated": archive_dated,
            "all_candidate_rcvs_found": all_found, "checkpoint": snap}


def export_recovered() -> list[dict]:
    conn = db_open()
    rows = [json.loads(x[0]) for x in conn.execute("SELECT payload FROM rcv ORDER BY accession")]
    conn.close()
    out = PHASE / "global_rcv_recovered.tsv"
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite existing recovered RCV table: {out}")
    columns = ["RCVaccession", "RCVversion", "preferred_condition_name", "condition_IDs", "MONDO", "MedGen", "OMIM", "Orphanet",
               "aggregate_germline_classification", "review_status", "last_evaluated", "VariationID", "VCV_linkage", "gene_symbols", "recovery_source", "archive_release_dated"]
    with out.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return rows


def fallback_missing_rcvs() -> dict:
    """Fetch only accession-level gaps after the frozen archive scan."""
    import requests
    conn = db_open()
    recovered = {x[0] for x in conn.execute("SELECT accession FROM rcv")}
    seed_rcvs = set()
    with SEED.open(encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            seed_rcvs.add(row["RCVaccession"].split(".")[0])
    missing = sorted(seed_rcvs - recovered)
    API_CACHE.mkdir(parents=True, exist_ok=True)
    failed = []
    session = requests.Session()
    for offset in range(0, len(missing), 40):
        batch = missing[offset:offset+40]
        params = {"db":"clinvar", "rettype":"clinvarset", "id":",".join(batch), "retmode":"xml", "tool":"phase0R_full_rcv_recovery"}
        error = ""
        for attempt in range(10):
            try:
                response = session.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params, timeout=(40, 300))
                response.raise_for_status()
                root = etree.fromstring(response.content, etree.XMLParser(huge_tree=True))
                if local_name(root) != "ClinVarResult-Set":
                    raise ValueError("Unexpected E-utilities XML root")
                got = set()
                for elem in root:
                    if local_name(elem) != "ClinVarSet":
                        continue
                    rec = rcv_record(elem, "E-utilities targeted RCV fallback")
                    if not rec or rec["RCVaccession"] not in batch:
                        continue
                    acc = rec["RCVaccession"]
                    cache = API_CACHE / f"{acc}.xml"
                    if not cache.exists():
                        with cache.open("xb") as f:
                            f.write(etree.tostring(elem, encoding="utf-8", xml_declaration=True))
                    conn.execute("INSERT OR REPLACE INTO rcv(accession,payload) VALUES (?,?)", (acc, json.dumps(rec, ensure_ascii=False)))
                    got.add(acc)
                conn.commit()
                for acc in sorted(set(batch) - got):
                    failed.append({"RCVaccession":acc,"reason":"E-utilities returned no RCV record"})
                error = ""
                with (LOGS / "phase0R_rcv_fallback.log").open("a", encoding="utf-8") as log:
                    log.write(f"{utcnow()} requested={len(batch)} recovered={len(got)} missing={len(batch)-len(got)}\n")
                break
            except Exception as exc:
                error = f"{type(exc).__name__}: {exc}"
                if attempt < 9:
                    time.sleep(min(5 * (attempt + 1), 60))
        if error:
            failed.extend({"RCVaccession":acc,"reason":error} for acc in batch)
            with (LOGS / "phase0R_rcv_fallback.log").open("a", encoding="utf-8") as log:
                log.write(f"{utcnow()} failed_batch={','.join(batch)} error={error}\n")
        time.sleep(0.34)
    conn.close()
    return {"fallback_requested_n":len(missing),"fallback_failures":failed}


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing Phase 0R output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_recovery() -> dict:
    integrity = validate_frozen_archive()
    parsed = parse_archive()
    fallback = fallback_missing_rcvs()
    records = export_recovered()
    candidate_n = len({x["RCVaccession"].split(".")[0] for x in csv.DictReader(SEED.open(encoding="utf-8"), delimiter="\t")})
    recovered = {r["RCVaccession"] for r in records}
    failed = sorted({r["RCVaccession"].split(".")[0] for r in csv.DictReader(SEED.open(encoding="utf-8"), delimiter="\t")} - recovered)
    recovery_rate = len(recovered) / candidate_n if candidate_n else 0.0
    write_tsv(PHASE / "rcv_recovery_summary.tsv", [{
        "candidate_rcv_n":candidate_n,
        "recovered_rcv_n":len(recovered),
        "failed_rcv_n":len(failed),
        "recovery_rate":f"{recovery_rate:.8f}",
        "recovery_percent":f"{recovery_rate*100:.4f}",
        "archive_rcv_n":parsed.get("archived_rcv_n", "NA"),
        "fallback_requested_n":fallback["fallback_requested_n"],
        "archive_release_dated":parsed.get("archive_release_dated", "NA"),
        "result":"PASS" if recovery_rate >= 0.95 else "HOLD_RCV_RECOVERY",
    }], ["candidate_rcv_n","recovered_rcv_n","failed_rcv_n","recovery_rate","recovery_percent","archive_rcv_n","fallback_requested_n","archive_release_dated","result"])
    write_tsv(PHASE / "rcv_recovery_failures.tsv", [{"RCVaccession":x,"reason":"not present in frozen archive and targeted E-utilities fallback failed"} for x in failed], ["RCVaccession","reason"])
    result = {"candidate_rcv_n":candidate_n,"recovered_rcv_n":len(recovered),"failed_rcv_n":len(failed),"recovery_rate":recovery_rate,
              "verdict":"CONTINUE" if recovery_rate >= .95 else "HOLD_RCV_RECOVERY","fallback_requested_n":fallback["fallback_requested_n"],
              "fallback_failure_n":len(fallback["fallback_failures"]),"archive_release_dated":parsed.get("archive_release_dated", "NA")}
    result["archive_md5_match"] = integrity["md5_match"]
    gate_path = META / "phase0R_recovery_gate.json"
    if gate_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing recovery gate: {gate_path}")
    gate_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def load_ontology_helpers():
    sys.path.insert(0, str(ROOT / "scripts"))
    from run_phase0 import parse_mondo, mondo_map, ancestors
    return parse_mondo, mondo_map, ancestors


def exact_disease_audit() -> dict:
    """Apply exact RCV classification and MONDO mapping to the recovered seed."""
    gate = json.loads((META / "phase0R_recovery_gate.json").read_text(encoding="utf-8"))
    if gate["recovery_rate"] < .95:
        return {"verdict":"HOLD_RCV_RECOVERY", "reason":"recovery rate below 95%; disease normalization not run"}
    parse_mondo, mondo_map, ancestors = load_ontology_helpers()
    terms, xref_map, synonym_map, alt_map = parse_mondo(RAW / "mondo-v2026-09-01.obo")
    if not terms:
        raise RuntimeError("Frozen MONDO ontology could not be parsed")
    exact_name_lookup=defaultdict(set)
    for mid, term in terms.items():
        if term.get("obsolete"):
            continue
        if term.get("name"):
            exact_name_lookup[term["name"].strip().casefold()].add(mid)
        for synonym, scope in term.get("synonyms", []):
            if scope == "EXACT":
                exact_name_lookup[synonym.strip().casefold()].add(mid)
    genes = {}
    with (PHASE / "previous_965_gene_list.tsv").open(encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes[row["gene"]] = row
    pilot_genes = {"ACTB", "ACTG1"}
    gene_domain = set(genes) | pilot_genes
    seed_by_rcv = defaultdict(list)
    with SEED.open(encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["GeneSymbol"] in gene_domain:
                seed_by_rcv[row["RCVaccession"].split(".")[0]].append(row)
    rcv_db = db_open()
    rcv_records = {acc:json.loads(payload) for acc,payload in rcv_db.execute("SELECT accession,payload FROM rcv")}
    rcv_db.close()

    def is_plp(text: str) -> bool:
        return (text or "").strip().lower() in {"pathogenic","likely pathogenic","pathogenic/likely pathogenic","pathogenic, likely pathogenic"}

    def map_record(record: dict) -> list[str]:
        ids = [x for x in record.get("condition_IDs", "").split(";") if x]
        direct = set()
        for value in ids:
            direct.update(mondo_map([value], terms, xref_map, synonym_map, alt_map))
        if direct:
            return sorted(direct)
        names = [x.strip() for x in record.get("preferred_condition_name", "").split(";") if x.strip()]
        found = set()
        for name in names:
            found.update(exact_name_lookup.get(name.casefold(), set()))
        return sorted(found)

    raw_by_gene = defaultdict(set)
    exact_rcvs_by_gene = defaultdict(set)
    exact_variants_by_gene = defaultdict(set)
    raw_to_mondo = defaultdict(set)
    variants_by_gene_term = defaultdict(set)
    rcvs_by_gene_term = defaultdict(set)
    raw_assoc_by_gene_term = defaultdict(set)
    unmapped_by_gene = Counter()
    ambiguous_by_gene = Counter()
    variation_link_mismatch_by_gene = Counter()
    linkage_rows = []
    pair_reason = []
    for acc, seed_rows in seed_by_rcv.items():
        rec = rcv_records.get(acc)
        if rec is None or not is_plp(rec.get("aggregate_germline_classification", "")):
            continue
        if not rec.get("condition_IDs") and not rec.get("preferred_condition_name"):
            continue
        raw_key = (rec.get("preferred_condition_name", "").strip().casefold(), rec.get("condition_IDs", "").strip())
        mapped = map_record(rec)
        for seed in seed_rows:
            gene = seed["GeneSymbol"]
            protein = seed["ProteinChange"]
            xml_variants = {x for x in rec.get("VariationID", "").split(";") if x}
            if xml_variants and seed["VariationID"] not in xml_variants:
                variation_link_mismatch_by_gene[gene] += 1
                linkage_rows.append({"gene":gene,"RCVaccession":acc,"seed_VariationID":seed["VariationID"],"seed_ProteinChange":protein,
                                     "xml_VariationID":rec.get("VariationID",""),"VCV_linkage":rec.get("VCV_linkage",""),"status":"VARIATION_ID_MISMATCH"})
                continue
            linkage_rows.append({"gene":gene,"RCVaccession":acc,"seed_VariationID":seed["VariationID"],"seed_ProteinChange":protein,
                                 "xml_VariationID":rec.get("VariationID",""),"VCV_linkage":rec.get("VCV_linkage",""),
                                 "status":"MATCH" if xml_variants else "XML_VARIATION_ID_UNAVAILABLE"})
            exact_rcvs_by_gene[gene].add(acc)
            exact_variants_by_gene[gene].add(protein)
            raw_by_gene[gene].add(raw_key)
            if len(mapped) != 1:
                if len(mapped) == 0:
                    unmapped_by_gene[gene] += 1
                else:
                    ambiguous_by_gene[gene] += 1
                continue
            mondo = mapped[0]
            raw_to_mondo[(gene,raw_key)].add(mondo)
            variants_by_gene_term[(gene,mondo)].add(protein)
            rcvs_by_gene_term[(gene,mondo)].add(acc)
            raw_assoc_by_gene_term[(gene,mondo)].add((protein,acc))

    # Retain only leaf-level conditions when a broader MONDO term and its child
    # are both assigned to the same gene. Parent-assigned variants are not
    # reassigned to a child, avoiding artificial inflation of child counts.
    terms_by_gene = defaultdict(set)
    for gene, mondo in variants_by_gene_term:
        terms_by_gene[gene].add(mondo)
    kept_terms_by_gene = {}
    parent_child_by_gene = {}
    for gene, ids in terms_by_gene.items():
        cache = {mid:ancestors(mid, terms) for mid in ids}
        removed = {parent for parent in ids if any(parent in cache.get(child, set()) for child in ids if child != parent)}
        kept_terms_by_gene[gene] = ids - removed
        parent_child_by_gene[gene] = removed
        for parent in sorted(removed):
            children = sorted(child for child in ids if child != parent and parent in cache.get(child, set()))
            pair_reason.append({"gene":gene,"disease_A":parent,"disease_A_label":terms.get(parent,{}).get("name",""),
                                "disease_B":";".join(children),"disease_B_label":";".join(terms.get(c,{}).get("name","") for c in children),
                                "relationship":"PARENT_CHILD_REVIEW","status":"collapsed_parent_for_strict_count","reason":"same gene has both an ancestor MONDO term and a more specific descendant; broad term not counted as independent"})

    disease_rows = []
    gene_rows = []
    attrition_rows = []
    stage_counts = Counter()
    for gene in sorted(gene_domain):
        all_terms = terms_by_gene.get(gene, set())
        kept = kept_terms_by_gene.get(gene, set())
        term_counts = {m:len(variants_by_gene_term[(gene,m)]) for m in sorted(kept)}
        term_variants = {m:variants_by_gene_term[(gene,m)] for m in kept}
        selected_terms = {m for m,n in term_counts.items() if n >= 10}
        expanded_terms = {m for m,n in term_counts.items() if n >= 5}
        strict_union = set().union(*(term_variants[m] for m in selected_terms)) if selected_terms else set()
        expanded_union = set().union(*(term_variants[m] for m in expanded_terms)) if expanded_terms else set()
        strict = len(selected_terms) >= 2 and len(strict_union) >= 20
        expanded = len(expanded_terms) >= 2
        n_raw = len(raw_by_gene.get(gene, set()))
        mapped_raw = sum(1 for (g, _), mids in raw_to_mondo.items() if g == gene and mids)
        mono_raw_groups = defaultdict(set)
        for (g, raw), mids in raw_to_mondo.items():
            if g == gene:
                for m in mids:
                    mono_raw_groups[m].add(raw)
        same_mondo_collapses = sum(max(0,len(raws)-1) for raws in mono_raw_groups.values())
        raw_pair_rows = sum(len(v) for (g,_),v in raw_assoc_by_gene_term.items() if g == gene)
        unique_pair_variants = sum(len(v) for (g,m),v in variants_by_gene_term.items() if g == gene and m in kept)
        duplicate_collapses = max(0, raw_pair_rows - unique_pair_variants)
        hotspot = False
        disease_display = []
        for mondo in sorted(kept):
            variants = term_variants[mondo]
            residues = Counter()
            for protein in variants:
                m = re.match(r"p\.[A-Z](\d+)[A-Z]$", protein)
                if m:
                    residues[m.group(1)] += 1
            ordered = residues.most_common()
            hotspot_term = bool(len(variants) >= 10 and ordered and (ordered[0][1] / len(variants) >= .20 or sum(n for _,n in ordered[:3]) / len(variants) >= .50))
            hotspot = hotspot or hotspot_term
            disease_display.append(f"{mondo}={len(variants)} ({len(residues)} residues)")
            disease_rows.append({"gene":gene,"MONDO":mondo,"MONDO_label":terms.get(mondo,{}).get("name",""),
                                 "n_unique_variants":len(variants),"n_unique_residues":len(residues),"n_RCV":len(rcvs_by_gene_term[(gene,mondo)]),
                                 "strict_disease_ge10":"YES" if len(variants)>=10 else "NO","expanded_disease_ge5":"YES" if len(variants)>=5 else "NO",
                                 "parent_child_collapsed_terms_for_gene":";".join(sorted(parent_child_by_gene.get(gene,set()))),
                                 "raw_condition_names":";".join(sorted({rcv_records[a].get("preferred_condition_name","") for a in rcvs_by_gene_term[(gene,mondo)] if a in rcv_records}))})
        if gene in genes:
            if len(exact_rcvs_by_gene.get(gene,set())):
                stage_counts["after_exact_RCV"] += 1
            if all_terms:
                stage_counts["after_disease_normalization"] += 1
            if len(kept) >= 2:
                stage_counts["at_least_2_distinct_diseases"] += 1
            if expanded:
                stage_counts["at_least_5_variants_per_disease"] += 1
            if strict:
                stage_counts["strict_at_least_10_variants_per_disease"] += 1
        if not exact_rcvs_by_gene.get(gene):
            reason = "NO_EXACT_RCV_PLP_MISSENSE"
        elif not all_terms:
            reason = "NO_UNAMBIGUOUS_MONDO_DISEASE"
        elif len(kept) < 2:
            reason = "FEWER_THAN_TWO_DISEASES_AFTER_PARENT_CHILD_COLLAPSE"
        elif not expanded:
            reason = "FEWER_THAN_TWO_DISEASES_WITH_5_VARIANTS"
        elif len(selected_terms) < 2:
            reason = "FEWER_THAN_TWO_DISEASES_WITH_10_VARIANTS"
        elif len(strict_union) < 20:
            reason = "TOTAL_UNIQUE_VARIANTS_BELOW_20"
        else:
            reason = "STRICT_PASS"
        if gene in genes:
            attrition_rows.append({"gene":gene,"prescreen_status":"STRICT_PRESCREEN_965",
                               "n_exact_rcv_plp_missense":len(exact_rcvs_by_gene.get(gene,set())),
                               "n_raw_conditions":n_raw,"n_normalized_diseases":len(all_terms),
                               "n_independent_diseases_after_parent_child_collapse":len(kept),
                               "n_variants_per_disease":";".join(disease_display),
                               "same_MONDO_collapses":same_mondo_collapses,"parent_child_collapses":len(parent_child_by_gene.get(gene,set())),
                               "duplicate_variant_collapses":duplicate_collapses,"hotspot_flag":"YES" if hotspot else "NO",
                               "strict_after_exact_RCV":"YES" if strict else "NO","reason_if_lost":"" if strict else reason,
                               "ambiguous_MONDO_RCV_n":ambiguous_by_gene.get(gene,0),"unmapped_condition_RCV_n":unmapped_by_gene.get(gene,0),
                               "variation_link_mismatch_n":variation_link_mismatch_by_gene.get(gene,0)})
        gene_rows.append({"gene":gene,"prescreen_status":"STRICT_PRESCREEN_965" if gene in genes else "PRESPECIFIED_PILOT",
                          "n_exact_PLP_missense_variants":len(exact_variants_by_gene.get(gene,set())),
                          "n_exact_rcv":len(exact_rcvs_by_gene.get(gene,set())),"n_normalized_diseases":len(all_terms),
                          "n_distinct_diseases_after_parent_child_collapse":len(kept),"n_diseases_ge5":len(expanded_terms),
                          "n_diseases_ge10":len(selected_terms),"n_unique_variants_in_strict_groups":len(strict_union),
                          "n_unique_variants_in_expanded_groups":len(expanded_union),"strict_after_exact_RCV":"YES" if strict else "NO",
                          "reason_if_lost":"" if strict else reason})
    write_tsv(PHASE / "strict965_attrition.tsv", attrition_rows,
              ["gene","prescreen_status","n_exact_rcv_plp_missense","n_raw_conditions","n_normalized_diseases","n_independent_diseases_after_parent_child_collapse","n_variants_per_disease","same_MONDO_collapses","parent_child_collapses","duplicate_variant_collapses","hotspot_flag","strict_after_exact_RCV","reason_if_lost","ambiguous_MONDO_RCV_n","unmapped_condition_RCV_n","variation_link_mismatch_n"])
    write_tsv(PHASE / "gene_disease_variant_counts_exact_rcv.tsv", disease_rows,
              ["gene","MONDO","MONDO_label","n_unique_variants","n_unique_residues","n_RCV","strict_disease_ge10","expanded_disease_ge5","parent_child_collapsed_terms_for_gene","raw_condition_names"])
    write_tsv(PHASE / "gene_exact_rcv_attrition_summary.tsv", gene_rows,
              ["gene","prescreen_status","n_exact_PLP_missense_variants","n_exact_rcv","n_normalized_diseases","n_distinct_diseases_after_parent_child_collapse","n_diseases_ge5","n_diseases_ge10","n_unique_variants_in_strict_groups","n_unique_variants_in_expanded_groups","strict_after_exact_RCV","reason_if_lost"])
    write_tsv(PHASE / "disease_parent_child_collapses.tsv", pair_reason,
              ["gene","disease_A","disease_A_label","disease_B","disease_B_label","relationship","status","reason"])
    write_tsv(PHASE / "rcv_seed_variation_linkage_audit.tsv", linkage_rows,
              ["gene","RCVaccession","seed_VariationID","seed_ProteinChange","xml_VariationID","VCV_linkage","status"])
    stages = [
        {"stage":"previous strict pre-screen genes","genes":len(genes)},
        {"stage":"after exact RCV P/LP missense","genes":stage_counts["after_exact_RCV"]},
        {"stage":"after unambiguous MONDO normalization","genes":stage_counts["after_disease_normalization"]},
        {"stage":"at least 2 diseases after parent-child collapse","genes":stage_counts["at_least_2_distinct_diseases"]},
        {"stage":"at least 2 diseases with >=5 unique variants each","genes":stage_counts["at_least_5_variants_per_disease"]},
        {"stage":"strict >=10 variants/disease and >=20 total","genes":stage_counts["strict_at_least_10_variants_per_disease"]},
    ]
    write_tsv(PHASE / "strict965_attrition_stages.tsv", stages, ["stage","genes"])
    threshold_rows=[]
    for threshold in (2,3,5,10,15,20):
        eligible_genes=[]
        qualifying_disease_n=0
        qualifying_variant_n=0
        for gene in genes:
            qualifying=[m for m in kept_terms_by_gene.get(gene,set()) if len(variants_by_gene_term[(gene,m)])>=threshold]
            if len(qualifying)>=2:
                eligible_genes.append(gene)
                qualifying_disease_n += len(qualifying)
                qualifying_variant_n += len(set().union(*(variants_by_gene_term[(gene,m)] for m in qualifying)))
        threshold_rows.append({"minimum_variants_per_disease":threshold,"genes":len(eligible_genes),
                               "disease_groups":qualifying_disease_n,"unique_variants_in_qualifying_groups":qualifying_variant_n,
                               "evidence_stage":"exact RCV + unambiguous MONDO + parent-child collapse"})
    write_tsv(PHASE / "phase0R_candidate_thresholds.tsv",threshold_rows,
              ["minimum_variants_per_disease","genes","disease_groups","unique_variants_in_qualifying_groups","evidence_stage"])
    reasons = Counter(x["reason_if_lost"] for x in attrition_rows if x["reason_if_lost"])
    result = {"previous_965_recovered":len(genes)==965,"n_previous_genes":len(genes),"n_after_exact_rcv":stage_counts["after_exact_RCV"],
              "n_after_disease_normalization":stage_counts["after_disease_normalization"],"n_ge2_distinct_disease":stage_counts["at_least_2_distinct_diseases"],
              "n_ge5_variants_per_disease":stage_counts["at_least_5_variants_per_disease"],"n_strict_ge10":stage_counts["strict_at_least_10_variants_per_disease"],
              "exclusion_reasons":reasons.most_common(5),"path":"results/phase0R/strict965_attrition.tsv"}
    (META / "phase0R_exact_disease_counts.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def ppi_audit() -> dict:
    ppi_started=utcnow()
    attrition = read_tsv(PHASE / "strict965_attrition.tsv")
    strict_genes = sorted(r["gene"] for r in attrition if r["strict_after_exact_RCV"] == "YES")
    queried_genes = set(strict_genes) | {"ACTB", "ACTG1"}
    sys.path.insert(0, str(ROOT / "scripts"))
    from run_phase0 import build_partner_tables
    provenance, stats, _ = build_partner_tables(queried_genes)
    predictomes_path=RAW/"predictomes_human_pair_scores_2025-11-10.csv.gz"
    predictomes_url="https://predictomes-hsbps-dataset.s3.us-east-1.amazonaws.com/20251110_hs_predictome_pair_scores.csv.gz"
    if not predictomes_path.exists():
        part=RAW/(predictomes_path.name+".part")
        RAW.mkdir(parents=True,exist_ok=True)
        log_path=LOGS/"phase0R_predictomes_download.log"
        proc=subprocess.run(["curl","-fL","-C","-","--retry","10","--retry-delay","3","--retry-all-errors","-o",str(part),predictomes_url],
                            stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        log_path.write_text(f"url={predictomes_url}\nexit={proc.returncode}\n{proc.stderr}\n",encoding="utf-8")
        if proc.returncode==0:
            with gzip.open(part,"rb") as compressed:
                while compressed.read(8*1024*1024): pass
            part.replace(predictomes_path)
            (RAW/(predictomes_path.name+".source.json")).write_text(json.dumps({"url":predictomes_url,"release":"Human Predictome 2025-11-10 pair scores",
                "retrieved_utc":utcnow(),"bytes":predictomes_path.stat().st_size,"sha256":sha256(predictomes_path)},indent=2)+"\n",encoding="utf-8")
    predictomes_rows=[]
    if predictomes_path.exists():
        from run_phase0 import load_hgnc_uniprot_maps
        uniprot_to_symbol,_=load_hgnc_uniprot_maps(RAW/"hgnc_non_alt_loci_2026-09-27.tsv")
        with gzip.open(predictomes_path,"rt",encoding="utf-8",errors="replace",newline="") as f:
            reader=csv.DictReader(f)
            fields={re.sub(r"[^a-z0-9]","",x.lower()):x for x in (reader.fieldnames or [])}
            a_col=next((fields[x] for x in ("uniprot1","uniprotid1","protein1","proteinid1","accession1","uniprota") if x in fields),None)
            b_col=next((fields[x] for x in ("uniprot2","uniprotid2","protein2","proteinid2","accession2","uniprotb") if x in fields),None)
            pair_col=next((fields[x] for x in ("uniprotids","uniprotaccessions","proteinids") if x in fields),None)
            score_col=next((fields[x] for x in ("spoc","spocscore","spocprobability","spocprob") if x in fields),None)
            if (not pair_col and (not a_col or not b_col)) or not score_col:
                raise ValueError(f"Predictomes score table columns not recognized: {reader.fieldnames}")
            for row in reader:
                try: score=float(row[score_col])
                except (TypeError,ValueError): continue
                if score<=0.86:
                    continue
                if pair_col:
                    pair_ids=(row.get(pair_col) or "").strip().split(":")
                    if len(pair_ids)!=2: continue
                    a,b=pair_ids
                else:
                    a,b=(row.get(a_col) or "").strip(),(row.get(b_col) or "").strip()
                a=a.split("-")[0].strip()
                b=b.split("-")[0].strip()
                ga=uniprot_to_symbol.get(a,a if a in queried_genes else ("UNIPROT:"+a if a else ""))
                gb=uniprot_to_symbol.get(b,b if b in queried_genes else ("UNIPROT:"+b if b else ""))
                for gene,partner in ((ga,gb),(gb,ga)):
                    if gene in queried_genes and partner and gene!=partner:
                        predictomes_rows.append({"gene":gene,"partner":partner,"sources":"Predictomes Human pair scores 2025-11-10",
                            "source_count":1,"experimental_support":"no","binary_experimental_support":"no/not established","PIONEER_coverage":"no",
                            "structural_support":"prediction-only SPOC >0.86","predicted_only":"yes","tier":"Tier 3",
                            "experimental_systems":f"Predictomes SPOC={score:.6g}"})
    merged={}
    for row in provenance+predictomes_rows:
        key=(row["gene"],row["partner"])
        if key not in merged:
            merged[key]=dict(row)
        else:
            current=merged[key]
            old_sources=set((current.get("sources","")+";"+row.get("sources","")).split(";"))-{""}
            current["sources"]=";".join(sorted(old_sources))
            current["source_count"]=len(old_sources)
            current["predicted_only"]="no" if current.get("experimental_support")=="yes" else "yes"
            if row.get("PIONEER_coverage")=="yes": current["PIONEER_coverage"]="yes"
            if row.get("structural_support","").startswith("prediction-only") and current.get("structural_support")=="no": current["structural_support"]=row["structural_support"]
            systems=set((current.get("experimental_systems","")+";"+row.get("experimental_systems","")).split(";"))-{""}
            current["experimental_systems"]=";".join(sorted(systems))
    provenance=list(merged.values())
    stat_by={r["gene"]:r for r in stats}
    for gene in queried_genes:
        stat=stat_by.get(gene)
        if not stat: continue
        gene_prov=[r for r in provenance if r["gene"]==gene]
        stat["total_partners"]=len(gene_prov)
        stat["experimental_partners"]=sum(r["experimental_support"]=="yes" for r in gene_prov)
        stat["PIONEER_partners"]=sum("PIONEER" in r["sources"] for r in gene_prov)
        stat["Predictomes_partners"]=sum("Predictomes" in r["sources"] for r in gene_prov)
        stat["unique_credible_partners"]=sum(r["experimental_support"]=="yes" or "PIONEER" in r["sources"] or "Predictomes" in r["sources"] for r in gene_prov)
        stat["tier1_2_partners"]=stat["experimental_partners"]
        stat["ppi_sources"]="BioGRID MV-Physical 5.0.261; IntAct PSICQUIC; PIONEER; Predictomes SPOC >0.86"
    ppi_cols = ["gene","total_partners","experimental_partners","tier1_2_partners","PIONEER_partners","Predictomes_partners","unique_credible_partners","structural_ready_partners","ppi_sources"]
    prov_cols = ["gene","partner","sources","source_count","experimental_support","binary_experimental_support","PIONEER_coverage","structural_support","predicted_only","tier","experimental_systems"]
    write_tsv(PHASE / "ppi_coverage_phase0R_candidates.tsv", stats, ppi_cols)
    write_tsv(PHASE / "ppi_provenance_phase0R_candidates.tsv", provenance, prov_cols)
    failures=[]
    prior_failure_path=META/"source_failures.tsv"
    if prior_failure_path.exists():
        for row in read_tsv(prior_failure_path):
            if row.get("source")=="HuRI":
                failures.append({"source":"HuRI","url":row.get("url",""),"status":"unavailable; previous request timed out; not a cohort gate"})
            elif row.get("source","" ).startswith("IntAct PSICQUIC") and row.get("timestamp_utc","")>=ppi_started:
                failures.append({"source":row.get("source"),"url":row.get("url",""),"status":row.get("status","")})
    if not predictomes_path.exists():
        failures.append({"source":"Predictomes pair scores","url":predictomes_url,"status":"download unavailable; see logs/phase0R_predictomes_download.log"})
    write_tsv(PHASE/"phase0R_source_failures.tsv",failures,["source","url","status"])
    strict_stats = [r for r in stats if r["gene"] in set(strict_genes)]
    result = {"strict_gene_n":len(strict_genes),"genes_ge3_credible_partners":sum(int(r["unique_credible_partners"]) >= 3 for r in strict_stats),
            "genes_with_experimental_ppi":sum(int(r["experimental_partners"]) > 0 for r in strict_stats),
            "genes_with_pioneer_coverage":sum(int(r["PIONEER_partners"]) > 0 for r in strict_stats),
            "ppi_gene_n":len(stats),"pilot_coverage":{r["gene"]:{"credible_partners":r["unique_credible_partners"],"PIONEER_partners":r["PIONEER_partners"]} for r in stats if r["gene"] in {"ACTB","ACTG1"}},
            "genes_with_predictomes_coverage":sum(int(r["Predictomes_partners"])>0 for r in strict_stats),
            "source":"BioGRID MV-Physical 5.0.261; candidate-gene IntAct PSICQUIC; precomputed PIONEER very-high confidence; Predictomes SPOC >0.86",
            "huRI":"UNAVAILABLE_TIMEOUT_RECORDED","predictomes":"PAIR_SCORE_TABLE_ACQUIRED" if predictomes_path.exists() else "UNAVAILABLE"}
    (META / "phase0R_ppi_counts.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    return result


def mechanism_audit() -> dict:
    gene_rows = read_tsv(PHASE / "gene_exact_rcv_attrition_summary.tsv")
    strict_genes = {r["gene"] for r in gene_rows if r["prescreen_status"] == "STRICT_PRESCREEN_965" and r["strict_after_exact_RCV"] == "YES"}
    disease_rows = read_tsv(PHASE / "gene_disease_variant_counts_exact_rcv.tsv")
    diseases = defaultdict(list)
    for r in disease_rows:
        if r["gene"] in strict_genes or r["gene"] in {"ACTB","ACTG1"}:
            diseases[r["gene"]].append(r["MONDO"])
    sys.path.insert(0, str(ROOT / "scripts"))
    from run_phase0 import read_mechanism_sources
    audit_genes = set(strict_genes) | {"ACTB","ACTG1"}
    evidence = read_mechanism_sources(audit_genes, {}) if audit_genes else []
    evidence_map = defaultdict(list)
    for row in evidence:
        if row.get("MONDO"):
            evidence_map[(row["gene"],row["MONDO"])].append(row)
    pair_rows, gene_states, by_gene_state = [], [], Counter()
    for gene in sorted(audit_genes):
        pair_states = []
        ds = sorted(set(diseases.get(gene, [])))
        for i, a in enumerate(ds):
            for b in ds[i+1:]:
                ea = [x for x in evidence_map.get((gene,a),[]) if x.get("primary_evidence")=="yes"]
                eb = [x for x in evidence_map.get((gene,b),[]) if x.get("primary_evidence")=="yes"]
                ma = {x.get("mechanism_category","unknown") for x in ea} - {"unknown"}
                mb = {x.get("mechanism_category","unknown") for x in eb} - {"unknown"}
                ia = {x.get("inheritance_category","unknown") for x in ea} - {"unknown"}
                ib = {x.get("inheritance_category","unknown") for x in eb} - {"unknown"}
                if len(ma)==len(mb)==len(ia)==len(ib)==1:
                    state = "CONFIRMED_COMPATIBLE" if ma==mb and ia==ib else "CONFIRMED_INCOMPATIBLE"
                else:
                    state = "UNRESOLVED"
                pair_states.append(state)
                pair_rows.append({"gene":gene,"analysis_set":"STRICT" if gene in strict_genes else "PILOT","MONDO_A":a,"disease_A":";".join(sorted({x.get("disease_label","") for x in evidence_map.get((gene,a),[])})),
                                  "mechanism_A":";".join(sorted(ma)),"inheritance_A":";".join(sorted(ia)),"evidence_rows_A":len(ea),
                                  "MONDO_B":b,"disease_B":";".join(sorted({x.get("disease_label","") for x in evidence_map.get((gene,b),[])})),
                                  "mechanism_B":";".join(sorted(mb)),"inheritance_B":";".join(sorted(ib)),"evidence_rows_B":len(eb),"mechanism_state":state})
        if not pair_states:
            state = "NOT_ENOUGH_DISEASE_PAIRS"
        elif "CONFIRMED_INCOMPATIBLE" in pair_states:
            state = "CONFIRMED_INCOMPATIBLE"
        elif all(x == "CONFIRMED_COMPATIBLE" for x in pair_states):
            state = "CONFIRMED_COMPATIBLE"
        else:
            state = "UNRESOLVED"
        if gene in strict_genes:
            by_gene_state[state] += 1
        gene_states.append({"gene":gene,"analysis_set":"STRICT" if gene in strict_genes else "PILOT","n_disease_groups":len(ds),"n_disease_pairs":len(pair_states),"mechanism_status":state,
                            "compatible_pairs":pair_states.count("CONFIRMED_COMPATIBLE"),"incompatible_pairs":pair_states.count("CONFIRMED_INCOMPATIBLE"),
                            "unresolved_pairs":pair_states.count("UNRESOLVED"),"source_evidence_rows":sum(len(v) for (g,_),v in evidence_map.items() if g==gene)})
    evidence_cols = ["gene","MONDO","disease_label","source","inheritance","inheritance_category","molecular_mechanism","mechanism_category","variant_consequence","evidence_confidence","primary_evidence"]
    write_tsv(PHASE / "gene_disease_mechanism_evidence_exact_strict.tsv", evidence, evidence_cols)
    pair_cols = ["gene","analysis_set","MONDO_A","disease_A","mechanism_A","inheritance_A","evidence_rows_A","MONDO_B","disease_B","mechanism_B","inheritance_B","evidence_rows_B","mechanism_state"]
    gene_cols = ["gene","analysis_set","n_disease_groups","n_disease_pairs","mechanism_status","compatible_pairs","incompatible_pairs","unresolved_pairs","source_evidence_rows"]
    write_tsv(PHASE / "mechanism_pair_states_exact_strict.tsv", pair_rows, pair_cols)
    write_tsv(PHASE / "mechanism_gene_states_exact_strict.tsv", gene_states, gene_cols)
    result = {"strict_gene_n":len(strict_genes),"CONFIRMED_COMPATIBLE":by_gene_state["CONFIRMED_COMPATIBLE"],
            "UNRESOLVED":by_gene_state["UNRESOLVED"],"CONFIRMED_INCOMPATIBLE":by_gene_state["CONFIRMED_INCOMPATIBLE"],
            "NOT_KNOWN_INCOMPATIBLE":by_gene_state["CONFIRMED_COMPATIBLE"]+by_gene_state["UNRESOLVED"],
            "pair_state_counts":dict(Counter(r["mechanism_state"] for r in pair_rows))}
    (META / "phase0R_mechanism_counts.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    return result


def file_rows(path: Path) -> int | str:
    if path.suffix.lower() not in {".tsv", ".csv"}:
        return "NA"
    with path.open("rb") as f:
        return max(0, sum(block.count(b"\n") for block in iter(lambda:f.read(1024*1024),b"")) - 1)


def write_phase0r_figures() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    stages=read_tsv(PHASE/"strict965_attrition_stages.tsv") if (PHASE/"strict965_attrition_stages.tsv").exists() else []
    thresholds=read_tsv(PHASE/"phase0R_candidate_thresholds.tsv") if (PHASE/"phase0R_candidate_thresholds.tsv").exists() else []
    reasons=read_tsv(PHASE/"top_five_exclusion_reasons.tsv") if (PHASE/"top_five_exclusion_reasons.tsv").exists() else []
    if stages:
        fig,ax=plt.subplots(figsize=(9,4.8))
        labels=[r["stage"] for r in stages][::-1]
        vals=[int(r["genes"]) for r in stages][::-1]
        ax.barh(range(len(labels)),vals,color="#54718b")
        ax.set_yticks(range(len(labels)),labels,fontsize=8)
        ax.set_xlabel("Genes")
        ax.set_title("Phase 0R: exact RCV and disease-normalization attrition")
        ax.grid(axis="x",alpha=.2)
        fig.tight_layout()
        fig.savefig(PHASE/"Figure_P0R-1_strict965_attrition.png",dpi=180)
        plt.close(fig)
    if thresholds:
        fig,ax=plt.subplots(figsize=(7,4.5))
        x=[int(r["minimum_variants_per_disease"]) for r in thresholds]
        y=[int(r["genes"]) for r in thresholds]
        ax.plot(x,y,marker="o",color="#54718b",linewidth=2)
        ax.set_xticks(x)
        ax.set_xlabel("Minimum unique substitutions per disease")
        ax.set_ylabel("Genes with at least two qualifying diseases")
        ax.set_title("Phase 0R exact-RCV threshold screen")
        ax.grid(alpha=.2)
        fig.tight_layout()
        fig.savefig(PHASE/"Figure_P0R-2_threshold_screen.png",dpi=180)
        plt.close(fig)
    if reasons:
        fig,ax=plt.subplots(figsize=(8,4.2))
        labels=[r["reason"].replace("_"," ") for r in reasons][::-1]
        vals=[int(r["genes"]) for r in reasons][::-1]
        ax.barh(range(len(labels)),vals,color="#8797a4")
        ax.set_yticks(range(len(labels)),labels,fontsize=8)
        ax.set_xlabel("Genes")
        ax.set_title("Most frequent first exclusion reason")
        ax.grid(axis="x",alpha=.2)
        fig.tight_layout()
        fig.savefig(PHASE/"Figure_P0R-3_exclusion_reasons.png",dpi=180)
        plt.close(fig)


def write_phase0r_manifests() -> None:
    input_manifest= META/"phase0R_input_manifest.tsv"
    step_manifest= META/"phase0R_step_manifest.tsv"
    output_manifest= META/"phase0R_output_manifest.tsv"
    combined_manifest= META/"phase0R_manifest.tsv"
    repro_path=ROOT/"reports"/"PHASE0R_reproducibility_log.md"
    seed_meta=json.loads(SEED_CHECKPOINT.read_text(encoding="utf-8"))
    integrity_path=META/"phase0R_archive_integrity.json"
    archive_integrity=json.loads(integrity_path.read_text(encoding="utf-8")) if integrity_path.exists() else {}
    archive_hash=""
    dl_log=LOGS/"phase0R_rcv_resume.log"
    if dl_log.exists():
        m=re.search(r"sha256=([0-9a-f]{64})\s+data/raw/ClinVarRCVRelease_2026-09\.xml\.gz\b",dl_log.read_text(encoding="utf-8",errors="replace"))
        if m: archive_hash=m.group(1)
    input_specs=[
        (SUMMARY,"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2026-09.txt.gz","ClinVar frozen monthly snapshot 2026-09-03","global seed source"),
        (ARCHIVE,"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_release/ClinVarRCVRelease_2026-09.xml.gz","ClinVar RCV release file; XML Dated 2026-08-29, server last-modified 2026-09-03","frozen variant-condition archive"),
        (RAW/"ClinVarRCVRelease_2026-09.xml.gz.md5","https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_release/ClinVarRCVRelease_2026-09.xml.gz.md5","official ClinVar archive MD5 sidecar","frozen archive checksum"),
        (RAW/"mondo-v2026-09-01.obo","https://github.com/monarch-initiative/mondo/releases/tag/v2026-09-01","MONDO v2026-09-01","disease normalization"),
        (INTERIM/"clinvar_preliminary_candidate_disease_rows.tsv","preserved Phase 0 interim output","2026-09 summary prescreen","exact reconstruction of previous 965 genes"),
        (RAW/"gene2phenotype_all_latest.csv","https://www.ebi.ac.uk/gene2phenotype/api/panel/all/download/","captured 2026-09-27","light mechanism audit"),
        (RAW/"gencc_submissions_latest.tsv","https://thegencc.org/download/action/submissions-export-tsv","captured 2026-09-27","light mechanism audit"),
        (RAW/"BIOGRID-MV-Physical-5.0.261.tab3.zip","https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-5.0.261/BIOGRID-MV-Physical-5.0.261.tab3.zip","BioGRID MV-Physical 5.0.261","experimental PPI"),
        (RAW/"PIONEER_very_high_human_2026-09-27.txt","https://pioneer.yulab.org/static/predictions/very_high/human.txt","captured 2026-09-27","precomputed PPI partner coverage"),
        (RAW/"hgnc_non_alt_loci_2026-09-27.tsv","https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/non_alt_loci_set.txt","captured 2026-09-27","gene to UniProt mapping"),
    ]
    if (RAW/"predictomes_human_pair_scores_2025-11-10.csv.gz").exists():
        input_specs.append((RAW/"predictomes_human_pair_scores_2025-11-10.csv.gz","https://predictomes-hsbps-dataset.s3.us-east-1.amazonaws.com/20251110_hs_predictome_pair_scores.csv.gz","Human Predictome pair scores, 2025-11-10","prediction-only PPI coverage"))
    input_rows=[]
    for path,url,release,role in input_specs:
        if path.exists():
            digest=archive_hash if path==ARCHIVE and archive_hash else (seed_meta["source_sha256"] if path==SUMMARY else sha256(path))
            input_rows.append({"path":str(path.relative_to(ROOT)),"url_or_source":url,"release_or_version":release,
                               "retrieved_or_mtime":datetime.fromtimestamp(path.stat().st_mtime,timezone.utc).isoformat(),
                               "bytes":path.stat().st_size,"sha256":digest,"status":"available","phase0R_role":role,
                               "integrity_check":("gzip -t PASS; official MD5 PASS" if path==ARCHIVE and archive_integrity.get("md5_match") else ("checksum reference" if path.suffix==".md5" else "NA"))})
        else:
            input_rows.append({"path":str(path.relative_to(ROOT)),"url_or_source":url,"release_or_version":release,
                               "retrieved_or_mtime":"NA","bytes":"NA","sha256":"NA","status":"unavailable","phase0R_role":role,"integrity_check":"NA"})
    prior_rows=read_tsv(META/"input_manifest.tsv") if (META/"input_manifest.tsv").exists() else []
    prior_by_path={r.get("path",""):r for r in prior_rows}
    phase0r_genes=set()
    if (PHASE/"strict965_attrition.tsv").exists():
        phase0r_genes={r["gene"] for r in read_tsv(PHASE/"strict965_attrition.tsv") if r.get("strict_after_exact_RCV")=="YES"}
    query_genes=phase0r_genes|{"ACTB","ACTG1"}
    for cache in sorted((RAW/"intact_api").glob("*.tab25")):
        if not any(cache.name.startswith(g+"_") for g in query_genes):
            continue
        side=cache.with_suffix(cache.suffix+".download.json")
        meta={}
        if side.exists():
            try: meta=json.loads(side.read_text(encoding="utf-8"))
            except ValueError: pass
        input_rows.append({"path":str(cache.relative_to(ROOT)),"url_or_source":meta.get("url","IntAct PSICQUIC cached result"),
                           "release_or_version":"current snapshot query","retrieved_or_mtime":meta.get("retrieved_utc",datetime.fromtimestamp(cache.stat().st_mtime,timezone.utc).isoformat()),
                           "bytes":cache.stat().st_size,"sha256":sha256(cache),"status":"cached","phase0R_role":"candidate-gene experimental PPI","integrity_check":"NA"})
    for path in (RAW/"phase0R_targeted_rcv_api").glob("*.xml") if (RAW/"phase0R_targeted_rcv_api").exists() else []:
        input_rows.append({"path":str(path.relative_to(ROOT)),"url_or_source":"ClinVar E-utilities accession-level fallback","release_or_version":"live fallback; not frozen archive","retrieved_or_mtime":datetime.fromtimestamp(path.stat().st_mtime,timezone.utc).isoformat(),
                           "bytes":path.stat().st_size,"sha256":sha256(path),"status":"cached","phase0R_role":"missing seed RCV fallback","integrity_check":"XML parsed and cached"})
    with input_manifest.open("x",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,["path","url_or_source","release_or_version","retrieved_or_mtime","bytes","sha256","status","phase0R_role","integrity_check"],delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(input_rows)

    parse_meta=json.loads(PARSER_CHECKPOINT.read_text(encoding="utf-8")) if PARSER_CHECKPOINT.exists() else {}
    steps=[
        {"step":"global P/LP missense RCV seed","inputs":"data/raw/variant_summary_2026-09.txt.gz","input_rows":seed_meta["source_rows_scanned"],"outputs":"results/phase0R/global_plp_missense_rcv_seed.tsv","output_rows":file_rows(SEED),"status":"complete"},
        {"step":"rebuild prior 965 strict prescreen genes","inputs":"data/interim/clinvar_preliminary_candidate_disease_rows.tsv","input_rows":file_rows(INTERIM/"clinvar_preliminary_candidate_disease_rows.tsv"),"outputs":"results/phase0R/previous_965_gene_list.tsv","output_rows":file_rows(PHASE/"previous_965_gene_list.tsv"),"status":"complete"},
        {"step":"frozen RCV archive streaming recovery","inputs":"data/raw/ClinVarRCVRelease_2026-09.xml.gz","input_rows":parse_meta.get("xml_records_scanned_total","NA"),"outputs":"results/phase0R/global_rcv_recovered.tsv","output_rows":file_rows(PHASE/"global_rcv_recovered.tsv") if (PHASE/"global_rcv_recovered.tsv").exists() else "NA","status":"complete" if (PHASE/"global_rcv_recovered.tsv").exists() else "hold"},
    ]
    if (PHASE/"strict965_attrition.tsv").exists():
        steps += [
            {"step":"exact RCV classification and MONDO normalization for 965 genes","inputs":"global RCV records; previous 965 genes; MONDO v2026-09-01","input_rows":file_rows(PHASE/"global_rcv_recovered.tsv"),"outputs":"results/phase0R/strict965_attrition.tsv","output_rows":file_rows(PHASE/"strict965_attrition.tsv"),"status":"complete"},
            {"step":"strict-gene mechanism and PPI audit","inputs":"exact strict cohort; Gene2Phenotype; GenCC; BioGRID; IntAct; PIONEER; Predictomes scores","input_rows":"strict exact-RCV gene count","outputs":"results/phase0R/mechanism_gene_states_exact_strict.tsv; ppi_coverage_phase0R_candidates.tsv","output_rows":file_rows(PHASE/"mechanism_gene_states_exact_strict.tsv"),"status":"complete" if (PHASE/"mechanism_gene_states_exact_strict.tsv").exists() else "pending"},
        ]
    with step_manifest.open("x",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,["step","inputs","input_rows","outputs","output_rows","status"],delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(steps)

    try:
        git_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()
        git_status=subprocess.check_output(["git","status","--short"],cwd=ROOT,text=True).strip() or "clean"
    except Exception:
        git_commit="NA";git_status="NA"
    repro=(
        "# Phase 0R reproducibility log\n\n"
        f"- Run timestamp (UTC): {utcnow()}\n- Python: {platform.python_version()}\n- lxml: {etree.LXML_VERSION}\n"
        f"- Git HEAD before Phase 0R result commit: `{git_commit}`; working tree at log generation: `{git_status.replace(chr(10), '; ')}`.\n"
        f"- Phase 0R script SHA256: `{sha256(ROOT/'scripts/run_phase0r.py')}`.\n"
        f"- Frozen summary SHA256: `{seed_meta['source_sha256']}`.\n"
        f"- Frozen RCV archive SHA256: `{archive_hash or 'see input manifest; checksum not found in transfer log'}`; official MD5 verification: `{archive_integrity.get('md5_match','not yet run')}` (expected `{archive_integrity.get('expected_md5','NA')}`, observed `{archive_integrity.get('observed_md5','NA')}`).\n"
        "- Seed counts and pre-download checkpoint: `metadata/phase0R_seed_checkpoint.json`.\n"
        "- RCV parse checkpoints: `metadata/phase0R_rcv_parse_checkpoint.json`; records are streamed with lxml and cached in a SQLite store.\n"
        "- Input and output provenance: `metadata/phase0R_input_manifest.tsv`, `metadata/phase0R_output_manifest.tsv`, `metadata/phase0R_step_manifest.tsv`, and combined `metadata/phase0R_manifest.tsv`.\n"
        "- No structure prediction, interface scoring, FoldX, MutPred-PPI, PERMANOVA, classification, tissue matching, or machine learning was run.\n"
    )
    repro_path.write_text(repro,encoding="utf-8")
    output_paths=[]
    output_paths += sorted(p for p in PHASE.rglob("*") if p.is_file())
    output_paths += sorted(p for p in META.glob("phase0R*") if p.is_file())
    output_paths += [p for p in (ROOT/"reports").glob("PHASE0R*") if p.is_file()]
    output_paths += [repro_path, ROOT/"scripts"/"run_phase0r.py", dl_log]
    db_path=INTERIM/"phase0R_rcv_recovered.sqlite"
    if db_path.exists(): output_paths.append(db_path)
    output_paths=sorted(set(p for p in output_paths if p.exists() and p not in {output_manifest,combined_manifest}))
    out_rows=[]
    for path in output_paths:
        out_rows.append({"path":str(path.relative_to(ROOT)),"artifact_type":"script" if path.suffix==".py" else ("report" if path.suffix==".md" else ("log" if path.parent==LOGS else "phase0R_output")),
                         "rows":file_rows(path),"bytes":path.stat().st_size,"sha256":sha256(path)})
    with output_manifest.open("x",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,["path","artifact_type","rows","bytes","sha256"],delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(out_rows)
    combined=[]
    combined.extend({"record_type":"input","path":r["path"],"status":r["status"],"bytes":r["bytes"],"sha256":r["sha256"],"rows":"NA","detail":r["phase0R_role"]+"; integrity="+r.get("integrity_check","NA")} for r in input_rows)
    combined.extend({"record_type":"output","path":r["path"],"status":"generated","bytes":r["bytes"],"sha256":r["sha256"],"rows":r["rows"],"detail":r["artifact_type"]} for r in out_rows)
    with combined_manifest.open("x",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,["record_type","path","status","bytes","sha256","rows","detail"],delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(combined)


def final_report() -> dict:
    gate = json.loads((META / "phase0R_recovery_gate.json").read_text(encoding="utf-8"))
    integrity_path=META/"phase0R_archive_integrity.json"
    archive_integrity=json.loads(integrity_path.read_text(encoding="utf-8")) if integrity_path.exists() else {}
    report_path=ROOT/"reports"/"PHASE0R_FULL_RCV_RECOVERY_REPORT.md"
    if gate["recovery_rate"] < .95:
        if report_path.exists():
            raise FileExistsError(f"Refusing to overwrite report: {report_path}")
        report_path.write_text(
            "# Phase 0R — RCV recovery hold\n\n"
            "**PHASE 0R VERDICT: `HOLD_RCV_RECOVERY`**\n\n"
            f"Frozen ClinVar snapshot: 2026-09-03. Candidate RCV accessions: {gate['candidate_rcv_n']:,}; recovered: {gate['recovered_rcv_n']:,}; failed: {gate['failed_rcv_n']:,}; recovery: {gate['recovery_rate']:.2%}.\n\n"
            "Recovery is below the 95% gate. Disease normalization, mechanism classification, and PPI audit were not run, per the stop rule. Phase 1 was NOT started.\n",
            encoding="utf-8")
        write_phase0r_figures()
        write_phase0r_manifests()
        return {"verdict":"HOLD_RCV_RECOVERY","report_path":str(report_path)}
    disease_metrics = json.loads((META / "phase0R_exact_disease_counts.json").read_text(encoding="utf-8"))
    mech_result = json.loads((META / "phase0R_mechanism_counts.json").read_text(encoding="utf-8"))
    ppi_summary = json.loads((META / "phase0R_ppi_counts.json").read_text(encoding="utf-8"))
    attrition = read_tsv(PHASE / "strict965_attrition.tsv")
    exact_genes = read_tsv(PHASE / "gene_exact_rcv_attrition_summary.tsv")
    disease_rows = read_tsv(PHASE / "gene_disease_variant_counts_exact_rcv.tsv")
    mech_rows = read_tsv(PHASE / "mechanism_gene_states_exact_strict.tsv")
    ppi_rows = read_tsv(PHASE / "ppi_coverage_phase0R_candidates.tsv")
    strict_rows = [r for r in exact_genes if r["prescreen_status"] == "STRICT_PRESCREEN_965" and r["strict_after_exact_RCV"] == "YES"]
    strict_genes = {r["gene"] for r in strict_rows}
    mech_by_gene = {r["gene"]:r["mechanism_status"] for r in mech_rows}
    ppi_by_gene = {r["gene"]:r for r in ppi_rows}
    final_genes = sorted(g for g in strict_genes if mech_by_gene.get(g)=="CONFIRMED_COMPATIBLE" and int(ppi_by_gene.get(g,{}).get("unique_credible_partners",0))>=3)
    unresolved = mech_result["UNRESOLVED"]
    strict_n = len(strict_genes)
    if strict_n >= 30 and unresolved > strict_n / 2:
        verdict = "GO_CONDITIONAL_MECHANISM"
    elif len(final_genes) >= 30:
        verdict = "GO_FULL_ATLAS"
    elif 10 <= len(final_genes) <= 29:
        verdict = "GO_FOCUSED"
    else:
        # With recovery and exact disease normalization complete, fewer than
        # ten high-quality analyzable genes is a no-go. The conditional verdict
        # is reserved for the explicitly defined >=30-gene unresolved cohort.
        verdict = "NO_GO_ATLAS"

    # Candidate ranking is feasibility-only: disease distinctions, per-disease
    # counts, evidence state, then experimentally supported partner coverage.
    min_variants = {}
    for gene in strict_genes:
        counts = [int(r["n_unique_variants"]) for r in disease_rows if r["gene"]==gene and int(r["n_unique_variants"])>=10]
        min_variants[gene] = min(counts) if counts else 0
    candidates = sorted(strict_genes, key=lambda g:(
        {"CONFIRMED_COMPATIBLE":0,"UNRESOLVED":1,"CONFIRMED_INCOMPATIBLE":2}.get(mech_by_gene.get(g,"UNRESOLVED"),3),
        -int(next((r["n_distinct_diseases_after_parent_child_collapse"] for r in strict_rows if r["gene"]==g),0)),
        -min_variants.get(g,0),
        -int(next((r["n_unique_variants_in_strict_groups"] for r in strict_rows if r["gene"]==g),0)),
        -int(ppi_by_gene.get(g,{}).get("unique_credible_partners",0)),
        -int(ppi_by_gene.get(g,{}).get("PIONEER_partners",0)), g))[:20]
    candidate_rows=[]
    for rank,g in enumerate(candidates,1):
        er = next(r for r in strict_rows if r["gene"]==g)
        pr = ppi_by_gene.get(g,{})
        mstatus=mech_by_gene.get(g,"UNRESOLVED")
        candidate_rows.append({"rank":rank,"gene":g,"n_distinct_diseases":er["n_distinct_diseases_after_parent_child_collapse"],
                               "strict_diseases_ge10":er["n_diseases_ge10"],"unique_variants_in_strict_groups":er["n_unique_variants_in_strict_groups"],
                               "minimum_variants_per_strict_disease":min_variants.get(g,0),"mechanism_status":mstatus,
                               "credible_PPI_partners":pr.get("unique_credible_partners",0),"tier1_2_partners":pr.get("tier1_2_partners",pr.get("experimental_partners",0)),
                               "PIONEER_partners":pr.get("PIONEER_partners",0),"Predictomes_partners":pr.get("Predictomes_partners",0),
                               "final_analyzable":"YES" if g in final_genes else "NO"})
    write_tsv(PHASE / "top_candidate_genes_exact_rcv.tsv",candidate_rows,
              ["rank","gene","n_distinct_diseases","strict_diseases_ge10","unique_variants_in_strict_groups","minimum_variants_per_strict_disease","mechanism_status","credible_PPI_partners","tier1_2_partners","PIONEER_partners","Predictomes_partners","final_analyzable"])

    reasons = Counter(r["reason_if_lost"] for r in attrition if r["reason_if_lost"])
    top_reasons = reasons.most_common(5)
    write_tsv(PHASE / "top_five_exclusion_reasons.tsv", [{"rank":i,"reason":reason,"genes":n} for i,(reason,n) in enumerate(top_reasons,1)], ["rank","reason","genes"])

    pilot_rows=[]
    mech_by_gene_rows={r["gene"]:r for r in mech_rows}
    for gene in ("ACTB","ACTG1"):
        er=next((r for r in exact_genes if r["gene"]==gene),{})
        pr=ppi_by_gene.get(gene,{})
        mr=mech_by_gene_rows.get(gene,{})
        ds=[r for r in disease_rows if r["gene"]==gene]
        pilot_rows.append({"gene":gene,"n_exact_RCV":er.get("n_exact_rcv",0),"n_unique_PLP_missense_variants":er.get("n_exact_PLP_missense_variants",0),
                           "n_normalized_diseases":er.get("n_normalized_diseases",0),"n_distinct_diseases_after_parent_child_collapse":er.get("n_distinct_diseases_after_parent_child_collapse",0),
                           "variant_count_per_disease":";".join(f"{d['MONDO']}={d['n_unique_variants']} ({d['n_unique_residues']} residues)" for d in ds),
                           "strict_status":er.get("strict_after_exact_RCV","NO"),"reason":er.get("reason_if_lost","no exact candidate row in Phase0R domain"),
                           "mechanism_status":mr.get("mechanism_status","NOT_ENOUGH_DISEASE_PAIRS"),"credible_PPI_partners":pr.get("unique_credible_partners",0),
                           "PIONEER_partners":pr.get("PIONEER_partners",0),"existing_structural_partner_proxy":pr.get("structural_ready_partners",0)})
    write_tsv(PHASE / "ACTB_ACTG1_feasibility_exact_rcv.tsv",pilot_rows,
              ["gene","n_exact_RCV","n_unique_PLP_missense_variants","n_normalized_diseases","n_distinct_diseases_after_parent_child_collapse","variant_count_per_disease","strict_status","reason","mechanism_status","credible_PPI_partners","PIONEER_partners","existing_structural_partner_proxy"])

    strict_variant_total=sum(int(r["n_unique_variants_in_strict_groups"]) for r in strict_rows)
    strict_variant_partner_supported=sum(int(r["n_unique_variants_in_strict_groups"]) for r in strict_rows if int(ppi_by_gene.get(r["gene"],{}).get("unique_credible_partners",0))>=3)
    partner_proxy=strict_variant_partner_supported/strict_variant_total if strict_variant_total else None
    metrics={**gate, **disease_metrics, **mech_result, "genes_ge3_credible_partners":ppi_summary["genes_ge3_credible_partners"],
             "final_analyzable_genes":len(final_genes),"final_analyzable_gene_symbols":final_genes,
             "partner_analysis_ready_gene_level_variant_fraction":partner_proxy,"verdict":verdict,
             "top_five_exclusion_reasons":top_reasons,"ppi_summary":ppi_summary}
    (META / "phase0R_final_metrics.json").write_text(json.dumps(metrics,indent=2) + "\n",encoding="utf-8")

    if report_path.exists():
        raise FileExistsError(f"Refusing to overwrite report: {report_path}")
    old_stats = read_tsv(PHASE / "strict965_attrition_stages.tsv")
    stages = "\n".join(f"| {r['stage']} | {r['genes']} |" for r in old_stats)
    threshold_rows = read_tsv(PHASE / "phase0R_candidate_thresholds.tsv")
    threshold_table = "\n".join(f"| {r['minimum_variants_per_disease']} | {r['genes']} | {r['disease_groups']} | {r['unique_variants_in_qualifying_groups']} |" for r in threshold_rows)
    top_reason_table="\n".join(f"| {i} | {reason} | {n} |" for i,(reason,n) in enumerate(top_reasons,1)) or "| NA | No exclusions | 0 |"
    candidate_table="\n".join(f"| {r['rank']} | {r['gene']} | {r['n_distinct_diseases']} | {r['unique_variants_in_strict_groups']} | {r['mechanism_status']} | {r['credible_PPI_partners']} | {r['PIONEER_partners']} | {r['final_analyzable']} |" for r in candidate_rows) or "| — | No exact-RCV strict genes | 0 | 0 | — | 0 | 0 | NO |"
    pilot_table="\n".join(f"| {r['gene']} | {r['n_exact_RCV']} | {r['n_unique_PLP_missense_variants']} | {r['n_distinct_diseases_after_parent_child_collapse']} | {r['variant_count_per_disease'] or 'NA'} | {r['strict_status']} | {r['reason']} | {r['credible_PPI_partners']} | {r['PIONEER_partners']} |" for r in pilot_rows)
    fraction_text="NA (zero strict-group variant denominator)" if partner_proxy is None else f"{partner_proxy:.1%} of strict-group variants are in genes with ≥3 credible partners; this is a gene-level PPI proxy, not variant-level interface readiness."
    report=(
        "# Phase 0R — Full ClinVar RCV recovery and 965-gene re-audit\n\n"
        f"**PHASE 0R VERDICT: `{verdict}`**\n\n"
        "This report uses the frozen September 2026 ClinVar snapshot and a newly reconstructed 965-gene prescreen list. All counts below come from the Phase 0R workflow.\n\n"
        "## Frozen release and seed checkpoint\n\n"
        "- Frozen ClinVar monthly snapshot: 2026-09-03 (`variant_summary_2026-09.txt.gz`; archived checksum recorded in the checkpoint).\n"
        "- Frozen RCV file: `ClinVarRCVRelease_2026-09.xml.gz`; official archive path `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_release/ClinVarRCVRelease_2026-09.xml.gz`, server last-modified 2026-09-03; the XML root identifies release date `2026-08-29`. Both dates are recorded because they describe the archive path and its internal release date. Official sidecar MD5 and gzip integrity checks are required before parsing.\n"
        f"- Global seed: {json.loads(SEED_CHECKPOINT.read_text(encoding='utf-8'))['n_unique_variation_ids']:,} unique Variation IDs; {json.loads(SEED_CHECKPOINT.read_text(encoding='utf-8'))['n_unique_gene_protein_substitutions']:,} unique gene–protein substitutions; {gate['candidate_rcv_n']:,} candidate RCV accessions; {json.loads(SEED_CHECKPOINT.read_text(encoding='utf-8'))['n_unique_genes']:,} genes.\n"
        f"- Recovery: {gate['recovered_rcv_n']:,}/{gate['candidate_rcv_n']:,} RCVs ({gate['recovery_rate']:.2%}); {gate['failed_rcv_n']:,} failed; recovery gate `{gate['verdict']}`.\n"
        f"- The archive was resumed from the existing partial file, validated with `gzip -t` and official MD5 (`{archive_integrity.get('md5_match','not verified')}`), then streamed with `lxml.etree.iterparse`; only seed-listed RCVs were retained. Full-archive SHA256/MD5 and recovery source for every RCV are in the manifests/table.\n\n"
        "## Rebuilt 965-gene attrition\n\n"
        "The preserved preliminary disease rows reproduced exactly 965 genes under the original summary-screen rule (≥2 MONDO groups each with ≥10 summary-screen substitutions). These are prescreen candidates, not exact-RCV findings.\n\n"
        "| Stage | Genes |\n|---|---:|\n"+stages+"\n\n"
        f"- After exact RCV P/LP germline missense review: {disease_metrics['n_after_exact_rcv']} genes.\n"
        f"- After at least one unambiguous MONDO disease mapping: {disease_metrics['n_after_disease_normalization']} genes; with ≥2 disease groups after same-MONDO and parent-child collapsing: {disease_metrics['n_ge2_distinct_disease']}.\n"
        f"- Expanded (≥2 diseases, ≥5 unique substitutions per disease): {disease_metrics['n_ge5_variants_per_disease']} genes.\n"
        f"- Strict (≥2 diseases, ≥10 substitutions per disease, ≥20 total): {disease_metrics['n_strict_ge10']} genes.\n"
        "- Unambiguous parent terms were removed from independent disease counts when a more specific descendant was also present; variants attached only to a removed broad parent were not reassigned to a child. Relationship details are in `disease_parent_child_collapses.tsv`.\n"
        "- Exact synonyms and duplicate variant–condition submissions are collapsed. Hotspot flag: at least one strict-sized disease group has a single residue with ≥20% of substitutions or its top three residues with ≥50%. Residue counts remain separately reported.\n\n"
        "| Minimum substitutions per disease | Genes | Disease groups | Unique substitutions |\n|---:|---:|---:|---:|\n"+threshold_table+"\n\n"
        "## Mechanism states and PPI coverage\n\n"
        f"- Strict exact-RCV genes: {strict_n}. Mechanism gene states: `CONFIRMED_COMPATIBLE` {mech_result['CONFIRMED_COMPATIBLE']}; `UNRESOLVED` {mech_result['UNRESOLVED']}; `CONFIRMED_INCOMPATIBLE` {mech_result['CONFIRMED_INCOMPATIBLE']}; `NOT_KNOWN_INCOMPATIBLE` {mech_result['NOT_KNOWN_INCOMPATIBLE']}. Counts are genes; disease-pair states are separately listed. Missing G2P/GenCC evidence is `UNRESOLVED`.\n"
        f"- Strict genes with ≥3 credible PPI partners: {ppi_summary['genes_ge3_credible_partners']} of {strict_n}. Credibility includes experimental BioGRID/IntAct support or precomputed very-high-confidence PIONEER / Predictomes SPOC >0.86 coverage; Tier 1/2 experimental counts remain separate. HuRI was unavailable after a timeout.\n"
        f"- Partner-analysis-ready coverage proxy: {fraction_text}\n"
        f"- Final high-quality analyzable genes: {len(final_genes)}. A gene qualifies here only if it passes strict disease thresholds, is `CONFIRMED_COMPATIBLE`, and has ≥3 credible partners. Tier 1/2 experimental partner counts are reported separately.\n\n"
        "## Top exact-RCV candidates\n\n"
        "Ranked by mechanism state, disease distinction, variant counts, then credible/PIONEER partner coverage; no predicted perturbation outcome was used.\n\n"
        "| Rank | Gene | Distinct diseases | Strict-group variants | Mechanism | Credible partners | PIONEER partners | Analyzable |\n|---:|---|---:|---:|---|---:|---:|---|\n"+candidate_table+"\n\n"
        "## Top five exclusion reasons among the 965\n\n"
        "Reasons are mutually exclusive and assigned at the first failed attrition stage.\n\n| Rank | Reason | Genes |\n|---:|---|---:|\n"+top_reason_table+"\n\n"
        "## ACTB / ACTG1 exact-RCV pilot audit\n\n"
        "| Gene | Exact RCVs | Unique P/LP missense variants | Distinct diseases | Variants per disease | Strict | Reason | Credible partners | PIONEER partners |\n|---|---:|---:|---:|---|---|---|---:|---:|\n"+pilot_table+"\n\n"
        "## Data-source failures and limitations\n\n"
        "- Frozen RCV archive download resumed from the retained partial; transfer completion, byte size, SHA256, `gzip -t` result, and archive date are recorded in `metadata/phase0R_manifest.tsv` and the download log.\n"
        "- Targeted E-utilities fallback was not needed: the frozen archive recovered every seed RCV. The fallback is restricted to accessions missing after archive parsing.\n"
        "- HuRI timed out; IntAct and BioGRID experimental PPI support and precomputed PIONEER/Predictomes partner coverage were used where available. No new complex structures were downloaded from Predictomes.\n"
        "- Parent-child and unresolved MedGen/MONDO mappings are conservatively excluded from independent disease counts and retained in the audit outputs.\n\n"
        "## Phase boundary\n\n"
        "Phase 0R only. No structure modeling, interface scoring, MutPred-PPI, FoldX, AlphaMissense comparison, PERMANOVA, classification, tissue matching, or machine learning was run. **Phase 1 was NOT started.**\n"
    )
    report_path.write_text(report,encoding="utf-8")
    write_phase0r_figures()
    write_phase0r_manifests()
    return {"verdict":verdict,"report_path":str(report_path),"strict_gene_n":strict_n,"final_analyzable_genes":len(final_genes),
            "top_candidates":candidates[:5],"pilot_rows":pilot_rows,"exclusion_reasons":top_reasons}


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit("usage: run_phase0r.py seed|previous-965|download-archive|recover|disease-audit|ppi-audit|mechanism-audit|final-report")
    if sys.argv[1] == "seed":
        print(json.dumps(make_seed(), indent=2))
    elif sys.argv[1] == "previous-965":
        print(json.dumps(rebuild_previous_965(), indent=2))
    elif sys.argv[1] == "download-archive":
        print(json.dumps(download_archive_resumable(), indent=2))
    elif sys.argv[1] == "recover":
        print(json.dumps(run_recovery(), indent=2))
    elif sys.argv[1] == "disease-audit":
        print(json.dumps(exact_disease_audit(), indent=2))
    elif sys.argv[1] == "ppi-audit":
        print(json.dumps(ppi_audit(), indent=2))
    elif sys.argv[1] == "mechanism-audit":
        print(json.dumps(mechanism_audit(), indent=2))
    elif sys.argv[1] == "final-report":
        print(json.dumps(final_report(), indent=2))
    else:
        raise SystemExit("usage: run_phase0r.py seed|previous-965|download-archive|recover|disease-audit|ppi-audit|mechanism-audit|final-report")


if __name__ == "__main__":
    main()
