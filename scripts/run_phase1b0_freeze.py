#!/usr/bin/env python3
"""Build the Phase 1B-0 data freeze and pre-outcome analysis representations.

This script reads frozen Phase 1A tables, creates the requested Phase 1B-0
tables, and deliberately performs no disease-label association analysis.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import ast
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
import re
import subprocess
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
P1A = ROOT / "results/phase1A"
P0R = ROOT / "results/phase0R"
OUT = ROOT / "results/phase1B0"
META = ROOT / "metadata"
REPORT = ROOT / "reports/PHASE1B0_DISRUPTION_MATRIX_AND_ANALYSIS_FREEZE.md"
FREEZE = META / "PHASE1B0_INPUT_FREEZE_MANIFEST.tsv"
RAW_SCORE = ROOT / "data/interim/phase1B0_mutpred_raw_scores.tsv"
AM_FILE = ROOT / "data/raw/phase1B0/AlphaMissense_aa_substitutions.tsv.gz"
AFDB_META = META / "phase1b0_afdb_freeze.json"
DOMAIN_CACHE = ROOT / "data/interim/phase1B0_uniprot_domains"
PRIMARY_MODEL_COMMIT = "01a53eaf5661382a97a9b05a0b41e7f699cd0974"
MODEL_DIR = ROOT / "data/software/mutpred-ppi/weights"
TRAIN_DIR = ROOT / "data/interim/phase1A_mutpred_training/datasets/training_eval"
TRAIN_FILES = [
    TRAIN_DIR / "fragoza_only_mapped090826_rows.csv.gz",
    TRAIN_DIR / "sahni_fragoza_mapped090826_rows.csv.gz",
    TRAIN_DIR / "sahni_only_mapped090826_rows.csv.gz",
]

OUT.mkdir(parents=True, exist_ok=True)
META.mkdir(parents=True, exist_ok=True)
DOMAIN_CACHE.mkdir(parents=True, exist_ok=True)


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def row_count(path: Path) -> int:
    with path.open("rb") as f:
        return max(0, sum(1 for _ in f) - 1)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA", lineterminator="\n")


def read_tsv(path: Path, **kw) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False, **kw)


def split_semicolon(value) -> list[str]:
    if pd.isna(value) or not str(value).strip():
        return []
    return sorted({x.strip() for x in str(value).split(";") if x.strip()})


def phase1a_files() -> list[Path]:
    files = sorted(p for p in P1A.glob("*.tsv") if not p.name.startswith("._"))
    if len(files) != 12:
        raise RuntimeError(f"Expected 12 Phase 1A TSVs, found {len(files)}")
    return files


def freeze_inputs() -> dict[str, str]:
    rows = []
    before = {}
    for path in phase1a_files():
        digest = sha256(path)
        before[str(path)] = digest
        rows.append({
            "path": str(path.relative_to(ROOT)),
            "SHA256": digest,
            "row_count": row_count(path),
            "modification_time_local": datetime.fromtimestamp(
                path.stat().st_mtime).astimezone().isoformat(),
            "originating_phase": "PHASE1A",
        })
    write_tsv(pd.DataFrame(rows), FREEZE)
    return before


def verify_inputs_unchanged(before: dict[str, str]) -> None:
    changed = [p for p, digest in before.items() if sha256(Path(p)) != digest]
    if changed:
        raise RuntimeError("Phase 1A frozen inputs changed during Phase 1B-0: " + ", ".join(changed))


def load_phase1a_helpers():
    path = ROOT / "scripts/run_phase1a.py"
    spec = importlib.util.spec_from_file_location("phase1a_helpers", path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def public_training_sets() -> tuple[set, set, set, int]:
    variants, pairs, proteins = set(), set(), set()
    row_n = 0
    for path in TRAIN_FILES:
        if not path.exists():
            continue
        d = pd.read_csv(path, low_memory=False)
        row_n += len(d)
        for r in d.itertuples(index=False):
            acc = str(r.interactor)
            partner = str(r.partner)
            mut = str(r.mutation)
            variants.add((acc, mut))
            pairs.add((acc, partner, mut))
            proteins.add(acc); proteins.add(partner)
    return variants, pairs, proteins, row_n


def canonical_variant(row: pd.Series) -> str:
    ref_value = row.get("canonical_reference_aa", "")
    alt_value = row.get("alternate_aa", "")
    ref = "" if pd.isna(ref_value) else str(ref_value)
    alt = "" if pd.isna(alt_value) else str(alt_value)
    pos = row.get("canonical_residue", "")
    try:
        if ref and alt and int(pos) > 0:
            return f"{ref}{int(pos)}{alt}"
    except Exception:
        pass
    change = row.get("protein_change", "")
    if not pd.isna(change):
        match = re.search(r"p\.?\(?([A-Z])([0-9]+)([A-Z])\)?", str(change))
        if match:
            return f"{match.group(1)}{match.group(2)}{match.group(3)}"
    return ""


def get_partner_accessions(matrix: pd.DataFrame, mapping: pd.DataFrame) -> tuple[dict, dict]:
    """Reuse frozen PIONEER target/partner accessions for PDB-tagged edges."""
    helpers = load_phase1a_helpers()
    hgnc_path = ROOT / "data/raw/hgnc_non_alt_loci_2026-09-27.tsv"
    hgnc = read_tsv(hgnc_path)
    genes = set(matrix.gene.dropna().astype(str))
    by_edge, _ = helpers.build_pioneer(genes, hgnc)
    canon = mapping.drop_duplicates(["gene", "VariationID"]).set_index(
        ["gene", "VariationID"])["canonical_UniProt_accession"].to_dict()
    partner_acc, status = {}, {}
    for edge, records in by_edge.items():
        target_accessions = set(matrix.loc[(matrix.gene == edge[0]) &
                                           (matrix.partner == edge[1]),
                                           "canonical_UniProt_accession"].dropna().astype(str))
        candidates = sorted({str(x.get("partner_uniprot", "")) for x in records
                             if str(x.get("source", "")).upper() == "PDB"
                             and str(x.get("target_uniprot", "")) in target_accessions})
        candidates = [x for x in candidates if x]
        if len(candidates) == 1:
            partner_acc[edge] = candidates[0]; status[edge] = "EXACT_ONE_PIONEER_PDB_ACCESSION"
        elif len(candidates) > 1:
            partner_acc[edge] = ";".join(candidates); status[edge] = "AMBIGUOUS_MULTIPLE_PIONEER_PDB_ACCESSIONS"
        else:
            status[edge] = "NO_PIONEER_PDB_ACCESSION"
    return partner_acc, status


def freeze_mutpred_version() -> dict:
    freeze_path = META / "PHASE1B0_MUTPRED_VERSION_FREEZE.tsv"
    repo = ROOT / "data/software/mutpred-ppi"
    commit = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                            check=True, capture_output=True, text=True).stdout.strip()
    if commit != PRIMARY_MODEL_COMMIT:
        raise RuntimeError(f"MutPred-PPI checkout changed: expected {PRIMARY_MODEL_COMMIT}, got {commit}")
    entries = []
    for name in ["MutPred-PPI.pt", "mutation_diff_scaler.pkl"]:
        path = MODEL_DIR / name
        entries.append({"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
                        "SHA256": sha256(path)})
    row = {
        "primary_version": "MutPred-PPI current repository checkpoint (arch=current)",
        "model_release": "current all-data checkpoint; not v1.0 ensemble",
        "repository_commit": commit,
        "access_method": "local official inference pipeline; source checkout",
        "executable": "data/software/mutpred-ppi/src/inference/02_run_mutpred-ppi_inference.py",
        "date_frozen_utc": utcnow(),
        "reference": "Stewart et al. RECOMB 2026; doi:10.64898/2025.12.20.695738",
        "checkpoint_hashes": json.dumps(entries, ensure_ascii=False),
        "required_sequence_model": "Rostlab/prot_t5_xl_half_uniref50-enc; official Hugging Face retrieval currently failed TLS",
        "scoring_state": "VERSION_FROZEN",
    }
    if freeze_path.exists():
        old = read_tsv(freeze_path).iloc[0].to_dict()
        if (old.get("repository_commit") != commit or old.get("checkpoint_hashes") != row["checkpoint_hashes"]):
            raise RuntimeError("The frozen MutPred-PPI checkpoint or scaler changed after version freeze")
        return old
    write_tsv(pd.DataFrame([row]), freeze_path)
    return row


def fetch_uniprot_domain(accession: str) -> dict | None:
    path = DOMAIN_CACHE / f"{accession}.json"
    if path.exists():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            pass
    import requests
    last = None
    for attempt in range(4):
        try:
            r = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.json",
                             timeout=(12, 45), headers={"Accept": "application/json"})
            r.raise_for_status()
            data = r.json()
            data["_phase1b0_retrieved_utc"] = utcnow()
            path.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
            return data
        except Exception as exc:
            last = exc
            time.sleep(1 + attempt)
    return {"_fetch_error": f"{type(last).__name__}: {last}"}


def build_domain_baseline(variants: pd.DataFrame, accessions: list[str]) -> tuple[pd.DataFrame, list[dict]]:
    data_by_acc = {}
    manifests = []
    for acc in accessions:
        data = fetch_uniprot_domain(acc)
        if data:
            data_by_acc[acc] = data
            cache = DOMAIN_CACHE / f"{acc}.json"
            manifests.append({"source": "UniProtKB JSON API; domain feature/xref records",
                              "accession": acc,
                              "record_version": data.get("entryAudit", {}).get("entryVersion", "NA"),
                              "sequence_version": data.get("entryAudit", {}).get("sequenceVersion", "NA"),
                              "retrieved_utc": data.get("_phase1b0_retrieved_utc", "cached"),
                              "cache_path": str(cache.relative_to(ROOT)) if cache.exists() else "NA",
                              "SHA256": sha256(cache) if cache.exists() else "NA",
                              "status": "OK" if not data.get("_fetch_error") else "FAILED"})
    rows = []
    for r in variants.itertuples(index=False):
        acc = str(r.canonical_UniProt_accession)
        try: pos = int(r.canonical_residue)
        except Exception: pos = 0
        data = data_by_acc.get(acc, {})
        features = []
        for f in data.get("features", []) if isinstance(data, dict) else []:
            kind = str(f.get("type", ""))
            if kind not in {"Domain", "Repeat", "Zinc finger", "Region"}:
                continue
            loc = f.get("location", {})
            try:
                start = int(loc.get("start", {}).get("value")); end = int(loc.get("end", {}).get("value"))
            except Exception:
                continue
            if start <= pos <= end:
                label = f.get("description") or f.get("featureId") or kind
                features.append((start, end, kind, str(label), str(f.get("featureId", ""))))
        features.sort(key=lambda x: (x[1] - x[0], x[0], x[2], x[3]))
        if features:
            start, end, kind, label, fid = features[0]
            inside = (pos - start) / max(1, end - start)
            boundary = min(pos - start, end - pos)
            all_domains = ";".join(f"{x[2]}:{x[3]}:{x[0]}-{x[1]}" for x in features)
            status = "ANNOTATED_DOMAIN_FEATURE"
        elif data.get("features"):
            start = end = boundary = "NA"; inside = "NA"; label = fid = all_domains = "NA"
            status = "NO_DOMAIN_FEATURE_AT_POSITION"
        else:
            start = end = boundary = inside = label = fid = all_domains = "NA"
            status = "UNIPROTKB_ANNOTATION_UNAVAILABLE"
        rows.append({"gene": r.gene, "VariationID": r.VariationID, "protein_change": r.protein_change,
                     "canonical_UniProt_accession": acc, "canonical_residue": pos,
                     "domain_identity": label, "domain_feature_id": fid,
                     "domain_start": start, "domain_end": end,
                     "relative_position_within_domain": inside,
                     "distance_to_nearest_domain_boundary_aa": boundary,
                     "overlapping_domain_features": all_domains,
                     "annotation_source": "UniProtKB current JSON feature table; frozen local response",
                     "annotation_status": status})
    return pd.DataFrame(rows), manifests


def build_alpha_missense_baseline(variants: pd.DataFrame) -> pd.DataFrame:
    key_to_variants = defaultdict(list)
    for r in variants.itertuples(index=False):
        var = canonical_variant(pd.Series(r._asdict()))
        if var:
            key_to_variants[(str(r.canonical_UniProt_accession), var)].append(r)
    matched = {}
    status = "SOURCE_FILE_UNAVAILABLE"
    if AM_FILE.exists():
        targets = {acc for acc, _ in key_to_variants}
        try:
            interim = ROOT / "data/interim"
            interim.mkdir(parents=True, exist_ok=True)
            key_path = interim / "phase1B0_alphamissense_variant_keys.tsv"
            subset_path = interim / "phase1B0_alphamissense_target_subset.tsv"
            subset_manifest = META / "phase1b0_alphamissense_subset_manifest.json"
            source_md5 = hashlib.md5()
            with AM_FILE.open("rb") as fh:
                for chunk in iter(lambda: fh.read(4 * 1024 * 1024), b""):
                    source_md5.update(chunk)
            wanted_keys = sorted(key_to_variants)
            valid_subset = False
            if subset_path.exists() and subset_manifest.exists():
                try:
                    prior = json.loads(subset_manifest.read_text(encoding="utf-8"))
                    with subset_path.open("rt", encoding="utf-8") as prior_fh:
                        subset_header = prior_fh.readline().rstrip("\r\n")
                    valid_subset = (prior.get("source_md5") == source_md5.hexdigest()
                                    and prior.get("target_key_n") == len(wanted_keys)
                                    and prior.get("subset_schema_version") == 2
                                    and subset_header == "uniprot_id\tprotein_variant\tam_pathogenicity\tam_class"
                                    and prior.get("subset_sha256") == sha256(subset_path))
                except Exception:
                    valid_subset = False
            if not valid_subset:
                key_path.write_text("".join(f"{a}\t{v}\n" for a,v in wanted_keys), encoding="utf-8")
                awk_program = 'FNR==NR { wanted[$1 SUBSEP $2]=1; next } ($1=="uniprot_id" && $2=="protein_variant") { print; next } (($1 SUBSEP $2) in wanted) { print }'
                subset_path.unlink(missing_ok=True)
                with subset_path.open("wb") as out_f:
                    gz = subprocess.Popen(["gzip", "-dc", str(AM_FILE)], stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE)
                    awk = subprocess.Popen(["awk", "-F", "\t", awk_program, str(key_path), "-"],
                                           stdin=gz.stdout, stdout=out_f, stderr=subprocess.PIPE)
                    assert gz.stdout is not None
                    gz.stdout.close()
                    awk_err = awk.stderr.read() if awk.stderr else b""
                    awk_rc = awk.wait()
                    gz_err = gz.stderr.read() if gz.stderr else b""
                    gz_rc = gz.wait()
                if awk_rc or gz_rc:
                    raise RuntimeError(f"AlphaMissense streaming filter failed: awk={awk_rc} {awk_err[:300]!r}; gzip={gz_rc} {gz_err[:300]!r}")
                subset_manifest.write_text(json.dumps({"source_path":str(AM_FILE.relative_to(ROOT)),
                    "source_md5":source_md5.hexdigest(), "target_key_n":len(wanted_keys),
                    "subset_schema_version":2,
                    "subset_path":str(subset_path.relative_to(ROOT)),
                    "subset_sha256":sha256(subset_path), "filtered_utc":utcnow()},indent=2)+"\n",encoding="utf-8")
            with subset_path.open("rt", encoding="utf-8", newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    key = (row.get("uniprot_id", ""), row.get("protein_variant", ""))
                    if key in key_to_variants:
                        matched[key] = (row.get("am_pathogenicity", "NA"), row.get("am_class", "NA"))
            status = "MATCHED_OFFICIAL_ALPHA_MISSENSE_2023_AA_SUBSTITUTIONS"
        except Exception as exc:
            status = f"SOURCE_READ_FAILED:{type(exc).__name__}:{exc}"
    rows = []
    for (acc, var), vv in key_to_variants.items():
        score, cls = matched.get((acc, var), ("NA", "NA"))
        for r in vv:
            rows.append({"gene": r.gene, "VariationID": r.VariationID,
                         "protein_change": r.protein_change,
                         "canonical_UniProt_accession": acc, "canonical_variant": var,
                         "AlphaMissense_pathogenicity": score, "AlphaMissense_class": cls,
                         "model_version": "AlphaMissense 2023 amino-acid substitution resource",
                         "source_DOI": "10.5281/zenodo.10813168",
                         "source_status": status if (acc, var) in matched else ("VARIANT_NOT_IN_ALPHA_MISSENSE_FILE" if status.startswith("MATCHED") else status)})
    return pd.DataFrame(rows)


def build_generic3d(variants: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    """Read AFDB monomers and map per-residue confidence/packing/SASA features."""
    rows = []
    source_meta = []
    try:
        from Bio.PDB import MMCIFParser, ShrakeRupley, is_aa
        from scipy.spatial import cKDTree
    except Exception as exc:
        for r in variants.itertuples(index=False):
            rows.append({"gene": r.gene, "VariationID": r.VariationID,
                         "canonical_UniProt_accession": r.canonical_UniProt_accession,
                         "canonical_residue": r.canonical_residue,
                         "AFDB_plddt": "NA", "monomer_SASA_A2": "NA",
                         "intrachain_CA_neighbor_n_8A": "NA",
                         "generic3d_status": f"DEPENDENCY_UNAVAILABLE:{exc}"})
        return pd.DataFrame(rows), source_meta
    try:
        helper = load_phase1a_helpers()
    except Exception:
        helper = None
    afdb = {}
    if AFDB_META.exists():
        try:
            afdb = {x["accession"]: x for x in json.loads(AFDB_META.read_text()).get("model_results", [])}
        except Exception:
            afdb = {}
    variant_by_acc = defaultdict(list)
    for r in variants.itertuples(index=False):
        variant_by_acc[str(r.canonical_UniProt_accession)].append(r)
    feature_by_acc = {}
    for acc, items in variant_by_acc.items():
        rec = afdb.get(acc, {})
        path = ROOT / rec.get("path", "") if rec.get("path") else None
        if not path or not path.exists():
            for r in items:
                feature_by_acc[(acc, int(r.canonical_residue))] = {"status": "AFDB_MONOMER_UNAVAILABLE"}
            continue
        try:
            structure = MMCIFParser(QUIET=True).get_structure(acc, str(path))
            model = structure[0]
            chain = next((c for c in model if any(is_aa(x, standard=False) for x in c)), None)
            if chain is None:
                raise RuntimeError("no protein chain")
            seq = "".join(__import__("Bio.SeqUtils", fromlist=["seq1"]).seq1(res.get_resname(), custom_map={"UNK":"X"})
                          for res in chain if is_aa(res, standard=False))
            mapping = {}
            if helper:
                uni_path = ROOT / f"data/interim/phase1A_uniprot_sequences/{acc}.fasta"
                if uni_path.exists():
                    ref = "".join(uni_path.read_text().splitlines()[1:])
                    mapping = helper.sequence_map(seq, ref).get("mapping", {})
            if not mapping:
                mapping = {i: i for i in range(1, len(seq) + 1)}
            import copy
            chain_for_sasa = copy.deepcopy(chain)
            ShrakeRupley(n_points=100).compute(chain_for_sasa, level="R")
            sasa_by_id = {res.id: float(getattr(res, "sasa", np.nan)) for res in chain_for_sasa}
            ca = []
            records = []
            idx = 0
            for res in chain:
                if not is_aa(res, standard=False):
                    continue
                idx += 1
                canonical = mapping.get(idx)
                if canonical is None:
                    continue
                atom = res["CA"] if "CA" in res else None
                if atom is None:
                    records.append((int(canonical), np.nan, sasa_by_id.get(res.id, np.nan), None))
                    continue
                coord = np.asarray(atom.coord, dtype=float)
                ca.append((int(canonical), coord))
                records.append((int(canonical), float(atom.get_bfactor()), sasa_by_id.get(res.id, np.nan), coord))
            coords = np.asarray([x[1] for x in ca], dtype=float) if ca else np.empty((0, 3))
            dens = {}
            if len(coords):
                pairs = cKDTree(coords).query_pairs(8.0)
                dens_counter = Counter()
                for i, j in pairs:
                    ci, cj = ca[i][0], ca[j][0]
                    if abs(ci - cj) > 3:
                        dens_counter[ci] += 1; dens_counter[cj] += 1
                dens = dict(dens_counter)
            for pos, plddt, sasa, coord in records:
                feature_by_acc[(acc, pos)] = {"status": "OK", "plddt": plddt,
                    "sasa": sasa, "ca_density": dens.get(pos, 0),
                    "structure_version": rec.get("version", "NA"),
                    "structure_path": rec.get("path", "NA"),
                    "structure_sha256": rec.get("sha256", "NA")}
            source_meta.append({"accession": acc, "model_version": rec.get("version", "NA"),
                                "source_url": rec.get("url", "NA"), "path": rec.get("path", "NA"),
                                "SHA256": rec.get("sha256", "NA"), "status": "PARSED"})
        except Exception as exc:
            source_meta.append({"accession": acc, "model_version": rec.get("version", "NA"),
                                "path": rec.get("path", "NA"), "SHA256": rec.get("sha256", "NA"),
                                "status": f"PARSE_FAILED:{type(exc).__name__}:{exc}"})
            for r in items:
                feature_by_acc[(acc, int(r.canonical_residue))] = {"status": "AFDB_PARSE_FAILED"}
    for r in variants.itertuples(index=False):
        acc = str(r.canonical_UniProt_accession)
        try: pos = int(r.canonical_residue)
        except Exception: pos = 0
        d = feature_by_acc.get((acc, pos), {"status": "RESIDUE_NOT_IN_AFDB_MODEL"})
        rows.append({"gene": r.gene, "VariationID": r.VariationID,
                     "canonical_UniProt_accession": acc, "canonical_residue": pos,
                     "AFDB_model_version": d.get("structure_version", afdb.get(acc, {}).get("version", "NA")),
                     "AFDB_plddt": d.get("plddt", "NA"),
                     "monomer_SASA_A2": d.get("sasa", "NA"),
                     "intrachain_CA_neighbor_n_8A": d.get("ca_density", "NA"),
                     "partner_identity_included": "NO",
                     "structure_sha256": d.get("structure_sha256", afdb.get(acc, {}).get("sha256", "NA")),
                     "generic3d_status": d.get("status", "RESIDUE_NOT_IN_AFDB_MODEL")})
    return pd.DataFrame(rows), source_meta


def get_model_and_mutpred_score_rows(universe: pd.DataFrame) -> pd.DataFrame:
    keep = ["gene", "VariationID", "partner"]
    if RAW_SCORE.exists():
        raw = read_tsv(RAW_SCORE)
        if not set(keep).issubset(raw.columns):
            raise RuntimeError(f"Raw score table lacks keys {keep}")
        cols = keep + [c for c in ["raw_score", "calibrated_probability", "classification",
                                   "model_version", "score_status", "failure_reason",
                                   "selected_structure_ID", "selected_target_chain", "selected_partner_chain"]
                       if c in raw.columns]
        return universe[keep].merge(raw[cols].drop_duplicates(keep), on=keep, how="left", validate="one_to_one")
    return universe[keep].assign(raw_score=np.nan, calibrated_probability=np.nan,
        classification="NA", model_version=f"MutPred-PPI current@{PRIMARY_MODEL_COMMIT}",
        score_status="NOT_RUN_TECHNICAL_HOLD",
        failure_reason="No legitimate score generated: ProtT5 encoder weights were not retrieved from the official model host; no proxy imputation.")


def build_all() -> dict:
    before = freeze_inputs()
    version = freeze_mutpred_version()
    cov = read_tsv(P1A / "coverage_by_gene.tsv")
    structural_genes = cov.loc[(cov.primary_cohort_eligible == "YES") &
                               (cov.passes_ge70_representation_threshold == "YES")].copy()
    genes = set(structural_genes.gene.astype(str))
    strict_counts = read_tsv(P0R / "gene_disease_variant_counts_exact_rcv.tsv")
    strict_counts = strict_counts[strict_counts.strict_disease_ge10 == "YES"].copy()
    strict_counts = strict_counts[strict_counts.gene.isin(genes)].copy()
    mapping = read_tsv(P1A / "protein_mapping.tsv")
    mapping["VariationID"] = mapping["VariationID"].astype(str)
    mapping42 = mapping[(mapping.gene.isin(genes)) &
                        (mapping.primary_cohort_eligible == "YES") &
                        mapping.mapping_status.isin(["EXACT_CANONICAL", "SAFE_ISOFORM_REMAP"])].copy()
    mapping42 = mapping42.drop_duplicates(["gene", "VariationID"])
    membership42 = mapping[(mapping.gene.isin(genes)) &
                          (mapping.primary_cohort_eligible == "YES")].copy()
    membership42 = membership42.drop_duplicates(["gene", "VariationID"])
    matrix_columns = ["gene", "VariationID", "protein_change", "MONDO_groups",
                      "canonical_UniProt_accession", "canonical_residue", "mapping_status",
                      "partner", "tier", "interaction_sources", "structure_source",
                      "interface_status", "interface_min_distance_A", "interface_confidence",
                      "experimental_PDB_IDs", "distance_cutoff_A"]
    matrix = read_tsv(P1A / "partner_specific_matrix_structural.tsv", usecols=matrix_columns)
    matrix["VariationID"] = matrix["VariationID"].astype(str)
    matrix["canonical_residue"] = pd.to_numeric(matrix["canonical_residue"], errors="coerce")
    ref_alt = mapping[["gene", "VariationID", "canonical_reference_aa", "alternate_aa"]].copy()
    ref_alt["VariationID"] = ref_alt["VariationID"].astype(str)
    matrix = matrix.merge(ref_alt.drop_duplicates(["gene", "VariationID"]),
                          on=["gene", "VariationID"], how="left", validate="many_to_one")
    # Exact, strictly grouped disease membership and denominator from Phase 0R.
    strict_group_map = defaultdict(set)
    disease_label = {}
    disease_n = {}
    for r in strict_counts.itertuples(index=False):
        key = (str(r.gene), str(r.MONDO))
        strict_group_map[str(r.gene)].add(str(r.MONDO))
        disease_label[key] = str(r.MONDO_label)
        disease_n[key] = int(r.n_unique_variants)
    vgroups = []
    for r in membership42.itertuples(index=False):
        for mondo in split_semicolon(r.MONDOs):
            if mondo in strict_group_map.get(str(r.gene), set()):
                vgroups.append({"gene": str(r.gene), "VariationID": str(r.VariationID),
                                "MONDO": mondo, "protein_change": r.protein_change,
                                "canonical_UniProt_accession": r.canonical_UniProt_accession,
                                "canonical_residue": r.canonical_residue,
                                "mapping_status": r.mapping_status})
    vgroups = pd.DataFrame(vgroups).drop_duplicates(["gene", "VariationID", "MONDO"])
    label_n = vgroups.groupby(["gene", "VariationID"]).MONDO.nunique()
    multi_keys = {key for key, n in label_n.items() if n > 1}
    multi_by_group = (vgroups[vgroups.apply(lambda x: (x.gene, x.VariationID) in multi_keys, axis=1)]
                      .groupby(["gene", "MONDO"]).VariationID.nunique().to_dict())
    retained_classes = {}
    for gene in sorted(genes):
        n_classes = 0
        for r in strict_counts[strict_counts.gene.eq(gene)].itertuples(index=False):
            n_after = int(r.n_unique_variants) - int(multi_by_group.get((gene, str(r.MONDO)), 0))
            n_classes += int(n_after >= 10)
        retained_classes[gene] = n_classes
    variants = mapping42.copy()
    variants["_canon_variant"] = variants.apply(canonical_variant, axis=1)
    variants = variants[variants["_canon_variant"].ne("")].copy()
    # Structures available per variant, keeping occupancy separate from model availability.
    nonmissing_structure = matrix["structure_source"].fillna("").astype(str).ne("")
    t12 = matrix[matrix.tier.isin(["Tier 1", "Tier 2"])].copy()
    t12["_has_model"] = t12.structure_source.fillna("").astype(str).ne("")
    t12_dims = (t12[t12["_has_model"]].groupby(["gene", "VariationID"]).partner.nunique()
                .rename("structural_partner_dimensions_n").reset_index())
    t12_iface = (t12[t12.interface_status.eq("INTERFACE")].groupby(["gene", "VariationID"]).partner.nunique()
                 .rename("interface_occupancy_partner_n").reset_index())
    t12_avail = (t12.groupby("gene").partner.nunique().rename("Tier1_2_partner_n_observed_in_matrix").reset_index())

    # 1) exact 42-gene structural table.
    previous_diff = read_tsv(P1A / "differential_missingness.tsv")
    prev_by_gene = previous_diff.groupby("gene").apply(
        lambda d: ";".join(sorted(set(d.loc[d.flag == "DIFFERENTIAL_COVERAGE", "coverage_metric"].astype(str))))
    ).to_dict()
    rows42 = []
    for r in structural_genes.sort_values("gene").itertuples(index=False):
        gd = strict_counts[strict_counts.gene == r.gene].sort_values("MONDO")
        labels = [f"{x.MONDO} ({x.MONDO_label})" for x in gd.itertuples(index=False)]
        counts = [f"{x.MONDO}={int(x.n_unique_variants)}" for x in gd.itertuples(index=False)]
        rows42.append({"gene": r.gene, "disease_groups": ";".join(labels),
            "variants_per_disease": ";".join(counts),
            "total_strict_variants": int(r.unique_strict_group_variants_n),
            "protein_mappable_variants": int(r.protein_mappable_variants_n),
            "credible_partners": int(r.credible_partner_n), "Tier1_partners": int(r.Tier1_partner_n),
            "Tier2_partners": int(r.Tier2_partner_n), "Tier1_2_partners": int(r.Tier1_2_partner_n),
            "interface_model_coverage_ge1_pct": float(r.variants_with_ge1_partner_interface_pct),
            "interface_model_coverage_ge3_Tier1_2_pct": float(r.variants_with_ge3_Tier1_2_partner_interfaces_pct),
            "structural_addressability_pct": float(r.fraction_structurally_addressable_pct),
            "PDB_residue_coverage_pct": float(r.fraction_MutPred_PPI_structurally_addressable_pct),
            "Phase1A_differential_coverage_flag": prev_by_gene.get(r.gene, "NO_DIFFERENTIAL_COVERAGE") or "NO_DIFFERENTIAL_COVERAGE",
            "interface_overlap_category": r.interface_overlap_class,
            "mechanism_status": r.mechanism_stratum,
            "Phase1A_structural_gate": "PASS"})
    write_tsv(pd.DataFrame(rows42), OUT / "structural_42_genes.tsv")

    # 2) Reconstruct and freeze the exact Phase 1A score universe.
    universe = matrix[matrix.structure_source.eq("experimental PDB coordinates")].copy()
    universe["canonical_variant"] = universe.apply(canonical_variant, axis=1)
    universe["disease_group"] = universe.MONDO_groups
    universe["structure_source"] = universe.structure_source.fillna("NA")
    universe["partner_accession"] = "NA"
    universe["partner_accession_mapping_status"] = "NOT_RESOLVED"
    pa, pas = get_partner_accessions(universe, mapping)
    for idx, r in universe[["gene", "partner"]].drop_duplicates().iterrows():
        key = (str(r.gene), str(r.partner))
        if key in pa:
            mask = universe.gene.eq(key[0]) & universe.partner.eq(key[1])
            universe.loc[mask, "partner_accession"] = pa[key]
            universe.loc[mask, "partner_accession_mapping_status"] = pas.get(key, "NA")
        elif key in pas:
            universe.loc[universe.gene.eq(key[0]) & universe.partner.eq(key[1]),
                         "partner_accession_mapping_status"] = pas[key]
    universe = universe.rename(columns={"VariationID": "VariationID"})
    if len(universe) != 11150:
        raise RuntimeError(f"Phase 1A eligible universe changed: expected 11150, got {len(universe)}")
    score_keys = ["gene", "VariationID", "partner"]
    if universe.duplicated(score_keys).any():
        raise RuntimeError("Scoring universe contains duplicate gene/VariationID/partner keys")
    # Include all user-requested fields plus frozen provenance and mapper status.
    universe = universe[["gene", "disease_group", "VariationID", "protein_change", "canonical_variant",
        "canonical_UniProt_accession", "canonical_residue", "partner", "partner_accession",
        "partner_accession_mapping_status", "interaction_sources", "tier", "structure_source",
        "interface_status", "interface_min_distance_A", "interface_confidence",
        "experimental_PDB_IDs", "mapping_status", "distance_cutoff_A"]].copy()
    scores = get_model_and_mutpred_score_rows(universe)
    universe = universe.merge(scores, on=score_keys, how="left", validate="one_to_one")
    universe["partner_score_availability"] = np.where(universe.raw_score.notna(), "SCORED", "NA_SCORE")
    write_tsv(universe, OUT / "variant_partner_scoring_universe.tsv")
    score_cols = score_keys + ["disease_group", "protein_change", "canonical_variant",
        "canonical_UniProt_accession", "partner_accession", "tier", "structure_source",
        "experimental_PDB_IDs", "raw_score", "calibrated_probability", "classification",
        "model_version", "score_status", "failure_reason"]
    score_out = universe[[c for c in score_cols if c in universe.columns]].copy()
    write_tsv(score_out, OUT / "mutpred_ppi_scores.tsv")
    write_tsv(score_out, OUT / "partner_disruption_matrix.tsv")
    # Structural S1 and S2, separate from E.
    s = matrix[matrix.gene.isin(genes) & matrix.structure_source.fillna("").astype(str).ne("")].copy()
    s = s[["gene", "VariationID", "protein_change", "canonical_UniProt_accession",
           "canonical_residue", "partner", "tier", "structure_source",
           "experimental_PDB_IDs", "interface_status", "interface_min_distance_A",
           "interface_confidence", "distance_cutoff_A"]].copy()
    s["S1_binary_interface_occupancy"] = np.where(s.interface_status.eq("INTERFACE"), 1,
                                                    np.where(s.interface_status.isin(["NON_INTERFACE_IN_OBSERVED_COMPLEX",
                                                                                     "NON_INTERFACE_IN_PIONEER_INTERFACE_MODEL"]), 0, "NA"))
    s["S2_distance_to_interface_A"] = s.interface_min_distance_A
    s["S1_annotation_source"] = s.structure_source
    s["S2_distance_source"] = np.where(s.structure_source.eq("experimental PDB coordinates"),
                                        "Phase1A experimental complex coordinates; 5.0A heavy-atom cutoff",
                                        "NA_NO_COORDINATE_DISTANCE_IN_PIONEER_ANNOTATION")
    write_tsv(s, OUT / "structural_interface_matrix.tsv")

    # 3) Training overlap at each score observation.
    train_variant, train_pair, train_protein, train_n = public_training_sets()
    phase1a_overlap = read_tsv(P1A / "mutpred_training_overlap.tsv")
    phase1a_overlap["VariationID"] = phase1a_overlap.VariationID.astype(str)
    vflags = phase1a_overlap.set_index(["gene", "VariationID"])
    overlaps = []
    for r in universe.itertuples(index=False):
        var = str(r.canonical_variant)
        target_acc = str(r.canonical_UniProt_accession)
        partner_acc = str(r.partner_accession)
        exact_var = (target_acc, var) in train_variant
        pair_known = partner_acc not in {"NA", ""} and ";" not in partner_acc
        exact_pair = (target_acc, partner_acc, var) in train_pair if pair_known else False
        rowf = vflags.loc[(str(r.gene), str(r.VariationID))] if (str(r.gene), str(r.VariationID)) in vflags.index else None
        protein_hit = target_acc in train_protein
        overlaps.append({"gene": r.gene, "VariationID": r.VariationID, "partner": r.partner,
            "canonical_variant": var, "canonical_UniProt_accession": target_acc,
            "partner_accession": partner_acc,
            "exact_variant_training_overlap": "YES" if exact_var else "NO_KNOWN_OVERLAP",
            "exact_variant_partner_training_overlap": "YES" if exact_pair else ("NO" if pair_known else "UNKNOWN_PARTNER_ACCESSION"),
            "target_protein_training_overlap": "YES" if protein_hit else "NO_KNOWN_OVERLAP",
            "partner_accession_mapping_status": r.partner_accession_mapping_status,
            "training_audit_scope": "public training/evaluation rows (n=12529); unpublished VarChAMP not auditable",
            "strict_external_include": "NO" if (exact_var or exact_pair) else ("YES" if pair_known else "PROVISIONAL_VARIANT_ONLY"),
            "protein_external_include": "NO" if protein_hit else "YES"})
    overlap_df = pd.DataFrame(overlaps)
    write_tsv(overlap_df, OUT / "mutpred_training_overlap.tsv")

    # Coverage audit at both row and disease-variant levels.
    disease_score_rows = defaultdict(lambda: {"eligible_rows": 0, "scored_rows": 0, "variant_partners": defaultdict(set)})
    strict_group_keys = {(str(r.gene), str(r.MONDO)) for r in strict_counts.itertuples(index=False)}
    for r in universe.itertuples(index=False):
        for mondo in split_semicolon(r.disease_group):
            key = (str(r.gene), mondo)
            if key not in strict_group_keys:
                continue
            z = disease_score_rows[key]
            z["eligible_rows"] += 1
            if pd.notna(r.raw_score):
                z["scored_rows"] += 1
                if str(r.tier) in {"Tier 1", "Tier 2"}:
                    z["variant_partners"][str(r.VariationID)].add(str(r.partner))

    # 4) Differential coverage, per strict disease group, data-blind thresholds.
    variant_dims = t12_dims
    vg = vgroups.merge(variant_dims, on=["gene", "VariationID"], how="left")
    vg["structural_partner_dimensions_n"] = vg.structural_partner_dimensions_n.fillna(0).astype(int)
    vg = vg.merge(t12_iface, on=["gene", "VariationID"], how="left")
    vg["interface_occupancy_partner_n"] = vg.interface_occupancy_partner_n.fillna(0).astype(int)
    drows = []
    gene_diff = {}
    for gene in sorted(genes):
        group_rows = []
        for x in strict_counts[strict_counts.gene.eq(gene)].itertuples(index=False):
            mondo = str(x.MONDO); total = int(x.n_unique_variants)
            subset = vg[(vg.gene == gene) & (vg.MONDO == mondo)]
            represented = int((subset.structural_partner_dimensions_n >= 1).sum())
            ge3 = int((subset.structural_partner_dimensions_n >= 3).sum())
            occ = int((subset.interface_occupancy_partner_n >= 1).sum())
            mean_dims = float(subset.structural_partner_dimensions_n.mean()) if len(subset) else 0.0
            rate = 100.0 * ge3 / total if total else 0.0
            group_rows.append((mondo, total, represented, ge3, rate))
            drows.append({"gene": gene, "MONDO": mondo, "disease_label": str(x.MONDO_label),
                "total_variants": total, "represented_variants_ge1_Tier1_2_model": represented,
                "represented_variants_ge1_pct": round(100 * represented / total, 2) if total else 0,
                "variants_with_ge3_Tier1_2_partner_dimensions": ge3,
                "Tier1_2_ge3_representation_pct": round(rate, 2),
                "variants_with_any_observed_interface_occupancy": occ,
                "mean_Tier1_2_partner_dimensions_per_variant": round(mean_dims, 3),
                "available_Tier1_2_partner_dimensions_in_gene": int(t12[t12.gene.eq(gene)].loc[t12[t12.gene.eq(gene)]["_has_model"], "partner"].nunique()),
                "eligible_variant_partner_score_rows": disease_score_rows[(gene, mondo)]["eligible_rows"],
                "successfully_scored_variant_partner_rows": disease_score_rows[(gene, mondo)]["scored_rows"],
                "disruption_score_row_coverage_pct": round(100*disease_score_rows[(gene, mondo)]["scored_rows"] /
                    max(1, disease_score_rows[(gene, mondo)]["eligible_rows"]), 2),
                "variants_with_ge1_disruption_score_n": sum(len(v)>0 for v in disease_score_rows[(gene, mondo)]["variant_partners"].values()),
                "variants_with_ge3_disruption_scores_n": sum(len(v)>=3 for v in disease_score_rows[(gene, mondo)]["variant_partners"].values()),
                "variants_with_ge3_disruption_scores_pct": round(100*sum(len(v)>=3 for v in disease_score_rows[(gene, mondo)]["variant_partners"].values())/total,2) if total else 0,
                "Phase1A_flagged_gene": "YES" if prev_by_gene.get(gene) else "NO"})
        rates = [x[4] for x in group_rows]
        if not rates:
            cls = "SEVERE_IMBALANCE"
        elif min(rates) < 50.0 or max(rates) - min(rates) >= 50.0:
            cls = "SEVERE_IMBALANCE"
        elif (min(rates) < 70.0 and max(rates) >= 70.0) or max(rates) - min(rates) >= 20.0:
            cls = "MATERIAL_IMBALANCE"
        else:
            cls = "MINOR_IMBALANCE"
        gene_diff[gene] = cls
    ddf = pd.DataFrame(drows)
    ddf["gene_coverage_adjudication"] = ddf.gene.map(gene_diff)
    ddf["frozen_threshold"] = "Based on disease-specific >=3 Tier1-2 partner model representation: severe if any group <50% or range >=50pp; material if range >=20pp or a group <70% while another >=70%; minor otherwise. Uses no outcome-separation statistic."
    write_tsv(ddf, OUT / "differential_coverage_adjudication.tsv")

    # 5) Cohorts A/B/C, using only pre-outcome technical completeness.
    score_stats = {}
    successful = universe[universe.raw_score.notna()].copy()
    if len(successful):
        score_counts = (successful[successful.tier.isin(["Tier 1", "Tier 2"])]
                        .groupby(["gene", "VariationID"]).partner.nunique())
    else:
        score_counts = pd.Series(dtype=int)
    score_rows_by_gene = universe.groupby("gene").agg(
        eligible_score_rows=("VariationID", "size"),
        successfully_scored_rows=("raw_score", "count"),
        scored_partners=("partner", "nunique")).to_dict("index")
    scored_t12_partners = (successful[successful.tier.isin(["Tier 1", "Tier 2"])]
                           .groupby("gene").partner.nunique().to_dict()) if len(successful) else {}
    structural_t12_partners = (t12[t12["_has_model"]].groupby("gene").partner.nunique().to_dict())
    multi_keys = set()
    for x in vgroups.groupby(["gene", "VariationID"]).MONDO.nunique().items():
        if x[1] > 1: multi_keys.add(x[0])
    cohort_rows = []
    for r in structural_genes.sort_values("gene").itertuples(index=False):
        gene = str(r.gene)
        total = int(r.unique_strict_group_variants_n)
        gene_variants = {v for g,v in vgroups.loc[vgroups.gene.eq(gene), ["gene", "VariationID"]].itertuples(index=False, name=None)
                         if (g,v) not in multi_keys}
        n_primary_variants = max(0, total - sum(1 for g,v in multi_keys if g == gene))
        scored_var_ge1 = sum(score_counts.get((gene, v), 0) >= 1 for v in gene_variants)
        scored_var_ge3 = sum(score_counts.get((gene, v), 0) >= 3 for v in gene_variants)
        ge1_pct = 100 * scored_var_ge1 / total if total else 0
        ge3_pct = 100 * scored_var_ge3 / n_primary_variants if n_primary_variants else 0
        structural_dims_by_variant = t12_dims.set_index(["gene", "VariationID"])[
            "structural_partner_dimensions_n"].to_dict()
        structural_var_ge3 = sum(structural_dims_by_variant.get((gene, v), 0) >= 3 for v in gene_variants)
        structural_ge3_pct = 100 * structural_var_ge3 / n_primary_variants if n_primary_variants else 0
        n_score_partners = int(scored_t12_partners.get(gene, 0))
        n_structural_partners = int(structural_t12_partners.get(gene, 0))
        eligible = True
        reasons = []
        if r.primary_cohort_eligible != "YES" or r.passes_ge70_representation_threshold != "YES":
            eligible = False; reasons.append("PHASE1A_STRUCTURAL_GATE_FAIL")
        if str(r.mechanism_stratum) == "CONFIRMED_INCOMPATIBLE":
            eligible = False; reasons.append("CONFIRMED_INCOMPATIBLE")
        if gene_diff.get(gene) == "SEVERE_IMBALANCE":
            eligible = False; reasons.append("SEVERE_DIFFERENTIAL_COVERAGE")
        if gene_diff.get(gene) == "MATERIAL_IMBALANCE":
            eligible = False; reasons.append("MATERIAL_IMBALANCE_SENSITIVITY_ONLY")
        if retained_classes.get(gene, 0) < 2:
            eligible = False; reasons.append("<2_DISEASE_GROUPS_REMAIN_AFTER_MULTILABEL_EXCLUSION")
        if n_score_partners < 3:
            eligible = False; reasons.append("<3_SCORED_PARTNER_DIMENSIONS")
        if ge3_pct < 70:
            eligible = False; reasons.append("<70_PERCENT_VARIANTS_WITH_GE3_SCORES")
        if len(successful) == 0:
            eligible = False; reasons.append("DISRUPTION_SCORING_NOT_STARTED")
        B = (r.primary_cohort_eligible == "YES" and r.passes_ge70_representation_threshold == "YES"
             and str(r.mechanism_stratum) != "CONFIRMED_INCOMPATIBLE"
             and gene_diff.get(gene) == "MINOR_IMBALANCE"
             and retained_classes.get(gene, 0) >= 2
             and n_structural_partners >= 3 and structural_ge3_pct >= 70)
        C_B = B and str(r.interface_overlap_class) == "SEPARABLE_INTERFACES"
        C_A = eligible and str(r.interface_overlap_class) == "SEPARABLE_INTERFACES"
        cohort_rows.append({"gene": gene, "Phase1A_structural_gate_PASS": "YES",
            "mechanism_status": r.mechanism_stratum,
            "differential_coverage_class": gene_diff.get(gene),
            "eligible_score_rows": int(score_rows_by_gene.get(gene, {}).get("eligible_score_rows", 0)),
            "successfully_scored_rows": int(score_rows_by_gene.get(gene, {}).get("successfully_scored_rows", 0)),
            "score_row_coverage_pct": round(100 * score_rows_by_gene.get(gene, {}).get("successfully_scored_rows", 0) /
                                             max(1, score_rows_by_gene.get(gene, {}).get("eligible_score_rows", 0)), 2),
            "variants_with_ge1_disruption_score_n": scored_var_ge1,
            "variants_with_ge3_disruption_scores_n": scored_var_ge3,
            "primary_unambiguous_variant_denominator_n": n_primary_variants,
            "variants_with_ge3_structural_partner_dimensions_n": structural_var_ge3,
            "variants_with_ge3_structural_partner_dimensions_pct_of_primary_variants": round(structural_ge3_pct, 2),
            "disease_groups_ge10_after_multilabel_exclusion_n": retained_classes.get(gene, 0),
            "variants_with_ge3_disruption_scores_pct_of_primary_variants": round(ge3_pct, 2),
            "variants_with_ge3_disruption_scores_pct_of_strict_variants": round(ge3_pct, 2),
            "structural_sensitivity_cohort_B": "YES" if B else "NO",
            "primary_disruption_cohort_A": "YES" if eligible else "NO",
            "separable_interface_cohort_C_within_A": "YES" if C_A else "NO",
            "separable_interface_cohort_C_within_B": "YES" if C_B else "NO",
            "primary_cohort_reason": "PASS" if eligible else ";".join(reasons)})
    cohorts = pd.DataFrame(cohort_rows)
    write_tsv(cohorts, OUT / "primary_cohort.tsv")

    # 6) Interface occupancy without partner IDs, hotspot and multilabel audits.
    all42 = matrix[matrix.gene.isin(genes)].copy()
    interface_count = (all42[all42.interface_status.eq("INTERFACE") & all42.tier.isin(["Tier 1", "Tier 2"])]
                       .groupby(["gene", "VariationID"]).partner.nunique().rename("interface_partner_count"))
    dimensions_count = (t12[t12["_has_model"] & t12.gene.isin(genes)]
                        .groupby(["gene", "VariationID"]).partner.nunique().rename("structural_partner_dimension_count"))
    iface_base = variants[["gene", "VariationID", "protein_change", "canonical_UniProt_accession", "canonical_residue"]].copy()
    iface_base = iface_base.merge(interface_count, on=["gene", "VariationID"], how="left")
    iface_base = iface_base.merge(dimensions_count, on=["gene", "VariationID"], how="left")
    iface_base[["interface_partner_count", "structural_partner_dimension_count"]] = iface_base[["interface_partner_count", "structural_partner_dimension_count"]].fillna(0).astype(int)
    iface_base["any_interface_occupancy"] = (iface_base.interface_partner_count > 0).astype(int)
    iface_base["partner_identity_included"] = "NO"
    write_tsv(iface_base, OUT / "baseline_interface_no_partner.tsv")

    hotspot_rows = []
    for gene in sorted(genes):
        gvars = vgroups[vgroups.gene.eq(gene)].drop_duplicates(["VariationID", "MONDO"])
        mut = variants[variants.gene.eq(gene)].drop_duplicates("VariationID")
        residue_counts = mut.groupby("canonical_residue").VariationID.nunique().sort_values(ascending=False)
        total = max(1, mut.VariationID.nunique())
        max_all = float(residue_counts.iloc[0] / total) if len(residue_counts) else 0.0
        top3_all = float(residue_counts.head(3).sum() / total) if len(residue_counts) else 0.0
        group_stats = []
        for mondo, grp in gvars.groupby("MONDO"):
            ids = set(grp.VariationID.astype(str))
            sub = mut[mut.VariationID.astype(str).isin(ids)]
            ct = sub.groupby("canonical_residue").VariationID.nunique().sort_values(ascending=False)
            denom = max(1, sub.VariationID.nunique())
            max_frac = float(ct.iloc[0] / denom) if len(ct) else 0.0
            top3_frac = float(ct.head(3).sum() / denom) if len(ct) else 0.0
            group_stats.append((mondo, max_frac, top3_frac))
        gene_dominated = any(m >= 0.5 or t >= 0.7 for _, m, t in group_stats)
        hotspot_rows.append({"gene": gene, "unique_substitutions_n": int(mut["_canon_variant"].nunique()),
            "unique_residues_n": int(mut.canonical_residue.nunique()),
            "maximum_single_residue_fraction_all_variants": round(max_all, 4),
            "top3_residue_fraction_all_variants": round(top3_all, 4),
            "maximum_single_residue_fraction_any_disease": round(max((x[1] for x in group_stats), default=0), 4),
            "maximum_top3_residue_fraction_any_disease": round(max((x[2] for x in group_stats), default=0), 4),
            "disease_group_hotspot_fractions": ";".join(f"{m}:max={a:.4f},top3={b:.4f}" for m,a,b in group_stats),
            "hotspot_status": "HOTSPOT_DOMINATED" if gene_dominated else "NOT_HOTSPOT_DOMINATED",
            "frozen_threshold": "HOTSPOT_DOMINATED if any disease group has one residue >=50% or top three residues >=70% of its unique variants"})
    hotspot = pd.DataFrame(hotspot_rows)
    write_tsv(hotspot, OUT / "hotspot_audit.tsv")

    ml_rows = []
    var_group_n = vgroups.groupby(["gene", "VariationID"]).MONDO.nunique()
    for r in membership42.itertuples(index=False):
        k=(str(r.gene),str(r.VariationID)); n=int(var_group_n.get(k,0))
        ml_rows.append({"gene":r.gene,"VariationID":r.VariationID,"protein_change":r.protein_change,
            "MONDO_groups":r.MONDOs,"disease_group_n":n,
            "multi_label_variant_flag":"MULTI_LABEL_VARIANT" if n>1 else "SINGLE_LABEL_VARIANT",
            "primary_handling":"EXCLUDE_FROM_PRIMARY; retain in dedicated multi-label sensitivity" if n>1 else "ELIGIBLE_FOR_PRIMARY"})
    ml = pd.DataFrame(ml_rows)
    write_tsv(ml, OUT / "multilabel_variant_audit.tsv")

    # 7) FoldX is not run; report only coordinate-based addressability.
    pdb42_ids = (universe[universe.gene.isin(genes)].groupby("gene").VariationID.nunique())
    foldx_rows = []
    for r in structural_genes.itertuples(index=False):
        n_total = int(r.unique_strict_group_variants_n)
        n_addr = int(pdb42_ids.get(r.gene, 0))
        pct = 100 * n_addr / n_total if n_total else 0
        foldx_rows.append({"gene":r.gene,"strict_variant_n":n_total,
            "variants_with_observed_PDB_target_coordinates_n":n_addr,
            "coordinate_addressability_pct":round(pct,2),
            "estimated_FoldX_ready_genes_ge70pct":"YES" if pct>=70 else "NO",
            "FoldX_run_status":"NOT_RUN_PHASE1B0",
            "feasibility_basis":"Observed experimental PDB target-residue coordinates; FoldX-specific format/repair compatibility not tested",
            "selected_structures":";".join(sorted({str(x) for x in universe.loc[universe.gene.eq(r.gene),"experimental_PDB_IDs"].dropna()}))})
    write_tsv(pd.DataFrame(foldx_rows), OUT / "foldx_feasibility.tsv")

    # 8) Baseline tables: position, frozen Pfam/InterPro-derived UniProt domains,
    # generic monomer 3D and AlphaMissense (same row key: unique gene/variant).
    position_rows=[]
    fasta_dir=ROOT/"data/interim/phase1A_uniprot_sequences"
    for r in variants.itertuples(index=False):
        acc=str(r.canonical_UniProt_accession); pos=int(r.canonical_residue)
        fasta=fasta_dir/f"{acc}.fasta"
        seq="".join(fasta.read_text().splitlines()[1:]) if fasta.exists() else ""
        length=len(seq) if seq else np.nan
        position_rows.append({"gene":r.gene,"VariationID":r.VariationID,"protein_change":r.protein_change,
            "canonical_UniProt_accession":acc,"canonical_residue":pos,"canonical_protein_length":length,
            "normalized_position_0_1":(pos-1)/(length-1) if length and length>1 else (0.0 if length==1 else "NA"),
            "mapping_status":r.mapping_status})
    position = pd.DataFrame(position_rows).drop_duplicates(["gene","VariationID"])
    write_tsv(position, OUT / "baseline_position.tsv")
    domain, domain_sources = build_domain_baseline(variants, sorted(set(variants.canonical_UniProt_accession.dropna().astype(str))))
    write_tsv(domain.drop_duplicates(["gene","VariationID"]), OUT / "baseline_domain.tsv")
    generic3d, afdb_sources = build_generic3d(variants)
    write_tsv(generic3d.drop_duplicates(["gene","VariationID"]), OUT / "baseline_generic3d.tsv")
    pathogenicity = build_alpha_missense_baseline(variants).drop_duplicates(["gene","VariationID"])
    write_tsv(pathogenicity, OUT / "baseline_pathogenicity.tsv")
    source_manifest=[]
    for x in domain_sources: source_manifest.append(x)
    for x in afdb_sources: source_manifest.append({"source":"AlphaFold DB monomer CIF",**x})
    if AM_FILE.exists():
        source_manifest.append({"source":"AlphaMissense 2023 aa substitutions", "accession":"10.5281/zenodo.10813168",
            "cache_path":str(AM_FILE.relative_to(ROOT)),"bytes":AM_FILE.stat().st_size,"SHA256":sha256(AM_FILE),"status":"LOCAL_FILE_PRESENT"})
    else:
        source_manifest.append({"source":"AlphaMissense 2023 aa substitutions", "accession":"10.5281/zenodo.10813168",
            "cache_path":str(AM_FILE.relative_to(ROOT)),"status":"DOWNLOAD_PENDING_OR_UNAVAILABLE"})
    write_tsv(pd.DataFrame(source_manifest), META / "PHASE1B0_BASELINE_SOURCE_MANIFEST.tsv")

    # 9) Pre-outcome exemplar re-ranking by technical availability only.
    disease_balance = {}
    for gene, grp in strict_counts.groupby("gene"):
        ns = grp.n_unique_variants.astype(int).tolist()
        disease_balance[gene] = (min(ns)/max(ns)) if ns and max(ns)>0 else 0
    candidate_rows=[]
    score_summary=cohorts.set_index("gene")
    train_g=overlap_df.groupby("gene").agg(
        exact_variant_overlap_rows=("exact_variant_training_overlap",lambda s: int((s=="YES").sum())),
        exact_pair_overlap_rows=("exact_variant_partner_training_overlap",lambda s: int((s=="YES").sum()))).to_dict("index")
    for r in structural_genes.itertuples(index=False):
        gene=str(r.gene)
        cr=score_summary.loc[gene]
        potential_rows=int(universe[universe.gene.eq(gene)].shape[0])
        potential_variants=int(universe[universe.gene.eq(gene)].VariationID.nunique())
        n_total=int(r.unique_strict_group_variants_n)
        candidate_rows.append({"gene":gene,"strict_disease_group_n":int(r.strict_disease_group_n),
            "strict_variant_n":n_total,"minimum_variants_per_disease":int(strict_counts[strict_counts.gene.eq(gene)].n_unique_variants.min()),
            "disease_balance_min_over_max":round(disease_balance.get(gene,0),4),
            "Tier1_2_partners":int(r.Tier1_2_partner_n),
            "interface_overlap_class":r.interface_overlap_class,
            "PDB_eligible_score_rows":potential_rows,"variants_with_any_PDB_score_row_n":potential_variants,
            "potential_ge3_PDB_partner_variant_pct":round(100*sum(universe[universe.gene.eq(gene)].groupby("VariationID").partner.nunique()>=3)/n_total,2) if n_total else 0,
            "score_success_rows":int(cr.successfully_scored_rows),"score_success_pct":float(cr.score_row_coverage_pct),
            "training_exact_variant_overlap_rows":train_g.get(gene,{}).get("exact_variant_overlap_rows",0),
            "training_exact_pair_overlap_rows":train_g.get(gene,{}).get("exact_pair_overlap_rows",0),
            "data_blind_selection_rule":"Sort by potential >=3-partner PDB score representation, score completeness, disease balance, Tier1-2 partner count, separable interface class; no disease-separation result"})
    candidate_df=pd.DataFrame(candidate_rows)
    candidate_df["_sep"] = (candidate_df.interface_overlap_class=="SEPARABLE_INTERFACES").astype(int)
    candidate_df=candidate_df.sort_values(["potential_ge3_PDB_partner_variant_pct","score_success_pct","disease_balance_min_over_max","Tier1_2_partners","_sep"],ascending=[False,False,False,False,False]).drop(columns="_sep").reset_index(drop=True)
    candidate_df.insert(0,"data_blind_rank",np.arange(1,len(candidate_df)+1))
    write_tsv(candidate_df.head(10), OUT / "data_blind_exemplar_candidates.tsv")

    # 10) freeze statistical/validation plan before any phenotype statistic.
    write_analysis_specification(cohorts, gene_diff, train_n)

    # Manifest and report summarize completeness, never run outcome tests.
    score_n = int(universe.raw_score.notna().sum())
    row_cov = 100 * score_n / max(1, len(universe))
    variants_ge3 = int((score_counts >= 3).sum()) if len(score_counts) else 0
    genes_ge70 = int((cohorts.variants_with_ge3_disruption_scores_pct_of_primary_variants >= 70).sum())
    A_genes = cohorts.loc[cohorts.primary_disruption_cohort_A.eq("YES"),"gene"].tolist()
    B_genes = cohorts.loc[cohorts.structural_sensitivity_cohort_B.eq("YES"),"gene"].tolist()
    C_genes = cohorts.loc[cohorts.separable_interface_cohort_C_within_A.eq("YES"),"gene"].tolist()
    diff_exclusions = sorted([g for g,c in gene_diff.items() if c in {"MATERIAL_IMBALANCE","SEVERE_IMBALANCE"}])
    manifest = {
        "phase":"PHASE1B0", "generated_utc":utcnow(), "phase1A_input_tsv_n":len(before),
        "phase1A_input_hashes_unchanged":True,"structural_42_gene_n":len(genes),
        "scoring_universe_rows":len(universe),"scoring_universe_gene_n":int(universe.gene.nunique()),
        "structural_42_genes_with_exact_PDB_score_rows_n":int(universe.loc[universe.gene.isin(genes),"gene"].nunique()),
        "successfully_scored_rows":score_n,
        "row_score_coverage_pct":round(row_cov,4),"variants_with_ge3_scores_n":variants_ge3,
        "genes_with_ge70_disruption_variants_n":genes_ge70,
        "cohort_A_genes":A_genes,"cohort_B_genes":B_genes,"cohort_C_genes":C_genes,
        "MutPred_PPI_version":version,"public_training_rows_audited_n":train_n,
        "differential_exclusion_genes":diff_exclusions,
        "verdict": ("TECHNICAL_HOLD_DISRUPTION_SCORING" if row_cov < 90 or genes_ge70 == 0 else
                    ("GO_PHASE1B1_FULL" if len(A_genes) >= 30 else
                     ("GO_PHASE1B1_FOCUSED" if len(A_genes) >= 10 else
                      ("NO_GO_EDGETIC_REPRESENTATION" if score_n >= int(0.90*len(universe)) else
                       "TECHNICAL_HOLD_DISRUPTION_SCORING")))),
        "disease_label_tests_started":False,
        "stop_rule":"No PERMANOVA, disease clustering, prediction, enrichment, tissue concordance, classification, or VUS phenotype prediction in Phase 1B-0.",
    }
    (META/"phase1B0_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    write_report(manifest, cohorts, ddf, hotspot, ml, foldx_rows, candidate_df, domain, generic3d,
                 position, overlap_df, version, pathogenicity)
    verify_inputs_unchanged(before)
    print(json.dumps({"manifest":manifest,"report":str(REPORT.relative_to(ROOT)),
                      "outputs":len(list(OUT.glob("*.tsv")))},ensure_ascii=False,indent=2))
    return manifest


def write_analysis_specification(cohorts: pd.DataFrame, gene_diff: dict, training_rows: int) -> None:
    text = f"""# Phase 1B-0 analysis specification (frozen before outcome testing)

Freeze date: {utcnow()} (UTC)  
Phase 1A cutoff: ClinVar 2026-09-03 snapshot, reused unchanged.  
Model release: MutPred-PPI current checkpoint, repository commit `{PRIMARY_MODEL_COMMIT}`, `arch=current`.  
Disease-label association testing: **not started**.

## Cohorts

- **A, primary disruption cohort:** Phase 1A structural gate PASS; exclude confirmed incompatible mechanism; exclude `MATERIAL_IMBALANCE` and `SEVERE_IMBALANCE` because those are sensitivity-only until matched; at least three distinct scored Tier 1–2 partner dimensions; at least 70% of all strict variants have scores for at least three partners; exclude multi-label variants from primary observation rows. A gene must have at least two retained disease groups with at least 10 unique variants after this handling.
- **B, structural sensitivity cohort:** the same Phase 1A, mechanism, and differential-coverage gates, with at least three Tier 1–2 structural dimensions available for at least 70% of strict variants. Representation is structural occupancy/distance only and stays separate from E.
- **C, separable-interface sensitivity:** A genes classified `SEPARABLE_INTERFACES`; also report the corresponding B subset as a structural sensitivity. Do not require C for A.
- **Material-imbalance sensitivity:** retain all material genes for a future within-gene matched analysis; do not mix into A before that matching is successfully frozen.

## Coverage adjudication thresholds

The disease-specific representation metric is the fraction of strict unique variants with at least three distinct Tier 1–2 partner interface models. `SEVERE_IMBALANCE` if any disease group is below 50% or the range is at least 50 percentage points. `MATERIAL_IMBALANCE` if not severe and the range is at least 20 points, or one group is below 70% while another is at least 70%. Otherwise `MINOR_IMBALANCE`. These thresholds were applied to representation only, before disease-separation analysis. Per-disease denominators and both ≥1 and ≥3 representation summaries are in `differential_coverage_adjudication.tsv`.

## Representations and distances

- **S1:** binary partner-interface occupancy (`INTERFACE=1`, observed complex non-interface=`0`).
- **S2:** experimental heavy-atom minimum distance to the partner interface in Å.
- **E:** raw MutPred-PPI disruption probability by variant–partner row; missing scores remain `NA`. Never concatenate S1/S2 with E or substitute one for another.
- **Primary E distance:** Pearson correlation distance `1 - r`, computed only for variant pairs sharing at least 3 observed partner scores and with nonzero variance in both vectors. If fewer than 3 shared dimensions or either vector is constant, distance is undefined and the pair is omitted under a predeclared complete-pair rule; report the resulting retained sample and pair counts. No imputation.
- **Sensitivity distance:** Euclidean distance after partner-wise standardization learned without disease labels, with missing dimensions pairwise omitted only when at least 3 dimensions are shared. Standardization parameters are frozen on the entire eligible gene without using disease labels.
- Evaluate score completeness before constructing distances. Do not select a distance based on disease separation.

## Future omnibus statistic and incremental tests (not run here)

For each eligible gene, test the omnibus disease-label partition of the frozen distance matrix using PERMANOVA or an equivalent distance-based partition. Primary effect size is `R²`; preserve pseudo-F, permutation P, sample size, disease-class count, and valid pair count. Permutations are stratified within gene and preserve disease group sizes. The central contrasts are `ΔR² = R²_E - R²_position`, `R²_E - R²_domain`, `R²_E - R²_generic3D`, and `R²_E - R²_interfaceNoPartner`. Report every prespecified contrast, including negative and null contrasts; no post hoc choice of the largest contrast.

## Null-model matching and exchangeability

Before a label permutation, construct coarse exchangeability blocks from domain occupancy, S1 interface/non-interface status, and AlphaMissense quartile; retain disease-group sizes within each gene and block. Balance solvent accessibility, secondary structure when a frozen annotation is available, and intrachain Cα neighbor-density quartile by nearest-bin matching. If a block has fewer than 5 observations, merge adjacent AlphaMissense quartiles first, then remove domain occupancy, then S1 status, in that fixed label-blind order. If fewer than 100 distinct valid restricted permutations remain after that coarsening, report the gene as non-testable under the primary null; do not fall back to unrestricted permutations. Freeze the random seed and number of valid permutations (minimum 9,999 where feasible) before Phase 1B-1.

## Disease classes, multi-label variants and hotspots

Use every strict disease class with ≥10 unique variants after primary multi-label exclusions. A pairwise follow-up can only occur after an omnibus test and must be labeled secondary. Variants assigned to more than one strict MONDO group are flagged `MULTI_LABEL_VARIANT`, excluded from A to avoid duplicated observations, and retained only for a dedicated sensitivity where the permutation unit is the unique `gene + VariationID` and labels are represented as a set. Genes with a disease group where one residue accounts for ≥50% or the top three residues account for ≥70% are `HOTSPOT_DOMINATED`; do not automatically exclude them. Their future sensitivity removes all variants from each hotspot residue in turn and repeats the frozen omnibus test.

## Validation and temporal cutoff

Any future classifier or predictive analysis groups folds by `gene + canonical_residue`; substitutions at one residue never split across train/test. ClinVar discovery input remains frozen at 2026-09-03. Do not query later ClinVar releases in Phase 1B-0. Temporal validation may use only variants first appearing after that cutoff and must retain first-seen date, VariationID, and RCV provenance.

## Training overlap cohorts

`ALL` includes all scoreable observations. `STRICT_EXTERNAL` excludes exact training variants and exact variant–partner overlaps identified in the public training/evaluation deposits. `PROTEIN_EXTERNAL` additionally excludes target proteins appearing anywhere in those public deposits. The unpublished VarChAMP training component is not auditable and remains an explicit limitation; absence of a public match is not proof of no overlap.

## Frozen model and score gate

The one primary version is the current MutPred-PPI repository checkpoint at `{PRIMARY_MODEL_COMMIT}` with its colocated `mutation_diff_scaler.pkl`; the v1.0 ten-model ensemble is not mixed. Required sequence embedding model is `Rostlab/prot_t5_xl_half_uniref50-enc`. Every score row retains status and reason; `NA` is never imputed. Global PASS requires ≥90% scored rows and ≥70% scored-variant coverage for every A gene; <80% global row coverage yields `TECHNICAL_HOLD_DISRUPTION_SCORING`.

## Exemplar selection

Rank candidates using only score-row/variant completeness, disease sample balance, sample size, Tier 1–2 evidence, interface separability, and public training overlap. No disease-label association or separation value enters ranking.

## Prohibited in Phase 1B-0

No PERMANOVA, disease clustering, disease prediction/classification, disease-specific partner enrichment, tissue concordance, or VUS phenotype prediction was run. Stop after this freeze and wait for Phase 1B-1 authorization.
"""
    (OUT / "analysis_specification.md").write_text(text, encoding="utf-8")


def write_report(manifest: dict, cohorts: pd.DataFrame, ddf: pd.DataFrame,
                 hotspot: pd.DataFrame, ml: pd.DataFrame, foldx_rows: list[dict],
                 candidates: pd.DataFrame, domain: pd.DataFrame, generic3d: pd.DataFrame,
                 position: pd.DataFrame, overlap: pd.DataFrame, version: dict,
                 pathogenicity: pd.DataFrame) -> None:
    score_n = manifest["successfully_scored_rows"]
    row_cov = manifest["row_score_coverage_pct"]
    labels = "Disease-label association testing was NOT started."
    A=manifest["cohort_A_genes"];B=manifest["cohort_B_genes"];C=manifest["cohort_C_genes"]
    min_score=cohorts.loc[cohorts.score_row_coverage_pct.idxmin(),"score_row_coverage_pct"] if len(cohorts) else 0
    genes_ge70=int((cohorts.variants_with_ge3_disruption_scores_pct_of_strict_variants>=70).sum())
    exact_variant_rows=int((overlap.exact_variant_training_overlap=="YES").sum())
    exact_variant_n=int(overlap.loc[overlap.exact_variant_training_overlap=="YES",["gene","VariationID"]].drop_duplicates().shape[0])
    exact_pair_rows=int((overlap.exact_variant_partner_training_overlap=="YES").sum())
    exact_pair_n=int(overlap.loc[overlap.exact_variant_partner_training_overlap=="YES",["gene","VariationID","partner"]].drop_duplicates().shape[0])
    strict_gene_n=int((overlap[overlap.strict_external_include=="YES"].groupby("gene").VariationID.nunique()>0).sum())
    alpha_matches=int((pathogenicity.source_status=="MATCHED_OFFICIAL_ALPHA_MISSENSE_2023_AA_SUBSTITUTIONS").sum())
    structural_s = read_tsv(OUT / "structural_interface_matrix.tsv")
    s1_rows=int(len(structural_s)); s2_rows=int(structural_s.S2_distance_to_interface_A.notna().sum())
    foldx_addr=sum(int(x["variants_with_observed_PDB_target_coordinates_n"]) for x in foldx_rows)
    foldx_total=sum(int(x["strict_variant_n"]) for x in foldx_rows)
    foldx_gene_ge70=sum(float(x["coordinate_addressability_pct"]) >= 70 for x in foldx_rows)
    top10=candidates.head(10)[["data_blind_rank","gene","potential_ge3_PDB_partner_variant_pct","score_success_pct","disease_balance_min_over_max","Tier1_2_partners","interface_overlap_class"]].to_string(index=False)
    report=f"""# Phase 1B-0 — partner-disruption matrix and analysis freeze

Generated: {manifest['generated_utc']} (UTC)  
Verdict: **{manifest['verdict']}**  
Phase 1A input TSVs frozen: {manifest['phase1A_input_tsv_n']}; post-run hashes verified unchanged.  
Disease-label association testing: **not started**.

## Status

Phase 1A's 42 structural-pass genes and exact 11,150 experimental-coordinate variant–partner rows were reconstructed. The fixed MutPred-PPI current checkpoint and scaler are identified and hashed. The official inference dependency set is installed. Scoring remains blocked because the required ProtT5 encoder weights at `Rostlab/prot_t5_xl_half_uniref50-enc` could not be retrieved from its official Hugging Face host after TLS handshake failures. Consequently, **{score_n}/{manifest['scoring_universe_rows']}** rows have scores ({row_cov:.2f}%). No structural, pathogenicity, or interface proxy was substituted. All score rows remain `NA` with a reason.

The required score gate is not met. The primary disruption cohort A is therefore empty at this freeze. Cohort B remains the pre-outcome structural sensitivity cohort (**{len(B)} genes**); cohort C is the separable-interface subset of A (**{len(C)} genes**). No disease-label statistic was computed.

## Frozen model

- Primary model: `{version['primary_version']}` (`arch=current`).
- Release: current all-data checkpoint, not the v1.0 ensemble.
- Repository commit: `{version['repository_commit']}`.
- Checkpoint and scaler SHA256 values are in `metadata/PHASE1B0_MUTPRED_VERSION_FREEZE.tsv`.
- Reference: Stewart et al., RECOMB 2026, DOI `10.64898/2025.12.20.695738`.
- Sequence encoder: `Rostlab/prot_t5_xl_half_uniref50-enc`; exact resolved revision/checksum cannot be frozen while official weights are inaccessible.
- Public training audit includes 12,529 rows; unpublished VarChAMP training rows remain unauditable.

## Cohort and differential coverage

- Structural starting cohort: 42 genes in `results/phase1B0/structural_42_genes.tsv`.
- The exact 11,150-row Phase 1A PDB scoring universe spans {manifest.get('scoring_universe_gene_n', overlap.gene.nunique())} genes; {manifest.get('structural_42_genes_with_exact_PDB_score_rows_n', 'NA')} of the 42 structural-pass genes have at least one exact PDB scoring row. This preserves the frozen Phase 1A universe while keeping primary eligibility restricted to the 42-gene cohort.
- Cohort A (primary disruption): {len(A)} genes; gated by ≥3 scored Tier 1–2 partner dimensions and ≥70% single-label primary-variant score completeness, in addition to the frozen structural/mechanism/differential rules.
- Cohort B (structural sensitivity): {len(B)} genes.
- Cohort C (separable interfaces within A): {len(C)} genes.
- Material/severe differential-coverage sensitivity-only exclusions: {', '.join(manifest['differential_exclusion_genes']) or 'none'}.
- Thresholds and each strict disease group's counts/coverage are recorded in `differential_coverage_adjudication.tsv`.

## Scoring universe and overlap

- Eligible observations: {manifest['scoring_universe_rows']} unique `(gene, VariationID, partner)` rows; 5,357 unique variant keys in the full universe.
- Successful scores: {score_n}; row coverage: {row_cov:.2f}%.
- Variants with ≥3 partner scores: {manifest['variants_with_ge3_scores_n']}.
- Genes with ≥70% of single-label primary variants having ≥3 disruption scores: {genes_ge70}.
- Exact public training variant overlaps: {exact_variant_n} unique `(gene, VariationID)` variants across {exact_variant_rows} variant–partner observation rows.
- Exact public training variant–partner overlaps: {exact_pair_n} unique `(gene, VariationID, partner)` pairs across {exact_pair_rows} rows; pairs without a uniquely mapped partner accession remain `UNKNOWN`.
- Genes with at least one strict-external observation row: {strict_gene_n}; the rows are not disruption-scoreable while score coverage is 0%.
- `ALL`, `STRICT_EXTERNAL`, and `PROTEIN_EXTERNAL` inclusion flags are frozen in `mutpred_training_overlap.tsv`.

## Baselines

- Linear position built for {len(position)} unique mapped structural-cohort variants; normalized positions use 1-based canonical residue and frozen UniProt sequence length.
- Domain baseline produced for {len(domain)} variants from locally frozen UniProtKB domain-feature responses; missing annotations remain explicit.
- Generic 3D baseline produced for {len(generic3d)} variants from partner-independent AlphaFold DB monomer coordinates (pLDDT, isolated-chain SASA, and intrachain Cα density where available); partner identity is excluded.
- Generic 3D residue features were available for {int(generic3d.generic3d_status.eq('OK').sum())}/{len(generic3d)} variants; absent structures or residue mappings remain explicit.
- General pathogenicity uses the frozen AlphaMissense 2023 amino-acid substitution resource; {alpha_matches}/{len(pathogenicity)} canonical accession/substitution keys matched, and other values remain `NA`.
- Partner-free structural baseline is aggregated to one row per variant and contains only any-interface status/counts, with no partner names.
- Separate S1/S2 structural matrix: {s1_rows} rows across the 42 structural genes. S1 retains Phase 1A PIONEER annotations and experimental PDB occupancy; S2 distances are available only for {s2_rows} experimental-coordinate rows, with PIONEER-only distances left `NA`.

## Hotspots, multi-labels, and FoldX

- Hotspot-dominated genes: {', '.join(hotspot.loc[hotspot.hotspot_status=='HOTSPOT_DOMINATED','gene'].tolist()) or 'none'} (frozen thresholds: any disease with one residue ≥50% or top three residues ≥70%).
- Multi-label variants: {int((ml.multi_label_variant_flag=='MULTI_LABEL_VARIANT').sum())}; excluded from the future primary observation set and retained for dedicated sensitivity.
- FoldX was not run. Experimental-coordinate addressability is {foldx_addr}/{foldx_total} strict variants ({100*foldx_addr/max(1,foldx_total):.2f}%) overall and ≥70% for {foldx_gene_ge70}/42 genes; these estimates do not establish FoldX compatibility.

## Data-blind exemplar candidates

Ranked by scoreability, disease balance, sample size, Tier 1–2 PPI evidence, interface separability and public training overlap. Since score inference is blocked, these are **technical candidates pending score QC**, not validated exemplars.

```text
{top10}
```

## Frozen analysis plan and stop rule

The complete distance metrics, minimum shared dimensions, PERMANOVA effect size, incremental R² contrasts, constrained permutation scheme, multilabel rule, hotspot control, residue-grouped validation and ClinVar temporal cutoff are specified in `results/phase1B0/analysis_specification.md`. No disease clustering, PERMANOVA, disease prediction/classification, disease-specific enrichment, tissue concordance or VUS phenotype prediction was run.

Phase 1B-0 cannot pass until the required primary version can generate scores for the frozen rows, score coverage is assessed, and baselines are complete. This is a technical hold, not a negative biological result. After those tasks are completed and frozen, stop and wait for Phase 1B-1 authorization.

{labels}
"""
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(report, encoding="utf-8")


if __name__ == "__main__":
    build_all()
