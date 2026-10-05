#!/usr/bin/env python3
"""Record the Phase 1B-0R2 exact MutPred-PPI v1 technical recovery state.

This script only audits frozen inputs and writes blocked-scoring outputs. It
does not run model inference, rebuild cohorts, or test disease-label effects.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase1B0R2"
REPORT = ROOT / "reports" / "PHASE1B0R2_EXACT_HALF_MODEL_RECOVERY.md"
MODEL_DIR = ROOT / "data/software/models/prot_t5_xl_half_uniref50-enc"
EXPECTED_ROWS = 11_150
MODEL_ID = "Rostlab/prot_t5_xl_half_uniref50-enc"
MODEL_SHA256 = "7f51ba885541c7dc569d46b796af57cc7a2ba7945107dced4f19d1b5ec091157"
MUTPRED_VERSION = "v1.0.0"
MUTPRED_COMMIT = "a6e2f47e18a3af5dbbc4d8accce8ab177049c56c"
FAILURE = "PROTT5_TLS_FAILURE"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(frame: pd.DataFrame, name: str) -> None:
    frame.to_csv(OUT / name, sep="\t", index=False, na_rep="NA")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (ROOT / "reports").mkdir(parents=True, exist_ok=True)

    freeze = pd.read_csv(ROOT / "metadata/PHASE1B0_INPUT_FREEZE_MANIFEST.tsv", sep="\t")
    freeze_rows = []
    for row in freeze.itertuples(index=False):
        p = ROOT / row.path
        actual = sha256(p) if p.is_file() else "MISSING"
        freeze_rows.append({
            "path": row.path,
            "expected_sha256": row.SHA256,
            "actual_sha256": actual,
            "status": "PASS" if actual == row.SHA256 else "INPUT_FREEZE_VIOLATION",
        })
    freeze_audit = pd.DataFrame(freeze_rows)
    write_tsv(freeze_audit, "input_freeze_audit.tsv")
    if not freeze_audit.status.eq("PASS").all():
        raise SystemExit("INPUT_FREEZE_VIOLATION")

    universe = pd.read_csv(
        ROOT / "results/phase1B0/variant_partner_scoring_universe.tsv",
        sep="\t", dtype={"VariationID": str},
    )
    if len(universe) != EXPECTED_ROWS:
        raise SystemExit(f"SCORING_UNIVERSE_MISMATCH: observed {len(universe)}")
    overlap = pd.read_csv(
        ROOT / "results/phase1B0/mutpred_training_overlap.tsv",
        sep="\t", dtype={"VariationID": str},
    )
    keys = ["gene", "VariationID", "partner"]
    overlap = overlap.drop_duplicates(keys)
    overlap_cols = keys + [
        "exact_variant_training_overlap", "exact_variant_partner_training_overlap",
        "strict_external_include", "protein_external_include",
    ]
    frame = universe.merge(overlap[overlap_cols], on=keys, how="left", validate="one_to_one")
    for col in overlap_cols[3:]:
        frame[col] = frame[col].fillna("UNKNOWN")

    score = pd.DataFrame({
        "gene": frame.gene,
        "VariationID": frame.VariationID,
        "variant": frame.canonical_variant,
        "partner": frame.partner,
        "structure_id": "NA",
        "structure_candidates": frame.experimental_PDB_IDs,
        "mutated_chain": "NA",
        "partner_chain": "NA",
        "MutPred_PPI_version": MUTPRED_VERSION,
        "score": np.nan,
        "status": "NOT_SCORED_TECHNICAL_HOLD",
        "failure_reason": FAILURE,
        "failure_detail": (
            "Official Hugging Face HTTPS connection terminated during TLS handshake; "
            "the exact half-encoder files were not obtained and inference was not started."
        ),
        "disease_group": frame.disease_group,
        "exact_variant_training_overlap": frame.exact_variant_training_overlap,
        "exact_variant_partner_training_overlap": frame.exact_variant_partner_training_overlap,
        "strict_external_include": frame.strict_external_include,
        "protein_external_include": frame.protein_external_include,
    })
    write_tsv(score, "mutpred_ppi_v1_scores.tsv")
    write_tsv(score[[
        "gene", "VariationID", "variant", "partner", "structure_candidates",
        "MutPred_PPI_version", "status", "failure_reason", "failure_detail",
    ]], "failure_registry.tsv")

    variant = (
        score.groupby(["gene", "VariationID", "variant"], dropna=False)
        .agg(eligible_partner_rows=("partner", "nunique"), scored_partner_rows=("score", "count"))
        .reset_index()
    )
    variant["variants_with_ge1_partner_score"] = variant.scored_partner_rows.ge(1)
    variant["variants_with_ge3_partner_scores"] = variant.scored_partner_rows.ge(3)
    variant["status"] = "NOT_EVALUABLE_TECHNICAL_HOLD"
    write_tsv(variant, "coverage_by_variant.tsv")

    structural = pd.read_csv(ROOT / "results/phase1B0/structural_42_genes.tsv", sep="\t")
    protein_map = pd.read_csv(
        ROOT / "results/phase1A/protein_mapping.tsv", sep="\t", dtype={"VariationID": str}
    )
    mondo_count = protein_map["MONDOs"].fillna("").astype(str).str.count(r"MONDO:")
    primary = protein_map[
        protein_map.primary_cohort_eligible.astype(str).str.upper().eq("YES") & mondo_count.eq(1)
    ]
    primary_n = primary.groupby("gene").VariationID.nunique().to_dict()
    row_n = score.groupby("gene").size().to_dict()
    var_n = score.groupby("gene").VariationID.nunique().to_dict()
    gene_rows = []
    for gene in structural.gene.astype(str):
        den = int(primary_n.get(gene, 0))
        gene_rows.append({
            "gene": gene,
            "eligible_single_label_primary_variants": den,
            "input_variant_partner_rows": int(row_n.get(gene, 0)),
            "input_variants": int(var_n.get(gene, 0)),
            "scored_rows": 0,
            "variants_ge1_partner_score": 0,
            "variants_ge3_partner_scores": 0,
            "primary_variant_ge3_coverage_pct": np.nan if den else np.nan,
            "passes_70pct_disruption_representation": "NOT_EVALUABLE",
            "status": "TECHNICAL_HOLD_NO_VALID_MODEL_SCORES",
        })
    gene_coverage = pd.DataFrame(gene_rows)
    write_tsv(gene_coverage, "coverage_by_gene.tsv")

    prior_diff = pd.read_csv(
        ROOT / "results/phase1B0/differential_coverage_adjudication.tsv", sep="\t"
    )
    prior_flags = set(prior_diff.loc[
        prior_diff.gene_coverage_adjudication.astype(str).str.upper().isin(
            ["SEVERE_IMBALANCE", "MATERIAL_IMBALANCE"]
        ), "gene"
    ].astype(str))
    diff = pd.DataFrame({"gene": structural.gene.astype(str)})
    diff["pre_scoring_representation_imbalance_flag"] = diff.gene.map(
        lambda g: "YES" if g in prior_flags else "NO"
    )
    diff["differential_score_coverage_status"] = "PENDING_NO_VALID_MUTPRED_SCORES"
    diff["differential_coverage_exclusion_applied"] = "NO"
    diff["score_coverage_recomputed"] = "NO"
    diff["note"] = "Exact half-encoder dependency unavailable at TLS handshake; differential coverage not adjudicated."
    write_tsv(diff, "differential_score_coverage.tsv")

    write_tsv(pd.DataFrame([{
        "rebuild_status": "NOT_REBUILT_PREREQUISITE_NOT_MET",
        "reason": "Exact v1 example and scoring coverage gates require validated model scores.",
        "rebuilt_primary_gene_count": "NA",
        "strict_external_gene_count": "NA",
        "disease_label_association_testing": "NOT_STARTED",
    }]), "rebuilt_primary_cohort.tsv")
    write_tsv(pd.DataFrame([{
        "example": "official v1.0.0 release example",
        "reference_rows": 10,
        "rows_scored": 0,
        "max_abs_difference": "NA",
        "status": "NOT_RUN_EXACT_HALF_ENCODER_UNAVAILABLE",
    }]), "official_example_validation.tsv")
    for n in (20, 100):
        write_tsv(pd.DataFrame([{
            "pilot_size": n,
            "rows_scored": 0,
            "status": "NOT_RUN_PREREQUISITE_BLOCKED",
            "failure_reason": FAILURE,
            "selection": "Pilot not run because exact v1 example prerequisite could not be met.",
        }]), f"pilot{n}.tsv")
    write_tsv(pd.DataFrame([{
        "batch_id": "NONE",
        "rows_attempted_in_model": 0,
        "rows_scored": 0,
        "rows_blocked": EXPECTED_ROWS,
        "primary_failure_reason": FAILURE,
        "status": "BLOCKED_BEFORE_BATCHING",
    }]), "batch_manifest.tsv")
    write_tsv(pd.DataFrame(columns=[
        "sequence_id", "sequence_sha256", "model", "embedding_path_or_cache_key",
        "success", "length", "timestamp",
    ]), "prott5_embedding_manifest.tsv")

    expected_files = [
        "config.json", "pytorch_model.bin", "spiece.model",
        "tokenizer_config.json", "special_tokens_map.json",
    ]
    downloads = []
    for filename in expected_files:
        p = MODEL_DIR / filename
        nbytes = p.stat().st_size if p.is_file() else 0
        actual = sha256(p) if p.is_file() else "NOT_OBTAINED"
        downloads.append({
            "model_id": MODEL_ID,
            "filename": filename,
            "expected_sha256": MODEL_SHA256 if filename == "pytorch_model.bin" else "NOT_SPECIFIED",
            "expected_size": "approximately 2.42 GB" if filename == "pytorch_model.bin" else "NOT_SPECIFIED",
            "local_path": str(p),
            "bytes_obtained": nbytes,
            "actual_sha256": actual,
            "status": "OBTAINED" if p.is_file() and nbytes else "NOT_OBTAINED",
            "method_attempted": "official HF HTTPS resolve URL and official hf CLI metadata request",
            "failure": "TLS EOF before HTTP response; no metadata or artifact transfer",
        })
    write_tsv(pd.DataFrame(downloads), "prott5_download_manifest.tsv")

    tls = pd.DataFrame([
        {
            "test": "curl -Iv https://huggingface.co/",
            "transport": "direct; all HTTP(S)_PROXY/ALL_PROXY/NO_PROXY variants unset",
            "ssl_backend": "curl 8.7.1; SecureTransport; LibreSSL 3.3.6",
            "dns_observation": "huggingface.co resolved by dscacheutil to IPv4-mapped 198.18.0.32",
            "http_status": "NONE",
            "redirects": 0,
            "result": "TLS handshake ended with SSL_ERROR_SYSCALL before HTTP response",
        },
        {
            "test": "curl -L official half-encoder resolve/main/config.json",
            "transport": "direct; no proxy configured",
            "ssl_backend": "curl 8.7.1; SecureTransport; LibreSSL 3.3.6",
            "dns_observation": "huggingface.co resolved to 198.18.0.32",
            "http_status": "000",
            "redirects": 0,
            "result": "curl exit 35; TLS SSL_ERROR_SYSCALL; no response body or redirect",
        },
        {
            "test": "hf download exact model config.json",
            "transport": "official Hugging Face CLI; huggingface_hub 0.36.2; no proxy configured",
            "ssl_backend": "Python requests / urllib3 SSL",
            "dns_observation": "same huggingface.co hostname",
            "http_status": "NONE",
            "redirects": "NOT_REACHED",
            "result": "five HEAD retries ended in SSLEOFError before metadata; no Xet handoff",
        },
        {
            "test": "Xet fallback / direct curl artifact transfer",
            "transport": "not attempted",
            "ssl_backend": "N/A",
            "dns_observation": "metadata prerequisite failed",
            "http_status": "N/A",
            "redirects": "N/A",
            "result": "not applicable because metadata and direct HTTPS were unavailable",
        },
    ])
    write_tsv(tls, "tls_diagnostics.tsv")

    dependency = pd.DataFrame([
        ("MutPred-PPI version", MUTPRED_VERSION, MUTPRED_VERSION, "PASS"),
        ("MutPred-PPI source commit", MUTPRED_COMMIT, MUTPRED_COMMIT, "PASS"),
        ("ProtT5 model ID", MODEL_ID, MODEL_ID, "FROZEN; WEIGHTS UNAVAILABLE"),
        ("ProtT5 weight SHA256", MODEL_SHA256, "NOT_OBTAINED", "UNVERIFIED"),
        ("HF CLI", "official hf download", "available; metadata request failed TLS", "NETWORK_BLOCKED"),
        ("curl SSL backend", "TLS to official HF host", "SecureTransport / LibreSSL 3.3.6", "TLS_EOF"),
        ("proxy variables", "record credential-free", "all checked variants unset", "NO_PROXY_CONFIGURED"),
        ("offline path patch", "only after verified model retrieval", "not applied", "NOT_REQUIRED_YET"),
        ("Xet fallback", "only after metadata success and Xet failure", "not attempted", "PREREQUISITE_FAILED"),
        ("inference precision semantics", "official GPU fp16 / CPU fp32 behavior", "not exercised", "NOT_RUN"),
        ("authorized off-machine transfer", "permitted if available", "no secondary machine exposed to this session", "UNAVAILABLE_HERE"),
    ], columns=["component", "expected", "observed", "status"])
    write_tsv(dependency, "dependency_manifest.tsv")

    model_checks = ROOT / "results/phase1B0R/mutpred_v1_model_checksums.tsv"
    if model_checks.is_file():
        (OUT / "mutpred_v1_model_checksums.tsv").write_bytes(model_checks.read_bytes())

    variant_n = int(variant.shape[0])
    structure_n = int(structural.gene.astype(str).nunique())
    report = f"""# PHASE 1B-0R2 — Exact MutPred-PPI v1 Recovery

## Verdict

**TECHNICAL_HOLD_DISRUPTION_SCORING**

The exact v1.0 ProtT5 half encoder is frozen as `{MODEL_ID}`. Official Hugging Face TLS failed before HTTP metadata was returned, so the weights, official example validation, pilots, and production inference could not run. This is a technical hold; no biological NO-GO was assessed.

## Frozen inputs and exact software

- MutPred-PPI version: `{MUTPRED_VERSION}`; source commit: `{MUTPRED_COMMIT}`.
- Frozen scoring universe: {len(universe):,} variant–partner rows; required count: {EXPECTED_ROWS:,}.
- Frozen-input audit: {len(freeze_audit)} files checked, {int(freeze_audit.status.eq('PASS').sum())} matched.
- Expected ProtT5 weight: `pytorch_model.bin`, approximately 2.42 GB, expected SHA256 `{MODEL_SHA256}`.
- Local half-encoder weight: **not obtained**; local weight bytes: 0; observed SHA256: `NOT_OBTAINED`; expected hash match: **not verifiable**.
- No offline-path-only patch was applied. No model architecture, tokenizer, precision, or weight substitution was made. The earlier full-precision-size target was not used.

## Acquisition and network diagnostics

The direct `curl -Iv https://huggingface.co/` test resolved the host to `198.18.0.32` and ended during the TLS handshake with `SSL_ERROR_SYSCALL` (curl 8.7.1 using SecureTransport / LibreSSL 3.3.6). All upper- and lower-case proxy environment variables checked were unset. The official half-model `resolve/main/config.json` request returned curl HTTP code `000`, zero redirects, and exit code 35. The official `hf download` config metadata request retried five times and ended with Python `SSLEOFError` before response metadata.

No Xet transfer occurred because metadata was unavailable. Direct curl resume/retry was not appropriate because direct HTTPS itself failed. A secondary authorized workstation was not accessible through this session, so no off-machine transfer was performed. Diagnostics are in `results/phase1B0R2/tls_diagnostics.tsv`; the exact-file acquisition ledger is in `results/phase1B0R2/prott5_download_manifest.tsv`.

## Reproduction gate and scoring

- Official v1 example: **not run**; exact ProtT5 files unavailable. Maximum absolute error: not available.
- Pilot-20 and pilot-100: **not run**, as required before example validation.
- Frozen scoring rows: {len(score):,}; successfully scored: 0; row coverage: 0.00%.
- Variants with ≥1 partner score: 0/{variant_n:,}. Variants with ≥3 partner scores: 0/{variant_n:,}.
- Structural genes with observed ≥70% disruption coverage: 0/{structure_n}; the coverage gate is **not evaluable** without valid scores.
- Most common row failure: `{FAILURE}` ({len(score):,} rows).
- Differential score-coverage exclusions applied: 0. The {len(prior_flags)} pre-scoring representation-imbalance flags remain pending; no score-based exclusion decision was made.
- Primary and strict-external cohorts were not rebuilt. Their counts are `NA`, not zero; the old Phase 1B-0 zero-gene result is not carried forward.
- Disease labels were not used for association or outcome analysis.

## Stop point

Remain at `TECHNICAL_HOLD_DISRUPTION_SCORING`. Resume only after the exact official half-encoder files are obtained and the local weight SHA256 matches the frozen value, followed by successful official example validation within absolute tolerance 1e-4. Then continue with the frozen pilots, full scoring, coverage gates, and de novo cohort reconstruction. Do not enter Phase 1B-1 from this technical-hold state.

Disease-label association testing was NOT started.
"""
    REPORT.write_text(report, encoding="utf-8")

    print(f"verdict=TECHNICAL_HOLD_DISRUPTION_SCORING")
    print(f"frozen_rows={len(universe)} scored_rows=0 coverage_pct=0.00")
    print(f"variants={variant_n} structural_genes={structure_n} differential_flags_pending={len(prior_flags)}")
    print(f"input_freeze={len(freeze_audit)} PASS={int(freeze_audit.status.eq('PASS').sum())}")
    print(f"report={REPORT}")


if __name__ == "__main__":
    main()
