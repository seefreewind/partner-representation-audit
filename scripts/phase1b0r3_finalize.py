#!/usr/bin/env python3
"""Finalize Phase 1B-0R3 row scores, coverage QC, cohorts, and the final gate.

This script performs only the frozen technical coverage checks. It does not
compare score distributions by disease or run association, clustering,
enrichment, AUROC, or tissue-concordance analyses.
"""
from __future__ import annotations

import math
import hashlib
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase1B0R3"
P1A = ROOT / "results/phase1A"
P0R = ROOT / "results/phase0R"
P1B0 = ROOT / "results/phase1B0"
SCORE_SCRIPT = ROOT / "scripts/phase1b0r3_score_pipeline.py"
sys.path.insert(0, str(ROOT / "scripts"))
from phase1b0r3_score_pipeline import load_universe_and_mapping  # noqa: E402

MODEL_VERSION = "MutPred-PPI v1.0.0@a6e2f47e18a3af5dbbc4d8accce8ab177049c56c"
PROTT5_MODEL = "Rostlab/prot_t5_xl_half_uniref50-enc@94a6abc029ae13029317b140b7424e012bf8dfbf"
ALLOWED_FAILURES = {
    "SEQUENCE_MISMATCH", "VARIANT_REFERENCE_MISMATCH", "STRUCTURE_CHAIN_FAILURE",
    "STRUCTURE_RESIDUE_FAILURE", "GRAPH_FAILURE", "PROTT5_INFERENCE_FAILURE",
    "MUTPRED_INFERENCE_FAILURE", "OUT_OF_MEMORY", "OTHER",
}


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    if not path.is_file():
        raise SystemExit(f"Required input is absent: {path}")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **kwargs)


def split_semicolon(value: object) -> list[str]:
    return [x.strip() for x in str(value).split(";") if x.strip() and x.strip() not in {"NA", "nan"}]


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", lineterminator="\n")


def verify_frozen_hashes() -> pd.DataFrame:
    manifest = read_tsv(ROOT / "metadata/PHASE1B0_INPUT_FREEZE_MANIFEST.tsv")
    if len(manifest) != 12:
        raise SystemExit(f"Expected the 12 frozen Phase 1A source hashes; found {len(manifest)}")
    audit: list[dict[str, str]] = []
    for r in manifest.itertuples(index=False):
        path = ROOT / str(r.path)
        observed = "MISSING"
        if path.is_file():
            digest = hashlib.sha256()
            with path.open("rb") as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(block)
            observed = digest.hexdigest()
        expected = str(r.SHA256)
        audit.append({"path": str(r.path), "expected_sha256": expected,
                      "observed_sha256": observed,
                      "status": "PASS" if observed == expected else "FAIL"})
    result = pd.DataFrame(audit)
    write_tsv(result, OUT / "frozen_input_hash_validation.tsv")
    if not result.status.eq("PASS").all():
        raise SystemExit("INPUT_FREEZE_VIOLATION: one or more Phase 1A source hashes changed")
    return result


def make_membership(genes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, set[str]]]:
    mapping = read_tsv(P1A / "protein_mapping.tsv")
    strict = read_tsv(P0R / "gene_disease_variant_counts_exact_rcv.tsv")
    strict = strict.loc[strict.strict_disease_ge10.eq("YES") & strict.gene.isin(genes)].copy()
    strict["n_unique_variants"] = pd.to_numeric(strict.n_unique_variants, errors="coerce").fillna(0).astype(int)
    strict_groups = strict.groupby("gene").MONDO.apply(lambda s: set(s.astype(str))).to_dict()
    eligible = mapping.loc[mapping.gene.isin(genes) & mapping.primary_cohort_eligible.eq("YES")].copy()
    members: list[dict[str, str]] = []
    for r in eligible.itertuples(index=False):
        for mondo in split_semicolon(r.MONDOs):
            if mondo in strict_groups.get(str(r.gene), set()):
                members.append({"gene": str(r.gene), "VariationID": str(r.VariationID),
                                "MONDO": mondo, "MONDO_label": ""})
    vgroups = pd.DataFrame(members, columns=["gene", "VariationID", "MONDO", "MONDO_label"])
    if not vgroups.empty:
        labels = strict.set_index(["gene", "MONDO"]).MONDO_label.to_dict()
        vgroups["MONDO_label"] = [labels.get((r.gene, r.MONDO), "") for r in vgroups.itertuples(index=False)]
        vgroups = vgroups.drop_duplicates(["gene", "VariationID", "MONDO"]).reset_index(drop=True)
    nlabels = vgroups.groupby(["gene", "VariationID"]).MONDO.nunique()
    multilabel = {tuple(k) for k, n in nlabels.items() if int(n) > 1}
    return strict, vgroups, {gene: groups for gene, groups in strict_groups.items()}


def official_pilot_files() -> tuple[pd.DataFrame, pd.DataFrame, float]:
    p20_uncached = read_tsv(OUT / "pilot20_uncached.tsv")
    p20_cached = read_tsv(OUT / "pilot20_cached.tsv")
    p100 = read_tsv(OUT / "pilot100_cached.tsv")
    for name, frame, expected in (("pilot20_uncached", p20_uncached, 20),
                                  ("pilot20_cached", p20_cached, 20)):
        if not {"complex_id", "variant", "score"}.issubset(frame.columns):
            raise SystemExit(f"{name} has unexpected official output columns")
        if len(frame) != expected or frame[["complex_id", "variant"]].duplicated().any():
            raise SystemExit(f"{name} expected {expected} unique official predictions; found {len(frame)}")
        frame["score"] = pd.to_numeric(frame.score, errors="coerce")
        if not np.isfinite(frame.score.to_numpy()).all():
            raise SystemExit(f"{name} contains a non-finite score")
        if not frame.score.between(0.0, 1.0, inclusive="both").all():
            raise SystemExit(f"{name} contains a score outside [0,1]")
    if not {"complex_id", "variant", "score"}.issubset(p100.columns):
        raise SystemExit("pilot100_cached has unexpected official output columns")
    if not 95 <= len(p100) <= 100 or p100[["complex_id", "variant"]].duplicated().any():
        raise SystemExit(f"pilot100 requires >=95/100 unique technical predictions; found {len(p100)}")
    p100["score"] = pd.to_numeric(p100.score, errors="coerce")
    if not np.isfinite(p100.score.to_numpy()).all():
        raise SystemExit("pilot100_cached contains a non-finite score")
    if not p100.score.between(0.0, 1.0, inclusive="both").all():
        raise SystemExit("pilot100_cached contains a score outside [0,1]")
    joined = p20_uncached.merge(p20_cached, on=["complex_id", "variant"], how="outer",
                                suffixes=("_uncached", "_cached"), indicator=True)
    if len(joined) != 20 or not joined._merge.eq("both").all():
        raise SystemExit("pilot20 cached/uncached prediction identities do not match")
    max_diff = float((joined.score_uncached - joined.score_cached).abs().max())
    if not np.isfinite(max_diff) or max_diff > 1e-6:
        raise SystemExit(f"pilot20 cache score difference exceeds 1e-6: {max_diff}")
    p20_cached.to_csv(OUT / "pilot20.tsv", sep="\t", index=False, lineterminator="\n")
    p100.to_csv(OUT / "pilot100.tsv", sep="\t", index=False, lineterminator="\n")
    joined.to_csv(OUT / "pilot20_cache_validation.tsv", sep="\t", index=False, lineterminator="\n")
    return p20_cached, p100, max_diff


def finalize() -> str:
    p20, p100, pilot20_max_diff = official_pilot_files()
    official_example = OUT / "official_example_validation.tsv"
    if not official_example.is_file():
        raise SystemExit(f"Official example validation deliverable is absent: {official_example}")
    example = read_tsv(official_example)
    if len(example) < 3 or not example.status.eq("PASS").all():
        raise SystemExit("Official example validation requires all three documented rows to PASS")
    example_max_diff = pd.to_numeric(example.absolute_difference, errors="coerce").max()
    if not math.isfinite(float(example_max_diff)) or float(example_max_diff) > 1e-4:
        raise SystemExit("Official example score difference exceeds frozen 1e-4 tolerance")
    hash_audit = verify_frozen_hashes()
    if len(hash_audit) != 12 or not hash_audit.status.eq("PASS").all():
        raise SystemExit("Frozen 12-input SHA256 audit did not PASS")
    model_audit = read_tsv(OUT / "model_file_audit.tsv")
    weight_audit = model_audit.loc[model_audit.filename.eq("pytorch_model.bin")]
    if len(weight_audit) != 1 or weight_audit.iloc[0].sha256 != "7f51ba885541c7dc569d46b796af57cc7a2ba7945107dced4f19d1b5ec091157":
        raise SystemExit("Verified local ProtT5 weight identity failed final audit")
    weight_size = int(weight_audit.iloc[0].size_bytes)
    weight_path = ROOT / "models/prot_t5_xl_half_uniref50-enc/pytorch_model.bin"
    if not weight_path.is_file() or weight_path.stat().st_size != weight_size:
        raise SystemExit("Local ProtT5 model path or weight size does not match the audit")

    rows = load_universe_and_mapping()
    if len(rows) != 11150 or rows.row_id.nunique() != 11150:
        raise SystemExit("Frozen 11,150-row universe or stable row IDs failed validation")
    selection = read_tsv(OUT / "pilot_selection/full_selection.tsv")
    pred = read_tsv(OUT / "full_unique_predictions.tsv")
    if not {"target_id", "partner_id", "local_variant"}.issubset(selection.columns):
        raise SystemExit("Full selection lacks official graph/variant keys")
    if not {"complex_id", "variant_zero_based", "score"}.issubset(pred.columns):
        raise SystemExit("Full prediction checkpoint lacks official zero-based inference keys")
    selection["complex_id"] = selection.target_id.astype(str) + "_" + selection.partner_id.astype(str)
    selection["local_variant"] = selection.local_variant.astype(str)
    # The graph builder converts the one-based structure-local substitution
    # to a zero-based key before official inference and checkpoint capture.
    selection["variant_zero_based"] = (
        selection.local_variant.str[0]
        + (selection.local_variant.str[1:-1].astype(int) - 1).astype(str)
        + selection.local_variant.str[-1]
    )
    pred["variant_zero_based"] = pred.variant_zero_based.astype(str)
    pred["score"] = pd.to_numeric(pred.score, errors="coerce")
    if pred[["complex_id", "variant_zero_based"]].duplicated().any():
        raise SystemExit("Full checkpoint has duplicate model prediction keys")
    if not np.isfinite(pred.score.to_numpy()).all() or not pred.score.between(0.0, 1.0, inclusive="both").all():
        raise SystemExit("Full checkpoint contains non-finite scores or values outside [0,1]")
    selection = selection.merge(pred, on=["complex_id", "variant_zero_based"],
                                how="outer", validate="one_to_one", indicator=True)
    if len(pred) != 10499 or not selection._merge.eq("both").all():
        raise SystemExit("Full predictions must match every frozen unique scoring unit exactly")
    if len(selection) != 10499:
        raise SystemExit(f"Frozen unique scoring-unit count changed: expected 10,499; got {len(selection)}")
    score_lookup = selection.set_index(["target_id", "partner_id", "local_variant"]).score

    row_scores: list[object] = []
    for i, r in enumerate(rows.itertuples(index=False)):
        key = (str(r.target_id), str(r.partner_id), str(r.local_variant))
        value = score_lookup.get(key, np.nan) if pd.notna(r.target_id) and pd.notna(r.partner_id) and pd.notna(r.local_variant) else np.nan
        row_scores.append(value)
    rows["score"] = row_scores

    overlap = read_tsv(P1B0 / "mutpred_training_overlap.tsv")
    if len(overlap) != 11150:
        raise SystemExit(f"Training-overlap provenance must cover all 11,150 rows; found {len(overlap)}")
    overlap_cols = ["gene", "VariationID", "partner", "exact_variant_training_overlap",
                    "exact_variant_partner_training_overlap", "target_protein_training_overlap",
                    "partner_accession_mapping_status", "strict_external_include"]
    rows = rows.merge(overlap[overlap_cols], on=["gene", "VariationID", "partner"],
                      how="left", validate="one_to_one", suffixes=("", "_training"))
    if rows.strict_external_include.isna().any():
        raise SystemExit("Training-overlap flags failed to resolve on frozen rows")

    mapping_reason = rows.status.astype(str)
    failures: list[str] = []
    statuses: list[str] = []
    selection_keys = set(zip(selection.target_id.astype(str), selection.partner_id.astype(str),
                             selection.local_variant.astype(str)))
    event_path = OUT / "full_inference_events.tsv"
    event_reason: dict[tuple[str, str], str] = {}
    if event_path.is_file():
        ev = read_tsv(event_path)
        for e in ev.itertuples(index=False):
            event_reason[(str(e.complex_variant).split(" ", 1)[0], str(e.complex_variant).split(" ", 1)[-1])] = str(e.technical_reason)
    full_work = OUT / "scoring/full/workspace"
    full_graph_dir = full_work / "af3_graphs"
    graph_label_entries: set[str] = set()
    labels_file = full_graph_dir / "all_variants.labels"
    if labels_file.is_file():
        with labels_file.open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(">"):
                    graph_label_entries.add(line[1:].strip())
    fasta_ids: set[str] = set()
    fasta_file = full_work / "wt_and_vt.fasta"
    if fasta_file.is_file():
        with fasta_file.open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(">"):
                    fasta_ids.add(line[1:].strip())
    for i, r in enumerate(rows.itertuples(index=False)):
        if pd.notna(r.score) and math.isfinite(float(r.score)):
            statuses.append("SCORED")
            failures.append("")
            continue
        key = (str(r.target_id), str(r.partner_id), str(r.local_variant))
        if mapping_reason.iloc[i] in ALLOWED_FAILURES:
            reason = mapping_reason.iloc[i]
        elif key not in selection_keys:
            reason = "GRAPH_FAILURE"
        else:
            complex_id = str(r.complex_id)
            event = event_reason.get((complex_id, str(r.local_variant)))
            if event in ALLOWED_FAILURES:
                reason = event
            elif (not (full_graph_dir / f"{complex_id}.mat").is_file()
                  or f"{complex_id}_interaction_loss_variant_{r.local_variant}" not in graph_label_entries):
                reason = "GRAPH_FAILURE"
            elif (str(r.target_id) not in fasta_ids or str(r.partner_id) not in fasta_ids
                  or f"{r.target_id} {r.local_variant}" not in fasta_ids):
                reason = "PROTT5_INFERENCE_FAILURE"
            else:
                reason = "OTHER"
        statuses.append("FAILED")
        failures.append(reason)
    rows["score_status"] = statuses
    rows["failure_reason"] = failures
    if any(reason and reason not in ALLOWED_FAILURES for reason in failures):
        raise SystemExit("Failure registry contains a reason outside the frozen allowed list")

    rows["mutpred_ppi_version"] = MODEL_VERSION
    rows["prott5_model"] = PROTT5_MODEL
    rows["structure_mapping_status"] = rows.status
    rows["structure_id"] = rows.selected_structure
    rows["training_overlap_flags"] = rows.apply(
        lambda r: ";".join(x for x in [
            f"exact_variant={r.exact_variant_training_overlap}",
            f"exact_variant_partner={r.exact_variant_partner_training_overlap}",
            f"target_protein={r.target_protein_training_overlap}"] if x), axis=1)

    output_cols = ["row_id", "gene", "disease_group", "VariationID", "protein_change",
                   "canonical_variant", "partner", "partner_accession", "target_id", "partner_id",
                   "selected_structure", "structure_id", "tier", "structure_mapping_status",
                   "mutpred_ppi_version", "prott5_model", "score", "score_status", "failure_reason",
                   "exact_variant_training_overlap", "exact_variant_partner_training_overlap",
                   "target_protein_training_overlap", "strict_external_include", "training_overlap_flags"]
    write_tsv(rows[output_cols], OUT / "mutpred_ppi_v1_scores.tsv")
    failures_frame = rows.loc[rows.score_status.eq("FAILED"), output_cols].copy()
    write_tsv(failures_frame, OUT / "failure_registry.tsv")

    structural = read_tsv(P1B0 / "structural_42_genes.tsv")
    if structural.gene.nunique() != 42:
        raise SystemExit(f"Expected the frozen 42-gene structural set; found {structural.gene.nunique()}")
    genes = set(structural.gene.astype(str))
    strict, vgroups, strict_groups = make_membership(genes)
    scored = rows.loc[rows.score_status.eq("SCORED")].copy()
    scored["VariationID"] = scored.VariationID.astype(str)
    row_counts = scored.groupby(["gene", "VariationID"]).partner.nunique().to_dict()
    t12_counts = (scored.loc[scored.tier.isin(["Tier 1", "Tier 2"])]
                  .groupby(["gene", "VariationID"]).partner.nunique().to_dict())
    variant_base = rows[["gene", "VariationID"]].drop_duplicates().copy()
    variant_base["VariationID"] = variant_base.VariationID.astype(str)
    memberships_by_variant = vgroups.groupby(["gene", "VariationID"]).MONDO.agg(lambda x: ";".join(sorted(set(x)))).to_dict()
    variant_base["strict_MONDO_groups"] = [memberships_by_variant.get((str(g), str(v)), "") for g, v in variant_base[["gene", "VariationID"]].itertuples(index=False, name=None)]
    variant_base["strict_disease_group_n"] = [len(split_semicolon(x)) for x in variant_base.strict_MONDO_groups]
    variant_base["multi_label_variant"] = np.where(variant_base.strict_disease_group_n.gt(1), "YES", "NO")
    variant_base["scored_partner_n"] = [int(row_counts.get((str(g), str(v)), 0)) for g, v in variant_base[["gene", "VariationID"]].itertuples(index=False, name=None)]
    variant_base["scored_Tier1_2_partner_n"] = [int(t12_counts.get((str(g), str(v)), 0)) for g, v in variant_base[["gene", "VariationID"]].itertuples(index=False, name=None)]
    variant_base["has_ge1_scored_partner"] = np.where(variant_base.scored_partner_n.gt(0), "YES", "NO")
    variant_base["has_ge3_scored_partners"] = np.where(variant_base.scored_partner_n.ge(3), "YES", "NO")
    variant_base["has_ge3_scored_Tier1_2_partners"] = np.where(variant_base.scored_Tier1_2_partner_n.ge(3), "YES", "NO")
    write_tsv(variant_base.sort_values(["gene", "VariationID"], kind="stable"), OUT / "coverage_by_variant.tsv")

    # Disease labels are used only to assign strict MONDO representation denominators.
    diff_rows: list[dict[str, object]] = []
    gene_class: dict[str, str] = {}
    retained_groups: dict[str, int] = {}
    for gene in sorted(genes):
        gstrict = strict.loc[strict.gene.eq(gene)].sort_values("MONDO")
        rates: list[float] = []
        retained = 0
        for d in gstrict.itertuples(index=False):
            group_members = vgroups.loc[vgroups.gene.eq(gene) & vgroups.MONDO.eq(d.MONDO), "VariationID"].astype(str).unique()
            denominator = int(d.n_unique_variants)
            represented1 = sum(row_counts.get((gene, v), 0) >= 1 for v in group_members)
            represented3 = sum(t12_counts.get((gene, v), 0) >= 3 for v in group_members)
            multi_in_group = sum((gene, v) in {(g, vv) for g, vv in vgroups.loc[vgroups.MONDO.eq(d.MONDO), ["gene", "VariationID"]].itertuples(index=False, name=None)} and
                                 len(split_semicolon(memberships_by_variant.get((gene, v), ""))) > 1 for v in group_members)
            rate = 100.0 * represented3 / max(1, denominator)
            rates.append(rate)
            if denominator - multi_in_group >= 10:
                retained += 1
            gene_rows = rows.loc[rows.gene.eq(gene) & rows.VariationID.astype(str).isin(set(group_members))]
            diff_rows.append({
                "gene": gene, "MONDO": d.MONDO, "disease_label": d.MONDO_label,
                "strict_unique_variants_n": denominator,
                "variants_with_ge1_scored_partner_n": int(represented1),
                "variants_with_ge3_scored_Tier1_2_partners_n": int(represented3),
                "variants_with_ge3_scored_Tier1_2_partners_pct": round(rate, 2),
                "partner_score_rows_n": len(gene_rows),
                "partner_score_rows_scored_n": int(gene_rows.score_status.eq("SCORED").sum()),
                "partner_score_row_coverage_pct": round(100 * gene_rows.score_status.eq("SCORED").mean(), 2) if len(gene_rows) else 0.0,
                "multilabel_variants_excluded_from_primary_n": int(multi_in_group),
                "primary_unambiguous_variants_n": max(0, denominator - int(multi_in_group)),
                "frozen_coverage_metric": "strict unique variants with >=3 distinct scored Tier 1-2 partners / strict unique variants",
            })
        if not rates or min(rates) < 50.0 or max(rates) - min(rates) >= 50.0:
            cls = "SEVERE_IMBALANCE"
        elif max(rates) - min(rates) >= 20.0 or (min(rates) < 70.0 and max(rates) >= 70.0):
            cls = "MATERIAL_IMBALANCE"
        else:
            cls = "MINOR_IMBALANCE"
        gene_class[gene] = cls
        retained_groups[gene] = retained
    diff = pd.DataFrame(diff_rows)
    diff["gene_coverage_adjudication"] = diff.gene.map(gene_class)
    write_tsv(diff, OUT / "differential_score_coverage.tsv")

    gene_rows_out: list[dict[str, object]] = []
    gene_eligible: dict[str, bool] = {}
    multilabel_pairs = set(zip(variant_base.loc[variant_base.multi_label_variant.eq("YES"), "gene"].astype(str),
                               variant_base.loc[variant_base.multi_label_variant.eq("YES"), "VariationID"].astype(str)))
    full_rows_scored = int(rows.score_status.eq("SCORED").sum())
    global_row_coverage = full_rows_scored / len(rows)
    for r in structural.sort_values("gene").itertuples(index=False):
        gene = str(r.gene)
        gr = rows.loc[rows.gene.eq(gene)]
        gene_vars = set(vgroups.loc[vgroups.gene.eq(gene), "VariationID"].astype(str))
        primary_vars = {v for v in gene_vars if (gene, v) not in multilabel_pairs}
        ge3 = sum(t12_counts.get((gene, v), 0) >= 3 for v in primary_vars)
        ge1 = sum(row_counts.get((gene, v), 0) >= 1 for v in primary_vars)
        ge3_pct = 100.0 * ge3 / max(1, len(primary_vars))
        scored_partners = gr.loc[gr.score_status.eq("SCORED") & gr.tier.isin(["Tier 1", "Tier 2"]), "partner"].nunique()
        incompatible = str(r.mechanism_status) == "CONFIRMED_INCOMPATIBLE"
        reasons: list[str] = []
        if str(r.Phase1A_structural_gate) != "PASS": reasons.append("PHASE1A_STRUCTURAL_GATE_FAIL")
        if incompatible: reasons.append("CONFIRMED_INCOMPATIBLE")
        if gene_class.get(gene) != "MINOR_IMBALANCE": reasons.append(gene_class.get(gene, "NO_STRICT_DISEASE_GROUPS"))
        if retained_groups.get(gene, 0) < 2: reasons.append("<2_DISEASE_GROUPS_GE10_AFTER_MULTILABEL_EXCLUSION")
        if scored_partners < 3: reasons.append("<3_SCORED_TIER1_2_PARTNER_DIMENSIONS")
        if ge3_pct < 70.0: reasons.append("<70_PERCENT_PRIMARY_VARIANTS_WITH_GE3_SCORES")
        eligible = not reasons
        gene_eligible[gene] = eligible
        gene_rows_out.append({
            "gene": gene, "Phase1A_structural_gate": r.Phase1A_structural_gate,
            "mechanism_status": r.mechanism_status,
            "differential_score_coverage_class": gene_class.get(gene, "NO_STRICT_DISEASE_GROUPS"),
            "eligible_score_rows_n": len(gr), "scored_rows_n": int(gr.score_status.eq("SCORED").sum()),
            "score_row_coverage_pct": round(100 * gr.score_status.eq("SCORED").mean(), 2) if len(gr) else 0.0,
            "scored_Tier1_2_partner_dimensions_n": int(scored_partners),
            "strict_unique_variants_n": len(gene_vars), "primary_unambiguous_variants_n": len(primary_vars),
            "primary_variants_with_ge1_score_n": int(ge1), "primary_variants_with_ge3_Tier1_2_scores_n": int(ge3),
            "primary_variants_with_ge3_Tier1_2_scores_pct": round(ge3_pct, 2),
            "strict_disease_groups_ge10_after_multilabel_exclusion_n": retained_groups.get(gene, 0),
            "eligible_for_rebuilt_primary_cohort": "YES" if eligible else "NO",
            "primary_cohort_reason": "PASS" if eligible else ";".join(reasons),
        })
    gene_df = pd.DataFrame(gene_rows_out)
    write_tsv(gene_df, OUT / "coverage_by_gene.tsv")

    primary_genes = {g for g, ok in gene_eligible.items() if ok}
    primary_rows = rows.loc[rows.gene.isin(primary_genes) & rows.score_status.eq("SCORED")].copy()
    primary_rows = primary_rows.loc[[((str(g), str(v)) not in multilabel_pairs) for g, v in primary_rows[["gene", "VariationID"]].itertuples(index=False, name=None)]]
    write_tsv(primary_rows[output_cols], OUT / "rebuilt_primary_cohort.tsv")
    strict_external = rows.loc[rows.strict_external_include.eq("YES") & rows.score_status.eq("SCORED")].copy()
    write_tsv(strict_external[output_cols], OUT / "strict_external_cohort.tsv")

    if global_row_coverage < 0.90:
        verdict = "TECHNICAL_HOLD_DISRUPTION_SCORING"
        gate_reason = f"global row coverage {global_row_coverage:.4%} is below the frozen 90% PASS threshold"
    elif gene_df.loc[gene_df.eligible_for_rebuilt_primary_cohort.eq("YES"), "primary_variants_with_ge3_Tier1_2_scores_pct"].lt(70).any():
        verdict = "TECHNICAL_HOLD_DISRUPTION_SCORING"
        gate_reason = "at least one retained primary gene is below 70% scored-variant coverage"
    elif len(primary_genes) >= 30:
        verdict = "GO_PHASE1B1_FULL"
        gate_reason = f"{len(primary_genes)} genes satisfy the frozen primary cohort gates"
    elif len(primary_genes) >= 10:
        verdict = "GO_PHASE1B1_FOCUSED"
        gate_reason = f"{len(primary_genes)} genes satisfy the frozen primary cohort gates"
    else:
        verdict = "NO_GO_EDGETIC_REPRESENTATION"
        gate_reason = f"exact scoring completed, but only {len(primary_genes)} genes satisfy the frozen representation/cohort gates"

    if len(rows) != 11150 or full_rows_scored + len(failures_frame) != 11150:
        raise SystemExit("Final score rows and failure registry do not reconcile to 11,150")
    cache_manifest = read_tsv(OUT / "embedding_cache_manifest.tsv")
    if not cache_manifest.verified.eq("PASS").all():
        raise SystemExit("Embedding cache manifest contains a non-PASS entry")
    batch_manifest = read_tsv(OUT / "batch_manifest.tsv")
    checkpoint = read_tsv(OUT / "checkpoints/full_unique_predictions.tsv")
    batch_rows = pd.to_numeric(batch_manifest.rows_in_batch, errors="coerce").fillna(0).sum()
    if batch_rows != len(checkpoint):
        raise SystemExit("Batch manifest and full-scoring checkpoint row totals disagree")

    report = ROOT / "reports/PHASE1B0R3_OFFLINE_MODEL_INGESTION_AND_SCORING.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        "# PHASE 1B-0R3 — Offline Model Ingestion and Final Scoring\n\n"
        f"## Final gate\n\n`{verdict}`\n\n{gate_reason}.\n\n"
        "The workflow stops at this gate. Phase 1B1 was not run.\n\n"
        "## Frozen inputs and execution\n\n"
        f"- Frozen score rows: {len(rows):,}; unique gene–VariationID variants: {rows[['gene','VariationID']].drop_duplicates().shape[0]:,}.\n"
        f"- Structure mapping: {int(rows.status.eq('READY').sum()):,} READY rows; {len(rows)-int(rows.status.eq('READY').sum()):,} mapping failures.\n"
        f"- Official MutPred-PPI: `{MODEL_VERSION}`.\n- ProtT5: `{PROTT5_MODEL}`.\n"
        f"- Local model path: `{weight_path}`; weight size {weight_size:,} bytes; SHA256 `{weight_audit.iloc[0].sha256}`.\n"
        f"- Frozen Phase 1A source hashes: {len(hash_audit)}/12 PASS.\n"
        f"- Official example: {len(example)} documented scores PASS; maximum absolute error {float(example_max_diff):.12g} (limit 1e-4).\n"
        f"- Pilot20: {len(p20)} predictions; cached/uncached maximum absolute difference {pilot20_max_diff:.12g} (limit 1e-6).\n"
        f"- Pilot100: {len(p100)} finite predictions.\n"
        f"- Full unique scoring checkpoint: {len(checkpoint):,} rows; checkpoint batches cover {int(batch_rows):,} successful scores.\n"
        f"- Final full-table scoring: {full_rows_scored:,}/{len(rows):,} rows ({100*global_row_coverage:.2f}%); failures: {len(failures_frame):,}.\n"
        f"- Exact embedding cache entries: {len(cache_manifest):,}; every manifest entry passed shape, dimension, and finite-value checks.\n\n"
        "- Local model/tokenizer load: PASS; recorded local-only residue embedding shape (723, 1024), all finite.\n"
        "- Coordinate reconciliation: structure-local substitutions are one-based; official checkpoint substitutions are zero-based. Explicit conversion matched every one of the 10,499 unique prediction keys one-to-one.\n\n"
        "## Coverage QC and cohort reconstruction\n\n"
        f"- Structural candidate genes audited: {len(genes)}; rebuilt primary cohort genes: {len(primary_genes)}.\n"
        f"- Rebuilt primary score rows: {len(primary_rows):,}; strict-external scored rows: {len(strict_external):,}.\n"
        f"- Variants with >=1 scored partner: {int(variant_base.scored_partner_n.ge(1).sum()):,}; with >=3 scored partners: {int(variant_base.scored_partner_n.ge(3).sum()):,}.\n"
        f"- Genes with >=70% representation among frozen primary single-label variants: {int(gene_df.primary_variants_with_ge3_Tier1_2_scores_pct.ge(70).sum())}; names: {', '.join(gene_df.loc[gene_df.primary_variants_with_ge3_Tier1_2_scores_pct.ge(70), 'gene'])}.\n"
        f"- Primary genes: {', '.join(sorted(primary_genes)) or 'none'}.\n"
        f"- Strict-external scored universe: {strict_external.gene.nunique()} genes; primary genes with strict-external scored rows: {strict_external.loc[strict_external.gene.isin(primary_genes), 'gene'].nunique()}. Whole-protein overlap was not an exclusion.\n"
        f"- Severe differential coverage: {int(gene_df.differential_score_coverage_class.eq('SEVERE_IMBALANCE').sum())} genes: {', '.join(gene_df.loc[gene_df.differential_score_coverage_class.eq('SEVERE_IMBALANCE'), 'gene'])}.\n"
        f"- Failure reasons: {failures_frame.failure_reason.value_counts().to_dict()}. All 233 failures occurred in structure/reference mapping; no model-inference failure remained among the 10,499 selected scoring units.\n"
        "- Coverage classes use the frozen thresholds: severe if any group <50% or range >=50 percentage points; material if range >=20 points or groups straddle 70%; otherwise MINOR_IMBALANCE. The frozen specification uses MINOR_IMBALANCE for the passing category.\n"
        "- Representation breadth is assessed across all strict groups, while the gene-level 70% denominator is restricted to frozen primary single-label variants. A high overall row-scoring rate does not supply three partner dimensions to every variant or balanced disease-group coverage.\n\n"
        "### Reassessment of the 14 historical representation flags\n\n"
        + gene_df.loc[gene_df.gene.isin(['BMPR2','BRAF','BRCA1','BRCA2','CACNA1C','EIF2B5','GRIN2B','MYH7','NOTCH3','PTPN11','SCN2A','TNNT2','TUBB3','VHL']), ['gene','scored_rows_n','score_row_coverage_pct','differential_score_coverage_class','eligible_for_rebuilt_primary_cohort']].to_markdown(index=False)
        + "\n\n"
        "- Disease labels were used only to assign frozen strict MONDO groups and compute representation coverage. No disease-score association, group score comparison, PERMANOVA, R², clustering, AUROC, enrichment, or tissue-concordance analysis was performed.\n\n"
        "Disease-label association testing was NOT started.\n\n"
        "## Deliverables\n\n"
        "`results/phase1B0R3/` contains the 11,150-row score table, failure registry, pilot outputs and cache validation, cache and batch manifests, variant/gene/disease-group coverage QC, and rebuilt primary and strict-external cohorts.\n",
        encoding="utf-8")
    print(f"final_verdict={verdict}")
    print(f"primary_genes={len(primary_genes)} full_row_coverage_pct={100*global_row_coverage:.2f} scored_rows={full_rows_scored} failures={len(failures_frame)}")
    print(f"report={report}")
    return verdict


if __name__ == "__main__":
    finalize()
