#!/usr/bin/env python3
"""Disease-blind pilot selection and checkpointed official MutPred-PPI scoring.

This wrapper calls the frozen MutPred-PPI v1.0.0 graph and inference code. It
does not alter model inputs, weights, architecture, scaler, or score calculation.
The only runtime addition is an exact, lossless ProtT5 embedding cache keyed by
the verified model-weight SHA256 and the input-sequence SHA256.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase1B0R3"
UNIVERSE = ROOT / "results/phase1B0/variant_partner_scoring_universe.tsv"
MAPPING = OUT / "row_structure_mapping.tsv"
STRUCTURES = OUT / "structure_inputs/structures"
REPO = ROOT / "data/software/mutpred-ppi-v1.0.0"
SRC = REPO / "src"
GRAPH_BUILDER = SRC / "01_make_contact_graphs_and_fasta.py"
INFERENCE = SRC / "02_run_mutpred-ppi_inference.py"
MODEL_AUDIT = OUT / "model_file_audit.tsv"
WEIGHT_FILE = ROOT / "models/prot_t5_xl_half_uniref50-enc/pytorch_model.bin"
WEIGHT_SHA256 = "7f51ba885541c7dc569d46b796af57cc7a2ba7945107dced4f19d1b5ec091157"
PILOT_SEED = "phase1b0r3-disease-blind-pilot-v1"
SELECTION_DIR = OUT / "pilot_selection"
SCORING_DIR = OUT / "scoring"
CACHE_DIR = OUT / "embedding_cache"
CHECKPOINT_DIR = OUT / "checkpoints"
LOG_DIR = ROOT / "logs/phase1B0R3"
SELECTION_COLUMNS = [
    "row_id", "gene", "VariationID", "variant", "partner",
    "target_accession", "partner_accession", "target_id", "partner_id",
    "local_variant", "selected_structure", "tier", "selection_sha256",
]
SCORE_KEY = ["target_id", "partner_id", "local_variant"]
CACHE_COLUMNS = [
    "model_weight_sha256", "sequence_sha256", "cache_file", "n_residues",
    "embedding_dim", "dtype", "size_bytes", "cache_event", "verified",
]


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def check_model_identity() -> None:
    if not WEIGHT_FILE.is_file() or WEIGHT_FILE.stat().st_size == 0:
        raise SystemExit(f"Verified ProtT5 weight file is absent or empty: {WEIGHT_FILE}")
    if not MODEL_AUDIT.is_file():
        raise SystemExit(f"Missing verified model audit: {MODEL_AUDIT}")
    audit = pd.read_csv(MODEL_AUDIT, sep="\t", dtype=str)
    row = audit.loc[audit.filename.eq("pytorch_model.bin")]
    if len(row) != 1 or row.iloc[0].sha256 != WEIGHT_SHA256:
        raise SystemExit("ProtT5 model audit does not match the frozen weight SHA256")
    if int(row.iloc[0].size_bytes) != WEIGHT_FILE.stat().st_size:
        raise SystemExit("ProtT5 weight size differs from the verified model audit")


def load_universe_and_mapping() -> pd.DataFrame:
    if not MAPPING.is_file():
        raise SystemExit(f"Structure mapping is not complete: {MAPPING} is absent")
    universe = pd.read_csv(UNIVERSE, sep="\t", dtype={"VariationID": str})
    mapping = pd.read_csv(MAPPING, sep="\t", dtype={"VariationID": str})
    if len(universe) != 11150 or universe[["gene", "VariationID"]].drop_duplicates().shape[0] != 5357:
        raise SystemExit("INPUT_FREEZE_VIOLATION: frozen universe dimensions differ")
    if len(mapping) != 11150:
        raise SystemExit(f"Structure mapping must contain 11,150 rows; found {len(mapping)}")
    # Pair-level checkpoints are emitted in first-seen pair order. Occurrence
    # numbers disambiguate duplicate row identities without consulting labels.
    keys = ["gene", "VariationID", "partner", "canonical_variant",
            "canonical_UniProt_accession", "partner_accession"]
    map_keys = ["gene", "VariationID", "partner", "variant",
                "target_accession", "partner_accession"]
    u = universe.copy()
    m = mapping.copy()
    u["__row_order"] = np.arange(len(u), dtype=np.int64)
    u["__occ"] = u.groupby(keys, sort=False, dropna=False).cumcount()
    m["__occ"] = m.groupby(map_keys, sort=False, dropna=False).cumcount()
    m = m.rename(columns={"variant": "canonical_variant",
                         "target_accession": "canonical_UniProt_accession"})
    join_keys = ["gene", "VariationID", "partner", "canonical_variant",
                 "canonical_UniProt_accession", "partner_accession", "__occ"]
    merged = u.merge(m, how="left", on=join_keys, validate="one_to_one", suffixes=("", "_mapped"))
    if len(merged) != 11150 or merged.status.isna().any():
        missing = int(merged.status.isna().sum()) if "status" in merged else -1
        raise SystemExit(f"INPUT_FREEZE_VIOLATION: mapping join did not resolve all rows ({missing} missing)")
    merged["row_id"] = [f"PHASE1B0R3-{i+1:05d}" for i in merged["__row_order"]]
    return merged.sort_values("__row_order", kind="stable").reset_index(drop=True)


def stage_selection_path(stage: str) -> Path:
    return SELECTION_DIR / f"{stage}_selection.tsv"


def stage_variants_path(stage: str) -> Path:
    return SCORING_DIR / stage / "variants.tsv"


def write_variants(selection: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["target_id", "local_variant", "partner_id"]
    unique = selection[cols].drop_duplicates().sort_values(cols, kind="stable")
    unique.to_csv(path, sep="\t", index=False, header=False)


def prepare_selections() -> None:
    SELECTION_DIR.mkdir(parents=True, exist_ok=True)
    SCORING_DIR.mkdir(parents=True, exist_ok=True)
    mapped = load_universe_and_mapping()
    ready = mapped.loc[mapped.status.eq("READY") & mapped.target_id.notna() &
                       mapped.partner_id.notna() & mapped.local_variant.notna()].copy()
    if ready.empty:
        raise SystemExit("No READY rows are available for pilot selection")
    ready["selection_sha256"] = ready.apply(
        lambda r: sha256_text(
            f"{PILOT_SEED}|{r.gene}|{r.VariationID}|{r.partner_accession}|{r.canonical_variant}"
        ), axis=1,
    )
    unique = ready.sort_values(
        ["selection_sha256", "gene", "VariationID", "partner_accession"], kind="stable"
    ).drop_duplicates(SCORE_KEY, keep="first")
    unique = unique.rename(columns={"canonical_variant": "variant"})
    unique["target_accession"] = unique["canonical_UniProt_accession"]
    unique["tier"] = unique["tier"].astype(str)
    selected = unique[SELECTION_COLUMNS].reset_index(drop=True)
    for stage, n in (("pilot20", 20), ("pilot100", 100), ("full", len(selected))):
        frame = selected.head(n).copy()
        if stage != "full" and len(frame) < n:
            raise SystemExit(f"PILOT_SELECTION_INSUFFICIENT: requested {n}, only {len(frame)} unique READY rows")
        frame.to_csv(stage_selection_path(stage), sep="\t", index=False, na_rep="NA")
        write_variants(frame, stage_variants_path(stage))
        print(f"{stage}: selected_unique_scoring_units={len(frame)}")
    print(f"ready_rows={len(ready)} unique_ready_scoring_units={len(selected)}")
    print(f"selection_rule=SHA256({PILOT_SEED}|gene|VariationID|partner_accession|variant); disease labels excluded")


def stage_workspace(stage: str) -> Path:
    return SCORING_DIR / stage / "workspace"


def build_graphs(stage: str, n_jobs: int) -> None:
    selection_path = stage_selection_path(stage)
    variants_path = stage_variants_path(stage)
    if not selection_path.is_file() or not variants_path.is_file():
        raise SystemExit(f"Selection for {stage} is missing; run --phase select first")
    work = stage_workspace(stage)
    work.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"{stage}_graph_builder.log"
    cmd = [sys.executable, str(GRAPH_BUILDER), str(work), str(STRUCTURES), str(variants_path), str(n_jobs)]
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] command={cmd!r}\n")
        log.flush()
        subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    print(f"graph_builder=PASS stage={stage} log={log_path}")


def cache_paths(sequence: str) -> tuple[str, Path]:
    sequence_sha = sha256_text(sequence)
    return sequence_sha, CACHE_DIR / WEIGHT_SHA256 / f"{sequence_sha}.npy"


def append_cache_manifest(sequence_sha: str, path: Path, array: np.ndarray, event: str,
                          seen: set[str]) -> None:
    if sequence_sha in seen:
        return
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = OUT / "embedding_cache_manifest.tsv"
    new = not manifest.exists()
    with manifest.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CACHE_COLUMNS, delimiter="\t", lineterminator="\n")
        if new:
            writer.writeheader()
        writer.writerow({
            "model_weight_sha256": WEIGHT_SHA256,
            "sequence_sha256": sequence_sha,
            "cache_file": str(path.relative_to(ROOT)),
            "n_residues": int(array.shape[0]) if array.ndim else 0,
            "embedding_dim": int(array.shape[-1]) if array.ndim else 0,
            "dtype": str(array.dtype),
            "size_bytes": path.stat().st_size,
            "cache_event": event,
            "verified": "PASS" if array.ndim == 2 and array.shape[1] == 1024 and np.isfinite(array).all() else "FAIL",
        })
        handle.flush()
        os.fsync(handle.fileno())
    seen.add(sequence_sha)


def save_embedding(path: Path, array: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".npy.tmp")
    with tmp.open("wb") as handle:
        np.save(handle, array, allow_pickle=False)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)


def load_cached_embedding(path: Path) -> np.ndarray | None:
    try:
        arr = np.load(path, allow_pickle=False)
    except Exception:
        return None
    if arr.ndim != 2 or arr.shape[1] != 1024 or not np.isfinite(arr).all():
        return None
    return arr


def checkpoint_rows(path: Path) -> dict[tuple[str, str], float]:
    found: dict[tuple[str, str], float] = {}
    if not path.exists():
        return found
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                score = float(row["score"])
                if np.isfinite(score):
                    found[(row["complex_id"], row["variant_zero_based"])] = score
            except Exception:
                continue
    return found


def run_cached_inference(stage: str, mode: str) -> None:
    check_model_identity()
    work = stage_workspace(stage)
    graph_dir = work / "af3_graphs"
    fasta_path = work / "wt_and_vt.fasta"
    if not graph_dir.is_dir() or not fasta_path.is_file():
        raise SystemExit(f"Official graph/FASTA inputs are missing for {stage}; run --phase graph first")
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"{stage}_{mode}_inference.log"
    cache_manifest = OUT / "embedding_cache_manifest.tsv"
    seen_cache = set()
    if cache_manifest.exists():
        try:
            seen_cache = set(pd.read_csv(cache_manifest, sep="\t", dtype=str).sequence_sha256.dropna())
        except Exception:
            seen_cache = set()

    sys.path.insert(0, str(SRC))
    from utils import inference_utils as iu  # imported only after frozen source path is selected

    original_get_t5 = iu.get_t5_emb
    original_get_t5_model = iu.get_T5_model
    cache_stats = {"hits": 0, "misses": 0, "captured": 0}

    def cached_get_t5(key: str, device: Any, fasta_dict: dict, dataset_name: str,
                      t5_model: Any, t5_vocab: Any):
        sequence = fasta_dict.get(key)
        if sequence is None:
            return original_get_t5(key, device, fasta_dict, dataset_name, t5_model, t5_vocab)
        sequence_sha, path = cache_paths(str(sequence))
        if mode == "readwrite":
            arr = load_cached_embedding(path) if path.is_file() else None
            if arr is not None:
                cache_stats["hits"] += 1
                append_cache_manifest(sequence_sha, path, arr, "CACHE_HIT", seen_cache)
                return arr
            if cache_only_model_bypass:
                raise RuntimeError(f"Exact embedding cache miss during cache-only inference: {sequence_sha}")
        arr = original_get_t5(key, device, fasta_dict, dataset_name, t5_model, t5_vocab)
        if arr is None:
            cache_stats["misses"] += 1
            return None
        arr = np.asarray(arr)
        if arr.ndim != 2 or arr.shape[1] != 1024 or not np.isfinite(arr).all():
            raise ValueError(f"Invalid ProtT5 embedding for sequence {sequence_sha[:12]}: {arr.shape}")
        if not path.is_file() or load_cached_embedding(path) is None:
            save_embedding(path, arr)
        cache_stats["misses"] += 1
        event = "UNCACHED_CAPTURE" if mode == "capture" else "COMPUTED_ON_CACHE_MISS"
        append_cache_manifest(sequence_sha, path, arr, event, seen_cache)
        return arr

    iu.get_t5_emb = cached_get_t5
    cache_only_model_bypass = False
    if mode == "readwrite":
        # Reuse the verified embeddings without reloading the 2.4-GB encoder
        # when every FASTA sequence is already cached. A cache miss aborts the
        # affected official prediction; it can never silently use a placeholder.
        fasta_sequences = iu.get_dict_from_fasta(str(fasta_path))
        # The full FASTA can contain >10,000 mutant sequences. Avoid a random
        # filesystem probe (and full NPY read) for every sequence when the
        # manifest already shows that this is a partial cache. Cache writes
        # happen before their manifest rows, so an unmanifested file is only
        # treated as a miss and will be safely reused by cached_get_t5 later.
        verified_sequences: set[str] = set()
        if cache_manifest.is_file():
            try:
                manifest_rows = pd.read_csv(
                    cache_manifest, sep="\t", dtype=str,
                    usecols=["model_weight_sha256", "sequence_sha256", "verified"],
                )
                verified_sequences = set(manifest_rows.loc[
                    manifest_rows.model_weight_sha256.eq(WEIGHT_SHA256)
                    & manifest_rows.verified.eq("PASS"),
                    "sequence_sha256",
                ].dropna())
            except Exception:
                verified_sequences = set()
        fasta_hashes = {sha256_text(str(sequence)) for sequence in fasta_sequences.values()}
        all_cached = fasta_hashes.issubset(verified_sequences)
        if all_cached:
            # Only read every embedding when the manifest says the cache is
            # complete; this validates it before the encoder bypass is allowed.
            all_cached = all(
                load_cached_embedding(CACHE_DIR / WEIGHT_SHA256 / f"{sequence_sha}.npy") is not None
                for sequence_sha in fasta_hashes
            )
        if all_cached:
            iu.get_T5_model = lambda model_dir=None, device=None, **kwargs: (None, None)
            cache_only_model_bypass = True
    is_full = stage == "full"
    checkpoint = CHECKPOINT_DIR / "full_unique_predictions.tsv"
    batch_manifest = OUT / "batch_manifest.tsv"
    completed = checkpoint_rows(checkpoint) if is_full else {}
    if is_full:
        manifested = 0
        if batch_manifest.is_file() and batch_manifest.stat().st_size:
            prior_batches = pd.read_csv(batch_manifest, sep="\t")
            manifested = int(pd.to_numeric(prior_batches.rows_in_batch, errors="coerce").fillna(0).sum())
        recovery_gap = len(completed) - manifested
        if recovery_gap < 0:
            raise SystemExit(f"BATCH_MANIFEST_AHEAD_OF_CHECKPOINT: manifest={manifested}, checkpoint={len(completed)}")
        if recovery_gap:
            is_new_manifest = not batch_manifest.exists() or batch_manifest.stat().st_size == 0
            with batch_manifest.open("a", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                if is_new_manifest:
                    writer.writerow(["batch_id", "rows_in_batch", "cumulative_successes", "first_complex_variant", "last_complex_variant", "timestamp_local"])
                writer.writerow([f"RECOVERED-{len(completed):06d}", recovery_gap, len(completed),
                                 "VALIDATED_CHECKPOINT", "VALIDATED_CHECKPOINT",
                                 time.strftime("%Y-%m-%d %H:%M:%S")])
                handle.flush()
                os.fsync(handle.fileno())
    pending: list[tuple[str, str, float]] = []
    active = {"key": None}
    failures: list[tuple[str, str, str]] = []
    failure_path = OUT / f"{stage}_inference_events.tsv"

    def record_failure(complex_variant: str, reason: str, detail: str) -> None:
        failures.append((complex_variant, reason, detail))
        if not is_full:
            return
        is_new = not failure_path.exists() or failure_path.stat().st_size == 0
        with failure_path.open("a", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            if is_new:
                writer.writerow(["complex_variant", "technical_reason", "detail"])
            writer.writerow([complex_variant, reason, detail])
            handle.flush()
            os.fsync(handle.fileno())
    original_get_variants = iu.get_vts_from_wt
    original_get_complex = iu.get_complex_and_vt_emb
    original_model_predict = iu.model_predict

    if is_full:
        def filter_completed(labels_path: str):
            records = original_get_variants(labels_path)
            return {complex_id: [variant for variant in variants
                                 if (complex_id, variant) not in completed]
                    for complex_id, variants in records.items()}
        iu.get_vts_from_wt = filter_completed

        def capture_current(complex_and_vt: str, *args, **kwargs):
            active["key"] = complex_and_vt.split(" ", 1)
            try:
                result = original_get_complex(complex_and_vt, *args, **kwargs)
            except Exception as exc:
                record_failure(complex_and_vt, "PROTT5_INFERENCE_FAILURE", repr(exc))
                raise
            if result is None or not result or result[0] is None:
                record_failure(complex_and_vt, "PROTT5_INFERENCE_FAILURE", "embedding_or_variant_sequence_unavailable")
            return result
        iu.get_complex_and_vt_emb = capture_current

        def flush_batch() -> None:
            if not pending:
                return
            CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
            is_new = not checkpoint.exists() or checkpoint.stat().st_size == 0
            first_key, last_key = pending[0], pending[-1]
            with checkpoint.open("a", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                if is_new:
                    writer.writerow(["complex_id", "variant_zero_based", "score"])
                for row in pending:
                    writer.writerow(row)
                handle.flush()
                os.fsync(handle.fileno())
            is_new_manifest = not batch_manifest.exists() or batch_manifest.stat().st_size == 0
            with batch_manifest.open("a", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                if is_new_manifest:
                    writer.writerow(["batch_id", "rows_in_batch", "cumulative_successes", "first_complex_variant", "last_complex_variant", "timestamp_local"])
                cumulative = len(completed) + len(pending)
                writer.writerow([f"FULL-{cumulative:06d}", len(pending), cumulative,
                                 " ".join(first_key[:2]), " ".join(last_key[:2]),
                                 time.strftime("%Y-%m-%d %H:%M:%S")])
                handle.flush()
                os.fsync(handle.fileno())
            for complex_id, variant, score in pending:
                completed[(complex_id, variant)] = score
            pending.clear()

        def checkpoint_predict(*args, **kwargs):
            key = active["key"]
            try:
                prediction = original_model_predict(*args, **kwargs)
            except Exception as exc:
                label = "OUT_OF_MEMORY" if "memory" in str(exc).lower() else "MUTPRED_INFERENCE_FAILURE"
                record_failure(" ".join(key) if key else "UNKNOWN", label, repr(exc))
                raise
            if key and len(key) == 2:
                pending.append((key[0], key[1], float(prediction)))
                if len(pending) >= 25:
                    flush_batch()
            return prediction
        iu.model_predict = checkpoint_predict

    try:
        # The official inference function writes its result table but does not
        # create the directory. Keep both its output and temporary FASTA/HDF5
        # files inside this project-owned stage workspace.
        result_dir = work / "results"
        result_dir.mkdir(parents=True, exist_ok=True)
        temp_dir = OUT / "runtime_tmp"
        temp_dir.mkdir(parents=True, exist_ok=True)
        tempfile.tempdir = str(temp_dir)
        with log_path.open("a", encoding="utf-8") as log, open(os.devnull, "w") as devnull:
            stdout, stderr = sys.stdout, sys.stderr
            sys.stdout = log
            sys.stderr = log
            try:
                log.write(f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] stage={stage} mode={mode} checkpointed={is_full} completed_before={len(completed)}\n")
                log.write(f"cache_only_model_bypass={cache_only_model_bypass}\n")
                log.flush()
                iu.run_inference_on_dataset("cpu", str(work), str(graph_dir), str(fasta_path), str(result_dir))
            finally:
                sys.stdout, sys.stderr = stdout, stderr
        if is_full:
            flush_batch()
    finally:
        iu.get_t5_emb = original_get_t5
        iu.get_T5_model = original_get_t5_model
        iu.get_vts_from_wt = original_get_variants
        iu.get_complex_and_vt_emb = original_get_complex
        iu.model_predict = original_model_predict

    result_file = work / "results/MutPred-PPI_preds.tsv"
    if not result_file.is_file():
        raise SystemExit(f"Official inference did not produce its result file: {result_file}")
    saved = OUT / f"{stage}_{'uncached' if mode == 'capture' else 'cached'}.tsv"
    shutil.copy2(result_file, saved)
    if is_full:
        merged = checkpoint_rows(checkpoint)
        out_rows = [(cid, variant, score) for (cid, variant), score in sorted(merged.items())]
        pred_path = OUT / "full_unique_predictions.tsv"
        with pred_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["complex_id", "variant_zero_based", "score"])
            writer.writerows(out_rows)
    if failures and not is_full:
        with failure_path.open("a", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            if failure_path.stat().st_size == 0:
                writer.writerow(["complex_variant", "technical_reason", "detail"])
            writer.writerows(failures)
    print(f"inference=PASS stage={stage} mode={mode} predictions={sum(1 for _ in open(saved, encoding='utf-8'))-1} cache_hits={cache_stats['hits']} cache_computed={cache_stats['misses']} output={saved}")


def compare_pilot20() -> None:
    uncached = OUT / "pilot20_uncached.tsv"
    cached = OUT / "pilot20_cached.tsv"
    if not uncached.is_file() or not cached.is_file():
        raise SystemExit("Run pilot20 --phase inference with capture and readwrite modes first")
    a = pd.read_csv(uncached, sep="\t", dtype={"variant": str})
    b = pd.read_csv(cached, sep="\t", dtype={"variant": str})
    keys = ["complex_id", "variant"]
    joined = a.merge(b, on=keys, how="outer", suffixes=("_uncached", "_cached"), indicator=True)
    if len(joined) != 20 or not joined._merge.eq("both").all():
        raise SystemExit(f"PILOT20_CACHE_VALIDATION_FAIL: expected 20 matching predictions, found {len(joined)}")
    joined["abs_difference"] = (joined.score_uncached - joined.score_cached).abs()
    max_diff = float(joined.abs_difference.max())
    joined.to_csv(OUT / "pilot20_cache_validation.tsv", sep="\t", index=False)
    if not np.isfinite(joined[["score_uncached", "score_cached", "abs_difference"]].to_numpy()).all() or max_diff > 1e-6:
        raise SystemExit(f"PILOT20_CACHE_VALIDATION_FAIL: maximum absolute score difference={max_diff}")
    print(f"pilot20_cached_vs_uncached=PASS n=20 max_abs_difference={max_diff:.12g}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", required=True,
                        choices=["select", "graph", "inference", "compare-pilot20"])
    parser.add_argument("--stage", choices=["pilot20", "pilot100", "full"])
    parser.add_argument("--cache-mode", choices=["capture", "readwrite"])
    parser.add_argument("--n-jobs", type=int, default=4)
    args = parser.parse_args()
    if args.phase == "select":
        prepare_selections()
    elif args.phase == "graph":
        if not args.stage:
            parser.error("--stage is required with --phase graph")
        build_graphs(args.stage, args.n_jobs)
    elif args.phase == "inference":
        if not args.stage or not args.cache_mode:
            parser.error("--stage and --cache-mode are required with --phase inference")
        if args.stage != "pilot20" and args.cache_mode == "capture":
            parser.error("capture mode is reserved for the uncached pilot20 reference run")
        if args.stage == "full":
            CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
            lock_path = CHECKPOINT_DIR / "full_scoring.lock"
            with lock_path.open("w") as lock:
                try:
                    fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    raise SystemExit("FULL_SCORING_ALREADY_RUNNING: checkpoint lock is held")
                run_cached_inference(args.stage, args.cache_mode)
        else:
            run_cached_inference(args.stage, args.cache_mode)
    else:
        compare_pilot20()


if __name__ == "__main__":
    main()
