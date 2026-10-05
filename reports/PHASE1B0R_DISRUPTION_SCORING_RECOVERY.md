# PHASE 1B-0R — MutPred-PPI / ProtT5 Scoring Recovery

> Historical record superseded by [PHASE 1B-0R2 exact half-model recovery](PHASE1B0R2_EXACT_HALF_MODEL_RECOVERY.md). Its full ProtT5 target and blocker description reflect the earlier, superseded instruction and are not the current model specification or project status.

## Verdict

**TECHNICAL_HOLD_DISRUPTION_SCORING**

The required ProtT5 full PyTorch checkpoint was not obtained. No MutPred-PPI inference call, official example scoring, pilot, or production batch was started. The failure is technical and does not support a biological or structural NO-GO.

## Frozen-input audit

- Phase 1A freeze manifest: 12 files checked; 12 matched; no `INPUT_FREEZE_VIOLATION`.
- Frozen variant–partner scoring universe: 11,150 rows (required 11,150); unchanged.
- The frozen training-overlap flags were carried into the row-level recovery status file. No rows were removed.

## Version and model artifacts

- Primary MutPred-PPI version: **v1.0.0**.
- Exact source commit: `a6e2f47e18a3af5dbbc4d8accce8ab177049c56c` (official `v1.0.0` tag).
- Ten tagged v1.0 ensemble checkpoints and their scaler are present under `data/software/mutpred-ppi-v1.0.0/models/`; SHA256 records are in `results/phase1B0R/mutpred_v1_model_checksums.tsv`.
- Requested ProtT5 model: `Rostlab/prot_t5_xl_uniref50`. The official model page reports `pytorch_model.bin` SHA256 `c06141d42e93c12b5f6d95c701952156bd4079661f5f2d981d0d2036ba96dae9` and size about 11.3 GB ([Hugging Face model file](https://huggingface.co/Rostlab/prot_t5_xl_uniref50/blob/main/pytorch_model.bin)). No local weight file was obtained.
- Source expectation conflict recorded: the tagged v1.0 loader defaults to `Rostlab/prot_t5_xl_half_uniref50-enc`; this recovery followed the explicit full-model target and did not substitute the half encoder. The full model therefore requires end-to-end validation against the tagged example before scoring.

## Download and environment recovery

- Official `hf download` standard path was attempted for `config.json`, `spiece.model`, tokenizer files, and `pytorch_model.bin`.
- Failure occurred before metadata resolution: `huggingface.co` returned `SSLEOFError` / `SSL_ERROR_SYSCALL`; no Xet file transfer began. The external-DNS HTTPS diagnostic also failed before returning a usable official address. Xet fallback was not appropriate because the failure preceded Xet handoff.
- Download log: `logs/phase1B0R_prott5_download.log`.
- `hf --help` succeeds; the available CLI uses `huggingface_hub 0.36.2`. A dedicated Python 3.10.20 environment installation was started and interrupted during external-volume I/O before a verified environment was available. Exact state is in `dependency_manifest.tsv` and `logs/phase1B0R_environment_install.log`.
- Storage was sufficient for the requested download; the failure was endpoint connectivity, not disk capacity.

## Scoring and coverage status

- Successfully scored: **0 / 11,150**; row coverage **0.00%**.
- Variants with at least one partner score: **0**; with at least three: **0**.
- Structural genes with observed ≥70% disruption representation: **0** because no score was produced; this is a technical missingness result, not an intrinsic representation finding.
- Most common row status reason: `PROTT5_FAILURE` (11,150 rows). Each input row remains present in `mutpred_ppi_v1_scores.tsv` and `failure_registry.tsv` with an explicit status.
- Official example, pilot-20, pilot-100, and production batching were **not run** because the required model was unavailable.
- The 14 pre-scoring representation-imbalance genes remain unresolved. Differential score coverage was not reassessed without valid scores.
- Primary and strict-external cohorts were **not rebuilt**. The old Phase 1B-0 zero-gene cohort result was not carried forward as a new result.
- Coverage summaries are technical completeness outputs only; no disease-specific score comparison or disease-label association test was run.

## Required output inventory

`results/phase1B0R/` contains the requested manifests, row-level status table, failure registry, coverage tables, blocked validation/pilot records, and a non-rebuilt cohort status record. The ProtT5 embedding manifest is empty because no sequence was embedded.

## Stop point

Do not start Phase 1B-1. Resume Phase 1B-0R only after the exact full ProtT5 model can be retrieved from the official host, loaded locally, and the official v1.0 example passes the frozen numerical tolerance. Disease-label association testing was NOT started.
