# Phase 0R reproducibility log

- Run timestamp (UTC): 2026-09-28T08:32:18.375696+00:00
- Python: 3.9.6
- lxml: (6, 1, 1, 0)
- Git HEAD before Phase 0R result commit: `ee3a1341611becc3b92434c63c65f5e53b535b9a`; working tree at log generation: `?? metadata/phase0R_archive_integrity.json; ?? metadata/phase0R_exact_disease_counts.json; ?? metadata/phase0R_final_metrics.json; ?? metadata/phase0R_input_manifest.tsv; ?? metadata/phase0R_mechanism_counts.json; ?? metadata/phase0R_ppi_counts.json; ?? metadata/phase0R_rcv_download_checkpoint.json; ?? metadata/phase0R_rcv_parse_checkpoint.json; ?? metadata/phase0R_recovery_gate.json; ?? metadata/phase0R_seed_checkpoint.json; ?? metadata/phase0R_step_manifest.tsv; ?? reports/PHASE0R_FULL_RCV_RECOVERY_REPORT.md; ?? reports/PHASE0R_reproducibility_log.md; ?? results/phase0R/; ?? scripts/run_phase0r.py`.
- Phase 0R script SHA256: `fff4e93d9137dee5dd4f138e62f911ff99113d5e101b071aa736691585cfefb8`.
- Frozen summary SHA256: `186ad2838a138f0bc7a79b1b7b6dde5533b105175d293e1131652889529804c9`.
- Frozen RCV archive SHA256: `93117d24517b1c20d0c4a5e77c2d8f32d546351cb160ea33d5bf11ecff9a7ed3`; official MD5 verification: `True` (expected `0e63d1e8ebd44e486b5e923b6fad8bf0`, observed `0e63d1e8ebd44e486b5e923b6fad8bf0`).
- Seed counts and pre-download checkpoint: `metadata/phase0R_seed_checkpoint.json`.
- RCV parse checkpoints: `metadata/phase0R_rcv_parse_checkpoint.json`; records are streamed with lxml and cached in a SQLite store.
- Input and output provenance: `metadata/phase0R_input_manifest.tsv`, `metadata/phase0R_output_manifest.tsv`, `metadata/phase0R_step_manifest.tsv`, and combined `metadata/phase0R_manifest.tsv`.
- No structure prediction, interface scoring, FoldX, MutPred-PPI, PERMANOVA, classification, tissue matching, or machine learning was run.
