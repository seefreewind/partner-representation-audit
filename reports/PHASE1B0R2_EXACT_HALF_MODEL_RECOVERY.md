# PHASE 1B-0R2 — Exact MutPred-PPI v1 Recovery

## Verdict

**TECHNICAL_HOLD_DISRUPTION_SCORING**

The exact v1.0 ProtT5 half encoder is frozen as `Rostlab/prot_t5_xl_half_uniref50-enc`. Official Hugging Face TLS failed before HTTP metadata was returned, so the weights, official example validation, pilots, and production inference could not run. This is a technical hold; no biological NO-GO was assessed.

## Frozen inputs and exact software

- MutPred-PPI version: `v1.0.0`; source commit: `a6e2f47e18a3af5dbbc4d8accce8ab177049c56c`.
- Frozen scoring universe: 11,150 variant–partner rows; required count: 11,150.
- Frozen-input audit: 12 files checked, 12 matched.
- Expected ProtT5 weight: `pytorch_model.bin`, approximately 2.42 GB, expected SHA256 `7f51ba885541c7dc569d46b796af57cc7a2ba7945107dced4f19d1b5ec091157`.
- Local half-encoder weight: **not obtained**; local weight bytes: 0; observed SHA256: `NOT_OBTAINED`; expected hash match: **not verifiable**.
- No offline-path-only patch was applied. No model architecture, tokenizer, precision, or weight substitution was made. The earlier full-precision-size target was not used.

## Acquisition and network diagnostics

The direct `curl -Iv https://huggingface.co/` test resolved the host to `198.18.0.32` and ended during the TLS handshake with `SSL_ERROR_SYSCALL` (curl 8.7.1 using SecureTransport / LibreSSL 3.3.6). All upper- and lower-case proxy environment variables checked were unset. The official half-model `resolve/main/config.json` request returned curl HTTP code `000`, zero redirects, and exit code 35. The official `hf download` config metadata request retried five times and ended with Python `SSLEOFError` before response metadata.

No Xet transfer occurred because metadata was unavailable. Direct curl resume/retry was not appropriate because direct HTTPS itself failed. A secondary authorized workstation was not accessible through this session, so no off-machine transfer was performed. Diagnostics are in `results/phase1B0R2/tls_diagnostics.tsv`; the exact-file acquisition ledger is in `results/phase1B0R2/prott5_download_manifest.tsv`.

## Reproduction gate and scoring

- Official v1 example: **not run**; exact ProtT5 files unavailable. Maximum absolute error: not available.
- Pilot-20 and pilot-100: **not run**, as required before example validation.
- Frozen scoring rows: 11,150; successfully scored: 0; row coverage: 0.00%.
- Variants with ≥1 partner score: 0/5,357. Variants with ≥3 partner scores: 0/5,357.
- Structural genes with observed ≥70% disruption coverage: 0/42; the coverage gate is **not evaluable** without valid scores.
- Most common row failure: `PROTT5_TLS_FAILURE` (11,150 rows).
- Differential score-coverage exclusions applied: 0. The 14 pre-scoring representation-imbalance flags remain pending; no score-based exclusion decision was made.
- Primary and strict-external cohorts were not rebuilt. Their counts are `NA`, not zero; the old Phase 1B-0 zero-gene result is not carried forward.
- Disease labels were not used for association or outcome analysis.

## Stop point

Remain at `TECHNICAL_HOLD_DISRUPTION_SCORING`. Resume only after the exact official half-encoder files are obtained and the local weight SHA256 matches the frozen value, followed by successful official example validation within absolute tolerance 1e-4. Then continue with the frozen pilots, full scoring, coverage gates, and de novo cohort reconstruction. Do not enter Phase 1B-1 from this technical-hold state.

Disease-label association testing was NOT started.
