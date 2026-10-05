# PHASE 2S — Representation Coverage Benchmark

## Current gate

`TECHNICAL_HOLD_SEQUENCE_MODELS`

The frozen inputs, partner set, prespecified analysis and available structure-dependent baselines have been prepared. Valid sequence inference is blocked by official Hugging Face TLS failures for both required encoding models. No rescue outcome can be inferred from this technical hold. The original edgetic verdict remains `NO_GO_EDGETIC_REPRESENTATION`.

## Frozen benchmark scope

- U1: 42 exact Phase 1A structural-ready genes. U2: 132 of 137 strict genes after the five known-incompatible exclusions.
- Frozen credible gene–partner pairs: 15,942.
- Required distinct canonical accessions: 7,300; existing frozen sequences available: 267; the remainder are explicitly pending sequence retrieval. No ambiguous accession or variant mapping was forced.
- Every consumed frozen source and available cached sequence is listed with actual SHA256 in `metadata/PHASE2S_INPUT_FREEZE_MANIFEST.tsv` and its results copy.

## Methods and technical validation (Table S1)

SPPIDER-seq is pinned to the [manuscript release](https://github.com/aporollo-lab/SPPIDER-seq/tree/v2.0-manuscript), not the retrained v2.2.0 model. The [published article](https://pmc.ncbi.nlm.nih.gov/articles/PMC13142516/) identifies v2.0-manuscript as its archived release. X-PAIR is pinned to the [July 2026 PyPI release 0.1.0](https://pypi.org/project/xpair/0.1.0/), with its bundled pretrained checkpoint.

| method                        | version         | source_revision                          |   size_bytes | sha256                                                           | encoding_model               | technical_validation               |
|:------------------------------|:----------------|:-----------------------------------------|-------------:|:-----------------------------------------------------------------|:-----------------------------|:-----------------------------------|
| SPPIDER_SEQ_PARTNER_INTERFACE | v2.0-manuscript | 9c12e8c826ad4e47b00477cab56f6d58715d8f60 |     32813717 | de00e13dd1edbe7a4eb76e3cbf298ba65abe2ef44cd85d9ce12dff5b5e01c113 | facebook/esm2_t33_650M_UR50D | NOT_RUN_ENCODING_MODEL_TLS_FAILURE |
| SPPIDER_SEQ_PARTNER_INTERFACE | v2.0-manuscript | 9c12e8c826ad4e47b00477cab56f6d58715d8f60 |     32813717 | c359fcdb5751d94d3a1cdc01bb285cb06ff2aade604a992ecdc9e8642252cfe0 | facebook/esm2_t33_650M_UR50D | NOT_RUN_ENCODING_MODEL_TLS_FAILURE |
| XPAIR_PARTNER_INTERFACE       | 0.1.0           | PyPI release 2026-07-27                  |      5691097 | d8c50d03cd6107d2b411d7095a34243ec1b31e6b7cc1d590329c9328c1e8df72 | ElnaggarLab/ankh-large       | NOT_RUN_ENCODING_MODEL_TLS_FAILURE |
| XPAIR_PARTNER_INTERFACE       | 0.1.0           | PyPI release 2026-07-27                  |      5225342 | bbd8b3934ae26f15784a77f4f987c248365d3059c0f8475addf793135541144c | ElnaggarLab/ankh-large       | NOT_RUN_ENCODING_MODEL_TLS_FAILURE |

The two SPPIDER neural-head weights and X-PAIR checkpoint were acquired and hashed. SPPIDER requires `facebook/esm2_t33_650M_UR50D`; X-PAIR checkpoint hyperparameters specify `ankh_large`, requiring `ElnaggarLab/ankh-large`. Both official model API requests failed with `SSLEOFError` in Python requests and SSL handshake errors in curl. No matching local Hugging Face cache was found. Raw acquisition evidence is in `logs/phase2S/encoding_model_acquisition.log`.

SPPIDER_EXAMPLE_PASS and XPAIR_EXAMPLE_PASS are **not established**. Technical success rates are **NA (not run)**, and all sequence representation, rescue, consensus and structural-validation outcomes are **NA**. Empty raw-prediction tables indicate that inference was not run; they do not encode zero biological coverage.

SPPIDER_RESCUE_STATUS = TECHNICAL_HOLD; XPAIR_RESCUE_STATUS = TECHNICAL_HOLD. Neither model has an evaluable FULL, PARTIAL or NO_RESCUE result.

The reproducible environment targets Python 3.12.13, xpair 0.1.0, torch 2.7, pytorch-lightning 2.5.1, torchmetrics 1.7.1, esm 3.2.3 and transformers 4.48.1. An initial transformers 4.48.3 constraint conflicted with esm requirements; this constraint was corrected. The compatible retry then exited after a PyPI zstd request timed out; evidence is recorded in `logs/phase2S/environment_install_compatible.log`. This environment is incomplete. Installed package versions must be verified before any example or study inference.

SPPIDER official code processes the complete sequence through overlapping ESM windows (1024 tokens including BOS/EOS; stride 512), merging query logits and pooling across partner chunks. Receptor-centric query probabilities are prespecified as primary; peptide-centric output is secondary. X-PAIR checkpoint length limits are 50–2000 residues for each input; unsupported inputs must be explicit failures, with no truncation. PPI-existence output must not filter interface probabilities.

## Available coverage baseline (Table S2)

| method                          | universe         |   gene_n |   eligible_variant_partner_rows | available_rows   | row_coverage_pct   | variant_ge1_n   | variant_ge3_n   | genes_ge70_n   | balanced_genes_n   | status                   |   eligible_variant_n | variant_ge1_pct    | variant_ge3_pct    |
|:--------------------------------|:-----------------|---------:|--------------------------------:|:-----------------|:-------------------|:----------------|:----------------|:---------------|:-------------------|:-------------------------|---------------------:|:-------------------|:-------------------|
| MUTPRED_PPI_DISRUPTION          | U1_STRUCTURAL_42 |       42 |                         1164709 | 7640.0           | 0.6559578401128522 | 2250            | 1235            | 9              | 3                  | EVALUATED                |                 4133 | 54.43987418340189  | 29.881442051778368 |
| MUTPRED_PPI_DISRUPTION          | U2_STRICT_132    |      132 |                         1559241 | 10917.0          | 0.7001483414045679 | 5207            | 1235            | 9              | 3                  | EVALUATED                |                13572 | 38.36575302092543  | 9.099616858237548  |
| STRUCTURAL_BINARY_INTERFACE     | U1_STRUCTURAL_42 |       42 |                         1164709 | 53448.0          | 4.588957413396822  | 4132            | 4132            | NA             | NA                 | EVALUATED                |                 4133 | 99.97580450036293  | 99.97580450036293  |
| STRUCTURAL_BINARY_INTERFACE     | U2_STRICT_132    |      132 |                         1559241 | 67092.0          | 4.302862738986469  | 11538           | 6485            | NA             | NA                 | EVALUATED                |                13572 | 85.0132625994695   | 47.78219864426761  |
| STRUCTURAL_INTERFACE_DISTANCE   | U1_STRUCTURAL_42 |       42 |                         1164709 | 7836.0           | 0.6727860778958521 | 2363            | 1246            | NA             | NA                 | EVALUATED                |                 4133 | 57.17396564239051  | 30.14759254778611  |
| STRUCTURAL_INTERFACE_DISTANCE   | U2_STRICT_132    |      132 |                         1559241 | 11150.0          | 0.7150915092663674 | 5357            | 1246            | NA             | NA                 | EVALUATED                |                13572 | 39.470969643383434 | 9.1806660772178    |
| SPPIDER_SEQ_PARTNER_INTERFACE   | U1_STRUCTURAL_42 |       42 |                         1164709 | NA               | NA                 | NA              | NA              | NA             | NA                 | NOT_RUN_MODEL_DEPENDENCY |                 4133 | NA                 | NA                 |
| SPPIDER_SEQ_PARTNER_INTERFACE   | U2_STRICT_132    |      132 |                         1559241 | NA               | NA                 | NA              | NA              | NA             | NA                 | NOT_RUN_MODEL_DEPENDENCY |                13572 | NA                 | NA                 |
| XPAIR_PARTNER_INTERFACE         | U1_STRUCTURAL_42 |       42 |                         1164709 | NA               | NA                 | NA              | NA              | NA             | NA                 | NOT_RUN_MODEL_DEPENDENCY |                 4133 | NA                 | NA                 |
| XPAIR_PARTNER_INTERFACE         | U2_STRICT_132    |      132 |                         1559241 | NA               | NA                 | NA              | NA              | NA             | NA                 | NOT_RUN_MODEL_DEPENDENCY |                13572 | NA                 | NA                 |
| SEQUENCE_CONSENSUS_AVAILABILITY | U1_STRUCTURAL_42 |       42 |                         1164709 | NA               | NA                 | NA              | NA              | NA             | NA                 | NOT_RUN_MODEL_DEPENDENCY |                 4133 | NA                 | NA                 |
| SEQUENCE_CONSENSUS_AVAILABILITY | U2_STRICT_132    |      132 |                         1559241 | NA               | NA                 | NA              | NA              | NA             | NA                 | NOT_RUN_MODEL_DEPENDENCY |                13572 | NA                 | NA                 |
| CONSENSUS_HIGH_SECONDARY        | U1_STRUCTURAL_42 |       42 |                         1164709 | NA               | NA                 | NA              | NA              | NA             | NA                 | NOT_RUN_MODEL_DEPENDENCY |                 4133 | NA                 | NA                 |
| CONSENSUS_HIGH_SECONDARY        | U2_STRICT_132    |      132 |                         1559241 | NA               | NA                 | NA              | NA              | NA             | NA                 | NOT_RUN_MODEL_DEPENDENCY |                13572 | NA                 | NA                 |

The original technical scoring success was 10,917/11,150 rows (97.91%), representing 10,499 unique inference units. Original ≥70% disruption representation covered 9 of the 42 genes, and the rebuilt primary cohort contained FGFR2, MEN1 and TUBA1A. Those historical results remain unchanged.

The table above uses the frozen credible variant–partner universe as a common benchmark denominator. The original ≥70% gate is copied directly from the immutable legacy outputs. Audit found that legacy score-row membership and the full-cohort single-label registry produce different primary denominators. The initial derived cross-method C4/C5 results were withdrawn; they must not be used for rescue conclusions. New structural C4/C5 are NA pending exact membership reconciliation. Details are in `logs/phase2S/denominator_reconciliation.log`. A row lacking a structure remains unavailable for structural methods; this denominator differs from the already structure-selected 11,150-row production scoring denominator. Binary interface membership and finite interface distance have separate availability. Missing distances are not derived from binary labels. Matrices are saved as sparse long coordinates keyed by gene, variant and frozen partner.

## Coverage QC and planned comparisons

`analysis_specification.md` was saved before study-model predictions. The 70% gene threshold and strict disease-group balance thresholds are unchanged. The common sequence-arm partner universe includes all frozen credible PPI partners without requiring solved complexes, with evidence strata retained. Ordinary availability does not depend on a probability cutoff. Secondary CONSENSUS_HIGH requires both probabilities ≥0.5 and remains separate from consensus availability. No mutant-versus-WT deltas are permitted.

Only paired descriptive comparisons are planned; no formal paired hypothesis test was selected. Structural validity would use experimental/PDB labels, protein/complex-grouped resampling, and frozen interface thresholds. Partner specificity and a small deterministic scrambled-partner control remain unevaluated until both examples pass. Technical failure-predictor covariates are preserved, but no most-important predictor has been estimated.

## Deliverables and unavailable analyses

| Table | File | Status |
|---|---|---|
| S1 Method/version manifest | `method_manifest.tsv` | Actual acquired checkpoint hashes; encoding models unavailable |
| S2 Coverage by method | `method_comparison_summary.tsv` | Row/variant structural availability evaluated; new C4/C5 deferred; sequence arms NA |
| S3 Coverage by gene | `coverage_by_gene.tsv` | Same |
| S4 Disease-group balance | `coverage_by_disease_group.tsv` | Coverage QC only |
| S5 Rescued genes | `sequence_rescued_genes.tsv` | Not evaluable |
| S6 Technical failures | `technical_failures_by_method.tsv` | Model-level acquisition errors |
| S7 PDB validation | `structural_validation.tsv` | Not run |

All requested machine-readable output paths are present. Unexecuted sequence prediction outputs contain headers only, with method-level status explained in the manifest and failure table. Figures S1–S6 are deferred because the required cross-method sequence results are unavailable; no zero-valued rescue bars or fabricated validation curves were produced. This is a recoverable technical gate, not a completed coverage comparison.

## Recovery steps

Acquire the two exact encoding models from their official repositories, pin their revisions and verify actual file hashes. Complete and audit the pinned environment. Reconcile full-cohort primary membership against the exact frozen legacy rules before releasing any new C4/C5 endpoint. Obtain missing canonical sequences only for the frozen accession requests, saving new files under this workspace and recording their provenance. Run the official example gates first. Then infer each unique WT protein pair once, cache full-length outputs and extract frozen variant positions. Only after valid sequence inference may rescue counts, phenotype-balanced coverage, structural validation, specificity controls and the planned figures be calculated. Do not begin Phase 2P.

Original NO_GO_EDGETIC_REPRESENTATION was not altered.

Disease-label phenotype association testing was NOT started.
