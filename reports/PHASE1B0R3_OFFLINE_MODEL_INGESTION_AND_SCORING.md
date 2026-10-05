# PHASE 1B-0R3 — Offline Model Ingestion and Final Scoring

## Final gate

`NO_GO_EDGETIC_REPRESENTATION`

exact scoring completed, but only 3 genes satisfy the frozen representation/cohort gates.

The workflow stops at this gate. Phase 1B1 was not run.

## Frozen inputs and execution

- Frozen score rows: 11,150; unique gene–VariationID variants: 5,357.
- Structure mapping: 10,917 READY rows; 233 mapping failures.
- Official MutPred-PPI: `MutPred-PPI v1.0.0@a6e2f47e18a3af5dbbc4d8accce8ab177049c56c`.
- ProtT5: `Rostlab/prot_t5_xl_half_uniref50-enc@94a6abc029ae13029317b140b7424e012bf8dfbf`.
- Local model path: `/Volumes/硬盘2/个基因为何导致多种不同疾病——用预测互作界面解释等位异质性/models/prot_t5_xl_half_uniref50-enc/pytorch_model.bin`; weight size 2,416,373,051 bytes; SHA256 `7f51ba885541c7dc569d46b796af57cc7a2ba7945107dced4f19d1b5ec091157`.
- Frozen Phase 1A source hashes: 12/12 PASS.
- Official example: 3 documented scores PASS; maximum absolute error 3.8509464263e-06 (limit 1e-4).
- Pilot20: 20 predictions; cached/uncached maximum absolute difference 0 (limit 1e-6).
- Pilot100: 100 finite predictions.
- Full unique scoring checkpoint: 10,499 rows; checkpoint batches cover 10,499 successful scores.
- Final full-table scoring: 10,917/11,150 rows (97.91%); failures: 233.
- Exact embedding cache entries: 8,455; every manifest entry passed shape, dimension, and finite-value checks.

- Local model/tokenizer load: PASS; recorded local-only residue embedding shape (723, 1024), all finite.
- Coordinate reconciliation: structure-local substitutions are one-based; official checkpoint substitutions are zero-based. Explicit conversion matched every one of the 10,499 unique prediction keys one-to-one.

## Coverage QC and cohort reconstruction

- Structural candidate genes audited: 42; rebuilt primary cohort genes: 3.
- Rebuilt primary score rows: 573; strict-external scored rows: 10,861.
- Variants with >=1 scored partner: 5,207; with >=3 scored partners: 1,235.
- Genes with >=70% representation among frozen primary single-label variants: 9; names: EIF2B5, FGFR2, GNAO1, HBB, KCNQ1, MEN1, TP53, TUBA1A, VHL.
- Primary genes: FGFR2, MEN1, TUBA1A.
- Strict-external scored universe: 73 genes; primary genes with strict-external scored rows: 3. Whole-protein overlap was not an exclusion.
- Severe differential coverage: 29 genes: ACTA1, ATRX, BMPR2, BRCA2, CACNA1C, CFTR, COL3A1, DYNC1H1, EIF2B5, FGFR3, GRIN2B, LMNA, MLH1, MYH7, NF1, NOTCH3, PMP22, PTEN, PTPN11, RAF1, SCN2A, SDHB, SOS1, STAT3, STK11, STXBP1, TNNI3, TNNT2, TUBB3.
- Failure reasons: {'STRUCTURE_RESIDUE_FAILURE': 182, 'VARIANT_REFERENCE_MISMATCH': 27, 'STRUCTURE_CHAIN_FAILURE': 24}. All 233 failures occurred in structure/reference mapping; no model-inference failure remained among the 10,499 selected scoring units.
- Coverage classes use the frozen thresholds: severe if any group <50% or range >=50 percentage points; material if range >=20 points or groups straddle 70%; otherwise MINOR_IMBALANCE. The frozen specification uses MINOR_IMBALANCE for the passing category.
- Representation breadth is assessed across all strict groups, while the gene-level 70% denominator is restricted to frozen primary single-label variants. A high overall row-scoring rate does not supply three partner dimensions to every variant or balanced disease-group coverage.

### Reassessment of the 14 historical representation flags

| gene    |   scored_rows_n |   score_row_coverage_pct | differential_score_coverage_class   | eligible_for_rebuilt_primary_cohort   |
|:--------|----------------:|-------------------------:|:------------------------------------|:--------------------------------------|
| BMPR2   |               0 |                     0    | SEVERE_IMBALANCE                    | NO                                    |
| BRAF    |             153 |                    99.35 | MATERIAL_IMBALANCE                  | NO                                    |
| BRCA1   |             716 |                   100    | MATERIAL_IMBALANCE                  | NO                                    |
| BRCA2   |               2 |                   100    | SEVERE_IMBALANCE                    | NO                                    |
| CACNA1C |               0 |                     0    | SEVERE_IMBALANCE                    | NO                                    |
| EIF2B5  |             115 |                    83.33 | SEVERE_IMBALANCE                    | NO                                    |
| GRIN2B  |              78 |                   100    | SEVERE_IMBALANCE                    | NO                                    |
| MYH7    |             322 |                    83.2  | SEVERE_IMBALANCE                    | NO                                    |
| NOTCH3  |               0 |                     0    | SEVERE_IMBALANCE                    | NO                                    |
| PTPN11  |              87 |                   100    | SEVERE_IMBALANCE                    | NO                                    |
| SCN2A   |              96 |                   100    | SEVERE_IMBALANCE                    | NO                                    |
| TNNT2   |              33 |                   100    | SEVERE_IMBALANCE                    | NO                                    |
| TUBB3   |              52 |                   100    | SEVERE_IMBALANCE                    | NO                                    |
| VHL     |            1980 |                    99.8  | MATERIAL_IMBALANCE                  | NO                                    |

- Disease labels were used only to assign frozen strict MONDO groups and compute representation coverage. No disease-score association, group score comparison, PERMANOVA, R², clustering, AUROC, enrichment, or tissue-concordance analysis was performed.

Disease-label association testing was NOT started.

## Deliverables

`results/phase1B0R3/` contains the 11,150-row score table, failure registry, pilot outputs and cache validation, cache and batch manifests, variant/gene/disease-group coverage QC, and rebuilt primary and strict-external cohorts.

## Final delivery audit

Completed 2026-10-02 (Asia/Shanghai). All prescribed files exist and are nonempty; score/failure counts, unique prediction keys, batch totals and strict-external exclusions reconcile. All actual local model files were re-read and their SHA256s match the verified model audit. Evidence: `logs/phase1B0R3/final_delivery_audit.log`.
