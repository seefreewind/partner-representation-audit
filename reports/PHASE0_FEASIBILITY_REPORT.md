# Phase 0 feasibility report

**PHASE0_VERDICT = `NO_GO_ATLAS`**

ClinVar variant_summary from 2026-09 was used for the broad same-release screen. Exact RCV validation was restricted to the five highest-volume preliminary strict genes plus the prespecified ACTB/ACTG1 pilots. Candidate condition assignments in that seven-gene panel were checked using official RCV XML retrieved live through E-utilities. The full fixed RCV archive transfer was stopped after measured throughput implied a multi-hour download. Genome-wide exact RCV totals are unavailable; all-panel counts below are explicitly labeled preliminary. Parent-child and possible same-spectrum pairs remain in the manual review queue and are excluded from strict eligibility.

## 1. Executive verdict

NO_GO_ATLAS. This decision uses only the prespecified thresholds and available Phase 0 metadata. Phase 1 was not started.

## 2. Exact numbers

1. ClinVar P/LP missense variant–condition records: 14,703 among candidate RCV records before protein/MONDO filters; 5,381 pass the primary candidate filters. The same-month whole-release pre-screen found 290,618 potentially eligible summary rows; that count is not RCV-verified.
2. Unique variants: 145,300 unique gene–protein substitutions in the whole-release preliminary screen; 3,265 in filtered candidate RCV records; 0 partner-ready substitutions in the final high-quality set.
3. Protein-mapped variants: 14,703 candidate RCV P/LP germline missense records (100.0% of candidate P/LP missense records).
4. MONDO-mapped variants: 5,381 candidate RCV records (36.6% of protein-mapped candidate records).
5. Raw pleiotropic genes: 3,090 genes in the whole-release preliminary screen with at least two mapped MONDO disease groups; this is a broad-screen value, not RCV-confirmed.
6. Genes with ≥2 normalized distinct diseases: 7 RCV-confirmed within the seven-gene panel; 2,191 were flagged by the whole-release summary pre-screen. 87 additional panel pairs require manual review.
7. Expanded genes (≥5 variants/disease): 5 RCV-confirmed within the panel; 1,452 whole-release summary-screen candidates.
8. Strict genes (≥10 variants/disease and ≥20 total): 5 RCV-confirmed within the panel; 965 whole-release summary-screen candidates.
9. Mechanism-compatible strict genes: 0.
10. Genes with ≥3 credible PPI partners: 6.
11. Final high-quality analyzable genes: 0.
12. Final usable variants: 0.
13. Median variants/disease: 4.0.
14. Median credible partners/gene: 160.0.
15. Partner-analysis-ready coverage: NA (zero denominator) across 0 variants in mechanism-compatible strict groups.
16. ACTB pilot ready: NO.
17. ACTG1 pilot ready: NO.

## Threshold screen

The first six rows are the whole-release ClinVar summary pre-screen; the next six are the selected seven-gene RCV panel. They are separate evidence stages.

| Minimum variants/disease | Genes | Disease groups | Unique substitutions | Evidence stage |
|---:|---:|---:|---:|---|
| 2 | 2191 | 7774 | 94299 | whole-release ClinVar summary; preliminary screen |
| 3 | 1840 | 6219 | 89117 | whole-release ClinVar summary; preliminary screen |
| 5 | 1452 | 4745 | 83290 | whole-release ClinVar summary; preliminary screen |
| 10 | 965 | 3057 | 73762 | whole-release ClinVar summary; preliminary screen |
| 15 | 742 | 2305 | 67656 | whole-release ClinVar summary; preliminary screen |
| 20 | 573 | 1774 | 62249 | whole-release ClinVar summary; preliminary screen |
| 2 | 6 | 52 | 3255 | selected RCV-validation panel |
| 3 | 6 | 46 | 3250 | selected RCV-validation panel |
| 5 | 5 | 39 | 3236 | selected RCV-validation panel |
| 10 | 5 | 29 | 3202 | selected RCV-validation panel |
| 15 | 5 | 21 | 3166 | selected RCV-validation panel |
| 20 | 5 | 18 | 3146 | selected RCV-validation panel |

## 3. Top candidate genes

The RCV validation panel was selected using preliminary strict-group variant volume (top five) plus the prespecified ACTB/ACTG1 pilots. Within the panel, final rankings use disease distinction, sample size, mechanism consistency, experimental PPI and PIONEER coverage. Unresolved disease pairs remain flagged.

| Gene | Diseases | Unique variants | Strict diseases | Mechanism | Credible partners | PIONEER partners | Status |
|---|---:|---:|---:|---|---:|---:|---|
| BRCA1 | 14 | 1231 | 8 | uncertain/different | 351 | 257 | screened candidate |
| TTN | 20 | 342 | 7 | uncertain/different | 164 | 82 | screened candidate |
| BRCA2 | 19 | 1818 | 7 | uncertain/different | 113 | 27 | screened candidate |
| FBN1 | 16 | 1328 | 3 | uncertain/different | 21 | 16 | screened candidate |
| USH2A | 7 | 414 | 4 | uncertain/different | 2 | 3 | screened candidate |
| ACTB | 6 | 15 | 0 | uncertain/different | 449 | 165 | screened candidate |
| ACTG1 | 3 | 4 | 0 | uncertain/different | 160 | 47 | screened candidate |

## 4. Empirical feasibility and bottlenecks

- Gene–disease group variant counts: median 4.0; IQR 17.0; range 1–865. Median disease groups per candidate gene: 14.0. Median experimental partners per candidate gene: 160.0.
- Candidate genes with ≥3 experimental partners: 85.7%; with PIONEER coverage: 100.0%; with experimental PPI support: 100.0%.
- Existing-interface structural proxy: NA (zero denominator) of strict mechanism-compatible variants have a PIONEER/PDB interface reference with reference-residue concordance. This is an annotation proxy, not structure validation.
- Largest bottleneck: **data access**.
- No formal statistical power calculation was performed; this is an empirical feasibility audit.

## 5. Phase 1 recommendation

NO_GO_ATLAS. If a focused follow-up is approved, use 0 high-quality genes and 0 partner-ready substitutions, beginning from `clinvar_variant_condition_filtered.tsv`, `disease_pairs_relationship_audit.tsv`, `disease_pairs_needing_manual_review.tsv`, `gene_disease_mechanism_pair_audit.tsv`, `ppi_provenance.tsv`, and `variant_coverage_proxy.tsv`. Phase 1 was not run.

## Data snapshot and limits

- ClinVar: archived monthly snapshot 2026-09-03 (`variant_summary_2026-09.txt.gz`). The whole-release screen flagged 965 possible strict genes; only the seven-gene validation panel (ACTB;ACTG1;BRCA1;BRCA2;FBN1;TTN;USH2A) received exact RCV retrieval. The fixed RCV XML archive is 6.0 GB compressed; measured transfer speed implied a multi-hour download. Raw API XML responses for queried panel records are preserved.
- MONDO: version v2026-09-01.
- Gene2Phenotype: all-panel download dated 2026-09-27.
- GenCC: current export captured on the run date; verify file completeness using the manifest and failure log.
- PPI: BioGRID MV-Physical 5.0.261, candidate-focused IntAct PSICQUIC queries, and PIONEER very-high-confidence human predictions. HuRI access timed out. Predictomes was not queried and is not counted.
- Manual review: MONDO exact-equivalent terms were collapsed; parent-child and possible same-spectrum pairs remain in `disease_pairs_needing_manual_review.tsv`.
- Missing values are `NA`; the report does not substitute inferred numbers for unavailable data.

## Required feasibility tables and figures

Tables: `clinvar_cohort_attrition.tsv`, `candidate_gene_thresholds.tsv`, `top_candidate_genes.tsv`, `gene_disease_mechanism_audit.tsv`, `gene_disease_mechanism_pair_audit.tsv`, `ppi_coverage_by_gene.tsv`, `ppi_provenance.tsv`, `ACTB_ACTG1_feasibility.tsv`, `disease_pairs_relationship_audit.tsv`, and `disease_pairs_needing_manual_review.tsv`.

Figures: `Figure_P0-1_cohort_attrition.png` through `Figure_P0-5_feasibility_heatmap.png`.

## Data-source failures and review queue

RCV API failed batches: 0.

See `metadata/source_failures.tsv` and `results/phase0/rcv_api_failures.tsv` for exact records and errors. The full RCV archive was not completed; HuRI timed out; IntAct was queried per candidate UniProt accession.

## Stop confirmation

Phase 0 artifacts are generated. **Phase 1 was NOT started.** No structure modeling, interface scoring, FoldX, MutPred-PPI, PERMANOVA, or machine learning was run.
