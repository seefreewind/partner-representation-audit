# Phase 1B-0 — partner-disruption matrix and analysis freeze

Generated: 2026-09-28T13:46:24.403138+00:00 (UTC)  
Verdict: **TECHNICAL_HOLD_DISRUPTION_SCORING**  
Phase 1A input TSVs frozen: 12; post-run hashes verified unchanged.  
Disease-label association testing: **not started**.

## Status

Phase 1A's 42 structural-pass genes and exact 11,150 experimental-coordinate variant–partner rows were reconstructed. The fixed MutPred-PPI current checkpoint and scaler are identified and hashed. The official inference dependency set is installed. Scoring remains blocked because the required ProtT5 encoder weights at `Rostlab/prot_t5_xl_half_uniref50-enc` could not be retrieved from its official Hugging Face host after TLS handshake failures. Consequently, **0/11150** rows have scores (0.00%). No structural, pathogenicity, or interface proxy was substituted. All score rows remain `NA` with a reason.

The required score gate is not met. The primary disruption cohort A is therefore empty at this freeze. Cohort B remains the pre-outcome structural sensitivity cohort (**8 genes**); cohort C is the separable-interface subset of A (**0 genes**). No disease-label statistic was computed.

## Frozen model

- Primary model: `MutPred-PPI current repository checkpoint (arch=current)` (`arch=current`).
- Release: current all-data checkpoint, not the v1.0 ensemble.
- Repository commit: `01a53eaf5661382a97a9b05a0b41e7f699cd0974`.
- Checkpoint and scaler SHA256 values are in `metadata/PHASE1B0_MUTPRED_VERSION_FREEZE.tsv`.
- Reference: Stewart et al., RECOMB 2026, DOI `10.64898/2025.12.20.695738`.
- Sequence encoder: `Rostlab/prot_t5_xl_half_uniref50-enc`; exact resolved revision/checksum cannot be frozen while official weights are inaccessible.
- Public training audit includes 12,529 rows; unpublished VarChAMP training rows remain unauditable.

## Cohort and differential coverage

- Structural starting cohort: 42 genes in `results/phase1B0/structural_42_genes.tsv`.
- The exact 11,150-row Phase 1A PDB scoring universe spans 73 genes; 32 of the 42 structural-pass genes have at least one exact PDB scoring row. This preserves the frozen Phase 1A universe while keeping primary eligibility restricted to the 42-gene cohort.
- Cohort A (primary disruption): 0 genes; gated by ≥3 scored Tier 1–2 partner dimensions and ≥70% single-label primary-variant score completeness, in addition to the frozen structural/mechanism/differential rules.
- Cohort B (structural sensitivity): 8 genes.
- Cohort C (separable interfaces within A): 0 genes.
- Material/severe differential-coverage sensitivity-only exclusions: BMPR2, BRAF, BRCA1, BRCA2, CACNA1C, EIF2B5, GRIN2B, MYH7, NOTCH3, PTPN11, SCN2A, TNNT2, TUBB3, VHL.
- Thresholds and each strict disease group's counts/coverage are recorded in `differential_coverage_adjudication.tsv`.

## Scoring universe and overlap

- Eligible observations: 11150 unique `(gene, VariationID, partner)` rows; 5,357 unique variant keys in the full universe.
- Successful scores: 0; row coverage: 0.00%.
- Variants with ≥3 partner scores: 0.
- Genes with ≥70% of single-label primary variants having ≥3 disruption scores: 0.
- Exact public training variant overlaps: 28 unique `(gene, VariationID)` variants across 62 variant–partner observation rows.
- Exact public training variant–partner overlaps: 21 unique `(gene, VariationID, partner)` pairs across 21 rows; pairs without a uniquely mapped partner accession remain `UNKNOWN`.
- Genes with at least one strict-external observation row: 73; the rows are not disruption-scoreable while score coverage is 0%.
- `ALL`, `STRICT_EXTERNAL`, and `PROTEIN_EXTERNAL` inclusion flags are frozen in `mutpred_training_overlap.tsv`.

## Baselines

- Linear position built for 4132 unique mapped structural-cohort variants; normalized positions use 1-based canonical residue and frozen UniProt sequence length.
- Domain baseline produced for 4132 variants from locally frozen UniProtKB domain-feature responses; missing annotations remain explicit.
- Generic 3D baseline produced for 4132 variants from partner-independent AlphaFold DB monomer coordinates (pLDDT, isolated-chain SASA, and intrachain Cα density where available); partner identity is excluded.
- Generic 3D residue features were available for 3619/4132 variants; absent structures or residue mappings remain explicit.
- General pathogenicity uses the frozen AlphaMissense 2023 amino-acid substitution resource; 4030/4132 canonical accession/substitution keys matched, and other values remain `NA`.
- Partner-free structural baseline is aggregated to one row per variant and contains only any-interface status/counts, with no partner names.
- Separate S1/S2 structural matrix: 53448 rows across the 42 structural genes. S1 retains Phase 1A PIONEER annotations and experimental PDB occupancy; S2 distances are available only for 7836 experimental-coordinate rows, with PIONEER-only distances left `NA`.

## Hotspots, multi-labels, and FoldX

- Hotspot-dominated genes: none (frozen thresholds: any disease with one residue ≥50% or top three residues ≥70%).
- Multi-label variants: 1637; excluded from the future primary observation set and retained for dedicated sensitivity.
- FoldX was not run. Experimental-coordinate addressability is 2363/4133 strict variants (57.17%) overall and ≥70% for 22/42 genes; these estimates do not establish FoldX compatibility.

## Data-blind exemplar candidates

Ranked by scoreability, disease balance, sample size, Tier 1–2 PPI evidence, interface separability and public training overlap. Since score inference is blocked, these are **technical candidates pending score QC**, not validated exemplars.

```text
 data_blind_rank   gene  potential_ge3_PDB_partner_variant_pct  score_success_pct  disease_balance_min_over_max  Tier1_2_partners interface_overlap_class
               1   TP53                                 100.00                0.0                        0.8642               270   PARTIALLY_OVERLAPPING
               2 TUBA1A                                 100.00                0.0                        0.7826                28    SEPARABLE_INTERFACES
               3 EIF2B5                                 100.00                0.0                        0.7667                12    SEPARABLE_INTERFACES
               4    VHL                                  98.80                0.0                        0.9500                53   PARTIALLY_OVERLAPPING
               5   MEN1                                  95.07                0.0                        0.7339                10   PARTIALLY_OVERLAPPING
               6  KCNQ1                                  86.39                0.0                        0.3190                 5    SEPARABLE_INTERFACES
               7    HBB                                  85.96                0.0                        0.3750                 5   PARTIALLY_OVERLAPPING
               8  FGFR2                                  81.82                0.0                        0.7333                14   PARTIALLY_OVERLAPPING
               9  GNAO1                                  81.08                0.0                        0.3793                14   PARTIALLY_OVERLAPPING
              10   BRAF                                  70.31                0.0                        0.5957                29   PARTIALLY_OVERLAPPING
```

## Frozen analysis plan and stop rule

The complete distance metrics, minimum shared dimensions, PERMANOVA effect size, incremental R² contrasts, constrained permutation scheme, multilabel rule, hotspot control, residue-grouped validation and ClinVar temporal cutoff are specified in `results/phase1B0/analysis_specification.md`. No disease clustering, PERMANOVA, disease prediction/classification, disease-specific enrichment, tissue concordance or VUS phenotype prediction was run.

Phase 1B-0 cannot pass until the required primary version can generate scores for the frozen rows, score coverage is assessed, and baselines are complete. This is a technical hold, not a negative biological result. After those tasks are completed and frozen, stop and wait for Phase 1B-1 authorization.

Disease-label association testing was NOT started.
