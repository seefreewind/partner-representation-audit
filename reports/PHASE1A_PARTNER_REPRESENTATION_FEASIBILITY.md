# Phase 1A — partner-resolved representation feasibility

Run date (UTC): 2026-09-28T10:38:41.224161+00:00

## Verdict

**GO_PHASE1B_FULL**

The frozen Phase 0R baseline was reused without revision. The true Phase 1A intersection is **128 genes**: 137 strict genes, excluding 5 explicitly CONFIRMED_INCOMPATIBLE genes, followed by the ≥3 credible-partner criterion. Phase 1A is a representation and missingness audit only. No disease-label association, classifier, disease-specific enrichment, UMAP interpretation, or phenotype prediction was run.

## Cohort and mechanism

- Strict starting genes: 137; strict ClinVar disease groups: 353.
- Explicitly incompatible genes excluded from primary: 5 (BEST1, COL1A2, COL2A1, FH, OPA1).
- Primary cohort: 128 genes; 133 of all strict genes have ≥3 credible partners before mechanism exclusion.
- Partner-count bands among non-incompatible genes: {"exactly_2": 2, "3_to_5": 8, "6_to_10": 11, "over_10": 109, "under_2": 2}. All-strict bands: {"exactly_2": 2, "3_to_5": 8, "6_to_10": 12, "over_10": 113, "under_2": 2}.
- Pairwise mechanism states remain in cohort_intersection.tsv; UNKNOWN_COMPATIBILITY has not been recoded as compatible.

## PPI and structure reuse

Primary-cohort credible edges by frozen hierarchy: Tier 1 **199**, Tier 2 **2396**, Tier 3 **11378**, Tier 4 **1969**. Thus 89 genes have ≥3 Tier 1–2 partners. Tier 4 remains explicitly prediction-only.

RCSB pairwise queries covered 245 edges from the subset already tagged source=PDB in the frozen very-high-confidence PIONEER table; the coordinate audit reviewed up to three lexically ordered search hits per edge and used RCSB-listed biological assemblies. 314 unique PDB entry IDs across 418 biological assemblies yielded mapped target-partner contacts. Within each assembly, mapped partner-chain copies were combined and each target residue was assigned its nearest heavy-atom distance to any mapped partner copy. `PDB_target_chain_observation_n` records target-chain observations; assembly IDs are retained per edge. The remaining pairs outside that query set are marked unqueried, not as having no structure. Interface membership uses a frozen **5.0 Å heavy-atom** cutoff for experimental coordinates; PIONEER interface residues are separately labeled as very-high-confidence predictions. New complex modelling and FoldX were not run.

## Protein mapping and representation

- Unique ClinVar variants across primary genes: 13089; EXACT_CANONICAL + SAFE_ISOFORM_REMAP: 13076 (99.90%).
- Variants with ≥3 Tier 1–2 partner-specific interface representations: 4132 (31.57%).
- Variants with ≥1 observed experimental PDB residue mapping suitable for a future interaction-specific model input: 5357 (40.93%).
- MutPred-PPI scores generated: 0/11150 structurally addressable variant-partner rows (0.00%). The score matrix is header-only; structural addressability is not treated as a disruption score.
- Genes passing the ≥70% protein-mapping and ≥70% ≥3-partner-representation thresholds: 42.
- Genes below the protein-mapping threshold: 0; genes below ≥3 credible partners within the primary cohort: 0; genes below the Tier 1–2 structural representation threshold: 86.
- Differential coverage flagged by the frozen rule: 7 genes: BMPR2, CDKN2A, CYP21A2, EIF2B5, GJB2, LIPA, PAH.

## MutPred-PPI

The official MutPred-PPI repository was pinned at commit 01a53eaf5661382a97a9b05a0b41e7f699cd0974; the current repository checkpoint is the primary generation. Its published v1.0 ensemble is a separate sensitivity generation and scores are not comparable. The training overlap table audits recoverable public rows; the repository documents an unpublished VarChAMP component, so an absent match means no known overlap in the public audit, not proof of no overlap in all training data. MutPred-PPI inference was not run in this pass: the 11150 residue-addressable variant-partner rows remain candidates, and each requires a matching mmCIF-derived contact graph and ProtT5 embedding before the frozen current checkpoint can score it. The score matrix is header-only; no score is imputed from interface occupancy.

Training rows recovered: 12529 across fragoza_only_mapped090826, sahni_fragoza_mapped090826, sahni_only_mapped090826. Public training audit files: data/interim/phase1A_mutpred_training/datasets/training_eval/fragoza_only_mapped090826_rows.csv.gz, data/interim/phase1A_mutpred_training/datasets/training_eval/sahni_fragoza_mapped090826_rows.csv.gz, data/interim/phase1A_mutpred_training/datasets/training_eval/sahni_only_mapped090826_rows.csv.gz.

## Interface overlap and future exemplars

Overlap thresholds were frozen before any phenotype analysis: Jaccard ≥0.50 is highly overlapping; Jaccard ≤0.10 is separable if both interfaces have ≥5 mapped residues; intermediate values are partial. Gene-level labels require ≥3 comparisons and ≥80% of pairwise comparisons in the corresponding class. Separable genes: 19.

Top 10 data-blind candidates use only variant count/balance, Tier 1–2 PPI quality, mapping coverage, ≥3-partner structural coverage, and interface overlap:

1. **FGFR3** — 32 variants across 2 diseases; 305 Tier 1–2 partners; mapping 100.00%; ≥3 Tier 1–2 interface representation 100.00%; PARTIALLY_OVERLAPPING.
2. **TP53** — 94 variants across 2 diseases; 270 Tier 1–2 partners; mapping 100.00%; ≥3 Tier 1–2 interface representation 100.00%; PARTIALLY_OVERLAPPING.
3. **MLH1** — 152 variants across 4 diseases; 144 Tier 1–2 partners; mapping 100.00%; ≥3 Tier 1–2 interface representation 100.00%; PARTIALLY_OVERLAPPING.
4. **BRCA1** — 201 variants across 4 diseases; 136 Tier 1–2 partners; mapping 100.00%; ≥3 Tier 1–2 interface representation 100.00%; PARTIALLY_OVERLAPPING.
5. **SMAD3** — 34 variants across 2 diseases; 123 Tier 1–2 partners; mapping 100.00%; ≥3 Tier 1–2 interface representation 100.00%; PARTIALLY_OVERLAPPING.
6. **LMNA** — 167 variants across 5 diseases; 122 Tier 1–2 partners; mapping 100.00%; ≥3 Tier 1–2 interface representation 100.00%; SEPARABLE_INTERFACES.
7. **STAT3** — 68 variants across 2 diseases; 111 Tier 1–2 partners; mapping 100.00%; ≥3 Tier 1–2 interface representation 100.00%; PARTIALLY_OVERLAPPING.
8. **RAF1** — 36 variants across 2 diseases; 79 Tier 1–2 partners; mapping 100.00%; ≥3 Tier 1–2 interface representation 100.00%; SEPARABLE_INTERFACES.
9. **PMP22** — 29 variants across 2 diseases; 59 Tier 1–2 partners; mapping 100.00%; ≥3 Tier 1–2 interface representation 100.00%; HIGHLY_OVERLAPPING.
10. **PTEN** — 181 variants across 4 diseases; 54 Tier 1–2 partners; mapping 100.00%; ≥3 Tier 1–2 interface representation 100.00%; SEPARABLE_INTERFACES.

## ACTB and ACTG1 secondary cases

Phase 0R exact-RCV pilot records: [{"gene": "ACTB", "n_exact_RCV": "119", "n_unique_PLP_missense_variants": "83", "n_normalized_diseases": "4", "n_distinct_diseases_after_parent_child_collapse": "4", "variant_count_per_disease": "MONDO:0009470=4 (4 residues);MONDO:0011823=1 (1 residues);MONDO:0100433=3 (3 residues);MONDO:0700092=2 (2 residues)", "strict_status": "NO", "reason": "FEWER_THAN_TWO_DISEASES_WITH_5_VARIANTS", "mechanism_status": "UNRESOLVED", "credible_PPI_partners": "509", "PIONEER_partners": "165", "existing_structural_partner_proxy": "8"}, {"gene": "ACTG1", "n_exact_RCV": "87", "n_unique_PLP_missense_variants": "51", "n_normalized_diseases": "3", "n_distinct_diseases_after_parent_child_collapse": "3", "variant_count_per_disease": "MONDO:0018838=2 (2 residues);MONDO:0019497=1 (1 residues);MONDO:0019719=1 (1 residues)", "strict_status": "NO", "reason": "FEWER_THAN_TWO_DISEASES_WITH_5_VARIANTS", "mechanism_status": "UNRESOLVED", "credible_PPI_partners": "202", "PIONEER_partners": "47", "existing_structural_partner_proxy": "3"}]

These remain historical/biological comparison cases and were not forced into the strict primary cohort. Their low variant counts within individual disease groups do not support quantitative primary disease-label testing without further data.

## Methods, thresholds, provenance, and limitations

1. Exact cohort source: results/phase0R/gene_disease_variant_counts_exact_rcv.tsv, restricted to strict_disease_ge10=YES; mechanism state from the frozen Phase 0R gene/pair tables; PPI counts/provenance from the Phase 0R frozen PPI outputs.
2. Variant-to-RCV-to-MONDO linkage uses the Phase 0R seed and recovered RCV table, requiring exact gene, VariationID, RCV accession, and one of the strict MONDO groups. Transcript HGVS was parsed from the frozen ClinVar variant_summary_2026-09.txt.gz. RefSeq protein sequence was obtained from NCBI E-utilities and aligned to HGNC-associated UniProt reference sequences. Equal-best alignments that map a residue differently are not forced.
3. PIONEER snapshot: PIONEER_very_high_human_2026-09-27.txt. Experimental evidence: BioGRID MV-Physical 5.0.261 and cached IntAct PSICQUIC evidence recorded in Phase 0R. HuRI remained unavailable as recorded in Phase 0R.
4. RCSB Search API exact UniProt accession pair query, entry assembly metadata, and biological assembly mmCIF files are cached under data/interim/phase1A_rcsb_pair_cache/, data/interim/phase1A_rcsb_entry_cache/, and data/raw/phase1A_rcsb_structures/. A focused pair-query scope is stated per edge in structure_coverage.tsv. The first three lexically ordered hits per edge were reviewed to bound coordinate processing while retaining the full search result list.
5. Coverage hard flag: within a gene, flag when disease-group coverage spans ≥20 percentage points, or at least one group is <70% while another is ≥70%. Reported separately for protein mappability and ≥3 Tier 1–2 interface partners.
6. Interface definitions do not estimate functional effect. MutPred-PPI predictions are interaction-specific functional scores; none were substituted by structural occupancy.
7. Public source endpoints: [ClinVar monthly variant summary](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2026-09.txt.gz), [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/), [UniProt REST API](https://www.uniprot.org/help/api), [RCSB Search API](https://search.rcsb.org/), [PIONEER frozen human predictions](https://pioneer.yulab.org/static/predictions/very_high/human.txt), [MutPred-PPI source and model documentation](https://github.com/rosstewart/mutpred-ppi), and [public MutPred-PPI training-data deposit](https://doi.org/10.5281/zenodo.22774603).
8. The main limitation is incomplete structure search outside the PIONEER source=PDB subset and unavailable unpublished VarChAMP training rows. The RCSB no-hit category is therefore limited to queried edges. Phase 1B admission should preserve this uncertainty and the differential-missingness flags.

## Output files

All twelve required result tables are in results/phase1A/. Re-run with python3 scripts/run_phase1a.py. The frozen inputs' SHA-256 values and analysis date are recorded in metadata/phase1A_manifest.json.
