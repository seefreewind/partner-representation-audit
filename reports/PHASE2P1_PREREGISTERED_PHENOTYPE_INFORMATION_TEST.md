# PHASE 2P-1 — Preregistered phenotype information test

## Final decision

**NO_GLOBAL_PARTNER_IDENTITY_SUPPORT**

The fixed 10-gene, 547-residue analysis was executed. Median partner-identity ΔR² was 0.00179553, with gene-bootstrap 95% CI [-0.00410526, 0.00495717]. Eight of ten genes were positive. Condition A failed because the interval included zero; Condition B passed because 80% were positive. The combined preregistered support rule failed. This is a completed biological test, with all 10 primary genes evaluable.

The analysis does not support the planned global claim across this frozen cohort. RS1 had a positive X-PAIR partner-identity contrast with Δ-family q=0.0389961; its labeled disease association did not pass that separate FDR family. Individual contrasts do not replace the global support rule. The original NO_GO_EDGETIC_REPRESENTATION remains unchanged.

## Immutable inputs and execution

Protocol SHA256: `d657fab6b24f9aa757f9773bb1e9e0bd5d850e743eab6bdf819f126326ced9c3`. All 346 frozen file checks passed before association and again during final integrity verification: original 124 inputs, 79 Phase 2S delivery entries, and 143 Phase 2P-0 preregistration input entries. Hashes, membership, labels, partner identities and matrices were preserved. One canonical gene-residue is one observation; substitutions, RCVs, ClinVar submissions and review stars were not used as weights.

The preflight verified each gene's two group counts against the frozen counts, exact labeled/sorted raw score multisets, unique residue IDs and matched row identities. Disease groups below the prospective residue gate were not reinstated. Common partners were frozen before association. Median common partner count was 17. All new computations and deliverables reside in this workspace.

The statistical pipeline used Python with the installed NumPy, pandas, SciPy and scikit-learn versions recorded in execution_manifest.tsv. It completed automatically under the fixed specification. A pandas read-only mask error in secondary baseline processing was repaired using an explicit copy; saved primary/confirmatory outputs were reused. BLAS emitted floating-point status warnings despite finite numerical results. Saved null arrays were checked for finiteness, recomputed, and the observed Gram-matrix R² was checked against an independent coordinate sums-of-squares formula for every method and gene. Details appear in statistical_integrity_audit.tsv. No observed effect triggered a cohort or model change.

## Exact statistical procedure

For each gene, each named partner dimension was standardized by its population mean/SD across all frozen primary residues, without disease labels. Raw probabilities were sorted descending within residue, then separately standardized by rank column. Exactly constant columns were set to zero. All transformations are saved in standardization_parameters.tsv. Euclidean distances were analyzed by one-factor PERMANOVA: R²=SS_between/SS_total; pseudo-F uses group−1 and residue−group degrees of freedom.

Each gene used 10,000 unique alternative disease assignments, excluding the observed assignment, preserving exact group sizes. NumPy PCG64 and the frozen gene-specific SHA-derived seeds were used. The same assignments were shared by X-PAIR labeled/sorted and SPPIDER labeled/sorted. Both primary R² permutation arrays and paired ΔR² arrays are saved. One-sided Δ P=(1+count(null ≥ observed))/(B+1); labeled P uses pseudo-F with the same rule. Four BH families remain separate: X-PAIR labeled, X-PAIR Δ, SPPIDER labeled, SPPIDER Δ. No hypothesis families were pooled.

Global summaries weight each gene once. The median CI used 10,000 gene-level percentile bootstrap resamples, PCG64 seed 20261005. The one-sided exact binomial sign test gave P=0.0546875; zero effects would be omitted from the sign test but retained in the global 70% denominator. There were no zero effects. Q1/Q3 use linear interpolation.

## Primary X-PAIR results

Median labeled R²: 0.01306520; median sorted R²: 0.00848280. Median ΔR²: 0.00179553; IQR [0.00014122, 0.00403078]; range [-0.02229768, 0.01249303]. No X-PAIR labeled association had q<0.05. RS1 was the only Δ-family q<0.05 gene. The median paired difference is distinct from the difference between the two marginal medians.

| gene | n_total | n_group1 | n_group2 | common_partner_n | R2_XPAIR_LABELED | R2_XPAIR_SORTED | DELTA_R2_PARTNER_ID | p_DELTA | q_XPAIR_DELTA | q_XPAIR_LABELED |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| COL3A1 | 108 | 97 | 11 | 29 | 0.0420686 | 0.0507722 | -0.0087036 | 0.977702 | 0.977702 | 0.20298 |
| COL4A1 | 33 | 17 | 16 | 21 | 0.0164396 | 0.00912894 | 0.00731068 | 0.0871913 | 0.435956 | 0.957793 |
| COL4A3 | 122 | 112 | 10 | 6 | 0.02002 | 0.0195269 | 0.000493077 | 0.266073 | 0.687056 | 0.419958 |
| GUCY2D | 24 | 11 | 13 | 10 | 0.0629963 | 0.0852939 | -0.0222977 | 0.89861 | 0.977702 | 0.672933 |
| KCNQ2 | 77 | 17 | 60 | 16 | 0.00137907 | 0.00135513 | 2.39354e-05 | 0.549645 | 0.687056 | 0.994401 |
| LHCGR | 20 | 10 | 10 | 6 | 0.0130042 | 0.0104005 | 0.00260366 | 0.407059 | 0.687056 | 0.957793 |
| MEN1 | 46 | 29 | 17 | 167 | 0.00941971 | 0.00783666 | 0.00158305 | 0.532347 | 0.687056 | 0.957793 |
| MMUT | 54 | 11 | 43 | 18 | 0.0084162 | 0.00640819 | 0.00200801 | 0.490551 | 0.687056 | 0.957793 |
| RS1 | 31 | 14 | 17 | 4 | 0.0131262 | 0.000633213 | 0.012493 | 0.00389961 | 0.0389961 | 0.957793 |
| TUBA1A | 32 | 11 | 21 | 554 | 0.0104924 | 0.0059859 | 0.00450648 | 0.476252 | 0.687056 | 0.957793 |

## Confirmatory SPPIDER and concordance

Median SPPIDER ΔR²: 0.000015874. Six genes were positive; COL3A1 passed both its separately corrected labeled and Δ families. Four genes were positive in both methods. Across 10 gene effects, Spearman ρ=-0.745455; Pearson r=-0.581819. Direction/magnitude agreement was weak and negative across methods. SPPIDER outcomes did not revise the X-PAIR primary decision.

| gene | R2_SPPIDER_LABELED | R2_SPPIDER_SORTED | DELTA_R2_SPPIDER | p_DELTA | q_SPPIDER_DELTA | q_SPPIDER_LABELED |
| --- | --- | --- | --- | --- | --- | --- |
| COL3A1 | 0.0944661 | 0.0905216 | 0.00394446 | 0.00409959 | 0.0409959 | 0.0469953 |
| COL4A1 | 0.0312757 | 0.0319715 | -0.000695816 | 0.928307 | 0.969003 | 0.629062 |
| COL4A3 | 0.000874639 | 0.000873865 | 7.73199e-07 | 0.919608 | 0.969003 | 0.964104 |
| GUCY2D | 0.140092 | 0.138933 | 0.00115912 | 0.274273 | 0.685681 | 0.0619938 |
| KCNQ2 | 0.0177262 | 0.0177512 | -2.49904e-05 | 0.586541 | 0.969003 | 0.483952 |
| LHCGR | 0.113035 | 0.112434 | 0.000601458 | 0.241176 | 0.685681 | 0.447705 |
| MEN1 | 0.0168735 | 0.0168425 | 3.09738e-05 | 0.50355 | 0.969003 | 0.629062 |
| MMUT | 0.00158486 | 0.00076668 | 0.000818181 | 0.220778 | 0.685681 | 0.964104 |
| RS1 | 0.0934256 | 0.111407 | -0.0179814 | 0.969003 | 0.969003 | 0.229977 |
| TUBA1A | 0.0314126 | 0.0346895 | -0.00327687 | 0.687631 | 0.969003 | 0.582942 |

| Spearman | Pearson | same_sign_n | positive_both_n | XPAIR_positive_only_n | SPPIDER_positive_only_n | n_genes |
| --- | --- | --- | --- | --- | --- | --- |
| -0.745455 | -0.581819 | 4 | 4 | 4 | 2 | 10 |

## Generic baseline comparisons

Every valid contrast recomputed X-PAIR and its comparator on identical frozen complete-case rows, retaining ≥10 residues per group. No missing baseline was imputed. DOMAIN_PROXY denotes frozen UniProtKB feature/domain proxies; formal Pfam/InterPro baselines were unavailable. These are marginal descriptive contrasts and do not establish conditional independence.

| baseline | n_genes | median_baseline_R2 | median_row_matched_XPAIR_R2 | median_XPAIR_increment |
| --- | --- | --- | --- | --- |
| DOMAIN_PROXY | 2 | 0.0414639 | 0.0249125 | -0.0165513 |
| GENERIC3D | 4 | 0.014955 | 0.00995605 | -0.00621503 |
| GENERIC_INTERFACE | 5 | 0.0047717 | 0.00941971 | 0.00464801 |
| POSITION | 10 | 0.0132872 | 0.0130652 | -0.00515769 |

The largest median available baseline R² was DOMAIN_PROXY (0.04146386), evaluated in only 2 genes. Its median row-matched X-PAIR increment was -0.01655132. This descriptive ranking compares differing availability sets; it is not a common-denominator declaration of baseline superiority. Position was evaluable in all 10 genes. All 40 gene-baseline attempts, subset sizes and attrition are retained in Table P3.

## Restricted and distance sensitivities

Fixed strata combined protein-position quintile, frozen domain membership and the prespecified structural context. No strata were regenerated. GUCY2D exhaustively evaluated 1,763 alternative assignments (space 1,764); LHCGR evaluated 3,149 (space 3,150). Other genes used 10,000 unique alternatives. At nominal Δ P<0.05, unrestricted/restricted status agreed in 10/10 genes; RS1 remained the only nominal positive contrast. Labeled-association nominal status differed for COL3A1 and COL4A3; these disagreements are retained, without selecting the smaller P.

| gene | B | exhaustive | restricted_p_XPAIR_LABELED | restricted_p_DELTA | unrestricted_p_DELTA | nominal_delta_concordance |
| --- | --- | --- | --- | --- | --- | --- |
| COL3A1 | 10000 | False | 0.0524948 | 0.945605 | 0.977702 | True |
| COL4A1 | 10000 | False | 0.843616 | 0.105489 | 0.0871913 | True |
| COL4A3 | 10000 | False | 0.0256974 | 0.285271 | 0.266073 | True |
| GUCY2D | 1763 | True | 0.195011 | 0.915533 | 0.89861 | True |
| KCNQ2 | 10000 | False | 0.996 | 0.584342 | 0.549645 | True |
| LHCGR | 3149 | True | 0.776508 | 0.315873 | 0.407059 | True |
| MEN1 | 10000 | False | 0.621738 | 0.524948 | 0.532347 | True |
| MMUT | 10000 | False | 0.705229 | 0.451055 | 0.490551 | True |
| RS1 | 10000 | False | 0.543246 | 0.00239976 | 0.00389961 | True |
| TUBA1A | 10000 | False | 0.774023 | 0.463054 | 0.476252 | True |

Spearman correlation-distance sensitivity was evaluable for KCNQ2, MEN1 and TUBA1A; all three ΔR² values were negative. COL4A3, LHCGR and RS1 failed the ≥10 nonconstant-dimension/row gate. COL3A1, COL4A1, GUCY2D and MMUT had zero-total-distance sorted rank controls and were unevaluable. Identical rank vectors were assigned exact zero Spearman distance to remove floating-point artifacts from corrcoef. The earlier secondary artifacts and repair rationale were retained in diagnostic_archive and correlation_implementation_correction.json; primary Euclidean outcomes were unchanged. Negative eigenvalues and Lingoes constants are recorded for valid non-Euclidean matrices. No excluded row was replaced.

| gene | status | usable_dimensions | uncorrected_negative_eigenvalues | Lingoes_constant | sorted_status | sorted_Lingoes_constant | R2_labeled | R2_sorted | DELTA_R2 | p_labeled | p_DELTA |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| COL3A1 | UNEVALUABLE_ZERO_TOTAL_SORTED_DISTANCE | 29 | 82 | 0.451938 | DEGENERATE_ZERO_TOTAL_DISTANCE | 0 | 0.0108892 | NA | NA | NA | NA |
| COL4A1 | UNEVALUABLE_ZERO_TOTAL_SORTED_DISTANCE | 21 | 19 | 0.0523468 | DEGENERATE_ZERO_TOTAL_DISTANCE | 0 | 0.0260351 | NA | NA | NA | NA |
| COL4A3 | UNEVALUABLE_DIMENSION_OR_ROW_GATE | 6 | NA | NA | NONCONSTANT_DIMENSION_OR_ROW_GATE_FAIL | NA | NA | NA | NA | NA | NA |
| GUCY2D | UNEVALUABLE_ZERO_TOTAL_SORTED_DISTANCE | 10 | 15 | 1.175 | DEGENERATE_ZERO_TOTAL_DISTANCE | 0 | 0.0373946 | NA | NA | NA | NA |
| KCNQ2 | EVALUABLE_LINGOES | 16 | 64 | 0.584788 | EVALUABLE_LINGOES | 8.32933e-05 | 0.0113036 | 0.0118708 | -0.000567233 | 0.976102 | 0.673733 |
| LHCGR | UNEVALUABLE_DIMENSION_OR_ROW_GATE | 6 | NA | NA | NONCONSTANT_DIMENSION_OR_ROW_GATE_FAIL | NA | NA | NA | NA | NA | NA |
| MEN1 | EVALUABLE_LINGOES | 167 | 23 | 0.147122 | EVALUABLE_LINGOES | 1.39069e-28 | 0.0145223 | 0.0379085 | -0.0233862 | 0.933307 | 0.976702 |
| MMUT | UNEVALUABLE_ZERO_TOTAL_SORTED_DISTANCE | 18 | 37 | 0.779378 | DEGENERATE_ZERO_TOTAL_DISTANCE | 0 | 0.0164215 | NA | NA | NA | NA |
| RS1 | UNEVALUABLE_DIMENSION_OR_ROW_GATE | 4 | NA | NA | NONCONSTANT_DIMENSION_OR_ROW_GATE_FAIL | NA | NA | NA | NA | NA | NA |
| TUBA1A | EVALUABLE_LINGOES | 554 | 14 | 0.0638871 | EVALUABLE_LINGOES | 6.8813e-15 | 0.0203126 | 0.0331486 | -0.012836 | 0.907409 | 0.907109 |

Raw-probability Euclidean and sort-after-standardization Euclidean controls were also executed exactly as prespecified and saved separately; they did not replace the primary endpoint.

| gene | sensitivity | delta | p_delta |
| --- | --- | --- | --- |
| COL3A1 | RAW_PROBABILITY_EUCLIDEAN | -0.00591322 | 0.952005 |
| COL3A1 | SORT_AFTER_STANDARDIZATION_EUCLIDEAN | -0.00958769 | 0.975402 |
| COL4A1 | RAW_PROBABILITY_EUCLIDEAN | 0.00673586 | 0.135286 |
| COL4A1 | SORT_AFTER_STANDARDIZATION_EUCLIDEAN | 0.0136724 | 0.10139 |
| COL4A3 | RAW_PROBABILITY_EUCLIDEAN | 0.000995915 | 0.10139 |
| COL4A3 | SORT_AFTER_STANDARDIZATION_EUCLIDEAN | 0.000436473 | 0.449255 |
| GUCY2D | RAW_PROBABILITY_EUCLIDEAN | 0.0129981 | 0.188981 |
| GUCY2D | SORT_AFTER_STANDARDIZATION_EUCLIDEAN | -0.00420175 | 0.70353 |
| KCNQ2 | RAW_PROBABILITY_EUCLIDEAN | 0.000505934 | 0.481552 |
| KCNQ2 | SORT_AFTER_STANDARDIZATION_EUCLIDEAN | 0.000217908 | 0.564244 |
| LHCGR | RAW_PROBABILITY_EUCLIDEAN | 0.000661812 | 0.477952 |
| LHCGR | SORT_AFTER_STANDARDIZATION_EUCLIDEAN | 0.00438376 | 0.351365 |
| MEN1 | RAW_PROBABILITY_EUCLIDEAN | 0.00121458 | 0.568543 |
| MEN1 | SORT_AFTER_STANDARDIZATION_EUCLIDEAN | 0.00129891 | 0.566943 |
| MMUT | RAW_PROBABILITY_EUCLIDEAN | 0.00139634 | 0.526147 |
| MMUT | SORT_AFTER_STANDARDIZATION_EUCLIDEAN | 0.00212887 | 0.515648 |
| RS1 | RAW_PROBABILITY_EUCLIDEAN | 0.0126064 | 0.00339966 |
| RS1 | SORT_AFTER_STANDARDIZATION_EUCLIDEAN | 0.0130251 | 0.00329967 |
| TUBA1A | RAW_PROBABILITY_EUCLIDEAN | 0.00440455 | 0.483252 |
| TUBA1A | SORT_AFTER_STANDARDIZATION_EUCLIDEAN | 0.0044762 | 0.505849 |

## Same-residue discordance and information ceiling

All cross-disease, multilabel and ambiguous residue classes remained excluded from primary disease assignment. The frozen broader registry retains the 39 exclusive cross-disease residues and all ambiguous/multilabel cases. Cross-disease-linked flags also include residues in multilabel classes, avoiding undercounting. Per-gene counts and disease combinations remain published in the frozen cross_disease_residue_audit.tsv; the registry, exclusions and ceiling files are unchanged. A WT residue-context model necessarily ties substitutions mapped to that residue. This is the prospective residue-information limitation, not an optimized classifier.

The following context uses all phenotype-annotated residues per gene for discordance/ceiling and the frozen single-disease subset for observed ΔR². Ceilings are descriptive variant-mass quantities conditional on the retained ontology. ΔR² was not normalized by them; no association between ceiling and effect was tested.

| gene | phenotype_annotated_residues | cross_disease_linked_residues | DiscordantResidueRate | RESIDUE_INFORMATION_CEILING | DELTA_R2_PARTNER_ID |
| --- | --- | --- | --- | --- | --- |
| COL3A1 | 108 | 0 | 0 | 1 | -0.0087036 |
| COL4A1 | 33 | 0 | 0 | 1 | 0.00731068 |
| COL4A3 | 150 | 21 | 0.14 | 0.929348 | 0.000493077 |
| GUCY2D | 28 | 4 | 0.142857 | 0.926471 | -0.0222977 |
| KCNQ2 | 146 | 69 | 0.472603 | 0.792157 | 2.39354e-05 |
| LHCGR | 20 | 0 | 0 | 1 | 0.00260366 |
| MEN1 | 84 | 38 | 0.452381 | 0.774648 | 0.00158305 |
| MMUT | 81 | 27 | 0.333333 | 0.849515 | 0.00200801 |
| RS1 | 54 | 23 | 0.425926 | 0.821053 | 0.012493 |
| TUBA1A | 64 | 32 | 0.5 | 0.746988 | 0.00450648 |

## Final training-overlap audit

Before association, official SPPIDER dataset metadata and X-PAIR official repository retrieval were attempted. Both requests failed with proxy tunnel HTTP503. Exact release-specific query protein and query–partner training identity lists remained inaccessible. All 40 primary gene × method × scope audit rows are TRAINING_OVERLAP_UNKNOWN. The models' interaction/interface tasks did not use the present ClinVar disease-label endpoint; no reliable exact training-list verification was available. Structural/interface-model training overlap could not be excluded. This analysis is not described as external validation, fully held-out evaluation or independent model validation. Membership was not changed on this basis. Audit rows and preserved official request responses are in Table P5 and training_sources/.

## Secondary prediction

All 10 genes met the fixed sample-size gate. A nearest-centroid classifier was evaluated by repeated stratified 5-fold CV for 10 repeats, seed 20261006. Within a repeat each residue appeared in exactly one validation fold; paired representations used identical folds. Standardization used training rows only. Probabilities used normalized exp(−squared distance), temperature 1 and floor 1e−15. No tuning, feature selection or test-fold calibration was performed. The following values are medians across the 10 gene-specific mean-over-repeat metrics; repeats/folds are not treated as independent biological samples.

| representation | balanced_accuracy | macro_F1 | log_loss |
| --- | --- | --- | --- |
| SPPIDER_LABELED | 0.550375 | 0.485769 | 2.52661 |
| SPPIDER_SORTED | 0.535672 | 0.46201 | 2.57472 |
| XPAIR_LABELED | 0.471002 | 0.462677 | 2.45659 |
| XPAIR_SORTED | 0.427371 | 0.415643 | 3.07723 |

X-PAIR labeled median balanced accuracy was 0.4710 and sorted 0.4274; macro-F1 was 0.4627 and 0.4156. These descriptive predictive outcomes do not establish clinical classification utility or replace the global phenotype-information endpoint. All per-gene/per-repeat metrics, paired folds and out-of-fold probabilities are preserved in Table P6 and predictive_outputs/.

## SPPIDER-only extension

Seven of eight frozen starting genes were evaluable (702 residues in total); NF1 failed its prospective sample-size gate. Five of seven had positive ΔR²; no extension Δ or labeled hypothesis passed its separately corrected q<0.05 family. Extension genes, partners and residues remained separate from the primary 10-gene median, bootstrap, sign test and FDR.

| gene | status | n_residues | partners | delta | q_delta_extension | q_labeled_extension |
| --- | --- | --- | --- | --- | --- | --- |
| COL7A1 | COMPLETE | 98 | 9 | 7.01647e-05 | 0.568168 | 0.559944 |
| DYNC1H1 | COMPLETE | 40 | 234 | 0.00661182 | 0.568168 | 0.979802 |
| MYO7A | COMPLETE | 28 | 11 | 0.00364165 | 0.568168 | 0.891291 |
| NF1 | FROZEN_SAMPLE_SIZE_GATE_FAIL | NA | NA | NA | NA | NA |
| RYR1 | COMPLETE | 51 | 19 | 6.78038e-06 | 0.636516 | 0.559944 |
| SCN1A | COMPLETE | 374 | 9 | -0.000119862 | 0.681132 | 0.822768 |
| USH2A | COMPLETE | 76 | 5 | 0.00269131 | 0.568168 | 0.979802 |
| VWF | COMPLETE | 35 | 31 | -0.0024848 | 0.681132 | 0.559944 |

## Figures, tables and source data

Figures are exported as editable SVG and 600-dpi PNG, at 183-mm canvas width. PDF was omitted under the user's output preference. Table P1 is the primary table; P2 confirmatory; P3 all generic baselines; P4 restricted sensitivity; P5 training overlap; P6 all valid prediction metrics. Source data and raw compressed nulls are retained. Every final panel was inspected after rendering; QA and export-contract exceptions are recorded in figures/visual_QA.tsv and figure_export_audit.json.

### Figure P1

![Figure P1](../results/phase2P1/figures/FigureP1.png)

Observed partner-identity ΔR² for all 10 frozen genes (547 gene-residues). Each point is one gene; the vertical line marks zero. Effects use Euclidean PERMANOVA on separately standardized labeled and raw-sorted X-PAIR vectors. No per-gene confidence intervals were preregistered. The global median CI is calculated from 10,000 equal-weight gene bootstrap resamples and reported in the text. Source: xpair_primary_per_gene.tsv.

### Figure P2

![Figure P2](../results/phase2P1/figures/FigureP2.png)

Paired labeled (blue circles) and sorted (ochre squares) X-PAIR R² values for all 10 genes. Segments connect the identical within-gene residue cohort. Sorting raw probabilities removes partner assignment; rank-column scaling is performed separately without labels. Source: TableP1.tsv.

### Figure P3

![Figure P3](../results/phase2P1/figures/FigureP3.png)

X-PAIR versus SPPIDER partner-identity ΔR² for all 10 genes. Solid reference lines mark zero; the dashed line is equality. Every gene is labeled. Spearman and Pearson describe concordance, without determining the primary verdict. Labels are repositioned for clearance with leader lines. Source: cross_method_concordance.tsv and Tables P1–P2.

### Figure P4

![Figure P4](../results/phase2P1/figures/FigureP4.png)

a, Marginal R² for X-PAIR labeled/sorted and four generic baselines. The first two columns use the full frozen primary cohort; baseline columns use their frozen complete cases. These marginal cells are not an equal-row comparison when baseline data are incomplete. NA marks failure of the ≥10-residues-per-group complete-case gate, including absent annotations. b, Each point is X-PAIR R² minus baseline R² recomputed on the identical complete-case rows. Colors and shapes identify position, DOMAIN_PROXY, generic interface and generic 3D. No inferential CI is attached to these descriptive marginal contrasts. Table P3 and the figure source table expose every subset count and missingness status.

### Figure P5

![Figure P5](../results/phase2P1/figures/FigureP5.png)

a, Fraction of phenotype-annotated residues linked to ≥2 retained disease groups, including multilabel classes. b, Frozen variant-mass residue information ceiling. c, Primary X-PAIR ΔR². All 10 genes appear in the same order, on separate axes. These are descriptive context; ΔR² is not divided by the ceiling. Denominators in a/b include all frozen phenotype-annotated residues for a gene, whereas c uses safe single-disease primary residues. Source: residue_ceiling_context.tsv.

### Figure P6

![Figure P6](../results/phase2P1/figures/FigureP6.png)

Paired ΔR² permutation null for COL3A1, selected by the neutral alphabetical-first rule. Histogram contains all 10,000 unique, group-size-preserving alternative disease assignments; the observed assignment is excluded. The blue line marks observed ΔR²; the dashed gray line marks zero. The one-sided plus-one P counts null effects at least as large as the observed effect. Source: xpair_permutations/COL3A1.npz.

## Limits and stop status

The frozen analysis has only 10 genes, unequal group sizes, uneven generic-structure annotation availability and unknown interface-model training overlap. Baseline availability limits broad comparisons. Sorted correlation profiles can collapse to identical ranks; the resulting missing sensitivity is reported explicitly. WT partner-interface probabilities describe residue context, not variant-induced disruption. Neither the global result nor isolated contrasts identify causal partners or demonstrate an edgetic mechanism.

The preregistered primary support rule failed. Results, Tables P1–P6, Figures P1–P6, claim audit, raw numerical outputs and provenance are complete. The workflow stops at Phase 2P-1 for manual review. No tissue concordance, VUS prioritization, mechanistic partner interpretation or manuscript drafting was started.

Original NO_GO_EDGETIC_REPRESENTATION was not altered.

No post-hoc cohort, partner, distance, or endpoint selection was performed.
