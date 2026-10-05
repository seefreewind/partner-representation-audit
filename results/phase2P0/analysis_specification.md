# Phase 2P-0 frozen phenotype cohort and prospective statistical specification

Status: preregistered before any disease-label association calculation. Date: 2026-10-04 Asia/Shanghai. Frozen original biological verdict: NO_GO_EDGETIC_REPRESENTATION. Phase 2P-1 requires separate manual authorization even when this preparation gate is GO.

## Question, representations and methods

Does partner-specific interface context contain disease-label information beyond generic residue location, generic interface propensity and interface-score distributions without partner identity? WT sequence interface probabilities are residue-context features. They are not substitution effects, disruption probabilities or pathogenicity scores.

PRIMARY_METHOD = XPAIR_PARTNER_INTERFACE, exact X-PAIR 0.1.0 with frozen Ankh-large encoder/head. Choice is fixed by phenotype-blind Phase 2S structural validity. CONFIRMATORY_METHOD = SPPIDER_SEQ_PARTNER_INTERFACE, exact v2.0-manuscript/ESM2 receptor output. Do not fuse models or switch the primary method after association results.

## Frozen cohort and identity of every observation

Starting shared cohort: COL3A1, COL4A1, COL4A3, FGFR2, GUCY2D, KCNQ2, LHCGR, MEN1, MFN2, MMUT, PTEN, RHO, RPGR, RS1, TUBA1A. Secondary extension: COL7A1, DYNC1H1, MYO7A, NF1, RYR1, SCN1A, USH2A, VWF.

Final primary genes (10): COL3A1, COL4A1, COL4A3, GUCY2D, KCNQ2, LHCGR, MEN1, MMUT, RS1, TUBA1A. Final primary residues: 547. All residue IDs, canonical coordinates, exact retained MONDO labels, substitution provenance and RCVs are frozen in residue_phenotype_registry.tsv and per-gene matrices/*_residue_membership.tsv; retained groups and counts are in phenotype_group_residue_counts.tsv. The exact file contents and SHA256 hashes are part of this protocol through preregistration_input_manifest.tsv, SHA256 `7e809146f191b2c8a4404c56c3d8281315257ea5c1658247c30efd2af3076a42`. This is an immutable input list; alterations require a separately dated amendment before association.

Independent unit: gene × canonical WT residue. Do not independently count substitutions at one residue. Registry includes all mapped coordinate-bearing frozen pathogenic substitutions, including one unsafe mapping where present; unsafe mapping cannot enter the primary cohort. Canonical mapping must be EXACT_CANONICAL or SAFE_ISOFORM_REMAP and match the accession, coordinate and WT amino acid. Mapping/label exclusions are explicit in phenotype_exclusion_log.tsv.

Frozen retained disease groups are the strict ≥10-variant MONDO groups used in Phase 2S. No ontology update. Residue class precedence is MULTILABEL_VARIANT_RESIDUE, then AMBIGUOUS_DISEASE_RESIDUE, then CROSS_DISEASE_RESIDUE, then SINGLE_DISEASE_RESIDUE. Multilabel variants remain explicit even when the residue also has cross-disease links. A separate cross_disease_linked flag counts every residue with ≥2 retained disease groups regardless of exclusive class. This avoids undercounting same-residue discordance.

Only safe SINGLE_DISEASE_RESIDUE rows with both models and ≥3 common named partners can enter the primary sample-size gate. Retain disease groups with ≥10 such unique residues, without lowering this threshold; drop smaller groups prospectively and record them. A gene requires ≥2 qualifying groups and ≥3 common partners. Counts alone determine eligibility. No gene may be selected by future effect strength. The 8 extension genes are subject to the same residue/group thresholds, using SPPIDER availability only; their membership and partner lists are frozen separately and excluded from the X-PAIR global result.

## Same-residue discordance and ceiling

DiscordantResidueRate = residues with ≥2 retained disease labels / all residues with at least one retained disease label. Report counts, percentages and disease combinations; neither multilabel nor cross-disease rows are assigned to a convenient single label. The residue information ceiling is descriptive: each distinct frozen Variation ID contributes one unit of mass, divided equally among its retained labels if multilabel. Sum the largest disease-label mass at each residue and divide by total annotated variant mass. Report a macro residue ceiling as well. Unresolved labels remain unknown and do not count as false mapping errors. This ceiling is conditional on the frozen retained ontology and annotatable labels, not an empirical classification result. No sequence score is involved.

## Feature matrices and frozen controls

COMMON_PARTNERS are lexicographically ordered frozen credible partners with valid, hash-verified full WT tracks in BOTH models, selected without labels. Freeze exact IDs, track paths and hashes in common_partner_set_by_gene.tsv. X-PAIR uses the official six-decimal raw probability export, not the three-decimal presentation file; SPPIDER uses cached receptor probabilities. Same genes, residue rows and named partner columns are used for both methods. No >0.5 filtering applies to primary features.

XPAIR_UNLABELED_SCORE_PROFILE and SPPIDER_UNLABELED_SCORE_PROFILE sort each residue's raw probabilities descending. Raw sorting retains the score multiset, dimension count, mean and variance while removing named partner assignment. Per-residue generic controls are mean, maximum, population SD and fraction ≥0.5; they are secondary. Exact raw matrices and controls are frozen now; no disease-group score summary is calculated.

Position baseline: residue number and residue/protein-length, fixed without tuning. Domain baseline: formal Pfam/InterPro annotations are not present in the frozen baseline files. Do not fetch new annotations or treat missing as no-domain. Available frozen UniProtKB feature annotations are an explicitly named secondary domain proxy with identity, normalized within-domain position and boundary distance; formal Pfam/InterPro comparison is unavailable. Generic interface baseline removes named partner identity and uses any-interface, count of known touching interfaces and finite minimum distance, with unknown structure remaining NA. Generic 3D uses frozen pLDDT, monomer SASA and intrachain CA-neighbor count only; no structure regeneration. All missing masks and source identities are frozen in baseline manifests.

AlphaMissense is never a primary comparator. If already frozen scores are available in Phase 2P-1, a separate variant-level secondary analysis may use them, grouping all substitutions at a residue together. Otherwise report unavailable; no new downloads or primary-row mixing.

## Distances and transformations

For labeled X-PAIR, compute each named partner column's mean and population SD over ALL final phenotype-ready residues WITHIN GENE, without using disease assignments. Standardize columns; exactly zero-variance columns become zero and remain documented. Primary distance is Euclidean on these standardized vectors. The sorted control is separately standardized by rank-column mean/SD over the identical rows, using the same zero-variance rule; its raw prestandardization score multiset is preserved. Freeze this choice now. Confirmatory SPPIDER uses the identical procedure. Store means/SDs and effective varying dimensions in Phase 2P-1 before testing; never tune normalization against separation.

Sensitivity distances: raw-probability Euclidean labeled versus sorted; descending sorting AFTER named-column standardization is a separate sensitivity control, clearly distinguished from sorting raw probabilities. Spearman correlation distance (1 − Spearman) is permitted only for genes with ≥10 nonconstant partner columns and nonconstant within-row profiles for every included residue. Use exactly the frozen primary rows; if invalid, mark unevaluable rather than select rows or replace the primary distance. For non-Euclidean distances, report negative eigenvalues and prespecify Lingoes correction before PERMANOVA; do not choose corrections by P values.

## Primary per-gene test and paired endpoint

In Phase 2P-1 only: distance ~ retained disease group, ordinary one-factor PERMANOVA with an intercept and categorical disease. Report R², pseudo-F, permutation P, groups and unique residue counts. R² = between-group sum of squares / total sum of squares; pseudo-F uses df_between = groups−1 and df_within = residues−groups. A degenerate zero-total-distance matrix is unevaluable, with no cohort replacement.

Primary novelty endpoint: ΔR²_partnerID = R²_XPAIR_labeled − R²_XPAIR_sorted. R² and ΔR² are primary effect sizes; P values are secondary. No covariate-adjusted or baseline contrast may replace this endpoint after seeing results.

## Permutations and multiple testing

Per-gene seed = first 16 hexadecimal SHA256 characters of 'PHASE2P1|20261004|GENE', modulo 2^32; exact seeds, group sizes and multiset assignment spaces are in permutation_plan.tsv. Use NumPy PCG64. B=10,000 unique alternative label assignments sampled without replacement when the space is larger; exclude the observed assignment. If there are fewer alternatives, enumerate all and set B=space−1. Preserve each group size exactly. Use the SAME label assignment for labeled/sorted and X-PAIR/SPPIDER on every permutation. No label permutation has been generated in Phase 2P-0.

One-sided empirical P_delta = (1 + number of permuted ΔR² ≥ observed ΔR²)/(1+B), including ties; labeled-association P uses the same plus-one rule on pseudo-F. BH-FDR across eligible primary genes is separate for Family A (X-PAIR labeled association) and Family B (X-PAIR ΔR²); prespecified reporting threshold q≤0.05. SPPIDER families A/B are separate confirmatory families, never pooled with X-PAIR. Secondary analyses have separately labeled families and do not revise primary claims.

## Global primary summary and support rule

Do not pool residues across genes. Summarize the distribution of per-gene primary ΔR² by median, IQR (linear interpolation), positive/total count and median 95% CI from 10,000 gene-level bootstrap resamples, PCG64 seed 20261005, percentile 2.5/97.5 intervals. Report an exact one-sided binomial sign test P (positive direction; zero deltas omitted from the sign test but retained in total/70% denominator). Degenerate genes remain explicitly unevaluable; if fewer than 10 genes remain evaluable at testing, withhold the planned multi-gene support claim.

PARTNER_IDENTITY_SUPPORT requires BOTH median ΔR² >0 with bootstrap CI strictly above zero AND ≥70% of all phenotype-ready genes having positive ΔR². Individual FDR-positive genes can still be reported when this global rule fails. Missing/unevaluable genes do not count as positive and are not replaced. No causal or mechanistic claim follows.

## Confirmatory and secondary analyses

SPPIDER repeats identical labeled-versus-sorted tests with same genes/residues/partners/permutations. Report direction and magnitude agreement transparently. Its nonconfirmation does not automatically invalidate X-PAIR. SPPIDER_ONLY_EXTENSION runs only after the shared-model analysis is frozen/completed; freeze 8 starting genes, record sample-size attrition, and exclude all extension genes from the primary global estimate.

Secondary marginal contrasts: R²_XPAIR−R²_position, R²_XPAIR−R²_domain_proxy and R²_XPAIR−R²_generic_interface, with paired gene summaries. They are comparisons of marginal summaries, not proof of conditional independence. Incomplete baselines use only prospectively frozen complete-case residues satisfying the same ≥10-per-group/≥2-group gate; recompute both X-PAIR and baseline on exactly the same subset, report attrition, and never compare unequal-row R² values. Generic-3D comparison follows the same complete-case rule. No baseline imputation by disease. Position numbers within a gene are collinear with normalized position; use the normalized feature for distance without double-counting.

## Restricted confounding sensitivity

Fixed strata: normalized-position quintile bins [0,.2,.4,.6,.8,1], frozen domain membership (IN_DOMAIN/NO_DOMAIN/DOMAIN_UNKNOWN), and surface/buried bin (monomer SASA≤1 Å², >1 Å², or SURFACE_UNKNOWN). Position=1 enters final quintile. No disease-informed binning, new domain data or finer subdivisions. The exact strata per residue and available multiset label-assignment spaces are frozen. Exchange labels only within these fixed coarse combined strata. If fewer than 1,000 alternative assignments exist, mark RESTRICTED_PERMUTATION_UNINFORMATIVE; do not merge strata after seeing association results or claim a precise restricted P. Otherwise use up to 10,000 unique alternatives with same paired assignments and deterministic gene seed+1 modulo 2^32.

## Falsification and prediction

Same-residue cross-disease audit and information ceiling are mandatory limitations/negative controls. Keep cross-disease/multilabel rows out of primary association, but publish their counts and combinations. Do not relabel them to improve separation. A deterministic residue-only model necessarily ties their context; this mathematical upper bound is the prespecified falsification analysis, not an optimized predictor.

Prediction is secondary only and must be separately authorized with Phase 2P-1. Fixed nearest-centroid classifier, repeated stratified 5-fold CV for 10 repeats when every retained group has ≥10 residues; PCG64 seed 20261006, with all substitutions of one residue inseparable. Fit standardization using training folds only and apply it to validation; compare labeled versus sorted on identical folds. Balanced accuracy and macro-F1 are primary prediction summaries. For log loss use normalized exp(−squared distance to each training class centroid), temperature=1, stable log-sum-exp, probability floor=1e−15; no post-hoc calibration or tuning. No complex ML, random substitution-level splits or selection by CV performance.

## Training overlap and interpretation

training_overlap_audit.tsv distinguishes query and query–partner scopes for both models. Bundled packages contain no complete training identity lists; official dataset/repository fetch attempts returned proxy HTTP503, with logs under training_sources. All genes are TRAINING_OVERLAP_UNKNOWN, not NO_KNOWN_OVERLAP. Current SPPIDER dataset releases cannot automatically establish exact v2.0-manuscript training membership. Renew public-list audit before any Phase 2P-1 association, preserving identity lists and hashes; use known-overlap sensitivity if reliable lists become available. Overlap is not automatically an exclusion because interface training did not use ClinVar disease labels, but unknown is never treated as held-out.

Allowed interpretation: associated with, carries phenotype information, partitions allelic phenotypes, partner-specific interface context is informative. Prohibited: determines disease, causes phenotype via partner X, proves binding disruption or edgetic mechanism. Original NO_GO_EDGETIC_REPRESENTATION remains frozen.

## Prospective registry QC amendment

Before any association, registry n_substitutions was corrected to count distinct amino-acid protein changes and n_Variation_IDs separately records ClinVar IDs. Some nucleotide variants encode the same amino-acid substitution. No residue, gene, label or partner membership changed. The prior specification hash and reason are retained in prospective_QC_amendment.json; this specification supersedes that preliminary hash.

## Hard stop and amendments

Phase 2P-0 constructs counts, matrices, provenance, ceilings and protocol only. No PERMANOVA, disease R², ΔR², disease clustering, phenotype prediction, disease-specific score means, enrichment or tissue concordance is calculated. Record this protocol SHA256 and stop. Any later changes require a new immutable amendment with reason and hash BEFORE computing association; do not silently edit this frozen protocol.
