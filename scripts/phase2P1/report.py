from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'results/phase2P1';P=R/'results/phase2P0';D=json.loads((O/'reporting_summary.json').read_text())
def read(n):return pd.read_csv(O/n,sep='\t')
def table(df):
 d=df.copy()
 def f(v):
  if pd.isna(v):return 'NA'
  if isinstance(v,(float,np.floating)):return f'{v:.6g}'
  return str(v).replace('|','/').replace('\n',' ')
 return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(f(v) for v in row)+' |' for row in d.itertuples(index=False,name=None))
x=read('xpair_primary_per_gene.tsv');s=read('sppider_confirmatory_per_gene.tsv');corr=read('correlation_distance_sensitivity.tsv');rest=read('restricted_permutation_results.tsv');bas=read('baseline_summary.tsv');pred=read('predictive_global_descriptive_summary.tsv');ext=read('sppider_extension_results.tsv');ctx=read('residue_ceiling_context.tsv');cross=read('cross_method_concordance.tsv');scaling=read('standardization_parameters.tsv');rawsens=read('raw_and_poststandardization_sensitivity.tsv');bs=D['strongest_generic_baseline'];n_frozen=len(read('input_hash_verification.tsv'))
claim=read('final_claim_audit.tsv');claim['sppider_confirmation_summary']=f"median_delta={D['median_SPPIDER_delta']:.9g};positive={int(s.DELTA_R2_SPPIDER.gt(0).sum())}/10;delta_FDR_genes="+';'.join(s.loc[s.q_SPPIDER_DELTA.lt(.05),'gene']);claim['restricted_sensitivity_summary']=f"Delta nominal agreement={int(rest.nominal_delta_concordance.sum())}/10;RS1 restricted p={rest.loc[rest.gene.eq('RS1'),'restricted_p_DELTA'].iloc[0]:.9g}";claim['correlation_distance_summary']=f"{D['correlation_evaluable']}/10 evaluable; all three deltas negative; 4 zero-distance sorted controls;3 dimensionality exclusions";claim['genes_q_delta_lt_0.05']=';'.join(D['genes_q_XPAIR_DELTA_lt_0_05']);claim['genes_q_labeled_lt_0.05']=';'.join(D['genes_q_XPAIR_LABELED_lt_0_05']);claim.to_csv(O/'final_claim_audit.tsv',sep='\t',index=False)
legends={
'P1':'Observed partner-identity ΔR² for all 10 frozen genes (547 gene-residues). Each point is one gene; the vertical line marks zero. Effects use Euclidean PERMANOVA on separately standardized labeled and raw-sorted X-PAIR vectors. No per-gene confidence intervals were preregistered. The global median CI is calculated from 10,000 equal-weight gene bootstrap resamples and reported in the text. Source: xpair_primary_per_gene.tsv.',
'P2':'Paired labeled (blue circles) and sorted (ochre squares) X-PAIR R² values for all 10 genes. Segments connect the identical within-gene residue cohort. Sorting raw probabilities removes partner assignment; rank-column scaling is performed separately without labels. Source: TableP1.tsv.',
'P3':'X-PAIR versus SPPIDER partner-identity ΔR² for all 10 genes. Solid reference lines mark zero; the dashed line is equality. Every gene is labeled. Spearman and Pearson describe concordance, without determining the primary verdict. Labels are repositioned for clearance with leader lines. Source: cross_method_concordance.tsv and Tables P1–P2.',
'P4':'a, Marginal R² for X-PAIR labeled/sorted and four generic baselines. The first two columns use the full frozen primary cohort; baseline columns use their frozen complete cases. These marginal cells are not an equal-row comparison when baseline data are incomplete. NA marks failure of the ≥10-residues-per-group complete-case gate, including absent annotations. b, Each point is X-PAIR R² minus baseline R² recomputed on the identical complete-case rows. Colors and shapes identify position, DOMAIN_PROXY, generic interface and generic 3D. No inferential CI is attached to these descriptive marginal contrasts. Table P3 and the figure source table expose every subset count and missingness status.',
'P5':'a, Fraction of phenotype-annotated residues linked to ≥2 retained disease groups, including multilabel classes. b, Frozen variant-mass residue information ceiling. c, Primary X-PAIR ΔR². All 10 genes appear in the same order, on separate axes. These are descriptive context; ΔR² is not divided by the ceiling. Denominators in a/b include all frozen phenotype-annotated residues for a gene, whereas c uses safe single-disease primary residues. Source: residue_ceiling_context.tsv.',
'P6':'Paired ΔR² permutation null for COL3A1, selected by the neutral alphabetical-first rule. Histogram contains all 10,000 unique, group-size-preserving alternative disease assignments; the observed assignment is excluded. The blue line marks observed ΔR²; the dashed gray line marks zero. The one-sided plus-one P counts null effects at least as large as the observed effect. Source: xpair_permutations/COL3A1.npz.'}
(O/'figures'/'figure_legends.json').write_text(json.dumps(legends,indent=2)+'\n')
summary=f"""# PHASE 2P-1 — Preregistered phenotype information test

## Final decision

**{D['verdict']}**

The fixed 10-gene, 547-residue analysis was executed. Median partner-identity ΔR² was {D['median_delta_R2']:.8f}, with gene-bootstrap 95% CI [{D['bootstrap_CI_low']:.8f}, {D['bootstrap_CI_high']:.8f}]. Eight of ten genes were positive. Condition A failed because the interval included zero; Condition B passed because 80% were positive. The combined preregistered support rule failed. This is a completed biological test, with all 10 primary genes evaluable.

The analysis does not support the planned global claim across this frozen cohort. RS1 had a positive X-PAIR partner-identity contrast with Δ-family q=0.0389961; its labeled disease association did not pass that separate FDR family. Individual contrasts do not replace the global support rule. The original NO_GO_EDGETIC_REPRESENTATION remains unchanged.

## Immutable inputs and execution

Protocol SHA256: `{D['preregistration_SHA256']}`. All {n_frozen} frozen file checks passed before association and again during final integrity verification: original 124 inputs, 79 Phase 2S delivery entries, and 143 Phase 2P-0 preregistration input entries. Hashes, membership, labels, partner identities and matrices were preserved. One canonical gene-residue is one observation; substitutions, RCVs, ClinVar submissions and review stars were not used as weights.

The preflight verified each gene's two group counts against the frozen counts, exact labeled/sorted raw score multisets, unique residue IDs and matched row identities. Disease groups below the prospective residue gate were not reinstated. Common partners were frozen before association. Median common partner count was {D['median_common_partners']:.0f}. All new computations and deliverables reside in this workspace.

The statistical pipeline used Python with the installed NumPy, pandas, SciPy and scikit-learn versions recorded in execution_manifest.tsv. It completed automatically under the fixed specification. A pandas read-only mask error in secondary baseline processing was repaired using an explicit copy; saved primary/confirmatory outputs were reused. BLAS emitted floating-point status warnings despite finite numerical results. Saved null arrays were checked for finiteness, recomputed, and the observed Gram-matrix R² was checked against an independent coordinate sums-of-squares formula for every method and gene. Details appear in statistical_integrity_audit.tsv. No observed effect triggered a cohort or model change.

## Exact statistical procedure

For each gene, each named partner dimension was standardized by its population mean/SD across all frozen primary residues, without disease labels. Raw probabilities were sorted descending within residue, then separately standardized by rank column. Exactly constant columns were set to zero. All transformations are saved in standardization_parameters.tsv. Euclidean distances were analyzed by one-factor PERMANOVA: R²=SS_between/SS_total; pseudo-F uses group−1 and residue−group degrees of freedom.

Each gene used 10,000 unique alternative disease assignments, excluding the observed assignment, preserving exact group sizes. NumPy PCG64 and the frozen gene-specific SHA-derived seeds were used. The same assignments were shared by X-PAIR labeled/sorted and SPPIDER labeled/sorted. Both primary R² permutation arrays and paired ΔR² arrays are saved. One-sided Δ P=(1+count(null ≥ observed))/(B+1); labeled P uses pseudo-F with the same rule. Four BH families remain separate: X-PAIR labeled, X-PAIR Δ, SPPIDER labeled, SPPIDER Δ. No hypothesis families were pooled.

Global summaries weight each gene once. The median CI used 10,000 gene-level percentile bootstrap resamples, PCG64 seed 20261005. The one-sided exact binomial sign test gave P={D['sign_test_P']:.7f}; zero effects would be omitted from the sign test but retained in the global 70% denominator. There were no zero effects. Q1/Q3 use linear interpolation.

## Primary X-PAIR results

Median labeled R²: {D['median_XPAIR_labeled_R2']:.8f}; median sorted R²: {D['median_XPAIR_sorted_R2']:.8f}. Median ΔR²: {D['median_delta_R2']:.8f}; IQR [{D['Q1']:.8f}, {D['Q3']:.8f}]; range [{D['minimum']:.8f}, {D['maximum']:.8f}]. No X-PAIR labeled association had q<0.05. RS1 was the only Δ-family q<0.05 gene. The median paired difference is distinct from the difference between the two marginal medians.

{table(x[['gene','n_total','n_group1','n_group2','common_partner_n','R2_XPAIR_LABELED','R2_XPAIR_SORTED','DELTA_R2_PARTNER_ID','p_DELTA','q_XPAIR_DELTA','q_XPAIR_LABELED']])}

## Confirmatory SPPIDER and concordance

Median SPPIDER ΔR²: {D['median_SPPIDER_delta']:.9f}. Six genes were positive; COL3A1 passed both its separately corrected labeled and Δ families. Four genes were positive in both methods. Across 10 gene effects, Spearman ρ={D['cross_method_Spearman']:.6f}; Pearson r={cross.Pearson.iloc[0]:.6f}. Direction/magnitude agreement was weak and negative across methods. SPPIDER outcomes did not revise the X-PAIR primary decision.

{table(s[['gene','R2_SPPIDER_LABELED','R2_SPPIDER_SORTED','DELTA_R2_SPPIDER','p_DELTA','q_SPPIDER_DELTA','q_SPPIDER_LABELED']])}

{table(cross)}

## Generic baseline comparisons

Every valid contrast recomputed X-PAIR and its comparator on identical frozen complete-case rows, retaining ≥10 residues per group. No missing baseline was imputed. DOMAIN_PROXY denotes frozen UniProtKB feature/domain proxies; formal Pfam/InterPro baselines were unavailable. These are marginal descriptive contrasts and do not establish conditional independence.

{table(bas)}

The largest median available baseline R² was DOMAIN_PROXY ({bs['median_baseline_R2']:.8f}), evaluated in only {int(bs['n_genes'])} genes. Its median row-matched X-PAIR increment was {bs['median_XPAIR_increment']:.8f}. This descriptive ranking compares differing availability sets; it is not a common-denominator declaration of baseline superiority. Position was evaluable in all 10 genes. All 40 gene-baseline attempts, subset sizes and attrition are retained in Table P3.

## Restricted and distance sensitivities

Fixed strata combined protein-position quintile, frozen domain membership and the prespecified structural context. No strata were regenerated. GUCY2D exhaustively evaluated 1,763 alternative assignments (space 1,764); LHCGR evaluated 3,149 (space 3,150). Other genes used 10,000 unique alternatives. At nominal Δ P<0.05, unrestricted/restricted status agreed in 10/10 genes; RS1 remained the only nominal positive contrast. Labeled-association nominal status differed for COL3A1 and COL4A3; these disagreements are retained, without selecting the smaller P.

{table(rest[['gene','B','exhaustive','restricted_p_XPAIR_LABELED','restricted_p_DELTA','unrestricted_p_DELTA','nominal_delta_concordance']])}

Spearman correlation-distance sensitivity was evaluable for KCNQ2, MEN1 and TUBA1A; all three ΔR² values were negative. COL4A3, LHCGR and RS1 failed the ≥10 nonconstant-dimension/row gate. COL3A1, COL4A1, GUCY2D and MMUT had zero-total-distance sorted rank controls and were unevaluable. Identical rank vectors were assigned exact zero Spearman distance to remove floating-point artifacts from corrcoef. The earlier secondary artifacts and repair rationale were retained in diagnostic_archive and correlation_implementation_correction.json; primary Euclidean outcomes were unchanged. Negative eigenvalues and Lingoes constants are recorded for valid non-Euclidean matrices. No excluded row was replaced.

{table(corr)}

Raw-probability Euclidean and sort-after-standardization Euclidean controls were also executed exactly as prespecified and saved separately; they did not replace the primary endpoint.

{table(rawsens[['gene','sensitivity','delta','p_delta']])}

## Same-residue discordance and information ceiling

All cross-disease, multilabel and ambiguous residue classes remained excluded from primary disease assignment. The frozen broader registry retains the 39 exclusive cross-disease residues and all ambiguous/multilabel cases. Cross-disease-linked flags also include residues in multilabel classes, avoiding undercounting. Per-gene counts and disease combinations remain published in the frozen cross_disease_residue_audit.tsv; the registry, exclusions and ceiling files are unchanged. A WT residue-context model necessarily ties substitutions mapped to that residue. This is the prospective residue-information limitation, not an optimized classifier.

The following context uses all phenotype-annotated residues per gene for discordance/ceiling and the frozen single-disease subset for observed ΔR². Ceilings are descriptive variant-mass quantities conditional on the retained ontology. ΔR² was not normalized by them; no association between ceiling and effect was tested.

{table(ctx[['gene','phenotype_annotated_residues','cross_disease_linked_residues','DiscordantResidueRate','RESIDUE_INFORMATION_CEILING','DELTA_R2_PARTNER_ID']])}

## Final training-overlap audit

Before association, official SPPIDER dataset metadata and X-PAIR official repository retrieval were attempted. Both requests failed with proxy tunnel HTTP503. Exact release-specific query protein and query–partner training identity lists remained inaccessible. All 40 primary gene × method × scope audit rows are TRAINING_OVERLAP_UNKNOWN. The models' interaction/interface tasks did not use the present ClinVar disease-label endpoint; no reliable exact training-list verification was available. Structural/interface-model training overlap could not be excluded. This analysis is not described as external validation, fully held-out evaluation or independent model validation. Membership was not changed on this basis. Audit rows and preserved official request responses are in Table P5 and training_sources/.

## Secondary prediction

All 10 genes met the fixed sample-size gate. A nearest-centroid classifier was evaluated by repeated stratified 5-fold CV for 10 repeats, seed 20261006. Within a repeat each residue appeared in exactly one validation fold; paired representations used identical folds. Standardization used training rows only. Probabilities used normalized exp(−squared distance), temperature 1 and floor 1e−15. No tuning, feature selection or test-fold calibration was performed. The following values are medians across the 10 gene-specific mean-over-repeat metrics; repeats/folds are not treated as independent biological samples.

{table(pred)}

X-PAIR labeled median balanced accuracy was 0.4710 and sorted 0.4274; macro-F1 was 0.4627 and 0.4156. These descriptive predictive outcomes do not establish clinical classification utility or replace the global phenotype-information endpoint. All per-gene/per-repeat metrics, paired folds and out-of-fold probabilities are preserved in Table P6 and predictive_outputs/.

## SPPIDER-only extension

Seven of eight frozen starting genes were evaluable (702 residues in total); NF1 failed its prospective sample-size gate. Five of seven had positive ΔR²; no extension Δ or labeled hypothesis passed its separately corrected q<0.05 family. Extension genes, partners and residues remained separate from the primary 10-gene median, bootstrap, sign test and FDR.

{table(ext[['gene','status','n_residues','partners','delta','q_delta_extension','q_labeled_extension']])}

## Figures, tables and source data

Figures are exported as editable SVG and 600-dpi PNG, at 183-mm canvas width. PDF was omitted under the user's output preference. Table P1 is the primary table; P2 confirmatory; P3 all generic baselines; P4 restricted sensitivity; P5 training overlap; P6 all valid prediction metrics. Source data and raw compressed nulls are retained. Every final panel was inspected after rendering; QA and export-contract exceptions are recorded in figures/visual_QA.tsv and figure_export_audit.json.

"""
for n in range(1,7):
 summary+=f"### Figure P{n}\n\n![Figure P{n}](../results/phase2P1/figures/FigureP{n}.png)\n\n{legends['P'+str(n)]}\n\n"
summary+="""## Limits and stop status

The frozen analysis has only 10 genes, unequal group sizes, uneven generic-structure annotation availability and unknown interface-model training overlap. Baseline availability limits broad comparisons. Sorted correlation profiles can collapse to identical ranks; the resulting missing sensitivity is reported explicitly. WT partner-interface probabilities describe residue context, not variant-induced disruption. Neither the global result nor isolated contrasts identify causal partners or demonstrate an edgetic mechanism.

The preregistered primary support rule failed. Results, Tables P1–P6, Figures P1–P6, claim audit, raw numerical outputs and provenance are complete. The workflow stops at Phase 2P-1 for manual review. No tissue concordance, VUS prioritization, mechanistic partner interpretation or manuscript drafting was started.

Original NO_GO_EDGETIC_REPRESENTATION was not altered.

No post-hoc cohort, partner, distance, or endpoint selection was performed.
"""
report=R/'reports/PHASE2P1_PREREGISTERED_PHENOTYPE_INFORMATION_TEST.md';report.write_text(summary)
print('REPORT_WRITTEN',flush=True)
