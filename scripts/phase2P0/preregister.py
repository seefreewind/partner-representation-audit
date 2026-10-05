from pathlib import Path
import json,hashlib,math,datetime
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'results/phase2P0';M=O/'matrices'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def read(n):return pd.read_csv(O/n,sep='\t',keep_default_na=False)
s=json.loads((O/'cohort_summary.json').read_text());ready=read('phenotype_ready_genes.tsv');reg=read('residue_phenotype_registry.tsv');plan=read('permutation_plan.tsv');strata=[]
for col in ['restricted_assignment_space','restricted_status','restricted_B']:
 plan[col]=plan[col].astype('object') if col in plan else pd.Series(['']*len(plan),dtype='object')
for gene in s['phenotype_ready_genes']:
 rr=pd.read_csv(M/(gene+'_residue_membership.tsv'),sep='\t');pos=pd.read_csv(M/(gene+'_baseline_position.tsv'),sep='\t');dom=pd.read_csv(M/(gene+'_baseline_domain.tsv'),sep='\t',keep_default_na=False);geo=pd.read_csv(M/(gene+'_baseline_generic3d.tsv'),sep='\t',keep_default_na=False)
 d=rr.merge(pos[['residue_id','normalized_position']],on='residue_id',validate='one_to_one').merge(dom[['residue_id','domain_identity','annotation_status']],on='residue_id',validate='one_to_one').merge(geo[['residue_id','monomer_SASA_A2']],on='residue_id',validate='one_to_one');d['position_quintile']=np.minimum(4,np.floor(d.normalized_position*5).astype(int));d['domain_bin']=['NO_DOMAIN' if x=='NO_DOMAIN_FEATURE_AT_POSITION' else 'IN_DOMAIN' if v not in ['','NA'] else 'DOMAIN_UNKNOWN' for x,v in zip(d.annotation_status,d.domain_identity)];sasa=pd.to_numeric(d.monomer_SASA_A2,errors='coerce');d['surface_bin']=np.where(sasa.isna(),'SURFACE_UNKNOWN',np.where(sasa<=1,'BURIED_SASA_LE1','SURFACE_SASA_GT1'))
 def space(cols):
  n=1
  for _,part in d.groupby(cols,dropna=False):
   counts=part.MONDO_labels.value_counts();n*=math.factorial(len(part))//math.prod(math.factorial(int(v)) for v in counts)
  return n
 contexts=['position_quintile','domain_bin','surface_bin'];count=space(contexts);status='RESTRICTED_PERMUTATION_AVAILABLE' if count-1>=1000 else 'RESTRICTED_PERMUTATION_UNINFORMATIVE';mask=plan.gene.eq(gene);plan.loc[mask,'restricted_assignment_space']=str(count);plan.loc[mask,'restricted_status']=status;plan.loc[mask,'restricted_B']=min(10000,count-1) if count-1>=1000 else 0
 for row in d.itertuples():strata.append(dict(gene=gene,residue_id=row.residue_id,position_quintile=row.position_quintile,domain_bin=row.domain_bin,surface_bin=row.surface_bin,stratum=f'{row.position_quintile}|{row.domain_bin}|{row.surface_bin}',status=status))
plan.to_csv(O/'permutation_plan.tsv',sep='\t',index=False);pd.DataFrame(strata).to_csv(O/'restricted_permutation_strata.tsv',sep='\t',index=False)
# Extension membership and partner identity stay secondary and cannot enter primary summary.
pairs=pd.read_csv(R/'results/phase2S_full/pair_status_sppider.tsv',sep='\t');ex=[]
for row in ready.loc[ready.universe.eq('SECONDARY_SPPIDER_EXTENSION8')].itertuples():
 ps=pairs.loc[pairs.gene.eq(row.gene)&pairs.status.eq('SCORED')].sort_values('partner')
 for p in ps.itertuples():ex.append(dict(gene=p.gene,partner=p.partner,score_array_path=p.score_array_path,score_array_SHA256=p.score_array_SHA256,selection='FROZEN_SPPIDER_VALID_TRACKS_NO_LABELS'))
 labels=set(str(row.retained_MONDOs).split(';'));sub=reg.loc[reg.gene.eq(row.gene)&reg.phenotype_status.eq('SINGLE_DISEASE_RESIDUE')&reg.canonical_mapping_safe.eq(True)&reg.MONDO_labels.isin(labels)] if row.PHENOTYPE_READY else reg.iloc[:0];sub.to_csv(M/(row.gene+'_secondary_extension_membership.tsv'),sep='\t',index=False)
pd.DataFrame(ex).to_csv(O/'sppider_extension_partner_set.tsv',sep='\t',index=False)
# Validate row membership, WT locations, model/identity-free matrices and manifest hashes.
for gene in s['phenotype_ready_genes']:
 mem=pd.read_csv(M/(gene+'_residue_membership.tsv'),sep='\t');assert not mem.residue_id.duplicated().any();assert mem.MONDO_labels.value_counts().min()>=10 and mem.MONDO_labels.nunique()>=2
 for model,labeled in [('xpair','primary'),('sppider','confirmatory')]:
  a=np.load(M/(gene+'_'+model+'_'+labeled+'.npz'));b=np.load(M/(gene+'_'+model+'_unlabeled_control.npz'));assert np.array_equal(a['residue_ids'],b['residue_ids']) and np.array_equal(a['residue_ids'],mem.residue_id.astype(str));assert a['scores'].shape[1]>=3 and np.isfinite(a['scores']).all();assert np.array_equal(np.sort(a['scores'],axis=1)[:,::-1],b['scores'])
for f in O.glob('*manifest.tsv'):
 if f.name.startswith('._') or f.name=='preregistration_input_manifest.tsv':continue
 d=pd.read_csv(f,sep='\t',keep_default_na=False)
 if {'path','SHA256'}.issubset(d):
  for p,h in d[['path','SHA256']].itertuples(index=False,name=None):
   if p not in ['NA','']:assert sha(R/p)==h
inputs=[]
for p in sorted(O.rglob('*')):
 if p.is_file() and not p.name.startswith('._') and p.name not in ['analysis_specification.md','preregistration_input_manifest.tsv','final_gate.json','cohort_summary.json']:inputs.append(dict(path=str(p.relative_to(R)),SHA256=sha(p),bytes=p.stat().st_size))
pd.DataFrame(inputs).to_csv(O/'preregistration_input_manifest.tsv',sep='\t',index=False);ih=sha(O/'preregistration_input_manifest.tsv')
spec=f'''# Phase 2P-0 frozen phenotype cohort and prospective statistical specification

Status: preregistered before any disease-label association calculation. Date: 2026-10-04 Asia/Shanghai. Frozen original biological verdict: NO_GO_EDGETIC_REPRESENTATION. Phase 2P-1 requires separate manual authorization even when this preparation gate is GO.

## Question, representations and methods

Does partner-specific interface context contain disease-label information beyond generic residue location, generic interface propensity and interface-score distributions without partner identity? WT sequence interface probabilities are residue-context features. They are not substitution effects, disruption probabilities or pathogenicity scores.

PRIMARY_METHOD = XPAIR_PARTNER_INTERFACE, exact X-PAIR 0.1.0 with frozen Ankh-large encoder/head. Choice is fixed by phenotype-blind Phase 2S structural validity. CONFIRMATORY_METHOD = SPPIDER_SEQ_PARTNER_INTERFACE, exact v2.0-manuscript/ESM2 receptor output. Do not fuse models or switch the primary method after association results.

## Frozen cohort and identity of every observation

Starting shared cohort: COL3A1, COL4A1, COL4A3, FGFR2, GUCY2D, KCNQ2, LHCGR, MEN1, MFN2, MMUT, PTEN, RHO, RPGR, RS1, TUBA1A. Secondary extension: COL7A1, DYNC1H1, MYO7A, NF1, RYR1, SCN1A, USH2A, VWF.

Final primary genes ({s['phenotype_ready_n']}): {', '.join(s['phenotype_ready_genes'])}. Final primary residues: {s['total_primary_residues']}. All residue IDs, canonical coordinates, exact retained MONDO labels, substitution provenance and RCVs are frozen in residue_phenotype_registry.tsv and per-gene matrices/*_residue_membership.tsv; retained groups and counts are in phenotype_group_residue_counts.tsv. The exact file contents and SHA256 hashes are part of this protocol through preregistration_input_manifest.tsv, SHA256 `{ih}`. This is an immutable input list; alterations require a separately dated amendment before association.

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
'''
(O/'analysis_specification.md').write_text(spec);h=sha(O/'analysis_specification.md');s['analysis_specification_SHA256']=h;s['primary_training_overlap_unknown_genes']=s['phenotype_ready_n'];s['primary_training_overlap_known_genes']=0
(O/'cohort_summary.json').write_text(json.dumps(s,indent=2)+'\n')
report='# Phase 2P-0 phenotype cohort freeze and statistical preregistration\n\n## Final gate\n\n`'+s['verdict']+'`\n\nThe 15 shared sequence-balanced genes yield 10 phenotype-ready genes and 547 independent gene–residue observations. No disease-label association statistic was calculated. Phase 2P-1 was NOT started. Original NO_GO_EDGETIC_REPRESENTATION was not altered.\n\n## Cohort accounting\n\n'+pd.DataFrame([s]).drop(columns=['phenotype_ready_genes']).to_markdown(index=False)+'\n\n'+ready.to_markdown(index=False)+'\n\nPrimary genes: '+', '.join(s['phenotype_ready_genes'])+'. Five starting genes fail the unique-residue sample-size gate: FGFR2, MFN2, PTEN, RHO and RPGR. The extension is a separate secondary universe.\n\n## Registry, discordance and ceiling\n\nExclusive residue classes sum to 1,034 in the starting shared 15 genes: 680 single-disease, 39 cross-disease without contributing multilabel variants, and 315 multilabel/ambiguous. The discordance audit counts ALL residues linked to ≥2 disease labels, including multilabel-variant residues; its median is 33.33%. This definition differs intentionally from the mutually exclusive CROSS_DISEASE_RESIDUE class.\n\n'+read('cross_disease_residue_audit.tsv').to_markdown(index=False)+'\n\n'+read('residue_information_ceiling.tsv').to_markdown(index=False)+'\n\nCeilings use variant-label mass with one fractional vote per variant and no score input; they are theoretical, conditional on retained ontology labels. Unsafe mapping remains recorded and excluded from primary membership. Disease-group labels are copied from frozen strict MONDO groups; no ontology update or association-dependent selection occurs. Groups with fewer than 10 eligible unique single-disease residues are prospectively excluded.\n\n## Frozen features and controls\n\nFor primary ready genes, X-PAIR and SPPIDER matrices use identical residue IDs and alphabetically frozen named common partners. Raw sorting is verified to preserve each row score multiset exactly. Mean, max, population SD and fraction≥0.5 controls are stored without any disease-specific score summary. Exact probabilities are taken from existing WT tracks; no neural inference was rerun. The common-partner median is 17; median retained group size is 16.5 residues.\n\nPosition, domain proxy, generic interface and generic 3D are frozen with explicit missing masks. The available domain files are UniProtKB features, not Pfam/InterPro; the formal Pfam/InterPro baseline is unavailable, and the protocol separates the frozen proxy. Structural missingness is retained. The primary endpoint remains labeled-versus-sorted partner identity.\n\n## Preregistered statistical plan\n\nPrimary X-PAIR; confirmatory SPPIDER; primary Euclidean distance on within-gene column-standardized vectors; independent unit gene×residue. Primary endpoint is prospective ΔR² labeled−sorted, with paired label permutations preserving exact group sizes. B=10,000 or exhaustive smaller assignment space; separate BH-FDR families for labeled association and ΔR². Global summary resamples genes 10,000 times. Support requires positive median with CI excluding zero and ≥70% genes with positive ΔR². Sensitivities, constant-dimension handling, complete-case baseline comparisons, restricted permutation strata, secondary prediction rules and overlap interpretation are fully frozen in the specification. None of these tests has been executed.\n\n'+plan.to_markdown(index=False)+'\n\n## Training overlap\n\nOfficial SPPIDER repository identifies its training/benchmark dataset at [Hugging Face](https://huggingface.co/datasets/aporollo-lab/SPPIDER-seq-datasets); X-PAIR uses its [official GitLab repository](https://gitlab.lcqb.upmc.fr/srescalli/X-PAIR). Local distribution files do not supply complete training identity lists. Direct source requests returned proxy HTTP503 and are recorded. All 15 starting genes, including 10 primary-ready genes, remain TRAINING_OVERLAP_UNKNOWN for both query and pair scopes. No non-overlap or held-out claim is made. The protocol requires another list audit before Phase 2P-1 association.\n\n## Protocol and deliverable verification\n\nAnalysis specification: `results/phase2P0/analysis_specification.md`\n\nSHA256: `'+h+'`\n\nThe preregistration input manifest freezes exact memberships, disease labels, partner identities, raw/control matrices, baseline rows, restricted strata and training-source evidence by SHA256. All 124 original frozen inputs were rehashed and passed; both Phase 2S gene lists matched exactly. Matrix checks confirm no duplicate residue IDs, ≥2 retained groups each with ≥10 residues, ≥3 common dimensions, matched model rows and exact descending sorting. All required deliverables are present.\n\n## Stop\n\nNo disease-label association statistic was calculated.\n\nPhase 2P-1 was NOT started.\n\nOriginal NO_GO_EDGETIC_REPRESENTATION was not altered.\n'
(R/'reports/PHASE2P0_PHENOTYPE_COHORT_AND_PREREGISTRATION.md').write_text(report)
required='frozen_shared15_genes.tsv sppider_extension8_genes.tsv residue_phenotype_registry.tsv cross_disease_residue_audit.tsv residue_information_ceiling.tsv phenotype_group_residue_counts.tsv phenotype_ready_genes.tsv phenotype_exclusion_log.tsv common_partner_set_by_gene.tsv xpair_primary_matrix_manifest.tsv sppider_confirmatory_matrix_manifest.tsv xpair_unlabeled_control_manifest.tsv sppider_unlabeled_control_manifest.tsv baseline_position_manifest.tsv baseline_domain_manifest.tsv baseline_generic_interface_manifest.tsv baseline_generic3d_manifest.tsv training_overlap_audit.tsv permutation_plan.tsv analysis_specification.md'.split();assert all((O/n).exists() and (O/n).stat().st_size>0 for n in required)
(O/'final_gate.json').write_text(json.dumps(dict(verdict=s['verdict'],analysis_specification_SHA256=h,ready_genes=s['phenotype_ready_genes'],required_deliverables_PASS=len(required),original_frozen_files_PASS=124,no_association_statistic_calculated=True,phase2P1_started=False,original_verdict='NO_GO_EDGETIC_REPRESENTATION',completed_at=datetime.datetime.now().astimezone().isoformat()),indent=2)+'\n');print('PHASE2P0_COMPLETED',s['verdict'],h)
