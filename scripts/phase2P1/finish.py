from pathlib import Path
import json,hashlib,datetime,shutil
import numpy as np,pandas as pd
from stats import *
R=Path(__file__).resolve().parents[2];O=R/'results/phase2P1';P=R/'results/phase2P0'
assert json.loads((O/'pipeline_state.json').read_text())['stage']=='COMPUTATIONS_COMPLETE'
def read(n):return pd.read_csv(O/n,sep='\t')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
 return h.hexdigest()
x=read('xpair_primary_per_gene.tsv');s=read('sppider_confirmatory_per_gene.tsv');g=read('xpair_global_summary.tsv').iloc[0];c=read('cross_method_concordance.tsv').iloc[0];b=read('baseline_comparisons.tsv');res=read('restricted_permutation_results.tsv');corr=read('correlation_distance_sensitivity.tsv');pred=read('predictive_sensitivity.tsv');ext=read('sppider_extension_results.tsv');audit=[]
def check(name,value,detail=''):
 audit.append(dict(check=name,PASS=bool(value),detail=str(detail)))
 if not value:raise RuntimeError(name+': '+str(detail))
check('EXACT_PREREGISTRATION_SHA256',sha(P/'analysis_specification.md')=='d657fab6b24f9aa757f9773bb1e9e0bd5d850e743eab6bdf819f126326ced9c3')
check('TEN_PRIMARY_GENES_547_RESIDUES',len(x)==10 and x.n_total.sum()==547 and x.gene.tolist()==s.gene.tolist())
for row in x.itertuples():
 gene=row.gene;assign=np.load(O/'paired_assignments'/(gene+'.npz'));y=assign['observed_labels'];perm=assign['label_assignments'];check(gene+'_UNIQUE_GROUP_PRESERVING_PERMUTATIONS',len(perm)==10000 and len({r.tobytes() for r in perm})==10000 and all(np.array_equal(np.bincount(r),np.bincount(y)) for r in perm) and not any(np.array_equal(r,y) for r in perm))
 for method,key in [('xpair','primary'),('sppider','confirmatory')]:
  raw=np.load(P/'matrices'/(gene+'_'+method+'_'+key+'.npz'))['scores'];sortedraw=np.load(P/'matrices'/(gene+'_'+method+'_unlabeled_control.npz'))['scores'];check(gene+'_'+method+'_EXACT_RAW_MULTISET',np.array_equal(np.sort(raw,axis=1)[:,::-1],sortedraw))
  z=standardize(raw)[0];zz=standardize(sortedraw)[0];obs,null=analyze(z,zz,y,perm);stored=np.load(O/(method+'_permutations')/(gene+'.npz'));check(gene+'_'+method+'_RAW_NULL_RECOMPUTATION',all(np.allclose(null[k],stored[k],rtol=1e-10,atol=1e-12) for k in null))
  # Independent sums-of-squares calculation from coordinates verifies the Gram form.
  centered=z-z.mean(axis=0);total=(centered**2).sum();between=sum((y==k).sum()*((z[y==k].mean(axis=0)-z.mean(axis=0))**2).sum() for k in np.unique(y));check(gene+'_'+method+'_PERMANOVA_SS_EQUIVALENCE',np.isclose(obs['labeled'],between/total,rtol=1e-10,atol=1e-12))
  rr=x.loc[x.gene.eq(gene)].iloc[0] if method=='xpair' else s.loc[s.gene.eq(gene)].iloc[0];prefix=method.upper();check(gene+'_'+method+'_OBSERVED_AND_P_MATCH',np.isclose(obs['labeled'],rr['R2_'+prefix+'_LABELED']) and np.isclose(obs['sorted'],rr['R2_'+prefix+'_SORTED']) and np.isclose(obs['p_delta'],rr.p_DELTA))
 check(gene+'_BH_DELTA_ALIAS',np.isclose(row.q_DELTA,row.q_XPAIR_DELTA))
for df,method in [(x,'XPAIR'),(s,'SPPIDER')]:
 check(method+'_SEPARATE_BH_FAMILIES',np.allclose(df['q_'+method+'_LABELED'],bh(df['p_'+method+'_LABELED'])) and np.allclose(df['q_'+method+'_DELTA'],bh(df.p_DELTA)))
gs=global_summary(x.DELTA_R2_PARTNER_ID.to_numpy());check('GLOBAL_RULE_EXACT',gs['global_support_rule_pass']==bool(g.global_support_rule_pass));check('BOOTSTRAP_INTERVAL_EXACT',np.isclose(gs['bootstrap_CI_low'],g.bootstrap_CI_low) and np.isclose(gs['bootstrap_CI_high'],g.bootstrap_CI_high))
fold=read('predictive_fold_manifest.tsv');check('PREDICTIVE_TEN_REPEATS_NO_RESIDUE_FOLD_DUPLICATION',len(fold)==5470 and not fold.duplicated(['gene','repeat','residue_id']).any() and set(fold.validation_fold)==set(range(5)))
check('SPPIDER_EXTENSION_SEPARATE',len(ext)==8 and set(ext.gene).isdisjoint(x.gene) and ext.status.eq('COMPLETE').sum()==7)
# Baseline summary uses available, valid complete-case marginal R²; no baseline selects the primary endpoint.
valid=b[b.status.eq('EVALUABLE')];summary=valid.groupby('baseline').agg(n_genes=('gene','count'),median_baseline_R2=('R2_BASELINE','median'),median_row_matched_XPAIR_R2=('R2_XPAIR_MATCHED','median'),median_XPAIR_increment=('DELTA_R2_BASELINE','median')).reset_index();summary.to_csv(O/'baseline_summary.tsv',sep='\t',index=False);strong=summary.sort_values(['median_baseline_R2','baseline'],ascending=[False,True]).iloc[0]
predsum=pred.groupby(['gene','representation'])[['balanced_accuracy','macro_F1','log_loss']].mean().reset_index();predsum.to_csv(O/'predictive_gene_summary.tsv',sep='\t',index=False)
pglobal=predsum.groupby('representation')[['balanced_accuracy','macro_F1','log_loss']].median().reset_index();pglobal.to_csv(O/'predictive_global_descriptive_summary.tsv',sep='\t',index=False)
# Preserve complete raw numerical output manifests, including extension/sensitivities and CV probabilities.
manifest=[]
for path in sorted(O.rglob('*.npz')):
 if path.name.startswith('._'):continue
 arr=np.load(path);manifest.append(dict(path=str(path.relative_to(R)),SHA256=sha(path),arrays=json.dumps({k:list(arr[k].shape) for k in arr.files})))
pd.DataFrame(manifest).to_csv(O/'numerical_output_manifest.tsv',sep='\t',index=False)
T=O/'tables';T.mkdir(exist_ok=True)
for idx,name in enumerate(['xpair_primary_per_gene.tsv','sppider_confirmatory_per_gene.tsv','baseline_comparisons.tsv','restricted_permutation_results.tsv','training_overlap_final_audit.tsv','predictive_sensitivity.tsv'],1):shutil.copyfile(O/name,T/('TableP'+str(idx)+'.tsv'))
# Final read-only hash pass over every frozen upstream input.
orig=read('input_hash_verification.tsv');print(orig.columns.tolist(),flush=True)
for row in orig.to_dict('records'):
 path=Path(row.get('path',row.get('file','')))
 if not path.is_absolute():path=R/path
 expected=row.get('expected_SHA256',row.get('expected_sha256',row.get('SHA256')))
 check('UPSTREAM_FINAL_'+str(path.relative_to(R)) if path.is_relative_to(R) else 'UPSTREAM_FINAL_EXTERNAL',sha(path)==expected)
pd.DataFrame(audit).to_csv(O/'statistical_integrity_audit.tsv',sep='\t',index=False)
# Saved details for the reporting and figure stage.
final=dict(verdict=g.verdict,preregistration_SHA256=sha(P/'analysis_specification.md'),primary_genes=10,primary_residues=547,median_common_partners=float(x.common_partner_n.median()),median_XPAIR_labeled_R2=float(x.R2_XPAIR_LABELED.median()),median_XPAIR_sorted_R2=float(x.R2_XPAIR_SORTED.median()),**gs,genes_q_XPAIR_DELTA_lt_0_05=x.loc[x.q_XPAIR_DELTA.lt(.05),'gene'].tolist(),genes_q_XPAIR_LABELED_lt_0_05=x.loc[x.q_XPAIR_LABELED.lt(.05),'gene'].tolist(),median_SPPIDER_delta=float(s.DELTA_R2_SPPIDER.median()),cross_method_Spearman=float(c.Spearman),positive_both=int(c.positive_both_n),restricted_nominal_concordant=int(res.nominal_delta_concordance.sum()),correlation_evaluable=int(corr.status.eq('EVALUABLE_LINGOES').sum()),strongest_generic_baseline=strong.to_dict(),training_overlap_status='TRAINING_OVERLAP_UNKNOWN',predictive_summary=pglobal.to_dict('records'),extension_complete_n=int(ext.status.eq('COMPLETE').sum()),extension_positive_n=int(ext.loc[ext.status.eq('COMPLETE'),'delta'].gt(0).sum()),extension_q_delta_genes=ext.loc[ext.q_delta_extension.lt(.05),'gene'].tolist(),original_verdict='NO_GO_EDGETIC_REPRESENTATION',post_hoc_selection=False)
(O/'reporting_summary.json').write_text(json.dumps(final,indent=2,default=lambda v:v.item())+'\n');print('STATISTICAL_AUDIT_COMPLETE',flush=True)
