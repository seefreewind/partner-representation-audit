from pathlib import Path
import sys,os,json,time,hashlib,shutil,urllib.request,datetime,fcntl
import numpy as np,pandas as pd
from scipy.stats import spearmanr,pearsonr
from scipy.special import logsumexp
from sklearn.metrics import balanced_accuracy_score,f1_score,log_loss
from stats import *
R=Path(__file__).resolve().parents[2];P=R/'results/phase2P0';O=R/'results/phase2P1';O.mkdir(exist_ok=True);lock=(O/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert json.loads((O/'preflight_PASS.json').read_text())['all_hashes_PASS'];assert hashlib.sha256((P/'analysis_specification.md').read_bytes()).hexdigest()=='d657fab6b24f9aa757f9773bb1e9e0bd5d850e743eab6bdf819f126326ced9c3'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(rows,n):pd.DataFrame(rows).to_csv(O/n,sep='\t',index=False,na_rep='NA')
def state(stage,**kw):(O/'pipeline_state.json').write_text(json.dumps(dict(stage=stage,pid=os.getpid(),updated_at=datetime.datetime.now().astimezone().isoformat(),**kw),indent=2)+'\n')
plan=pd.read_csv(P/'permutation_plan.tsv',sep='\t',keep_default_na=False).set_index('gene');genes=json.loads((P/'cohort_summary.json').read_text())['phenotype_ready_genes'];data={};allperm={};scaled={};scaling=[];execution_hashes=[]
if '--resume' not in sys.argv:
 # Final official training-source attempt precedes result interpretation.
 state('TRAINING_OVERLAP_AUDIT');T=O/'training_sources';T.mkdir(exist_ok=True);training=[]
 for model,url in [('sppider','https://huggingface.co/api/datasets/aporollo-lab/SPPIDER-seq-datasets'),('xpair','https://gitlab.lcqb.upmc.fr/srescalli/X-PAIR/-/raw/main/README.md')]:
  try:
   response=urllib.request.urlopen(url,timeout=20).read();(T/(model+'_official_response.txt')).write_bytes(response);reason='OFFICIAL_METADATA_RETRIEVED_WITHOUT_VERIFIED_EXACT_TRAINING_IDENTITIES'
  except Exception as e:(T/(model+'_official_response.txt')).write_text(type(e).__name__+': '+str(e));reason='OFFICIAL_REQUEST_FAILED:'+str(e)
  for gene in genes:
   for scope in ['QUERY_PROTEIN','QUERY_PARTNER_PAIR']:training.append(dict(gene=gene,method=model,scope=scope,status='TRAINING_OVERLAP_UNKNOWN',source=url,reason=reason,ClinVar_disease_endpoint_used_in_training='NO_EVIDENCE_OF_PRESENT_ENDPOINT_TRAINING;MODELS_TRAINED_ON_INTERACTION_INTERFACE_TASKS',membership_changed=False))
 write(training,'training_overlap_final_audit.tsv')
# Generate one deterministic disease-label assignment array per gene; reuse across all representations.
state('PREPARE_SCALED_MATRICES_AND_PAIRED_PERMUTATIONS');(O/'paired_assignments').mkdir(exist_ok=True)
for gene in genes:
 mem=pd.read_csv(P/'matrices'/(gene+'_residue_membership.tsv'),sep='\t');labels=sorted(mem.MONDO_labels.unique());y=np.array([labels.index(x) for x in mem.MONDO_labels],dtype=np.int16);data[gene]=(mem,y,labels)
 permutation=np.load(O/'paired_assignments'/(gene+'.npz'))['label_assignments'] if '--resume' in sys.argv else permutations(y,int(plan.loc[gene,'seed']),int(plan.loc[gene,'permutation_B']))[0];allperm[gene]=permutation;path=O/'paired_assignments'/(gene+'.npz');np.savez_compressed(path,label_assignments=permutation,observed_labels=y,residue_ids=mem.residue_id.to_numpy(dtype=str),labels=np.array(labels));execution_hashes.append(dict(path=str(path.relative_to(R)),SHA256=sha(path),kind='PAIRED_UNIQUE_LABEL_ASSIGNMENTS',B=len(permutation)))
 for method,key in [('xpair','primary'),('sppider','confirmatory')]:
  a=np.load(P/'matrices'/(gene+'_'+method+'_'+key+'.npz'));b=np.load(P/'matrices'/(gene+'_'+method+'_unlabeled_control.npz'));scaled[gene,method]={}
  for typ,x in [('labeled',a['scores']),('sorted',b['scores'])]:
   z,mu,sd=standardize(x);scaled[gene,method][typ]=z;scaling.append(dict(gene=gene,method=method,representation=typ,original_dimensions=x.shape[1],usable_dimensions=int((sd>0).sum()),constant_dimensions=int((sd==0).sum()),means=json.dumps(mu.tolist()),population_SDs=json.dumps(sd.tolist()),rule='WITHIN_GENE_ALL_FROZEN_RESIDUES_DISEASE_BLIND'))
  scaled[gene,method]['raw']=a['scores'];scaled[gene,method]['raw_sorted']=b['scores']
write(scaling,'standardization_parameters.tsv')
if '--resume' not in sys.argv:
 primary=[];confirm=[];multiple=[]
 for method,outdir,output in [('xpair','xpair_permutations',primary),('sppider','sppider_permutations',confirm)]:
  state('PRIMARY_XPAIR' if method=='xpair' else 'CONFIRMATORY_SPPIDER');(O/outdir).mkdir(exist_ok=True)
  for gene in genes:
   mem,y,labels=data[gene];z=scaled[gene,method];obs,null=analyze(z['labeled'],z['sorted'],y,allperm[gene]);assert all(np.isfinite(list(obs.values()))),'DEGENERATE_PREREGISTERED_TEST'
   path=O/outdir/(gene+'.npz');np.savez_compressed(path,**null);execution_hashes.append(dict(path=str(path.relative_to(R)),SHA256=sha(path),kind=method+'_NULL',B=len(allperm[gene])))
   row=dict(gene=gene,n_total=len(y),group1=labels[0],n_group1=int((y==0).sum()),group2=labels[1],n_group2=int((y==1).sum()),common_partner_n=z['labeled'].shape[1],usable_partner_dimensions=int((z['labeled'].std(axis=0)>0).sum()),B=len(allperm[gene]),status='PREREGISTERED_COMPLETE')
   if method=='xpair':row.update(R2_XPAIR_LABELED=obs['labeled'],pseudoF_XPAIR_LABELED=obs['pseudoF'],p_XPAIR_LABELED=obs['p_labeled'],R2_XPAIR_SORTED=obs['sorted'],DELTA_R2_PARTNER_ID=obs['delta'],p_DELTA=obs['p_delta'],delta_positive=obs['delta']>0)
   else:row.update(R2_SPPIDER_LABELED=obs['labeled'],pseudoF_SPPIDER_LABELED=obs['pseudoF'],p_SPPIDER_LABELED=obs['p_labeled'],R2_SPPIDER_SORTED=obs['sorted'],DELTA_R2_SPPIDER=obs['delta'],p_DELTA=obs['p_delta'],delta_positive=obs['delta']>0)
   output.append(row)
  df=pd.DataFrame(output);labelcol='p_'+method.upper()+'_LABELED';qcol='q_'+method.upper()+'_LABELED';df[qcol]=bh(df[labelcol]);df['q_DELTA' if method=='xpair' else 'q_SPPIDER_DELTA']=bh(df.p_DELTA)
  if method=='xpair':df['q_XPAIR_DELTA']=df['q_DELTA']
  for row in df.itertuples():
   for family,pv,qv in [('LABELED',getattr(row,labelcol),getattr(row,qcol)),('DELTA',row.p_DELTA,getattr(row,'q_DELTA' if method=='xpair' else 'q_SPPIDER_DELTA'))]:multiple.append(dict(method=method,gene=row.gene,family=family,p=pv,q=qv,threshold=.05))
  write(df,'xpair_primary_per_gene.tsv' if method=='xpair' else 'sppider_confirmatory_per_gene.tsv');write(df[['gene','B',labelcol,'p_DELTA']],'xpair_permutation_summary.tsv' if method=='xpair' else 'sppider_permutation_summary.tsv')
  if method=='xpair':
   gs=global_summary(df.DELTA_R2_PARTNER_ID.to_numpy());gs['verdict']='PARTNER_IDENTITY_SUPPORT' if gs['global_support_rule_pass'] else 'NO_GLOBAL_PARTNER_IDENTITY_SUPPORT';write([gs],'xpair_global_summary.tsv');np.savez_compressed(O/'gene_bootstrap_medians.npz',medians=np.median(df.DELTA_R2_PARTNER_ID.to_numpy()[np.random.Generator(np.random.PCG64(20261005)).integers(0,10,size=(10000,10))],axis=1))
else:
 multiple=pd.read_csv(O/'multiple_testing.tsv',sep='\t').to_dict('records');gs=pd.read_csv(O/'xpair_global_summary.tsv',sep='\t').iloc[0].to_dict()
write(multiple,'multiple_testing.tsv');xp=pd.read_csv(O/'xpair_primary_per_gene.tsv',sep='\t');sp=pd.read_csv(O/'sppider_confirmatory_per_gene.tsv',sep='\t');delta=xp.DELTA_R2_PARTNER_ID.to_numpy();sd=sp.DELTA_R2_SPPIDER.to_numpy();write([dict(Spearman=spearmanr(delta,sd).statistic,Pearson=pearsonr(delta,sd).statistic,same_sign_n=int((np.sign(delta)==np.sign(sd)).sum()),positive_both_n=int(((delta>0)&(sd>0)).sum()),XPAIR_positive_only_n=int(((delta>0)&(sd<=0)).sum()),SPPIDER_positive_only_n=int(((delta<=0)&(sd>0)).sum()),n_genes=10)],'cross_method_concordance.tsv')
# Complete-case generic baselines, reporting missingness and row-matched XPAIR results.
state('GENERIC_BASELINES');baseline=[]
for gene in genes:
 mem,y,labels=data[gene];raw=scaled[gene,'xpair']['raw']
 for name in ['position','domain','generic_interface','generic3d']:
  d=pd.read_csv(P/'matrices'/(gene+'_baseline_'+name+'.tsv'),sep='\t',keep_default_na=False);assert d.residue_id.tolist()==mem.residue_id.tolist();cols={'position':['normalized_position'],'domain':['relative_position_within_domain','distance_to_nearest_domain_boundary_aa'],'generic_interface':['any_interface','known_interfaces_touching_residue','minimum_interface_distance'],'generic3d':['AFDB_plddt','monomer_SASA_A2','intrachain_CA_neighbor_n_8A']}[name];num=d[cols].replace({'True':1.,'False':0.}).apply(pd.to_numeric,errors='coerce');mask=num.notna().all(axis=1).to_numpy(copy=True)
  if name=='domain':mask&=~d.domain_identity.isin(['','NA']).to_numpy()
  yy=y[mask];counts=np.bincount(yy,minlength=len(labels));row=dict(gene=gene,baseline='DOMAIN_PROXY' if name=='domain' else name.upper(),n_primary=len(y),n_complete_case=len(yy),excluded_missing=int((~mask).sum()),group_counts=json.dumps(counts.tolist()))
  if len(yy)==0 or (counts<10).any():row.update(status='UNEVALUABLE_COMPLETE_CASE_GROUP_GATE',R2_XPAIR_MATCHED=np.nan,R2_BASELINE=np.nan,DELTA_R2_BASELINE=np.nan)
  else:
   x=num.loc[mask].to_numpy(dtype=float)
   if name=='domain':x=np.column_stack([x,pd.get_dummies(d.loc[mask,'domain_identity'],dtype=float).to_numpy()])
   xb=standardize(x)[0];xr=standardize(raw[mask])[0];rxb=float(r2(gram(xr),yy)[0]);rb=float(r2(gram(xb),yy)[0]);row.update(status='EVALUABLE' if np.isfinite(rb) else 'DEGENERATE_BASELINE',R2_XPAIR_MATCHED=rxb,R2_BASELINE=rb,DELTA_R2_BASELINE=rxb-rb,transformation='EUCLIDEAN_DISEASE_BLIND_NUMERIC_COLUMN_SCALING;DOMAIN_ID_ONEHOT')
  baseline.append(row)
write(baseline,'baseline_comparisons.tsv')
# Exact frozen restricted strata, using unique or exhaustive alternatives.
state('RESTRICTED_SENSITIVITY');strata=pd.read_csv(P/'restricted_permutation_strata.tsv',sep='\t');restricted=[];(O/'restricted_permutations').mkdir(exist_ok=True)
for gene in genes:
 mem,y,labels=data[gene];st=strata.loc[strata.gene.eq(gene)].set_index('residue_id').loc[mem.residue_id,'stratum'].tolist();pp=plan.loc[gene];perm,space=permutations(y,(int(pp.seed)+1)%(2**32),int(pp.restricted_B),st);assert space==int(pp.restricted_assignment_space) and len(perm)==int(pp.restricted_B)
 z=scaled[gene,'xpair'];obs,null=analyze(z['labeled'],z['sorted'],y,perm);path=O/'restricted_permutations'/(gene+'.npz');np.savez_compressed(path,label_assignments=perm,**null);execution_hashes.append(dict(path=str(path.relative_to(R)),SHA256=sha(path),kind='RESTRICTED_NULL',B=len(perm)));restricted.append(dict(gene=gene,restricted_p_XPAIR_LABELED=obs['p_labeled'],restricted_p_DELTA=obs['p_delta'],restricted_status=pp.restricted_status,B=len(perm),assignment_space=str(space),exhaustive=space-1==len(perm),unrestricted_p_DELTA=float(xp.loc[xp.gene.eq(gene),'p_DELTA'].iloc[0]),nominal_delta_concordance=(obs['p_delta']<.05)==(float(xp.loc[xp.gene.eq(gene),'p_DELTA'].iloc[0])<.05)))
write(restricted,'restricted_permutation_results.tsv')
state('DISTANCE_SENSITIVITY');corr=[];rawsens=[];(O/'sensitivity_permutations').mkdir(exist_ok=True)
for gene in genes:
 mem,y,labels=data[gene];z=scaled[gene,'xpair'];G,info=correlation_gram(z['raw']);H,info2=correlation_gram(z['raw_sorted']);row=dict(gene=gene,**info,sorted_status=info2['status'])
 if G is not None and H is not None and info['status']=='EVALUABLE_LINGOES' and info2['status']=='EVALUABLE_LINGOES':
  a=float(r2(G,y)[0]);b=float(r2(H,y)[0]);na=r2(G,allperm[gene]);nb=r2(H,allperm[gene]);row.update(R2_labeled=a,R2_sorted=b,DELTA_R2=a-b,p_labeled=empirical(F(na,len(y),2),F(a,len(y),2)),p_DELTA=empirical(na-nb,a-b),sorted_Lingoes_constant=info2.get('Lingoes_constant',np.nan));path=O/'sensitivity_permutations'/(gene+'_correlation.npz');np.savez_compressed(path,R2_labeled_perm=na,R2_sorted_perm=nb,delta_R2_perm=na-nb)
 else:row['status']='UNEVALUABLE_ZERO_TOTAL_SORTED_DISTANCE' if info2['status']=='DEGENERATE_ZERO_TOTAL_DISTANCE' else 'UNEVALUABLE_DIMENSION_OR_ROW_GATE'
 corr.append(row)
 for name,a,b in [('RAW_PROBABILITY_EUCLIDEAN',z['raw'],z['raw_sorted']),('SORT_AFTER_STANDARDIZATION_EUCLIDEAN',z['labeled'],np.sort(z['labeled'],axis=1)[:,::-1])]:
  obs,null=analyze(a,b,y,allperm[gene]);rawsens.append(dict(gene=gene,sensitivity=name,**obs));np.savez_compressed(O/'sensitivity_permutations'/(gene+'_'+name+'.npz'),**null)
write(corr,'correlation_distance_sensitivity.tsv');write(rawsens,'raw_and_poststandardization_sensitivity.tsv')
ceil=pd.read_csv(P/'residue_information_ceiling.tsv',sep='\t');discord=pd.read_csv(P/'cross_disease_residue_audit.tsv',sep='\t');context=xp[['gene','DELTA_R2_PARTNER_ID']].merge(ceil,on='gene').merge(discord,on='gene');write(context,'residue_ceiling_context.tsv')
# Prediction only after all shared association results are finalized. Train-only standardization, paired folds.
state('PREDICTIVE_SECONDARY');pred=[];foldlog=[];rng=np.random.Generator(np.random.PCG64(20261006));(O/'predictive_outputs').mkdir(exist_ok=True)
for gene in genes:
 mem,y,labels=data[gene];n=len(y);classes=np.unique(y);outputs={name:[] for name in ['XPAIR_LABELED','XPAIR_SORTED','SPPIDER_LABELED','SPPIDER_SORTED']}
 if min(np.bincount(y))<10:pred.append(dict(gene=gene,status='SAMPLE_SIZE_UNEVALUABLE'));continue
 for repeat in range(10):
  fold=np.zeros(n,dtype=int)
  for c in classes:
   ix=rng.permutation(np.flatnonzero(y==c))
   for k,chunk in enumerate(np.array_split(ix,5)):fold[chunk]=k
  for rowid,k in zip(mem.residue_id,fold):foldlog.append(dict(gene=gene,repeat=repeat,residue_id=rowid,validation_fold=int(k)))
  for method,typ,name in [('xpair','raw','XPAIR_LABELED'),('xpair','raw_sorted','XPAIR_SORTED'),('sppider','raw','SPPIDER_LABELED'),('sppider','raw_sorted','SPPIDER_SORTED')]:
   x=scaled[gene,method][typ];probs=np.zeros((n,len(classes)))
   for k in range(5):
    tr=fold!=k;te=fold==k;_,mu,sd=standardize(x[tr]);ok=sd>0;zt=np.zeros_like(x,dtype=float);zt[:,ok]=(x[:,ok]-mu[ok])/sd[ok];centers=np.array([zt[tr&(y==c)].mean(axis=0) for c in classes]);distance=((zt[te,None,:]-centers[None,:,:])**2).sum(axis=2);logits=-distance;pp=np.exp(logits-logsumexp(logits,axis=1)[:,None]);pp=np.maximum(pp,1e-15);pp/=pp.sum(axis=1,keepdims=True);probs[te]=pp
   outputs[name].append(probs);guess=probs.argmax(axis=1);pred.append(dict(gene=gene,representation=name,repeat=repeat,status='SECONDARY_5FOLD_PAIRED_TRAIN_ONLY_SCALING',balanced_accuracy=balanced_accuracy_score(y,guess),macro_F1=f1_score(y,guess,average='macro'),log_loss=log_loss(y,probs,labels=classes)))
 for name,probs in outputs.items():np.savez_compressed(O/'predictive_outputs'/(gene+'_'+name+'.npz'),probabilities=np.array(probs),labels=y,residue_ids=mem.residue_id.to_numpy(dtype=str))
write(pred,'predictive_sensitivity.tsv');write(foldlog,'predictive_fold_manifest.tsv')
# Frozen extension cohort, separate model, partner set and FDR; never merged into primary.
state('SPPIDER_ONLY_EXTENSION');extension=[];ready=pd.read_csv(P/'phenotype_ready_genes.tsv',sep='\t');ep=pd.read_csv(P/'sppider_extension_partner_set.tsv',sep='\t');(O/'sppider_extension_permutations').mkdir(exist_ok=True)
for row in ready.loc[ready.universe.eq('SECONDARY_SPPIDER_EXTENSION8')].itertuples():
 gene=row.gene
 if not row.PHENOTYPE_READY:extension.append(dict(gene=gene,status='FROZEN_SAMPLE_SIZE_GATE_FAIL',analysis='SPPIDER_ONLY_EXTENSION'));continue
 mem=pd.read_csv(P/'matrices'/(gene+'_secondary_extension_membership.tsv'),sep='\t');labels=sorted(mem.MONDO_labels.unique());y=np.array([labels.index(x) for x in mem.MONDO_labels],dtype=np.int16);ps=ep.loc[ep.gene.eq(gene)].sort_values('partner');arrays=[]
 for p in ps.itertuples():
  path=R/p.score_array_path;assert sha(path)==p.score_array_SHA256;arrays.append(np.load(path)['receptor'][mem.canonical_residue.to_numpy(dtype=int)-1])
 raw=np.column_stack(arrays);a=standardize(raw)[0];b=standardize(np.sort(raw,axis=1)[:,::-1])[0];seed=int(hashlib.sha256(('PHASE2P1|20261004|'+gene).encode()).hexdigest()[:16],16)%(2**32);perm,space=permutations(y,seed,10000);obs,null=analyze(a,b,y,perm);extension.append(dict(gene=gene,analysis='SPPIDER_ONLY_EXTENSION',status='COMPLETE',n_residues=len(y),group_sizes=json.dumps(np.bincount(y).tolist()),partners=len(ps),seed=seed,B=len(perm),**obs));np.savez_compressed(O/'sppider_extension_permutations'/(gene+'.npz'),label_assignments=perm,**null)
ed=pd.DataFrame(extension);mask=ed.status.eq('COMPLETE');ed.loc[mask,'q_labeled_extension']=bh(ed.loc[mask,'p_labeled']);ed.loc[mask,'q_delta_extension']=bh(ed.loc[mask,'p_delta']);write(ed,'sppider_extension_results.tsv')
write(execution_hashes,'permutation_array_manifest.tsv');claim=dict(gs);claim.update(genes_q_delta_lt_0_05=';'.join(xp.loc[xp.q_DELTA.lt(.05),'gene']),genes_q_labeled_lt_0_05=';'.join(xp.loc[xp.q_XPAIR_LABELED.lt(.05),'gene']),sppider_confirmation_summary='median_delta='+str(float(np.median(sd)))+';positive='+str(int((sd>0).sum())),restricted_sensitivity_summary='nominal_delta_concordant='+str(sum(x['nominal_delta_concordance'] for x in restricted))+'/10',training_overlap_status='TRAINING_OVERLAP_UNKNOWN',original_verdict='NO_GO_EDGETIC_REPRESENTATION',post_hoc_selection=False);write([claim],'final_claim_audit.tsv');state('COMPUTATIONS_COMPLETE',verdict=gs['verdict']);print('PHASE2P1_COMPUTATIONS_COMPLETE',flush=True)
