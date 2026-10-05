from common import *
import numpy as np,json,re,resource
from statsmodels.discrete.discrete_model import Logit
from statsmodels.tools.tools import add_constant
u=rd('results/phase2S/frozen_gene_universes.tsv');p=rd('results/phase2S_full/pair_plan.tsv');cov=pd.read_csv(O/'coverage_by_gene_all_methods.tsv',sep='\t');a=rd('results/phase1A/coverage_by_gene.tsv').set_index('gene');variants=rd('results/phase2S/frozen_variants.tsv');seqs=rd('results/phase2SR/canonical_sequence_recovery.tsv').set_index('requested_accession');features=[]
for g in sorted(u.gene):
    pp=p.loc[p.gene.eq(g)];v=variants.loc[variants.gene.eq(g)];acc=pp.query_accession.iloc[0] if len(pp) else next(iter(v.canonical_UniProt_accession),'');length=int(seqs.loc[acc,'sequence_length']) if acc in seqs.index else np.nan
    ar=a.loc[g];partners=pd.to_numeric(pp.partner_length,errors='coerce');features.append(dict(gene=g,query_length=length,median_partner_length=partners.loc[partners>0].median(),credible_partner_n=len(pp),Tier1_2_partner_count=int(ar.Tier1_2_partner_n),PDB_variant_mapping_fraction=float(ar.variants_with_observed_PDB_residue_coordinates_n)/max(1,int(ar.unique_strict_group_variants_n)),structurally_mappable_variant_fraction=float(ar.fraction_structurally_addressable_pct)/100,disease_variant_count=int(ar.unique_strict_group_variants_n),disease_group_count=int(ar.strict_disease_group_n),canonical_mapping_completeness=float(ar.protein_mapping_rate_pct)/100,query_XPAIR_length_supported=bool(50<=length<=2000) if np.isfinite(length) else False,disorder_fraction='NA_NOT_PREEXISTING'))
f=pd.DataFrame(features);wr(f,'failure_predictor_features.tsv');models=[]
for method,c in cov.groupby('method'):
    d=c.merge(f,on='gene',validate='one_to_one');d['C4_failure']=~d.C4_ge70;d['balance_failure']=~d.balanced_PASS
    for outcome in ['C4_failure','balance_failure']:
        for name in ['query_length','median_partner_length','credible_partner_n','Tier1_2_partner_count','PDB_variant_mapping_fraction','structurally_mappable_variant_fraction','disease_variant_count','disease_group_count','canonical_mapping_completeness']:
            part=d[[name,outcome]].dropna();y=part[outcome].astype(int);x=part[name];neg=x.loc[y.eq(0)];pos=x.loc[y.eq(1)];std=x.std(ddof=0)
            row=dict(method=method,outcome=outcome,analysis='UNIVARIABLE_DESCRIPTIVE',feature=name,n=len(part),failure_n=int(y.sum()),failure_median=pos.median(),pass_median=neg.median(),standardized_mean_difference=(pos.mean()-neg.mean())/std if std>0 and len(pos) and len(neg) else np.nan,status='DESCRIPTIVE_NO_CAUSAL_INTERPRETATION')
            models.append(row)
        chosen=['query_length','credible_partner_n','canonical_mapping_completeness'];part=d[chosen+[outcome]].dropna();y=part[outcome].astype(int);minor=min(int(y.sum()),len(y)-int(y.sum()));k=min(3,minor//10)
        use=[name for name in chosen if part[name].std(ddof=0)>0][:k]
        if use and len(y)>=30:
            x=(part[use]-part[use].mean())/part[use].std(ddof=0)
            try:
                fit=Logit(y,add_constant(x)).fit(disp=False,maxiter=200)
                if not fit.mle_retvals['converged'] or np.linalg.cond(add_constant(x))>1e6 or (fit.bse>10).any():raise ValueError('Separation/nonidentifiability')
                for name in use:models.append(dict(method=method,outcome=outcome,analysis='MULTIVARIABLE_LOGISTIC',feature=name,n=len(y),failure_n=int(y.sum()),OR_per_SD=np.exp(fit.params[name]),CI95_low=np.exp(fit.conf_int().loc[name,0]),CI95_high=np.exp(fit.conf_int().loc[name,1]),status='IDENTIFIABLE_EP_VARIABLE_GE10'))
            except Exception as e:models.append(dict(method=method,outcome=outcome,analysis='MULTIVARIABLE_LOGISTIC',status='NOT_IDENTIFIABLE',reason=str(e)))
        else:models.append(dict(method=method,outcome=outcome,analysis='MULTIVARIABLE_LOGISTIC',status='INSUFFICIENT_EVENTS_OR_VARIATION',failure_n=minor))
wr(pd.DataFrame(models),'failure_predictor_models.tsv')
# Pair-level and variant-level technical attrition: descriptive, no independent-residue inference.
causes=[];stages=[];cost=[]
for method in ['sppider','xpair']:
    status=pd.read_csv(O/('pair_status_'+method+'.tsv'),sep='\t');sc=status.loc[status.status.eq('SCORED')]
    for reason,n in status.status.value_counts().items():causes.append(dict(method=method,unit='FROZEN_GENE_PARTNER_PAIR',cause=reason,n=n,denominator=len(status),fraction=n/len(status)))
    for name in ['query_length','partner_length']:
        for failed,sub in [(False,sc),(True,status.loc[~status.status.eq('SCORED')])]:stages.append(dict(method=method,unit='PAIR',outcome='failure' if failed else 'scored',feature=name,n=len(sub),median=pd.to_numeric(sub[name]).median()))
    vc=pd.read_csv(O/'coverage_by_variant_all_methods.tsv',sep='\t');fullmethod='SPPIDER_SEQ_PARTNER_INTERFACE' if method=='sppider' else 'XPAIR_PARTNER_INTERFACE';vc=vc.loc[vc.method.eq(fullmethod)].merge(f,on='gene')
    for failed,sub in [(False,vc.loc[vc.GE3]),(True,vc.loc[~vc.GE3])]:
        for name in ['query_length','credible_partner_n','PDB_variant_mapping_fraction','canonical_mapping_completeness']:stages.append(dict(method=method,unit='VARIANT',outcome='GE3_failure' if failed else 'GE3_pass',feature=name,n=len(sub),median=sub[name].median(),independent_unit='NO_INFERENTIAL_CI_VARIANTS_CLUSTERED_IN_GENE'))
    files={R/r.score_array_path for r in sc.itertuples(index=False)};size=sum(z.stat().st_size for z in files)
    embeddings=sum(z.stat().st_size for z in (O/'embedding_cache'/method).glob('*.pt') if not z.name.startswith('._'))
    wall=pd.to_numeric(sc['wall_seconds'],errors='coerce') if 'wall_seconds' in sc else pd.to_numeric(sc['allocated_batch_wall_seconds'],errors='coerce');cpu=pd.to_numeric(status['cpu_seconds'],errors='coerce').sum()/3600 if 'cpu_seconds' in status else np.nan
    cost.append(dict(method=method,frozen_pairs_total=len(status),unique_pairs_attempted=int(status[method+'_preflight'].eq('PENDING').sum()),pairs_scored=len(sc),CPU_hours=cpu,GPU_hours=0 if method=='sppider' else np.nan,peak_RAM_bytes=pd.to_numeric(status.peak_RAM_bytes).max() if 'peak_RAM_bytes' in status else np.nan,peak_VRAM_bytes=0 if method=='sppider' else np.nan,residue_cache_bytes=size,embedding_cache_bytes=embeddings,median_time_per_pair_seconds=wall.median(),median_time_per_query_residue_seconds=(wall/pd.to_numeric(sc.query_length)).median(),time_definition='wall measured per pair including embedding/cache IO' if method=='sppider' else 'amortized batch wall proxy; not isolated forward latency',unavailable_resource_reason='' if method=='sppider' else 'MPS GPU busy time/VRAM and child CPU require profiling; report measured batch wall separately, never fabricate'))
# Instrumented official batches provide actual CPU and forward timing; GPU busy
# profiling remains unavailable on MPS and is never replaced by wall elapsed.
measurefile=O/'xpair_resource_batches.jsonl'
if measurefile.exists():
    measurements={r['batch']:r for r in (json.loads(x) for x in measurefile.read_text().splitlines())}
    row=next(r for r in cost if r['method']=='xpair');row['CPU_hours']=sum(x['CPU_seconds'] for x in measurements.values())/3600;row['peak_RAM_bytes']=max(x['peak_RAM_bytes'] for x in measurements.values())
    durations=[t['forward_wall_seconds'] for x in measurements.values() for t in x['pair_timings']]
    row['median_measured_forward_seconds_per_pair']=np.median(durations) if durations else np.nan;row['batch_wall_hours']=sum(x['total_wall_seconds'] for x in measurements.values())/3600
    row['measurement_scope']='instrumented batches; recovery-pilot cache generation excluded'
# Original computational accounting is reused when available, otherwise NA.
for method in ['MUTPRED_PPI_DISRUPTION','STRUCTURAL_BINARY_INTERFACE','STRUCTURAL_INTERFACE_DISTANCE']:cost.append(dict(method=method,status='FROZEN_OUTPUT_REUSED_NO_RECOMPUTATION',CPU_hours=0,GPU_hours=0,original_generation_cost='NA_NOT_COMPARABLE_WITH_NEW_INFERENCE'))
wr(pd.DataFrame(causes),'failure_causes_by_method.tsv');wr(pd.DataFrame(stages),'failure_stage_descriptive_effects.tsv');wr(pd.DataFrame(cost),'computational_cost.tsv')
long=[];xs=pd.read_csv(O/'pair_status_xpair.tsv',sep='\t');sp=pd.read_csv(O/'pair_status_sppider.tsv',sep='\t');audit=[]
for g in sorted(u.gene):
    pp=p.loc[p.gene.eq(g)];ss=sp.loc[sp.gene.eq(g)];xx=xs.loc[xs.gene.eq(g)];length=f.loc[f.gene.eq(g),'query_length'].iloc[0];before=pp.sppider_preflight.eq('PENDING').sum();after=pp.xpair_preflight.eq('PENDING').sum()
    long.append(dict(gene=g,query_length=length,gt1022=length>1022,gt2000=length>2000,gt3000=length>3000,sppider_scored_partners=int(ss.status.eq('SCORED').sum()),sppider_remains_scoreable=ss.status.eq('SCORED').any()))
    audit.append(dict(gene=g,query_length=length,query_gt2000=length>2000,query_lt50=length<50,sequence_supported_pairs=before,length_supported_pairs=after,query_length_rejected_pairs=int(((pd.to_numeric(pp.query_length)<50)|(pd.to_numeric(pp.query_length)>2000)).loc[pp.sppider_preflight.eq('PENDING')].sum()),partner_length_rejected_pairs=int(((pd.to_numeric(pp.partner_length)<50)|(pd.to_numeric(pp.partner_length)>2000)).loc[pp.sppider_preflight.eq('PENDING')].sum()),both_length_rejected_pairs=int((((pd.to_numeric(pp.query_length)<50)|(pd.to_numeric(pp.query_length)>2000))&((pd.to_numeric(pp.partner_length)<50)|(pd.to_numeric(pp.partner_length)>2000))).loc[pp.sppider_preflight.eq('PENDING')].sum()),length_only_GE3_unrepresentable=before>=3 and after<3,length_only_completely_unavailable=before>0 and after==0))
wr(pd.DataFrame(long),'sppider_long_sequence_audit.tsv');wr(pd.DataFrame(audit),'xpair_length_constraint_by_gene.tsv');print('FAILURE_AND_COST_COMPLETE',flush=True)
