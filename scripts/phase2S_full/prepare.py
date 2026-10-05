from common import *
import numpy as np,shutil
O.mkdir(exist_ok=True);L.mkdir(exist_ok=True)
f=rd('metadata/PHASE2SR_INPUT_FREEZE_MANIFEST.tsv');checks=[]
for r in f.itertuples(index=False):
    actual=digest(R/r.path);checks.append({'path':r.path,'expected_SHA256':r.SHA256,'actual_SHA256':actual,'status':'PASS' if actual==r.SHA256 else 'INPUT_FREEZE_VIOLATION'})
wr(pd.DataFrame(checks),'frozen_input_hash_validation.tsv')
assert len(checks)==124 and all(x['status']=='PASS' for x in checks),'INPUT_FREEZE_VIOLATION'
u=rd('results/phase2S/frozen_gene_universes.tsv');p=rd('results/phase2S/frozen_partner_sets.tsv');m=rd('results/phase2SR/common_primary_membership.tsv')
assert len(u)==132 and u.U1_STRUCTURAL_42.eq('True').sum()==42 and len(p)==15942 and m.included_in_common_primary_denominator.eq('True').sum()==9978,'INPUT_FREEZE_VIOLATION'
assert not p[['gene','partner']].duplicated().any()
for src in ['results/phase2SR/common_primary_membership.tsv','results/phase2S/frozen_gene_universes.tsv','results/phase2S/frozen_partner_sets.tsv','results/phase2S/frozen_variants.tsv','results/phase2S/frozen_disease_groups.tsv']:
    dest=O/Path(src).name
    if dest.exists():assert digest(dest)==digest(R/src),'INPUT_FREEZE_VIOLATION'
    else:shutil.copy2(R/src,dest)
s=rd('results/phase2SR/canonical_sequence_recovery.tsv').set_index('requested_accession');assert len(s)==7300 and s.status.eq('EXACT').all()
plan=[]
for r in p.itertuples(index=False):
    a,b=r.query_accession,r.partner_accession;qs='';ps='';ql=pl=0
    state='PENDING'
    if not a or a not in s.index:state='QUERY_SEQUENCE_MISSING'
    elif not b or b not in s.index:state='PARTNER_SEQUENCE_MISSING'
    else:qs=s.loc[a,'sequence_sha256'];ps=s.loc[b,'sequence_sha256'];ql=int(s.loc[a,'sequence_length']);pl=int(s.loc[b,'sequence_length'])
    xp=state if state!='PENDING' else 'PENDING' if 50<=ql<=2000 and 50<=pl<=2000 else 'MODEL_LENGTH_UNSUPPORTED'
    plan.append({'gene':r.gene,'partner':r.partner,'query_accession':a,'partner_accession':b,'query_sha256':qs,'partner_sha256':ps,'query_length':ql,'partner_length':pl,'sppider_preflight':state,'xpair_preflight':xp,'tier':r.tier,'evidence_stratum':r.evidence_stratum})
plan=pd.DataFrame(plan).sort_values(['gene','partner']);wr(plan,'pair_plan.tsv')
# Predeclare all secondary analysis choices before full outcomes.
contract={'date':now(),'source_authorization':'/Users/zy/.codex/attachments/0b506bf4-13d6-4bc0-a869-9578417a770c/已粘贴的文本.txt','primary_membership_sha256':digest(O/'common_primary_membership.tsv'),'consensus_standardization':'Each model standardized globally over common-primary scored variant-partner units shared by both methods; population SD ddof=0; undefined if constant. Disease labels excluded.','structural_validation':{'label_source':'Frozen experimental biological-assembly contacts only; threshold <=5 Angstrom','prediction_threshold':0.5,'bootstrap_n':1000,'bootstrap_unit':'query protein','seed':20261003,'macro':'arithmetic mean of per-query metrics; one-class groups excluded from AUROC/AUPRC; explicit status','training_overlap':'Unknown for SPPIDER/XPAIR unless independently verified; validation is not held-out performance'},'partner_specificity':{'meaningful_dispersion':'Across-partner population residue SD >=0.05','constant_track':'correlation NA','partners':'all successful frozen partners'},'shuffle':{'max_genes':20,'selection':'first20 lexicographic genes with at least3 sequence-supported pairs per method; first3 lexicographic partners per gene','length_bins':[0,50,200,500,1022,2000,3000,100000],'replacement':'lexicographic frozen accession in same bin, excluding query and all its credible partners; deterministic hash offset seed20261003','no_phenotype_used':True},'coverage':{'GE3':3,'C4_threshold':0.70,'retained_groups_min':2,'retained_group_unambiguous_n_min':10,'severe_min':0.50,'severe_range':0.50,'material_range':0.20,'historical_legacy':'9/3 retained unchanged; disease group membership uses recovered legacy make_membership; gene denominator exactly common registry'},'runtime':{'SPPIDER':'CPU','XPAIR_encoder':'CPU','XPAIR_head':'MPS auto per official CLI','serial_methods':True},'stop':'After Phase2S final verdict; no Phase2P'}
cp=O/'analysis_contract.json'
if not cp.exists():cp.write_text(json.dumps(contract,indent=2)+'\n')
manifest=[]
for method,version,dirname,head in [('SPPIDER_SEQ_PARTNER_INTERFACE','v2.0-manuscript','esm2_t33_650M_UR50D','data/software/phase2S/SPPIDER-seq-v2.0-manuscript/models/crossattn_rec_run27_best.pt'),('XPAIR_PARTNER_INTERFACE','0.1.0','ankh-large','data/software/phase2S/xpair-0.1.0/xpair/default_model/multitask_xfair.ckpt')]:
    mf=json.loads((R/'external_models'/dirname/'OFFICIAL_FILE_MANIFEST.json').read_text());weight=next(x for x in mf['files'] if x['path'].endswith(('.bin','.safetensors')))
    manifest.append({'method':method,'version':version,'encoder':mf['repo_id'],'revision':mf['revision'],'encoder_SHA256':weight['sha256'],'head_SHA256':digest(R/head),'runtime':'/Users/zy/.codex/tools/phase2s_'+('sppider' if method.startswith('SPPIDER') else 'xpair'),'status':'OFFLINE_EXAMPLE_PASS','started_at':now(),'contract_SHA256':digest(cp)})
wr(pd.DataFrame(manifest),'run_manifest.tsv')
for method in ['sppider','xpair']:
    out=plan.copy();out['status']=out[method+'_preflight'];wr(out,'pair_status_'+method+'.tsv') if not (O/('pair_status_'+method+'.tsv')).exists() else None
# Dedicated method-unavailability audit, fixed without looking at scores.
x=plan.copy();x['query_length_rejected']=(x.query_length<50)|(x.query_length>2000);x['partner_length_rejected']=(x.partner_length<50)|(x.partner_length>2000);x['both_length_rejected']=x.query_length_rejected&x.partner_length_rejected
x.loc[x.xpair_preflight.isin(['QUERY_SEQUENCE_MISSING','PARTNER_SEQUENCE_MISSING']),['query_length_rejected','partner_length_rejected','both_length_rejected']]=False
wr(x,'xpair_length_constraint_audit.tsv')
print(json.dumps({'freeze':'124/124 PASS','U1':42,'U2':132,'pairs':15942,'common_primary':9978,'sppider_preflight':plan.sppider_preflight.value_counts().to_dict(),'xpair_preflight':plan.xpair_preflight.value_counts().to_dict()}),flush=True)
