from common import *
import os,time,fcntl,subprocess,resource,shutil
import numpy as np
lock=(O/'xpair.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
M='XPAIR_PARTNER_INTERFACE';run=rd('results/phase2S_full/run_manifest.tsv').set_index('method').loc[M];enc=run.encoder_SHA256
C=O/'residue_cache/xpair';E=O/'embedding_cache/xpair';B=O/'xpair_batches'
for p in [C,E,B]:p.mkdir(parents=True,exist_ok=True)
X='/Users/zy/.codex/tools/phase2s_xpair/bin/python';env=dict(os.environ,HF_HUB_OFFLINE='1',TRANSFORMERS_OFFLINE='1',TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD='1',XPAIR_ANKH_LOCAL_DIR=str(R/'external_models/ankh-large'),OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4')
plan=rd('results/phase2S_full/pair_plan.tsv');journal=O/'xpair_pair_checkpoint.jsonl';records={}
if journal.exists():
    for line in journal.read_text().splitlines():
        try:
            r=json.loads(line);records[(r['gene'],r['partner'])]=r
        except json.JSONDecodeError:pass
# Reuse exact canonical embeddings from the passed pilot; validate each tensor.
import torch
for acc in sorted(set(plan.query_accession)|set(plan.partner_accession)):
    if not acc:continue
    src=R/'results/phase2SR/xpair_pilot/embeddings'/(acc+'.pt');dst=E/(acc+'.pt')
    if src.exists() and not dst.exists():
        a=torch.load(src,map_location='cpu',weights_only=False);assert a.shape==(len(seq(acc)),1536) and torch.isfinite(a).all();shutil.copy2(src,dst)
if journal.exists() and journal.stat().st_size and not journal.read_bytes().endswith(b'\n'):
    with journal.open('ab') as f:f.write(b'\n')
started=time.monotonic()
def base(r):return dict(gene=r.gene,partner=r.partner,query_accession=r.query_accession,partner_accession=r.partner_accession,method=M,model_hash=run.head_SHA256,encoder_sha256=enc,query_sha256=r.query_sha256,partner_sha256=r.partner_sha256,status=r.xpair_preflight)
def save(d):
    records[(d['gene'],d['partner'])]=d
    with journal.open('a') as f:f.write(json.dumps(d)+'\n')
def publish():
    rows=[];cache=[]
    for r in plan.itertuples(index=False):
        d=records.get((r.gene,r.partner),base(r));rows.append({**r._asdict(),**d})
        if d['status']=='SCORED':cache.append({k:d.get(k,'') for k in ['gene','partner','query_accession','partner_accession','method','model_hash','encoder_sha256','query_sha256','partner_sha256','residue_count','score_array_path','score_array_SHA256','cache_key','status','precision']})
    wr(pd.DataFrame(rows),'pair_status_xpair.tsv');wr(pd.DataFrame(cache),'residue_cache_manifest_xpair.tsv');checkpoint('xpair','FULL_INFERENCE',processed=len(records),scored=sum(x['status']=='SCORED' for x in records.values()),total=len(plan),elapsed_seconds=time.monotonic()-started,pid=os.getpid())
for r in plan.itertuples(index=False):
    if r.xpair_preflight!='PENDING' and (r.gene,r.partner) not in records:save(base(r))
# Recovery pilot already contains raw query/partner tracks; include them in an input lookup.
rawpilot=pd.read_csv(R/'results/phase2SR/xpair_pilot/predictions_six_decimal.csv');pilot={(r.id1,r.id2):r for r in rawpilot.itertuples(index=False)}
def materialize(r,pred,elapsed,batch,cache_hit=False):
    d=base(r);a=np.fromstring(pred.interface1_probs,sep=' ');b=np.fromstring(pred.interface2_probs,sep=' ')
    if len(a)!=int(r.query_length) or len(b)!=int(r.partner_length):d.update(status='OUTPUT_LENGTH_MISMATCH')
    elif not np.isfinite(a).all() or not np.isfinite(b).all():d.update(status='NONFINITE_OUTPUT')
    elif not ((0<=a)&(a<=1)).all() or not ((0<=b)&(b<=1)).all():d.update(status='OTHER',error='PROBABILITY_OUT_OF_RANGE')
    else:
        key=hashlib.sha256(json.dumps(['XPAIR-0.1.0',enc,r.query_sha256,r.partner_sha256]).encode()).hexdigest();file=C/(key+'.npz');tmp=file.with_suffix('.tmp')
        with tmp.open('wb') as f:np.savez(f,query=a,partner=b,query_sha=r.query_sha256,partner_sha=r.partner_sha256)
        tmp.replace(file);d.update(status='SCORED',cache_key=key,cache_hit=cache_hit,residue_count=len(a),partner_residue_count=len(b),score_array_path=str(file.relative_to(R)),score_array_SHA256=digest(file),precision='official_six_decimal_export')
    d.update(batch=batch,allocated_batch_wall_seconds=elapsed,timestamp=now());save(d)
for r in plan.itertuples(index=False):
    old=records.get((r.gene,r.partner))
    if r.xpair_preflight=='PENDING' and not old and (r.query_accession,r.partner_accession) in pilot:materialize(r,pilot[(r.query_accession,r.partner_accession)],0,'RECOVERY_PILOT',True)
publish()
pending=[r for r in plan.itertuples(index=False) if r.xpair_preflight=='PENDING' and not (records.get((r.gene,r.partner),{}).get('status')=='SCORED' and (R/records[(r.gene,r.partner)]['score_array_path']).exists())]
# Batch grouping reuses WT embeddings across every future partner; no variant inference.
for start in range(0,len(pending),128):
    subset=pending[start:start+128];signature=hashlib.sha256(json.dumps([(r.gene,r.partner) for r in subset]).encode()).hexdigest()[:16];D=B/signature;D.mkdir(exist_ok=True)
    pd.DataFrame([{'id1':r.query_accession,'id2':r.partner_accession} for r in subset]).drop_duplicates().to_csv(D/'pairs.tsv',sep='\t',index=False)
    accessions=sorted({r.query_accession for r in subset}|{r.partner_accession for r in subset});(D/'sequences.fasta').write_text(''.join(f'>{a}\n{seq(a)}\n' for a in accessions))
    raw=D/'predictions_six_decimal.csv';t=time.monotonic()
    if not raw.exists():
        cmd=[X,'scripts/phase2S_full/xpair_instrumented.py','--predict_pairs_path',str(D/'pairs.tsv'),'--predict_fasta_path',str(D/'sequences.fasta'),'--predict_output_path',str(D/'predictions.csv'),'--emb_dir_prediction',str(E),'--num_workers','0','--batch_size_prediction','1','--interaction_threshold','0']
        with (L/('xpair_batch_'+signature+'.log')).open('a') as f:subprocess.run(cmd,cwd=R,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
    rawdf=pd.read_csv(raw);lookup={(r.id1,r.id2):r for r in rawdf.itertuples(index=False)}
    assert set(lookup)=={(r.query_accession,r.partner_accession) for r in subset},'Missing pair outputs; preserve batch and diagnose'
    elapsed=time.monotonic()-t
    if (D/'resource_measurements.json').exists():
        measure=json.loads((D/'resource_measurements.json').read_text());measure.update(batch=signature,pair_n=len(subset));
        with (O/'xpair_resource_batches.jsonl').open('a') as f:f.write(json.dumps(measure)+'\n')
    for r in subset:materialize(r,lookup[(r.query_accession,r.partner_accession)],elapsed/len(subset),signature)
    publish();print('XPAIR_FULL',len(records),len(plan),'batch',signature,flush=True)
publish();checkpoint('xpair','COMPLETE',processed=len(records),scored=sum(x['status']=='SCORED' for x in records.values()),total=len(plan),pid=os.getpid());print('XPAIR_FULL_COMPLETE',flush=True)
