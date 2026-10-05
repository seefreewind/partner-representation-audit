from common import *
import os,importlib.util,time,resource,gc,fcntl
os.environ['HF_HUB_OFFLINE']='1';os.environ['TRANSFORMERS_OFFLINE']='1'
import torch,numpy as np
lock=(O/'sppider.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
M='SPPIDER_SEQ_PARTNER_INTERFACE';run=rd('results/phase2S_full/run_manifest.tsv').set_index('method').loc[M];enc=run.encoder_SHA256
C=O/'residue_cache/sppider';E=O/'embedding_cache/sppider';C.mkdir(parents=True,exist_ok=True);E.mkdir(parents=True,exist_ok=True)
spec=importlib.util.spec_from_file_location('sppider',R/'data/software/phase2S/SPPIDER-seq-v2.0-manuscript/cli/sppider_seq.py');s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
args=s.build_parser().parse_args(['--query','unused','--partners','unused','--output-dir',str(C),'--esm-model-name',str(R/'external_models/esm2_t33_650M_UR50D'),'--receptor-model-path',str(R/'data/software/phase2S/SPPIDER-seq-v2.0-manuscript/models/crossattn_rec_run27_best.pt'),'--peptide-model-path',str(R/'data/software/phase2S/SPPIDER-seq-v2.0-manuscript/models/crossattn_pep_run07_best.pt'),'--device','cpu'])
s.setup_runtime(args);s.load_ppi_models();assert s.MAX_TOKENS==1024 and s.STRIDE==512 and s.EMBED_DIM==1280
heads=[digest(args.receptor_model_path),digest(args.peptide_model_path)]
cache_paths={str(f.relative_to(R)) for directory in [C,R/'results/phase2SR/inference_cache/sppider'] for f in directory.iterdir() if f.name.endswith('.npz') and not f.name.startswith('._')}
embedding_memory={}
def chunks(q):
    sh=hashlib.sha256(q.encode()).hexdigest();key=hashlib.sha256(json.dumps([enc,sh,1024,512]).encode()).hexdigest();file=E/(key+'.pt')
    if sh in embedding_memory:return embedding_memory[sh]
    if file.exists():d=torch.load(file,map_location='cpu',weights_only=False);items,n=d['chunks'],d['length'];assert d['sequence_sha256']==sh and d['encoder_sha256']==enc
    else:
        items,n=s.embed_sequence_chunks(q)
        tmp=file.with_suffix('.tmp');torch.save({'chunks':items,'length':n,'sequence_sha256':sh,'encoder_sha256':enc},tmp);tmp.replace(file)
    coverage=np.zeros(len(q),dtype=bool)
    for x in items:
        a,b=x['start_idx'],x['end_idx'];assert x['embedding'].shape==(b-a,1280) and torch.isfinite(x['embedding']).all();coverage[a:b]=True
    assert n==len(q) and coverage.all() and items[0]['start_idx']==0 and items[-1]['end_idx']==len(q)
    if len(embedding_memory)>=8:embedding_memory.clear();gc.collect()
    embedding_memory[sh]=items;return items
# Full disk tracks, no recomputation per variant. Reuse validated pilot caches.
plan=rd('results/phase2S_full/pair_plan.tsv');journal=O/'sppider_pair_checkpoint.jsonl';records={}
if journal.exists():
    for line in journal.read_text().splitlines():
        try:
            r=json.loads(line);records[(r['gene'],r['partner'])]=r
        except json.JSONDecodeError:pass
if journal.exists() and journal.stat().st_size and not journal.read_bytes().endswith(b'\n'):
    with journal.open('ab') as f:f.write(b'\n')
consecutive=0;started=time.monotonic()
def publish():
    rows=[];cache=[]
    for r in plan.to_dict('records'):
        d=records.get((r['gene'],r['partner']),dict(status=r['sppider_preflight']));rows.append({**r,**d})
        if d['status']=='SCORED':cache.append({k:d.get(k,'') for k in ['gene','partner','query_accession','partner_accession','method','model_hash','encoder_sha256','query_sha256','partner_sha256','residue_count','score_array_path','score_array_SHA256','cache_key','status','precision']})
    wr(pd.DataFrame(rows),'pair_status_sppider.tsv');wr(pd.DataFrame(cache),'residue_cache_manifest_sppider.tsv')
    checkpoint('sppider','FULL_INFERENCE',processed=len(records),scored=sum(x['status']=='SCORED' for x in records.values()),total=len(plan),elapsed_seconds=time.monotonic()-started,pid=os.getpid())
for i,r in enumerate(plan.itertuples(index=False)):
    old=records.get((r.gene,r.partner));
    if old and old['status']=='SCORED' and old['score_array_path'] in cache_paths:continue
    t=time.monotonic();cpu=time.process_time();d=dict(gene=r.gene,partner=r.partner,query_accession=r.query_accession,partner_accession=r.partner_accession,method=M,model_hash=run.head_SHA256,encoder_sha256=enc,query_sha256=r.query_sha256,partner_sha256=r.partner_sha256,status=r.sppider_preflight)
    if d['status']=='PENDING':
        try:
            q,p=seq(r.query_accession),seq(r.partner_accession);assert hashlib.sha256(q.encode()).hexdigest()==r.query_sha256 and hashlib.sha256(p.encode()).hexdigest()==r.partner_sha256
            key=hashlib.sha256(json.dumps(['v2.0-manuscript',enc,heads,r.query_sha256,r.partner_sha256]).encode()).hexdigest();file=C/(key+'.npz');pilot=R/'results/phase2SR/inference_cache/sppider'/(key+'.npz')
            hit=file.exists() or pilot.exists()
            if hit:
                if not file.exists():file=pilot
                z=np.load(file);rec,pep=z['receptor'],z['peptide'];assert str(z['query_sha'])==r.query_sha256 and str(z['partner_sha'])==r.partner_sha256
            else:
                qc,pc=chunks(q),chunks(p);out=[]
                for model in [s.receptor_site_model,s.peptide_site_model]:
                    with torch.no_grad():
                        pooled=model([x['embedding'].to(s.device) for x in qc],[x['embedding'].to(s.device) for x in pc]);raw=s.merge_chunk_logits_to_full(pooled,[x['start_idx'] for x in qc],len(q));out.append(torch.sigmoid(raw).cpu().numpy())
                rec,pep=out
            if len(rec)!=len(q) or len(pep)!=len(q):d['status']='OUTPUT_LENGTH_MISMATCH'
            elif not np.isfinite(rec).all() or not np.isfinite(pep).all():d['status']='NONFINITE_OUTPUT'
            elif not ((0<=rec)&(rec<=1)).all() or not ((0<=pep)&(pep<=1)).all():d['status']='OTHER';d['error']='PROBABILITY_OUT_OF_RANGE'
            else:
                if not hit:
                    tmp=file.with_suffix('.tmp')
                    with tmp.open('wb') as f:np.savez(f,receptor=rec,peptide=pep,query_sha=r.query_sha256,partner_sha=r.partner_sha256)
                    tmp.replace(file)
                d.update(status='SCORED',cache_key=key,cache_hit=hit,residue_count=len(rec),score_array_path=str(file.relative_to(R)),score_array_SHA256=digest(file),precision='float32')
        except Exception as ex:d.update(status='MODEL_FAILURE',error=repr(ex))
    d.update(wall_seconds=time.monotonic()-t,cpu_seconds=time.process_time()-cpu,peak_RAM_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,timestamp=now())
    records[(r.gene,r.partner)]=d
    with journal.open('a') as f:f.write(json.dumps(d)+'\n')
    consecutive=consecutive+1 if d['status'] in ['MODEL_FAILURE','OUTPUT_LENGTH_MISMATCH','NONFINITE_OUTPUT','OTHER'] else 0
    if (i+1)%10==0 or i==len(plan)-1:publish();print('SPPIDER_FULL',i+1,len(plan),d['status'],flush=True)
    if consecutive>=3:publish();raise RuntimeError('Three consecutive model failures; stop for diagnosis')
publish();checkpoint('sppider','COMPLETE',processed=len(records),scored=sum(x['status']=='SCORED' for x in records.values()),total=len(plan),pid=os.getpid());print('SPPIDER_FULL_COMPLETE',flush=True)
