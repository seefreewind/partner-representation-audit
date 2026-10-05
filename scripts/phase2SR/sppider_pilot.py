from pathlib import Path
import os,importlib.util,hashlib,json,gc
os.environ['HF_HUB_OFFLINE']='1';os.environ['TRANSFORMERS_OFFLINE']='1'
import numpy as np,pandas as pd,torch
R=Path(__file__).resolve().parents[2];P=R/'results/phase2SR';C=P/'inference_cache/sppider';C.mkdir(parents=True,exist_ok=True)
spec=importlib.util.spec_from_file_location('sppider',R/'data/software/phase2S/SPPIDER-seq-v2.0-manuscript/cli/sppider_seq.py');s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
args=s.build_parser().parse_args(['--query','unused','--partners','unused','--output-dir',str(C),'--esm-model-name',str(R/'external_models/esm2_t33_650M_UR50D'),'--receptor-model-path',str(R/'data/software/phase2S/SPPIDER-seq-v2.0-manuscript/models/crossattn_rec_run27_best.pt'),'--peptide-model-path',str(R/'data/software/phase2S/SPPIDER-seq-v2.0-manuscript/models/crossattn_pep_run07_best.pt'),'--device','cpu'])
s.setup_runtime(args);s.load_ppi_models()
assert s.MAX_TOKENS==1024 and s.STRIDE==512 and s.EMBED_DIM==1280
manifest=json.loads((R/'external_models/esm2_t33_650M_UR50D/OFFICIAL_FILE_MANIFEST.json').read_text())
encoder=next(x['sha256'] for x in manifest['files'] if x['path']=='model.safetensors')
def sha(data):return hashlib.sha256(data).hexdigest()
heads=[sha(Path(args.receptor_model_path).read_bytes()),sha(Path(args.peptide_model_path).read_bytes())]
def sequence(a):return ''.join(x.strip() for x in (R/'external_models/canonical_sequences'/f'{a}.fasta').read_text().splitlines() if not x.startswith('>'))
def chunks(seq):
    items,n=s.get_or_compute_embeddings(seq)
    assert n==len(seq) and items
    covered=np.zeros(n,dtype=bool)
    for d in items:
        a,b=d['start_idx'],d['end_idx'];assert d['embedding'].shape==(b-a,1280);covered[a:b]=True
    assert covered.all() and items[0]['start_idx']==0 and items[-1]['end_idx']==n
    return items
def predict(q,p):
    qc,pc=chunks(q),chunks(p)
    outputs=[]
    for model in [s.receptor_site_model,s.peptide_site_model]:
        with torch.no_grad():
            pooled=model([d['embedding'].to(s.device) for d in qc],[d['embedding'].to(s.device) for d in pc])
            raw=s.merge_chunk_logits_to_full(pooled,[d['start_idx'] for d in qc],len(q))
            assert raw.shape==(len(q),) and torch.isfinite(raw).all()
            probs=torch.sigmoid(raw).cpu().numpy()
        assert len(probs)==len(q) and np.isfinite(probs).all() and ((0<=probs)&(probs<=1)).all()
        outputs.append(probs)
    return outputs
pairs=pd.read_csv(P/'pilot_inputs/pilot_pairs.tsv',sep='\t',dtype=str,keep_default_na=False)
v=pd.read_csv(R/'results/phase2S/frozen_variants.tsv',sep='\t',dtype=str,keep_default_na=False)
rows=[];cache=[];window=[];variants=[]
if (P/'sppider_window_validation.tsv').exists():
    window=pd.read_csv(P/'sppider_window_validation.tsv',sep='\t').to_dict('records')
for i,r in enumerate(pairs.itertuples(index=False)):
    if not r.query_accession or not r.partner_accession:
        rows.append({'gene':r.gene,'partner':r.partner,'status':'UNRESOLVED_CANONICAL_ACCESSION'});continue
    q,p=sequence(r.query_accession),sequence(r.partner_accession)
    qh,ph=sha(q.encode()),sha(p.encode());key=sha(json.dumps(['v2.0-manuscript',encoder,heads,qh,ph]).encode());f=C/(key+'.npz')
    hit=f.exists()
    if hit:
        z=np.load(f);rec,pep=z['receptor'],z['peptide'];assert str(z['query_sha'])==qh and str(z['partner_sha'])==ph
    else:
        rec,pep=predict(q,p)
        if not window or (len(q)>1022 and not any(x['query_length']>1022 for x in window)):
            again=predict(q,p);assert np.array_equal(rec,again[0]) and np.array_equal(pep,again[1])
            window.append({'gene':r.gene,'partner':r.partner,'query_length':len(q),'context_length':len(p),'MAX_TOKENS':s.MAX_TOKENS,'STRIDE':s.STRIDE,'no_internal_gaps':True,'endpoint_alignment':True,'cache_repeat_exact':True,'status':'PASS'})
        with (C/(key+'.tmp')).open('wb') as stream:np.savez(stream,receptor=rec,peptide=pep,query_sha=qh,partner_sha=ph)
        (C/(key+'.tmp')).replace(f)
    assert len(rec)==len(q) and len(pep)==len(q)
    assert np.isfinite(rec).all() and np.isfinite(pep).all() and ((0<=rec)&(rec<=1)).all() and ((0<=pep)&(pep<=1)).all()
    rows.append({'gene':r.gene,'partner':r.partner,'query_length':len(q),'partner_length':len(p),'status':'SCORED','cache_hit':hit,'cache_key':key})
    cache.append({'method_version':'SPPIDER-v2.0-manuscript','encoder_model_SHA256':encoder,'query_sequence_SHA256':qh,'partner_sequence_SHA256':ph,'cache_path':str(f.relative_to(R)),'status':'VALIDATED'})
    for x in v.loc[v.gene.eq(r.gene)].itertuples(index=False):
        status='MAPPING_UNAVAILABLE';pos=None;score='NA'
        if x.primary_cohort_eligible=='YES' and x.canonical_residue.isdigit():
            pos=int(x.canonical_residue)
            if 1<=pos<=len(q) and q[pos-1]==x.canonical_reference_aa:status='SCORED';score=float(rec[pos-1])
            else:status='REFERENCE_OR_COORDINATE_MISMATCH'
        variants.append({'gene':r.gene,'VariationID':x.VariationID,'partner':r.partner,'canonical_residue':pos,'interface_probability':score,'status':status})
    pd.DataFrame(rows).to_csv(P/'pilot_predictions_sppider.tsv',sep='\t',index=False)
    pd.DataFrame(cache).to_csv(P/'sppider_cache_manifest.tsv',sep='\t',index=False)
    pd.DataFrame(window).to_csv(P/'sppider_window_validation.tsv',sep='\t',index=False)
    pd.DataFrame(variants).to_csv(P/'pilot_variant_scores_sppider.tsv',sep='\t',index=False)
    print('PILOT_SPPIDER',i+1,'/',len(pairs),r.gene,r.partner,len(q),len(p),flush=True)
    # RAM eviction changes only caching, not embeddings or predictions.
    if i%3==2:s.embedding_cache.clear();gc.collect()
print('PILOT_SPPIDER_COMPLETE',flush=True)
