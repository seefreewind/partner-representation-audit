"""Official encoder acquisition (local recovery authorized after network switch).
Resolve an official commit once, pin it, download official files and hash them.
"""
from pathlib import Path
import argparse, hashlib, json, os, time
os.environ.setdefault('HF_HOME',str(Path(__file__).resolve().parents[2]/'external_models/.hf_home'))
os.environ.setdefault('HF_XET_CACHE',str(Path(__file__).resolve().parents[2]/'external_models/.hf_home/xet'))
from huggingface_hub import HfApi, snapshot_download

p=argparse.ArgumentParser()
p.add_argument('--output',type=Path,required=True)
a=p.parse_args(); a.output.mkdir(parents=True,exist_ok=True)
for repo,directory in [('facebook/esm2_t33_650M_UR50D','esm2_t33_650M_UR50D'),('ElnaggarLab/ankh-large','ankh-large')]:
    target=a.output/directory; target.mkdir(exist_ok=True)
    pin=target/'OFFICIAL_REVISION.json'
    if pin.exists():
        record=json.loads(pin.read_text()); assert record['repo_id']==repo
    else:
        for attempt in range(4):
            try:
                info=HfApi().model_info(repo,files_metadata=True)
                break
            except Exception:
                if attempt==3:raise
                time.sleep(10*(attempt+1))
        names=[x.rfilename for x in info.siblings]
        # Choose one official weight format, preferring safetensors.
        weights=[x for x in names if x.endswith('.safetensors')]
        if not weights:weights=[x for x in names if x.endswith('.bin') and 'pytorch_model' in x]
        assert weights,'No official weights found'
        support=[x for x in names if x.endswith(('.json','.txt','.model'))]
        selected=sorted(set(weights+support))
        metadata={}
        for item in info.siblings:
            if item.rfilename not in selected:continue
            lfs=item.lfs
            metadata[item.rfilename]={'bytes':item.size,'sha256':getattr(lfs,'sha256',None) if lfs else None}
        record={'repo_id':repo,'revision':info.sha,'files':selected,'official_metadata':metadata}
        pin.write_text(json.dumps(record,indent=2)+'\n')
    snapshot_download(repo_id=repo,revision=record['revision'],local_dir=str(target),allow_patterns=record['files'],max_workers=2)
    rows=[]
    for name in record['files']:
        f=target/name; h=hashlib.sha256()
        with f.open('rb') as stream:
            for chunk in iter(lambda:stream.read(8*1024*1024),b''):h.update(chunk)
        expected=record.get('official_metadata',{}).get(name,{})
        if expected.get('bytes') is not None:assert f.stat().st_size==expected['bytes'],name
        if expected.get('sha256'):assert h.hexdigest()==expected['sha256'],name
        rows.append({'path':name,'bytes':f.stat().st_size,'sha256':h.hexdigest()})
    config=json.loads((target/'config.json').read_text())
    if directory=='esm2_t33_650M_UR50D':
        assert config['num_hidden_layers']==33 and config['hidden_size']==1280 and config['max_position_embeddings']==1026,config
    else:assert config['d_model']==1536,config
    (target/'OFFICIAL_FILE_MANIFEST.json').write_text(json.dumps({'repo_id':repo,'revision':record['revision'],'files':rows},indent=2)+'\n')
    print(repo,record['revision'],'COMPLETE',flush=True)
