from pathlib import Path
import hashlib,json,datetime
import pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'results/phase2S_full';L=R/'logs/phase2S_full'
def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def rd(p):return pd.read_csv(R/p,sep='\t',dtype=str,keep_default_na=False)
def wr(d,n):
    p=O/n;t=p.with_suffix(p.suffix+'.tmp');d.to_csv(t,sep='\t',index=False,na_rep='NA');t.replace(p)
def now():return datetime.datetime.now().astimezone().isoformat()
def seq(acc):return ''.join(x.strip() for x in (R/'external_models/canonical_sequences'/f'{acc}.fasta').read_text().splitlines() if not x.startswith('>'))
def checkpoint(method,stage,**kw):
    p=O/(method+'_state.json');t=p.with_suffix('.tmp');t.write_text(json.dumps(dict(method=method,stage=stage,updated_at=now(),**kw),indent=2)+'\n');t.replace(p)
