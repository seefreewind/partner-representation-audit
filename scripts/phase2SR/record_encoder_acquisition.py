"""Record verified official downloads without claiming offline inference PASS."""
from pathlib import Path
import json,hashlib
import pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'results/phase2SR'
summaries=[]
for method,directory,repo in [('sppider','esm2_t33_650M_UR50D','facebook/esm2_t33_650M_UR50D'),('xpair','ankh-large','ElnaggarLab/ankh-large')]:
    D=R/'external_models'/directory
    pin=json.loads((D/'OFFICIAL_REVISION.json').read_text());mf=json.loads((D/'OFFICIAL_FILE_MANIFEST.json').read_text())
    assert pin['repo_id']==mf['repo_id']==repo and pin['revision']==mf['revision']
    byname={r['path']:r for r in mf['files']};assert set(byname)==set(pin['files'])
    weight=next(x for x in byname if x.endswith(('.safetensors','.bin')))
    assert byname[weight]['sha256']==pin['official_metadata'][weight]['sha256']
    for name,r in byname.items():assert (D/name).stat().st_size==r['bytes']
    configsha=byname['config.json']['sha256']
    tok={k:v['sha256'] for k,v in byname.items() if k!=weight and k!='config.json'}
    rows=[{'repo_id':repo,'revision':mf['revision'],'file':n,'local_path':str(D/n),'byte_size':r['bytes'],'SHA256':r['sha256'],'weight_SHA256':byname[weight]['sha256'],'config_SHA256':configsha,'tokenizer_hashes':json.dumps(tok,sort_keys=True),'status':'ACQUIRED_OFFICIAL_HASH_VERIFIED','offline_load':'NOT_RUN','example':'NOT_RUN'} for n,r in byname.items()]
    pd.DataFrame(rows).to_csv(O/(method+'_encoder_manifest.tsv'),sep='\t',index=False)
    summaries.append({'repo':repo,'revision':mf['revision'],'weight_file':weight,'bytes':byname[weight]['bytes'],'sha256':byname[weight]['sha256']})
(O/'encoder_download_provenance.json').write_text(json.dumps({'user_authorized_local_download_after_network_switch':True,'date':'2026-10-02','huggingface_hub':'1.8.0','python_runtime':'/Users/zy/.hf-prott5-download-venv/bin/python','encoders':summaries},indent=2)+'\n')
report=R/'reports/PHASE2SR_OFFLINE_ENCODER_AND_DENOMINATOR_RECOVERY.md'
s=report.read_text()
s=s.replace('Both official encoders await external transfer. No Hugging Face acquisition was attempted in this recovery phase on this computer.','Both official encoders have now been downloaded locally and checked against official file sizes and LFS SHA256. On 2026-10-02 the user explicitly authorized local acquisition after switching the network, superseding the earlier external-transfer restriction.')
s+='\n## Local encoder acquisition completed — 2026-10-02\n\n'+pd.DataFrame(summaries).to_markdown(index=False)+'\n\nAll selected config and tokenizer files were also downloaded and SHA256 recorded. Exact encoder manifests are now populated. A transient proxy 503 during Ankh metadata acquisition was recovered without changing models or revisions. The download-only task is complete. **Offline model loading, official examples and the pilot remain unvalidated; this acquisition does not establish GO_RESUME_PHASE2S.** Independent environment installation/validation logs remain authoritative.\n'
report.write_text(s)
print(json.dumps(summaries,indent=2))
