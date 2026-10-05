from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,datetime,importlib.metadata as im
R=Path(__file__).resolve().parents[2];O=R/'results/phase2P1';P=R/'results/phase2P0'
EXPECTED='d657fab6b24f9aa757f9773bb1e9e0bd5d850e743eab6bdf819f126326ced9c3'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
assert sha(P/'analysis_specification.md')==EXPECTED,'PREREGISTRATION_HASH_VIOLATION'
audit=[]
for family,file in [('ORIGINAL124','metadata/PHASE2SR_INPUT_FREEZE_MANIFEST.tsv'),('PHASE2S','results/phase2S_full/delivery_manifest.tsv'),('PHASE2P0','results/phase2P0/preregistration_input_manifest.tsv')]:
 d=pd.read_csv(R/file,sep='\t');assert len(d)==124 if family=='ORIGINAL124' else True
 for row in d.itertuples():
  observed=sha(R/row.path);audit.append(dict(family=family,path=row.path,expected_SHA256=row.SHA256,observed_SHA256=observed,status='PASS' if observed==row.SHA256 else 'FAIL'))
pd.DataFrame(audit).to_csv(O/'input_hash_verification.tsv',sep='\t',index=False);assert all(d['status']=='PASS' for d in audit),'FROZEN_INPUT_VIOLATION'
counts={'COL3A1':[97,11],'COL4A1':[17,16],'COL4A3':[112,10],'GUCY2D':[13,11],'KCNQ2':[60,17],'LHCGR':[10,10],'MEN1':[29,17],'MMUT':[43,11],'RS1':[17,14],'TUBA1A':[21,11]};ready=pd.read_csv(P/'phenotype_ready_genes.tsv',sep='\t');assert set(ready.loc[ready.universe.eq('PRIMARY_SHARED15')&ready.PHENOTYPE_READY,'gene'])==set(counts)
for gene,ns in counts.items():
 mem=pd.read_csv(P/'matrices'/(gene+'_residue_membership.tsv'),sep='\t');assert sorted(mem.MONDO_labels.value_counts().tolist())==sorted(ns),'PHENOTYPE_MEMBERSHIP_MISMATCH';assert not mem.residue_id.duplicated().any()
 for method,key in [('xpair','primary'),('sppider','confirmatory')]:
  a=np.load(P/'matrices'/(gene+'_'+method+'_'+key+'.npz'));b=np.load(P/'matrices'/(gene+'_'+method+'_unlabeled_control.npz'));assert np.array_equal(a['residue_ids'],mem.residue_id.astype(str));assert np.array_equal(np.sort(a['scores'],axis=1)[:,::-1],b['scores']);assert np.array_equal(a['residue_ids'],b['residue_ids']);assert a['scores'].shape[1]>=3
assert sum(sum(x) for x in counts.values())==547
manifest=[dict(field='preregistration_SHA256',value=EXPECTED),dict(field='source_authorization',value='c7847cb1-b57e-494b-9a7d-a964f791f547'),dict(field='preflight_timestamp',value=datetime.datetime.now().astimezone().isoformat()),dict(field='primary_genes',value=';'.join(counts)),dict(field='primary_residues',value=547),dict(field='independent_unit',value='GENE_CANONICAL_WT_RESIDUE'),dict(field='original_verdict',value='NO_GO_EDGETIC_REPRESENTATION')]+[dict(field='version_'+x,value=im.version(x)) for x in ['numpy','pandas','scipy','scikit-learn']]
pd.DataFrame(manifest).to_csv(O/'execution_manifest.tsv',sep='\t',index=False);(O/'preflight_PASS.json').write_text(json.dumps(dict(specification_SHA256=EXPECTED,input_files=len(audit),all_hashes_PASS=True,primary_genes=10,primary_residues=547),indent=2)+'\n');print('PREFLIGHT_PASS',len(audit))
