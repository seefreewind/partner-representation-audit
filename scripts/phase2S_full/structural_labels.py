"""Phenotype-blind full observed residue labels from frozen experimental assemblies."""
from common import *
import importlib.util,numpy as np
spec=importlib.util.spec_from_file_location('phase1a',R/'scripts/run_phase1a.py');s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
p=rd('results/phase2S_full/pair_plan.tsv').set_index(['gene','partner']);coverage=rd('results/phase1A/structure_coverage.tsv');C=O/'experimental_residue_labels';C.mkdir(exist_ok=True);rows=[]
for r in coverage.loc[coverage.observed_experimental_complex_IDs.ne('')].itertuples(index=False):
    if (r.gene,r.partner) not in p.index:continue
    info=p.loc[(r.gene,r.partner)];a,b=info.query_accession,info.partner_accession
    if not a or not b:continue
    q,ctx=seq(a),seq(b)
    for assembly in r.observed_experimental_complex_IDs.split(';'):
        name=assembly.replace(':assembly','-assembly');cif=R/'data/interim/phase1A_parsed_cif'/(name+'.cif');key=hashlib.sha256(json.dumps([r.gene,r.partner,assembly,info.query_sha256,info.partner_sha256]).encode()).hexdigest();file=C/(key+'.npz')
        if not cif.exists():rows.append(dict(gene=r.gene,partner=r.partner,complex=assembly,status='FROZEN_CIF_MISSING'));continue
        if not file.exists():
            contacts,status=s.parse_cif_contacts(cif,a,b,q,ctx);dist={}
            for c in contacts:
                for pos,val in c['interface_distances'].items():dist[pos]=min(dist.get(pos,np.inf),val)
            valid=sorted(pos for pos,val in dist.items() if 1<=pos<=len(q) and np.isfinite(val))
            if not valid:rows.append(dict(gene=r.gene,partner=r.partner,complex=assembly,status=status));continue
            positions=np.array(valid);distances=np.array([dist[x] for x in valid]);tmp=file.with_suffix('.tmp')
            with tmp.open('wb') as f:np.savez(f,positions=positions,distance_A=distances,label=(distances<=5).astype(int))
            tmp.replace(file)
        z=np.load(file);rows.append(dict(gene=r.gene,partner=r.partner,complex=assembly,query_accession=a,partner_accession=b,observed_residues=len(z['label']),positive_residues=int(z['label'].sum()),status='EXPERIMENTAL_LABELS',label_path=str(file.relative_to(R)),label_SHA256=digest(file),cif_SHA256=digest(cif),training_overlap='NOT_ESTABLISHED_NON_INDEPENDENT_VALIDATION'))
    if len(rows)%20==0:wr(pd.DataFrame(rows),'experimental_validation_label_manifest.tsv')
wr(pd.DataFrame(rows),'experimental_validation_label_manifest.tsv');print('EXPERIMENTAL_LABELS_COMPLETE',len(rows),flush=True)
