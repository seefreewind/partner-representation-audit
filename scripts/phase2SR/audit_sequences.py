from pathlib import Path
import hashlib
import pandas as pd
R=Path(__file__).resolve().parents[2];P=R/'results/phase2SR';D=R/'external_models/canonical_sequences'
old=pd.read_csv(R/'results/phase2S/sequence_manifest.tsv',sep='\t',dtype=str,keep_default_na=False)
new=pd.read_csv(P/'canonical_sequence_recovery.tsv',sep='\t',dtype=str,keep_default_na=False)
assert len(new)==7300 and set(old.UniProt_accession)==set(new.requested_accession)
audit=[]
for r in old.loc[old.status.eq('AVAILABLE_FROZEN')].itertuples(index=False):
    current=new.loc[new.requested_accession.eq(r.UniProt_accession)].iloc[0]
    equal=current.sequence_sha256==r.sequence_sha256
    audit.append({'accession':r.UniProt_accession,'frozen_sha256':r.sequence_sha256,'retrieved_sha256':current.sequence_sha256,'exact_match':equal})
    if not equal:
        seq=''.join(x.strip() for x in (R/r.source).read_text().splitlines() if not x.startswith('>'))
        assert hashlib.sha256(seq.encode()).hexdigest()==r.sequence_sha256
        (D/f'{r.UniProt_accession}.fasta').write_text(f'>{r.UniProt_accession} source=frozen\n{seq}\n')
        mask=new.requested_accession.eq(r.UniProt_accession)
        new.loc[mask,['sequence_sha256','sequence_length','resolved_accession','status','source','replacement_reason']]=[r.sequence_sha256,len(seq),r.UniProt_accession,'EXACT',r.source,'FROZEN_SEQUENCE_RETAINED_RETRIEVAL_DIFFERENCE']
pd.DataFrame(audit).to_csv(P/'frozen_sequence_identity_audit.tsv',sep='\t',index=False)
new.to_csv(P/'canonical_sequence_recovery.tsv',sep='\t',index=False)
print('Frozen sequence identity audit',len(audit),'match',sum(x['exact_match'] for x in audit),'differences retained frozen',sum(not x['exact_match'] for x in audit))
