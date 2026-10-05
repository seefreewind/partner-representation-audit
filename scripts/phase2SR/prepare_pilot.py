from pathlib import Path
import pandas as pd,hashlib,json
R=Path(__file__).resolve().parents[2];P=R/'results/phase2SR';D=P/'pilot_inputs';D.mkdir(exist_ok=True)
genes=pd.read_csv(P/'pilot_gene_list.tsv',sep='\t');pp=pd.read_csv(R/'results/phase2S/frozen_partner_sets.tsv',sep='\t',dtype=str,keep_default_na=False)
rows=[]
for g in genes.gene:
    group=pp.loc[pp.gene.eq(g)].sort_values('partner').head(3)
    for r in group.itertuples(index=False):rows.append({'gene':g,'partner':r.partner,'query_accession':r.query_accession,'partner_accession':r.partner_accession,'tier':r.tier,'selection':'FIRST_3_FROZEN_PARTNERS_LEXICOGRAPHIC_BEFORE_PREDICTION'})
table=pd.DataFrame(rows);table.to_csv(D/'pilot_pairs.tsv',sep='\t',index=False)
seq=pd.read_csv(P/'canonical_sequence_recovery.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('requested_accession')
fasta=[];xrows=[];status=[]
used=set()
for r in table.itertuples(index=False):
    a,b=r.query_accession,r.partner_accession
    state='PENDING'
    if a not in seq.index or b not in seq.index:state='UNRESOLVED_CANONICAL_ACCESSION'
    elif not 50<=int(seq.loc[a,'sequence_length'])<=2000 or not 50<=int(seq.loc[b,'sequence_length'])<=2000:state='MODEL_LENGTH_UNSUPPORTED'
    else:
        xrows.append({'id1':a,'id2':b})
        for acc in [a,b]:
            if acc not in used:
                fasta.append((R/'external_models/canonical_sequences'/f'{acc}.fasta').read_text());used.add(acc)
    status.append({'gene':r.gene,'partner':r.partner,'status':state})
pd.DataFrame(xrows).drop_duplicates().to_csv(D/'xpair_pairs.tsv',sep='\t',index=False)
(D/'xpair_sequences.fasta').write_text(''.join(fasta))
pd.DataFrame(status).to_csv(D/'xpair_preflight.tsv',sep='\t',index=False)
(D/'selection_contract.json').write_text(json.dumps({'genes':genes.gene.tolist(),'pair_n':len(rows),'partners_per_gene':3,'purpose':'technical pilot only; full frozen partner universe remains unchanged','no_disease_labels_used':True,'selection_sha256':hashlib.sha256((D/'pilot_pairs.tsv').read_bytes()).hexdigest()},indent=2)+'\n')
print('Pilot fixed',len(genes),'genes',len(rows),'pairs',len(xrows),'X-PAIR valid pairs')
