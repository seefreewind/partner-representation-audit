from pathlib import Path
import pandas as pd,numpy as np,json,sys
from Bio import SeqIO
R=Path(__file__).resolve().parents[2];P=R/'results/phase2SR';method=sys.argv[1]
rows=[]
if method=='sppider':
    query=list(SeqIO.parse(R/'data/software/phase2S/SPPIDER-seq-v2.0-manuscript/examples/examples_1-2_input-query.fasta','fasta'))[0]
    files=list((P/'sppider_official_example').glob('*.txt'));files=[f for f in files if 'query_as_' in f.name and not f.name.startswith('._')]
    assert len(files)==6,len(files)
    for f in files:
        d=pd.read_csv(f,sep='\t',comment='#');p=d.Probability.to_numpy()
        assert len(d)==len(query.seq) and np.array_equal(d.Position.to_numpy(),np.arange(1,len(query.seq)+1))
        assert ''.join(d.AminoAcid)==str(query.seq)
        assert np.isfinite(p).all() and ((0<=p)&(p<=1)).all()
        rows.append({'test':'official_example','file':str(f.relative_to(R)),'length':len(d),'finite':True,'probability_range':True,'coordinate_alignment':True,'result':'PASS'})
else:
    D=P/'xpair_official_workflow';seq={x.id:str(x.seq) for x in SeqIO.parse(D/'sequences.fasta','fasta')}
    d=pd.read_csv(D/'predictions.csv');pairs=pd.read_csv(D/'pairs.tsv',sep='\t');assert len(d)==len(pairs)==3
    for r in d.itertuples(index=False):
        assert np.isfinite(r.interaction_probs) and 0<=r.interaction_probs<=1
        for number in [1,2]:
            probs=np.fromstring(getattr(r,'interface'+str(number)+'_probs'),sep=' ');sid=getattr(r,'id'+str(number))
            assert len(probs)==len(seq[sid]),(sid,len(probs),len(seq[sid]))
            assert np.isfinite(probs).all() and ((0<=probs)&(probs<=1)).all()
        rows.append({'test':'official_documented_workflow','id1':r.id1,'id2':r.id2,'interaction_output':True,'finite':True,'length_match':True,'result':'PASS'})
pd.DataFrame(rows).to_csv(P/(method+'_example_validation.tsv'),sep='\t',index=False)
print(method.upper()+'_EXAMPLE_PASS=YES',len(rows))
