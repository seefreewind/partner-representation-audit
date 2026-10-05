from common import *
import numpy as np
p=rd('results/phase2S_full/pair_plan.tsv');lengths={r.requested_accession:int(r.sequence_length) for r in rd('results/phase2SR/canonical_sequence_recovery.tsv').itertuples(index=False)};bins=[0,50,200,500,1022,2000,3000,100000];out=[]
for method in ['sppider','xpair']:
    available=p.loc[p[method+'_preflight'].eq('PENDING')];genes=sorted(g for g,x in available.groupby('gene') if x.partner.nunique()>=3)[:20]
    for g in genes:
        subset=available.loc[available.gene.eq(g)].sort_values('partner').head(3);forbidden=set(p.loc[p.gene.eq(g),'partner_accession'])|set(p.loc[p.gene.eq(g),'query_accession'])
        for r in subset.itertuples(index=False):
            n=lengths[r.partner_accession];binid=int(np.searchsorted(bins,n,side='right')-1);candidates=sorted(a for a,v in lengths.items() if a not in forbidden and bins[binid]<=v<bins[binid+1] and (method=='sppider' or 50<=v<=2000))
            if not candidates:out.append(dict(method=method,gene=g,partner=r.partner,status='NO_SAME_BIN_UNRELATED_CONTEXT'));continue
            seed=int(hashlib.sha256(f'20261003|{method}|{g}|{r.partner}'.encode()).hexdigest(),16);a=candidates[seed%len(candidates)]
            out.append(dict(method=method,gene=g,partner=r.partner,query_accession=r.query_accession,true_partner_accession=r.partner_accession,shuffled_partner_accession=a,true_partner_length=n,shuffled_length=lengths[a],length_bin_low=bins[binid],length_bin_high=bins[binid+1],status='PRESPECIFIED_NO_PHENOTYPE',selection='FIRST20_INPUT_SUPPORTED_GENES_FIRST3_FROZEN_PARTNERS'))
file=O/'shuffled_control_plan.tsv'
if file.exists():raise RuntimeError('Control selection already frozen; do not overwrite')
wr(pd.DataFrame(out),'shuffled_control_plan.tsv');print('SHUFFLE_PLAN_FROZEN',len(out),digest(file),flush=True)
