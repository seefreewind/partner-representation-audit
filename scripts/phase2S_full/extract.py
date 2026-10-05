from common import *
import numpy as np,csv
v=rd('results/phase2S/frozen_variants.tsv');m=rd('results/phase2SR/common_primary_membership.tsv');primary=set(zip(m.loc[m.included_in_common_primary_denominator.eq('True'),'gene'],m.loc[m.included_in_common_primary_denominator.eq('True'),'VariationID']))
sequences={};output={}
for method in ['sppider','xpair']:
    pairs=rd('results/phase2S_full/pair_status_'+method+'.tsv');assert not pairs.status.eq('PENDING').any()
    path=O/('variant_partner_scores_'+method+'.tsv');tmp=path.with_suffix('.tmp');rows_scored=0
    fields=['gene','VariationID','partner','query_accession','partner_accession','canonical_residue','common_primary','interface_probability','status','pair_status']
    with tmp.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
        for r in pairs.itertuples(index=False):
            arr=None;q=''
            if r.status=='SCORED':
                file=R/r.score_array_path;assert digest(file)==r.score_array_SHA256;z=np.load(file);arr=z['receptor'] if method=='sppider' else z['query'];assert len(arr)==int(r.query_length) and np.isfinite(arr).all()
                if r.query_accession not in sequences:sequences[r.query_accession]=seq(r.query_accession)
                q=sequences[r.query_accession];assert hashlib.sha256(q.encode()).hexdigest()==r.query_sha256
            for x in v.loc[v.gene.eq(r.gene)].itertuples(index=False):
                state=r.status if arr is None else 'CANONICAL_MAPPING_UNAVAILABLE';score='NA';pos=x.canonical_residue
                if arr is not None and x.primary_cohort_eligible=='YES' and pos.isdigit():
                    n=int(pos)
                    if x.canonical_UniProt_accession!=r.query_accession or not 1<=n<=len(q) or q[n-1]!=x.canonical_reference_aa:state='VARIANT_REFERENCE_MISMATCH'
                    else:state='SCORED';score=float(arr[n-1]);rows_scored+=1
                w.writerow(dict(gene=r.gene,VariationID=x.VariationID,partner=r.partner,query_accession=r.query_accession,partner_accession=r.partner_accession,canonical_residue=pos,common_primary=(x.gene,x.VariationID) in primary,interface_probability=score,status=state,pair_status=r.status))
    tmp.replace(path);output[method]=rows_scored
# Shared availability is determined exclusively by valid canonical WT extraction.
a=pd.read_csv(O/'variant_partner_scores_sppider.tsv',sep='\t',dtype=str,keep_default_na=False);b=pd.read_csv(O/'variant_partner_scores_xpair.tsv',sep='\t',dtype=str,keep_default_na=False)
keys=['gene','VariationID','partner'];d=a[keys+['interface_probability','status']].merge(b[keys+['interface_probability','status']],on=keys,suffixes=('_sppider','_xpair'),validate='one_to_one')
d['common_primary']=[(g,vid) in primary for g,vid in d[['gene','VariationID']].itertuples(index=False,name=None)]
both=d.status_sppider.eq('SCORED')&d.status_xpair.eq('SCORED');d['both_available']=both;d['consensus_high']=False;d['continuous_consensus']=np.nan;params={}
zs=[]
for method in ['sppider','xpair']:
    arr=pd.to_numeric(d.loc[both,'interface_probability_'+method]);standard=arr.loc[d.loc[both,'common_primary']];mu=standard.mean();sd=standard.std(ddof=0);params[method]={'mean':mu,'population_SD':sd,'standardization_rows':len(standard),'membership':'common_primary shared valid pairs','disease_labels_used':False};zs.append((arr-mu)/sd if sd>0 else arr*np.nan)
d.loc[both,'continuous_consensus']=(zs[0]+zs[1])/2;d.loc[both,'consensus_high']=(pd.to_numeric(d.loc[both,'interface_probability_sppider'])>=.5)&(pd.to_numeric(d.loc[both,'interface_probability_xpair'])>=.5)
wr(d,'sequence_consensus_scores.tsv');(O/'consensus_standardization_parameters.json').write_text(json.dumps(params,indent=2)+'\n');print('VARIANT_EXTRACTION_COMPLETE',output,flush=True)
