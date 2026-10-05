from pathlib import Path
import hashlib,json
import pandas as pd,numpy as np
from Bio import SeqIO
R=Path(__file__).resolve().parents[2];P=R/'results/phase2SR';D=P/'pilot_inputs'
pairs=pd.read_csv(D/'pilot_pairs.tsv',sep='\t',dtype=str,keep_default_na=False)
v=pd.read_csv(R/'results/phase2S/frozen_variants.tsv',sep='\t',dtype=str,keep_default_na=False)
seq={x.id:str(x.seq) for x in SeqIO.parse(D/'xpair_sequences.fasta','fasta')}
raw=pd.read_csv(P/'xpair_pilot/predictions_six_decimal.csv');again=pd.read_csv(P/'xpair_pilot_repeat/predictions_six_decimal.csv')
assert len(raw)==len(again)
raw=raw.sort_values(['id1','id2']).reset_index(drop=True);again=again.sort_values(['id1','id2']).reset_index(drop=True)
expected=pd.read_csv(D/'xpair_pairs.tsv',sep='\t',dtype=str)
assert len(raw)==len(expected) and set(zip(raw.id1,raw.id2))==set(zip(expected.id1,expected.id2))
assert raw.equals(again),'X-PAIR cache repeat changed predictions'
pred={};cache=[];encoder=json.loads((R/'external_models/ankh-large/OFFICIAL_FILE_MANIFEST.json').read_text())
enc=next(x['sha256'] for x in encoder['files'] if x['path']=='pytorch_model.bin')
for r in raw.itertuples(index=False):
    assert np.isfinite(r.interaction_probs) and 0<=r.interaction_probs<=1
    for n in [1,2]:
        sid=getattr(r,'id'+str(n));a=np.fromstring(getattr(r,'interface'+str(n)+'_probs'),sep=' ')
        assert len(a)==len(seq[sid]) and np.isfinite(a).all() and ((0<=a)&(a<=1)).all()
    probs=np.fromstring(r.interface1_probs,sep=' ');pred[(r.id1,r.id2)]=probs
    qh=hashlib.sha256(seq[r.id1].encode()).hexdigest();ph=hashlib.sha256(seq[r.id2].encode()).hexdigest()
    key=hashlib.sha256(json.dumps(['XPAIR-0.1.0',enc,qh,ph]).encode()).hexdigest()
    cache.append({'method_version':'XPAIR-0.1.0','encoder_model_SHA256':enc,'query_sequence_SHA256':qh,'partner_sequence_SHA256':ph,'cache_path':str((P/'xpair_pilot/predictions_six_decimal.csv').relative_to(R)),'cache_key':key,'status':'VALIDATED_RAW_TRACKS'})
pre=pd.read_csv(D/'xpair_preflight.tsv',sep='\t',dtype=str).set_index(['gene','partner'])
rows=[];variants=[]
for r in pairs.itertuples(index=False):
    status=pre.loc[(r.gene,r.partner),'status'];probs=pred.get((r.query_accession,r.partner_accession))
    if probs is not None:status='SCORED'
    assert status!='PENDING',(r.gene,r.partner)
    rows.append({'gene':r.gene,'partner':r.partner,'status':status,'cache_repeat_exact':probs is not None})
    for x in v.loc[v.gene.eq(r.gene)].itertuples(index=False):
        state=status if probs is None else 'MAPPING_UNAVAILABLE';score='NA';pos=None
        if probs is not None and x.primary_cohort_eligible=='YES' and x.canonical_residue.isdigit():
            pos=int(x.canonical_residue);q=seq[r.query_accession]
            if 1<=pos<=len(q) and q[pos-1]==x.canonical_reference_aa:state='SCORED';score=float(probs[pos-1])
            else:state='REFERENCE_OR_COORDINATE_MISMATCH'
        variants.append({'gene':r.gene,'VariationID':x.VariationID,'partner':r.partner,'canonical_residue':pos,'interface_probability':score,'status':state})
pd.DataFrame(rows).to_csv(P/'pilot_predictions_xpair.tsv',sep='\t',index=False)
pd.DataFrame(variants).to_csv(P/'pilot_variant_scores_xpair.tsv',sep='\t',index=False)
pd.DataFrame(cache).to_csv(P/'xpair_cache_manifest.tsv',sep='\t',index=False)
pd.concat([pd.read_csv(P/'sppider_cache_manifest.tsv',sep='\t'),pd.DataFrame(cache)],ignore_index=True).to_csv(P/'inference_cache_manifest.tsv',sep='\t',index=False)
summary=[]
for m in ['sppider','xpair']:
    d=pd.read_csv(P/('pilot_predictions_'+m+'.tsv'),sep='\t');sc=pd.read_csv(P/('pilot_variant_scores_'+m+'.tsv'),sep='\t',dtype=str,keep_default_na=False)
    assert not sc.status.eq('REFERENCE_OR_COORDINATE_MISMATCH').any(),m+' coordinate mismatch'
    counts=sc.loc[sc.status.eq('SCORED')].groupby(['gene','VariationID']).partner.nunique().to_dict()
    for r in v.loc[v.gene.isin(pairs.gene)].itertuples(index=False):summary.append({'method':m,'gene':r.gene,'VariationID':r.VariationID,'scored_partner_n':counts.get((r.gene,r.VariationID),0),'ge3':counts.get((r.gene,r.VariationID),0)>=3})
    print(m,'technical pairs',int(d.status.eq('SCORED').sum()),'/',len(d),flush=True)
pd.DataFrame(summary).to_csv(P/'pilot_variant_partner_count_validation.tsv',sep='\t',index=False)
window=pd.read_csv(P/'sppider_window_validation.tsv',sep='\t')
assert (window.query_length<1022).any() and (window.query_length>1022).any()
assert window.cache_repeat_exact.all() and window.no_internal_gaps.all() and window.endpoint_alignment.all()
print('PILOT_VALIDATION_PASS; no rescue analysis; full benchmark not started',flush=True)
