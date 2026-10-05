"""Preserve Phase2S, independently reproduce legacy gates, recover sequences."""
from pathlib import Path
import hashlib,json,datetime,time,csv,io
import pandas as pd
import requests
R=Path(__file__).resolve().parents[2]; O=R/'results/phase2SR'; O.mkdir(exist_ok=True)
def sha(f):
    h=hashlib.sha256()
    with f.open('rb') as s:
        for x in iter(lambda:s.read(8*1024*1024),b''):h.update(x)
    return h.hexdigest()
def rd(n):return pd.read_csv(R/n,sep='\t',dtype=str,keep_default_na=False)
def wr(d,n):d.to_csv(O/n,sep='\t',index=False)
def freeze():
    target=R/'metadata/PHASE2SR_INPUT_FREEZE_MANIFEST.tsv'
    if target.exists():return
    files=[]
    for directory in ['results/phase2S','data/software/phase2S/SPPIDER-seq-v2.0-manuscript','data/software/phase2S/xpair-0.1.0']:
        files.extend(p for p in (R/directory).rglob('*') if p.is_file() and '__pycache__' not in str(p))
    files.extend([R/'reports/PHASE2S_REPRESENTATION_COVERAGE_BENCHMARK.md',R/'reports/PHASE1B0R3_OFFLINE_MODEL_INGESTION_AND_SCORING.md',R/'scripts/phase1b0r3_finalize.py'])
    for name in ['mutpred_ppi_v1_scores.tsv','coverage_by_gene.tsv','differential_score_coverage.tsv']:files.append(R/'results/phase1B0R3'/name)
    rows=[]
    for f in sorted(set(files)):
        n='NA'
        if f.suffix=='.tsv':
            with f.open('rb') as s:n=max(0,sum(b.count(b'\n') for b in iter(lambda:s.read(8*1024*1024),b''))-1)
        rows.append({'path':str(f.relative_to(R)),'SHA256':sha(f),'row_count':n,'modification_timestamp':datetime.datetime.fromtimestamp(f.stat().st_mtime,datetime.timezone.utc).isoformat()})
    d=pd.DataFrame(rows);wr(d,'input_freeze_manifest.tsv');d.to_csv(target,sep='\t',index=False)
    print('FREEZE_PASS',len(rows),flush=True)
def denominator():
    import importlib.util
    spec=importlib.util.spec_from_file_location('legacy',R/'scripts/phase1b0r3_finalize.py');legacy=importlib.util.module_from_spec(spec);spec.loader.exec_module(legacy)
    u=rd('results/phase2S/frozen_gene_universes.tsv');u1=set(u.loc[u.U1_STRUCTURAL_42.eq('True'),'gene']);u2=set(u.gene)
    strict,vg,_=legacy.make_membership(u1)
    rows=rd('results/phase1B0R3/mutpred_ppi_v1_scores.tsv'); rows.VariationID=rows.VariationID.astype(str)
    memberships=vg.groupby(['gene','VariationID']).MONDO.agg(set).to_dict()
    scorekeys=set(zip(rows.gene,rows.VariationID))
    multi={k for k in scorekeys if len(memberships.get(k,set()))>1}
    scored=rows.loc[rows.score_status.eq('SCORED')&rows.tier.isin(['Tier 1','Tier 2'])]
    counts=scored.groupby(['gene','VariationID']).partner.nunique().to_dict()
    registry=[];out=[];groupaudit=[]
    for gene in sorted(u1):
        ids=set(vg.loc[vg.gene.eq(gene),'VariationID']);primary={v for v in ids if (gene,v) not in multi}
        for v in sorted(ids):
            for mondo in sorted(memberships[(gene,v)]):registry.append({'gene':gene,'VariationID':v,'MONDO_group':mondo,'single_label_status':'MULTILABEL' if len(memberships[(gene,v)])>1 else 'SINGLE_LABEL','included_in_legacy_denominator':v in primary,'reason':'EXCLUDED_SCORE_ROW_MULTILABEL' if (gene,v) in multi else 'INCLUDED_LEGACY_RULE','score_row_present':(gene,v) in scorekeys})
        rates=[];retained=0
        for d in strict.loc[strict.gene.eq(gene)].itertuples(index=False):
            members=set(vg.loc[vg.gene.eq(gene)&vg.MONDO.eq(d.MONDO),'VariationID'])
            n=sum(counts.get((gene,v),0)>=3 for v in members);den=int(d.n_unique_variants)
            nm=sum(len(memberships.get((gene,v),set()))>1 for v in members)
            rates.append(100*n/max(1,den));retained+=den-nm>=10
            groupaudit.append({'gene':gene,'MONDO':d.MONDO,'denominator':den,'ge3':n,'coverage_pct':round(rates[-1],2),'multilabel_n':nm})
        cls='SEVERE_IMBALANCE' if not rates or min(rates)<50 or max(rates)-min(rates)>=50 else 'MATERIAL_IMBALANCE' if max(rates)-min(rates)>=20 or min(rates)<70<=max(rates) else 'MINOR_IMBALANCE'
        n=sum(counts.get((gene,v),0)>=3 for v in primary);rate=100*n/max(1,len(primary))
        passgene=rate>=70 and cls=='MINOR_IMBALANCE' and retained>=2 and scored.loc[scored.gene.eq(gene),'partner'].nunique()>=3
        out.append({'gene':gene,'primary_denominator':len(primary),'ge3':n,'ge3_pct':round(rate,2),'coverage_class':cls,'retained_groups':retained,'balanced':passgene})
    old=rd('results/phase1B0R3/coverage_by_gene.tsv').set_index('gene')
    for x in out:
        o=old.loc[x['gene']]
        assert x['primary_denominator']==int(o.primary_unambiguous_variants_n),x
        assert x['ge3']==int(o.primary_variants_with_ge3_Tier1_2_scores_n),x
        assert x['ge3_pct']==float(o.primary_variants_with_ge3_Tier1_2_scores_pct),x
        assert x['coverage_class']==o.differential_score_coverage_class,x
        assert x['balanced']==(o.eligible_for_rebuilt_primary_cohort=='YES'),x
    names=sorted(x['gene'] for x in out if x['balanced']);assert names==['FGFR2','MEN1','TUBA1A'];assert sum(x['ge3_pct']>=70 for x in out)==9
    wr(pd.DataFrame(registry),'legacy_primary_membership_exact.tsv');wr(pd.DataFrame(out),'legacy_denominator_reproduction.tsv');wr(pd.DataFrame(groupaudit),'legacy_disease_denominator_reproduction.tsv')
    # Freeze the same algorithm independent of future method availability.
    _,fullvg,_=legacy.make_membership(u2); fullmem=fullvg.groupby(['gene','VariationID']).MONDO.agg(set).to_dict()
    common=[]
    for (g,v),labels in sorted(fullmem.items()):
        included=(g,v) not in multi if g in u1 else len(labels)==1
        common.append({'gene':g,'VariationID':v,'MONDO_groups':';'.join(sorted(labels)),'included_in_common_primary_denominator':included,'reason':'EXACT_LEGACY_U1' if g in u1 else 'FROZEN_MAPPING_ELIGIBLE_SINGLE_LABEL_U2_EXTENSION','mapping_eligible':'YES'})
    wr(pd.DataFrame(common),'common_primary_membership.tsv')
    pilot=sorted(u1)[:10]+sorted(u2-u1)[:10];wr(pd.DataFrame({'gene':pilot,'universe':['U1']*10+['U2_OUTSIDE_U1']*10,'selection':'LEXICOGRAPHIC_LABEL_BLIND'}),'pilot_gene_list.tsv')
    print('LEGACY_DENOMINATOR_EXACT_PASS',9,3,names,'COMMON',sum(x['included_in_common_primary_denominator'] for x in common),flush=True)
def sequences():
    src=rd('results/phase2S/sequence_manifest.tsv');accessions=list(src.UniProt_accession);D=R/'external_models/canonical_sequences';D.mkdir(exist_ok=True)
    recovered={};session=requests.Session()
    # Batch search returns primary AND secondary accessions, permitting only
    # explicit UniProt secondary->primary replacement metadata.
    for start in range(0,len(accessions),100):
        subset=accessions[start:start+100];f=D/f'batch_{start:05d}.json'
        if f.exists():payload=json.loads(f.read_text())
        else:
            query=' OR '.join('accession:'+a for a in subset)
            try:
                resp=session.get('https://rest.uniprot.org/uniprotkb/search',params={'query':query,'format':'json','size':500},timeout=(15,90));resp.raise_for_status();payload=resp.json()
                assert 'results' in payload
                f.write_text(json.dumps(payload))
            except Exception as e:
                print('SEQUENCE_BATCH_ACCESS_FAILURE',start,str(e),flush=True);break
        for entry in payload['results']:
            primary=entry['primaryAccession'];aliases={primary,*entry.get('secondaryAccessions',[])}
            seq=entry.get('sequence',{}).get('value','')
            if not seq:continue
            for a in set(subset)&aliases:
                recovered.setdefault(a,[]).append((primary,seq))
        print('SEQUENCE_BATCH',start,len(recovered),flush=True)
    rows=[]
    for a in accessions:
        entries=set(recovered.get(a,[]));status='NOT_FOUND';resolved='';seq='';reason='NO_AUTHORITATIVE_ENTRY_RETURNED_OR_ACCESS_PENDING'
        if len(entries)==1:
            resolved,seq=next(iter(entries));status='EXACT' if resolved==a else 'OBSOLETE_REPLACED';reason='' if status=='EXACT' else 'EXPLICIT_UNIPROT_SECONDARY_ACCESSION'
        elif len(entries)>1:status='AMBIGUOUS';reason='MULTIPLE_AUTHORITATIVE_REPLACEMENTS'
        if not entries:
            r=src.loc[src.UniProt_accession.eq(a)].iloc[0]
            if r.status=='AVAILABLE_FROZEN':
                seq=''.join(x.strip() for x in (R/r.source).read_text().splitlines() if not x.startswith('>'));assert hashlib.sha256(seq.encode()).hexdigest()==r.sequence_sha256
                status='EXACT';resolved=a;reason='FROZEN_UNIPROT_CANONICAL_SEQUENCE'
        rows.append({'requested_accession':a,'status':status,'resolved_accession':resolved,'sequence_length':len(seq) if seq else 'NA','sequence_sha256':hashlib.sha256(seq.encode()).hexdigest() if seq else 'NA','source':'FROZEN_UNIPROT_CACHE' if reason=='FROZEN_UNIPROT_CANONICAL_SEQUENCE' else 'https://rest.uniprot.org/uniprotkb/search','replacement_reason':reason})
        if seq:(D/f'{a}.fasta').write_text(f'>{a} resolved={resolved}\n{seq}\n')
    wr(pd.DataFrame(rows),'canonical_sequence_recovery.tsv');print('SEQUENCE_RECOVERED',sum(x['status'] in ['EXACT','OBSOLETE_REPLACED'] for x in rows),'/',len(rows),flush=True)
if __name__=='__main__':
    freeze();denominator();sequences()
