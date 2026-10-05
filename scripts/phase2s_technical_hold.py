#!/usr/bin/env python3
"""Materialize available benchmark baselines; keep unexecuted model outcomes NA."""
from pathlib import Path
import hashlib
import json
import platform
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'results/phase2S'

def rd(name): return pd.read_csv(P/name,sep='\t',dtype=str,keep_default_na=False)
def write(frame,name): frame.to_csv(P/name,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    assert (P/'input_freeze_manifest.tsv').exists(), 'Run phase2s_prepare.py first'
    variants=rd('frozen_variants.tsv'); partners=rd('frozen_partner_sets.tsv')
    universes=rd('frozen_gene_universes.tsv'); groups=rd('frozen_disease_groups.tsv')
    u1=set(universes.loc[universes.U1_STRUCTURAL_42.eq('True'),'gene']);u2=set(universes.gene)
    pk=set(zip(partners.gene,partners.partner)); vk=set(zip(variants.gene,variants.VariationID))
    eligible={g:int(n) for g,n in partners.groupby('gene').partner.nunique().items()}
    model_rows=[]
    source=ROOT/'data/software/phase2S'
    for method,version,commit,files in [
        ('SPPIDER_SEQ_PARTNER_INTERFACE','v2.0-manuscript','9c12e8c826ad4e47b00477cab56f6d58715d8f60',['SPPIDER-seq-v2.0-manuscript/models/crossattn_rec_run27_best.pt','SPPIDER-seq-v2.0-manuscript/models/crossattn_pep_run07_best.pt']),
        ('XPAIR_PARTNER_INTERFACE','0.1.0','PyPI release 2026-07-27',['xpair-0.1.0/xpair/default_model/multitask_xfair.ckpt','xpair-0.1.0-py3-none-any.whl']),
    ]:
        for name in files:
            f=source/name
            model_rows.append({'method':method,'version':version,'source_revision':commit,'file':str(f.relative_to(ROOT)),'size_bytes':f.stat().st_size,'sha256':digest(f),'checkpoint_status':'ACQUIRED','technical_validation':'NOT_RUN_ENCODING_MODEL_TLS_FAILURE','encoding_model':'facebook/esm2_t33_650M_UR50D' if method.startswith('SPPIDER') else 'ElnaggarLab/ankh-large','encoding_model_revision':'NA','encoding_model_sha256':'NA','device':'CPU','environment':'.venv-phase2S (installation incomplete)','length_rule':'official full-length chunking, 1024 tokens / stride512' if method.startswith('SPPIDER') else '50<=query_length,partner_length<=2000; no truncation','interface_threshold':'0.5 secondary only'})
    write(pd.DataFrame(model_rows),'method_manifest.tsv')
    records=[]
    counts={}
    a=pd.read_csv(ROOT/'results/phase1B0R3/mutpred_ppi_v1_scores.tsv',sep='\t',dtype=str,keep_default_na=False)
    a=a.loc[a.gene.isin(u2)]
    aa=a.loc[a.score_status.eq('SCORED')]
    aa=aa.loc[[ (g,v) in vk and (g,p) in pk for g,v,p in aa[['gene','VariationID','partner']].itertuples(index=False,name=None)]]
    counts['MUTPRED_PPI_DISRUPTION']=aa.groupby(['gene','VariationID']).partner.nunique().to_dict()
    # Stream the frozen structural table instead of reconstructing interfaces.
    path=ROOT/'results/phase1A/partner_specific_matrix_structural.tsv'
    structural_path=P/'representation_matrix_structural.tsv'
    first=True;bcounts={};dcounts={}
    for chunk in pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False,chunksize=100000):
        chunk=chunk.loc[chunk.gene.isin(u2)]
        chunk=chunk.loc[[ (g,v) in vk and (g,p) in pk for g,v,p in chunk[['gene','VariationID','partner']].itertuples(index=False,name=None)]]
        available=chunk.interface_status.isin(['INTERFACE','NON_INTERFACE_IN_OBSERVED_COMPLEX','NON_INTERFACE_IN_PIONEER_INTERFACE_MODEL'])
        chunk['binary_interface']=np.where(available,np.where(chunk.interface_status.eq('INTERFACE'),1,0),np.nan)
        chunk['binary_available']=available
        numeric=pd.to_numeric(chunk.interface_min_distance_A,errors='coerce')
        chunk['distance_available']=np.isfinite(numeric)
        chunk[['gene','VariationID','partner','tier','structure_source','interface_status','binary_interface','binary_available','interface_min_distance_A','distance_available']].to_csv(structural_path,sep='\t',index=False,mode='w' if first else 'a',header=first,na_rep='NA')
        first=False
        for flag,lookup in [('binary_available',bcounts),('distance_available',dcounts)]:
            values=chunk.loc[chunk[flag]].groupby(['gene','VariationID']).partner.nunique()
            for key,n in values.items():lookup[key]=lookup.get(key,0)+int(n)
    counts['STRUCTURAL_BINARY_INTERFACE']=bcounts;counts['STRUCTURAL_INTERFACE_DISTANCE']=dcounts
    methods=list(counts)+['SPPIDER_SEQ_PARTNER_INTERFACE','XPAIR_PARTNER_INTERFACE','SEQUENCE_CONSENSUS_AVAILABILITY','CONSENSUS_HIGH_SECONDARY']
    geneout=[];diseaseout=[];varout=[];summary=[]
    primary=variants.frozen_primary_single_label.eq('True')
    for method in methods:
        run=method in counts; lookup=counts.get(method,{})
        for r in variants.itertuples(index=False):
            n=lookup.get((r.gene,r.VariationID),0) if run else np.nan
            varout.append({'method':method,'gene':r.gene,'VariationID':r.VariationID,'primary_single_label':r.frozen_primary_single_label,'available_partner_n':n,'eligible_partner_n':eligible.get(r.gene,0),'coverage_ge1':n>=1 if run else 'NA','coverage_ge3':n>=3 if run else 'NA','status':'EVALUATED' if run else 'NOT_RUN_MODEL_DEPENDENCY'})
        for gene in sorted(u2):
            v=variants.loc[variants.gene.eq(gene)];vp=v.loc[v.frozen_primary_single_label.eq('True')]
            n3=sum(lookup.get((gene,x),0)>=3 for x in vp.VariationID) if run else np.nan
            rate=100*n3/len(vp) if run and len(vp) else np.nan
            rates=[];retained=0
            for r in groups.loc[groups.gene.eq(gene)].itertuples(index=False):
                members=v.loc[v.strict_MONDO_groups.str.split(';').map(lambda x:r.MONDO in x)]
                represented=sum(lookup.get((gene,x),0)>=3 for x in members.loc[members.primary_cohort_eligible.eq('YES'),'VariationID']) if run else np.nan
                denominator=int(r.n_unique_variants)
                drate=100*represented/denominator if run else np.nan
                rates.append(drate)
                unambig=denominator-int(members.multilabel.eq('True').sum())
                retained+=unambig>=10
                diseaseout.append({'method':method,'gene':gene,'MONDO':r.MONDO,'strict_variant_n':denominator,'available_ge3_n':represented,'coverage_ge3_pct':drate,'status':'EVALUATED' if run else 'NOT_RUN_MODEL_DEPENDENCY'})
            cls='NA'
            if run:
                cls='SEVERE_IMBALANCE' if not rates or min(rates)<50 or max(rates)-min(rates)>=50 else 'MATERIAL_IMBALANCE' if max(rates)-min(rates)>=20 or min(rates)<70<=max(rates) else 'MINOR_IMBALANCE'
            passing=bool(rate>=70 and cls=='MINOR_IMBALANCE' and retained>=2) if run else 'NA'
            geneout.append({'method':method,'gene':gene,'in_U1':gene in u1,'primary_single_label_variant_n':len(vp),'primary_ge3_variant_n':n3,'primary_ge3_coverage_pct':rate,'differential_coverage_class':cls,'balanced_primary_compatible':passing,'frozen_partner_n':eligible.get(gene,0),'disease_group_coverage_range_pp':max(rates)-min(rates) if run and rates else 'NA','status':'EVALUATED' if run else 'NOT_RUN_MODEL_DEPENDENCY'})
        for universe,genes in [('U1_STRUCTURAL_42',u1),('U2_STRICT_132',u2)]:
            gv=[x for x in geneout if x['method']==method and x['gene'] in genes]
            vs=variants.loc[variants.gene.isin(genes)]
            denominator=sum(eligible.get(r.gene,0) for r in vs.itertuples(index=False))
            available=sum(lookup.get((r.gene,r.VariationID),0) for r in vs.itertuples(index=False)) if run else np.nan
            summary.append({'method':method,'universe':universe,'gene_n':len(genes),'eligible_variant_partner_rows':denominator,'available_rows':available,'row_coverage_pct':100*available/denominator if run and denominator else 'NA','variant_ge1_n':sum(lookup.get((r.gene,r.VariationID),0)>=1 for r in vs.itertuples(index=False)) if run else 'NA','variant_ge3_n':sum(lookup.get((r.gene,r.VariationID),0)>=3 for r in vs.itertuples(index=False)) if run else 'NA','genes_ge70_n':sum(x['primary_ge3_coverage_pct']>=70 for x in gv) if run else 'NA','balanced_genes_n':sum(x['balanced_primary_compatible'] is True for x in gv) if run else 'NA','status':'EVALUATED' if run else 'NOT_RUN_MODEL_DEPENDENCY'})
    write(pd.DataFrame(varout),'coverage_by_variant.tsv');write(pd.DataFrame(geneout),'coverage_by_gene.tsv');write(pd.DataFrame(diseaseout),'coverage_by_disease_group.tsv');write(pd.DataFrame(geneout),'balanced_gene_status.tsv');write(pd.DataFrame(summary),'method_comparison_summary.tsv')
    # C2/C3 retain explicit denominators and fractions alongside counts.
    for row in summary:
        universe_genes=u1 if row['universe']=='U1_STRUCTURAL_42' else u2
        total=int(variants.gene.isin(universe_genes).sum())
        row['eligible_variant_n']=total
        for endpoint in ['variant_ge1','variant_ge3']:
            row[endpoint+'_pct']=100*row[endpoint+'_n']/total if row['status']=='EVALUATED' and total else 'NA'
    write(pd.DataFrame(summary),'method_comparison_summary.tsv')
    for method in ['sppider','xpair']:
        write(pd.DataFrame(columns=['query_sequence_sha256','partner_sequence_sha256','model_sha256','residue','interface_probability']),method+'_raw_predictions.tsv')
        write(pd.DataFrame(columns=['gene','VariationID','partner','interface_probability','status']),method+'_variant_partner_scores.tsv')
        write(pd.DataFrame(columns=['gene','VariationID','partner','interface_probability','status']),'representation_matrix_'+method+'.tsv')
    write(pd.DataFrame([{'method':m,'status':'NOT_RUN_ENCODING_MODEL_TLS_FAILURE','AUROC':'NA','AUPRC':'NA','MCC':'NA','bootstrap_CI':'NA'} for m in ['SPPIDER_SEQ_PARTNER_INTERFACE','XPAIR_PARTNER_INTERFACE']]),'structural_validation.tsv')
    write(pd.DataFrame([{'method':m,'status':'NOT_RUN_ENCODING_MODEL_TLS_FAILURE','partner_specificity':'NA','null_control':'NA'} for m in ['SPPIDER_SEQ_PARTNER_INTERFACE','XPAIR_PARTNER_INTERFACE']]),'partner_specificity_audit.tsv')
    write(pd.DataFrame([{'gene':g,'rescue_status':'NOT_EVALUABLE_MODEL_DEPENDENCY','sppider':'NA','xpair':'NA','both':'NA'} for g in sorted(u2)]),'sequence_rescued_genes.tsv')
    write(pd.DataFrame([{'method':m,'stage':'ENCODING_MODEL_ACQUISITION','reason':'OFFICIAL_HUGGINGFACE_TLS_FAILURE','affected_units':'MODEL_LEVEL','counts_as_representation_failure':'NO'} for m in ['SPPIDER_SEQ_PARTNER_INTERFACE','XPAIR_PARTNER_INTERFACE']]),'technical_failures_by_method.tsv')
    predictors=pd.read_csv(ROOT/'results/phase1A/coverage_by_gene.tsv',sep='\t',dtype=str)
    predictors=predictors.loc[predictors.gene.isin(u2)].copy();predictors['structural_42_retained']=predictors.gene.isin(u1)
    write(predictors[['gene','structural_42_retained','credible_partner_n','Tier1_2_partner_n','strict_disease_group_n','unique_strict_group_variants_n','protein_mapping_rate_pct','variants_with_observed_PDB_residue_coordinates_n']],'coverage_failure_predictors.tsv')
    print(json.dumps({'baseline_summary':summary[:6],'status':'TECHNICAL_HOLD_SEQUENCE_MODELS'},ensure_ascii=False),flush=True)

if __name__=='__main__':main()
