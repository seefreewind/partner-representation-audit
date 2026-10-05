from common import *
import numpy as np,importlib.util
METHODS=['MUTPRED_PPI_DISRUPTION','STRUCTURAL_BINARY_INTERFACE','STRUCTURAL_INTERFACE_DISTANCE','SPPIDER_SEQ_PARTNER_INTERFACE','XPAIR_PARTNER_INTERFACE','SEQUENCE_CONSENSUS_AVAILABILITY','CONSENSUS_HIGH_SECONDARY']
u=rd('results/phase2S/frozen_gene_universes.tsv');u1=set(u.loc[u.U1_STRUCTURAL_42.eq('True'),'gene']);u2=set(u.gene);p=rd('results/phase2S/frozen_partner_sets.tsv');pk=set(zip(p.gene,p.partner));partner_n=p.groupby('gene').partner.nunique().to_dict()
v=rd('results/phase2S/frozen_variants.tsv');m=rd('results/phase2SR/common_primary_membership.tsv');primary=set(zip(m.loc[m.included_in_common_primary_denominator.eq('True'),'gene'],m.loc[m.included_in_common_primary_denominator.eq('True'),'VariationID']))
spec=importlib.util.spec_from_file_location('legacy',R/'scripts/phase1b0r3_finalize.py');legacy=importlib.util.module_from_spec(spec);spec.loader.exec_module(legacy);strict,vg,_=legacy.make_membership(u2);vg.VariationID=vg.VariationID.astype(str)
multilabel=vg.groupby(['gene','VariationID']).MONDO.nunique().gt(1).to_dict()
wr(vg,'frozen_group_membership_exact.tsv')
def counts_for(method):
    if method=='MUTPRED_PPI_DISRUPTION':
        d=rd('results/phase1B0R3/mutpred_ppi_v1_scores.tsv');d=d.loc[d.score_status.eq('SCORED')]
    elif method.startswith('STRUCTURAL'):
        sets={};binary=method=='STRUCTURAL_BINARY_INTERFACE'
        for d in pd.read_csv(R/'results/phase1A/partner_specific_matrix_structural.tsv',sep='\t',dtype=str,keep_default_na=False,chunksize=100000):
            mask=d.interface_status.isin(['INTERFACE','NON_INTERFACE_IN_OBSERVED_COMPLEX','NON_INTERFACE_IN_PIONEER_INTERFACE_MODEL']) if binary else np.isfinite(pd.to_numeric(d.interface_min_distance_A,errors='coerce'))
            d=d.loc[mask & d.gene.isin(u2)]
            for g,vid,par in d[['gene','VariationID','partner']].itertuples(index=False,name=None):
                if (g,par) in pk:sets.setdefault((g,vid),set()).add(par)
        return {k:len(x) for k,x in sets.items()},sets
    else:
        name={'SPPIDER_SEQ_PARTNER_INTERFACE':'variant_partner_scores_sppider.tsv','XPAIR_PARTNER_INTERFACE':'variant_partner_scores_xpair.tsv','SEQUENCE_CONSENSUS_AVAILABILITY':'sequence_consensus_scores.tsv','CONSENSUS_HIGH_SECONDARY':'sequence_consensus_scores.tsv'}[method]
        if not (O/name).exists():return None,None
        sets={}
        for d in pd.read_csv(O/name,sep='\t',dtype=str,keep_default_na=False,chunksize=100000):
            flag=d.status.eq('SCORED') if 'status' in d else d.both_available.eq('True')
            if method=='CONSENSUS_HIGH_SECONDARY':flag&=d.consensus_high.eq('True')
            for g,vid,par in d.loc[flag,['gene','VariationID','partner']].itertuples(index=False,name=None):sets.setdefault((g,vid),set()).add(par)
        return {k:len(x) for k,x in sets.items()},sets
    d=d.loc[d.gene.isin(u2)];sets={}
    for g,vid,par in d[['gene','VariationID','partner']].itertuples(index=False,name=None):
        if (g,par) in pk:sets.setdefault((g,vid),set()).add(par)
    return {k:len(x) for k,x in sets.items()},sets
varrows=[];generows=[];diseaserows=[];summaries=[]
for method in METHODS:
    counts,sets=counts_for(method)
    if counts is None:continue
    for (g,vid) in sorted(primary):
        n=counts.get((g,vid),0);varrows.append({'method':method,'gene':g,'VariationID':vid,'common_primary':True,'partners_with_score':n,'GE1':n>=1,'GE3':n>=3,'frozen_partner_n':partner_n.get(g,0)})
    for g in sorted(u2):
        ids={vid for gg,vid in primary if gg==g};n1=sum(counts.get((g,x),0)>=1 for x in ids);n3=sum(counts.get((g,x),0)>=3 for x in ids);c4=n3/len(ids) if ids else np.nan;rates=[];retained=0
        for group in strict.loc[strict.gene.eq(g)].itertuples(index=False):
            members=set(vg.loc[vg.gene.eq(g)&vg.MONDO.eq(group.MONDO),'VariationID']);den=int(group.n_unique_variants);multi=sum(multilabel.get((g,x),False) for x in members);retain=den-multi>=10
            a=sum(counts.get((g,x),0)>=1 for x in members);b=sum(counts.get((g,x),0)>=3 for x in members);fraction=b/den if den else np.nan
            # Frozen legacy uses all strict >=10 groups for imbalance, while
            # retaining >=2 unambiguous groups is a separate eligibility gate.
            rates.append(fraction);retained+=retain
            diseaserows.append({'method':method,'gene':g,'MONDO':group.MONDO,'frozen_variant_denominator':den,'mapped_members_n':len(members),'multilabel_n':multi,'retained_unambiguous_ge10':retain,'GE1_n':a,'GE3_n':b,'GE3_fraction':fraction,'GE3_pct':100*fraction})
        rg=max(rates)-min(rates) if rates else np.nan
        cls='SEVERE_IMBALANCE' if not rates or min(rates)<.5 or rg>=.5 else 'MATERIAL_IMBALANCE' if rg>=.2 or min(rates)<.7<=max(rates) else 'MINOR_IMBALANCE'
        dimensions=len(set().union(*(sets.get((g,x),set()) for x in ids))) if ids else 0
        balanced=bool(c4>=.7 and cls=='MINOR_IMBALANCE' and retained>=2 and dimensions>=3)
        generows.append({'method':method,'gene':g,'U1':g in u1,'primary_denominator':len(ids),'GE1_n':n1,'GE3_n':n3,'C4_fraction':c4,'C4_pct':100*c4,'C4_ge70':bool(c4>=.7),'available_partners_n':dimensions,'balance_range_pp':100*rg,'balance_class':cls,'retained_disease_groups_n':retained,'balanced_PASS':balanced,'method_scoreable':dimensions>=1})
    for label,gs in [('U1',u1),('U2',u2)]:
        subset=[r for r in generows if r['method']==method and r['gene'] in gs];vs=[r for r in varrows if r['method']==method and r['gene'] in gs]
        den=sum(r['frozen_partner_n'] for r in vs);avail=sum(r['partners_with_score'] for r in vs)
        summaries.append({'method':method,'universe':label,'gene_n':len(gs),'variant_denominator':len(vs),'variant_partner_denominator':den,'available_variant_partner_rows':avail,'variant_pair_coverage_fraction':avail/den if den else np.nan,'variant_GE1':sum(r['GE1'] for r in vs),'variant_GE3':sum(r['GE3'] for r in vs),'method_scoreable_genes':sum(r['method_scoreable'] for r in subset),'genes_C4_ge70':sum(r['C4_ge70'] for r in subset),'balanced_genes':sum(r['balanced_PASS'] for r in subset),'status':'EVALUATED_COMMON_DENOMINATOR'})
wr(pd.DataFrame(varrows),'coverage_by_variant_all_methods.tsv');wr(pd.DataFrame(generows),'coverage_by_gene_all_methods.tsv');wr(pd.DataFrame(diseaserows),'coverage_by_disease_group_all_methods.tsv');wr(pd.DataFrame(generows),'balanced_gene_status_all_methods.tsv');wr(pd.DataFrame(summaries),'method_comparison_summary.tsv')
df=pd.DataFrame(generows);b={method:set(df.loc[df.method.eq(method)&df.balanced_PASS,'gene']) for method in METHODS};historical={'FGFR2','MEN1','TUBA1A'};rescued=[]
for g in sorted(u2):
    sp=g in b['SPPIDER_SEQ_PARTNER_INTERFACE'];xp=g in b['XPAIR_PARTNER_INTERFACE'];struct=any(g in b[x] for x in METHODS[:3])
    rescued.append({'gene':g,'historical_MutPred_balanced':g in historical,'SPPIDER_RESCUED':sp and g not in historical,'XPAIR_RESCUED':xp and g not in historical,'SPPIDER_only_rescued':sp and not xp and g not in historical,'XPAIR_only_rescued':xp and not sp and g not in historical,'BOTH_RESCUED':sp and xp,'SEQUENCE_ONLY_RESCUED':(sp or xp) and not struct,'evaluability':'FULL' if METHODS[3] in set(df.method) and METHODS[4] in set(df.method) else 'BASELINES_ONLY_SEQUENCE_PENDING'})
wr(pd.DataFrame(rescued),'rescued_gene_classification.tsv')
wr(df.loc[df.U1].pivot(index='gene',columns='method',values='C4_fraction').reset_index(),'U1_paired_C4_matrix.tsv')
print('COVERAGE_WRITTEN',sorted(df.method.unique()),'Historical9/3unchanged',flush=True)
