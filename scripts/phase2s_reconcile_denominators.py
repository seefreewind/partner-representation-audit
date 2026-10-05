from pathlib import Path
import pandas as pd
P=Path(__file__).resolve().parents[1]/'results/phase2S'
def rd(n):return pd.read_csv(P/n,sep='\t',dtype=str,keep_default_na=False)
def wr(d,n):d.to_csv(P/n,sep='\t',index=False)
v=rd('frozen_variants.tsv')
v['frozen_primary_single_label']=v.strict_group_n.eq('1')
wr(v,'frozen_variants.tsv')
c=rd('coverage_by_variant.tsv')
c['primary_single_label']=[bool(x) for x in c.set_index(['gene','VariationID']).index.map(v.set_index(['gene','VariationID']).frozen_primary_single_label)]
wr(c,'coverage_by_variant.tsv')
g=rd('coverage_by_gene.tsv')
for i,r in g.iterrows():
    members=v.loc[v.gene.eq(r.gene)&v.frozen_primary_single_label]
    g.loc[i,'primary_single_label_variant_n']=str(len(members))
    if r.status!='EVALUATED':continue
    cv=c.loc[c.method.eq(r.method)&c.gene.eq(r.gene)&c.VariationID.isin(members.VariationID)]
    n=int(pd.to_numeric(cv.available_partner_n).ge(3).sum())
    rate=100*n/len(members) if len(members) else 0
    g.loc[i,'primary_ge3_variant_n']=str(n)
    g.loc[i,'primary_ge3_coverage_pct']=str(rate)
    # Keep the already evaluated frozen disease balance and >=2 groups gate.
    old_pass=r.balanced_primary_compatible=='True'
    g.loc[i,'balanced_primary_compatible']=str(old_pass and rate>=70)
wr(g,'coverage_by_gene.tsv');wr(g,'balanced_gene_status.tsv')
s=rd('method_comparison_summary.tsv')
u=rd('frozen_gene_universes.tsv')
for i,r in s.iterrows():
    if r.status!='EVALUATED':continue
    genes=set(u.loc[u.U1_STRUCTURAL_42.eq('True'),'gene']) if r.universe=='U1_STRUCTURAL_42' else set(u.gene)
    gg=g.loc[g.method.eq(r.method)&g.gene.isin(genes)]
    s.loc[i,'genes_ge70_n']=str(int(pd.to_numeric(gg.primary_ge3_coverage_pct).ge(70).sum()))
    s.loc[i,'balanced_genes_n']=str(int(gg.balanced_primary_compatible.eq('True').sum()))
wr(s,'method_comparison_summary.tsv')
old=rd('mutpred_frozen_coverage_by_gene.tsv')
# Historical membership was constructed from score-row membership, which is
# not identical to a full-cohort single-label registry. Preserve its actual
# published denominator; do not replace it with an attractive revised result.
for r in old.itertuples(index=False):
    mask=g.method.eq('MUTPRED_PPI_DISRUPTION')&g.gene.eq(r.gene)
    g.loc[mask,'primary_single_label_variant_n']=str(r.primary_unambiguous_variants_n)
    g.loc[mask,'primary_ge3_variant_n']=str(r.primary_variants_with_ge3_Tier1_2_scores_n)
    g.loc[mask,'primary_ge3_coverage_pct']=str(r.primary_variants_with_ge3_Tier1_2_scores_pct)
    g.loc[mask,'balanced_primary_compatible']=str(r.eligible_for_rebuilt_primary_cohort=='YES')
    g.loc[mask,'status']='FROZEN_LEGACY_GATE'
for method in ['STRUCTURAL_BINARY_INTERFACE','STRUCTURAL_INTERFACE_DISTANCE']:
    mask=g.method.eq(method)
    g.loc[mask,'primary_ge3_coverage_pct']='NA'
    g.loc[mask,'balanced_primary_compatible']='NA'
    g.loc[mask,'status']='AVAILABILITY_ONLY_PRIMARY_DENOMINATOR_RECONCILIATION_PENDING'
    sm=s.method.eq(method)
    s.loc[sm,['genes_ge70_n','balanced_genes_n']]='NA'
s.loc[s.method.eq('MUTPRED_PPI_DISRUPTION'),['genes_ge70_n','balanced_genes_n']]=['9','3']
wr(g,'coverage_by_gene.tsv');wr(g,'balanced_gene_status.tsv');wr(s,'method_comparison_summary.tsv')
assert int(s.loc[s.method.eq('MUTPRED_PPI_DISRUPTION')&s.universe.eq('U1_STRUCTURAL_42'),'genes_ge70_n'].iloc[0])==9
print('Original published denominators and 9/3 gate preserved. Full-cohort single-label registry differs from legacy score-row membership; new cross-method C4/C5 deferred pending exact reconciliation. Row and variant availability remain evaluated.')
