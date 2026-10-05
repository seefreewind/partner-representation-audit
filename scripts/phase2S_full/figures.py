from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt,numpy as np
from matplotlib.colors import ListedColormap
F=O/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'svg.fonttype':'none','pdf.fonttype':42,'font.size':8,'axes.spines.right':False,'axes.spines.top':False,'legend.frameon':False})
methods=['MUTPRED_PPI_DISRUPTION','STRUCTURAL_BINARY_INTERFACE','STRUCTURAL_INTERFACE_DISTANCE','SPPIDER_SEQ_PARTNER_INTERFACE','XPAIR_PARTNER_INTERFACE'];display=['MutPred','Structural binary','Structural distance','SPPIDER','X-PAIR'];colors=['#777777','#A7B9C7','#728DA5','#67A798','#A790BA'];g=pd.read_csv(O/'coverage_by_gene_all_methods.tsv',sep='\t');s=pd.read_csv(O/'method_comparison_summary.tsv',sep='\t');s=s.loc[s.universe.eq('U2')].set_index('method');v=pd.read_csv(O/'structural_validation_metrics.tsv',sep='\t');features=pd.read_csv(O/'failure_predictor_features.tsv',sep='\t')
contract={'archetype':'quantitative grid','backend':'python','export':'SVG editable text and PNG at600dpi; no PDF','claims':{'S1':'Method attrition from132 genes to balanced availability','S2':'Alphabetical gene-level C4 and sequence balance classes','S3':'Within-gene paired U1 C4 comparisons','S4':'Historical MutPred3 and common-denominator balanced counts remain distinguished','S5':'Grouped experimental interface validity with1000 protein bootstraps; training overlap unknown','S6':'Technical input features and representation failure are compared descriptively'},'source_data':'phase2S_full TSVs, unchanged numbers','glyph_min_pt':6,'QA':'Final rendered visual review pending agent inspection'};(F/'figure_contract.json').write_text(json.dumps(contract,indent=2)+'\n')
def save(fig,name):
    fig.savefig(F/(name+'.svg'),bbox_inches='tight');fig.savefig(F/(name+'.png'),dpi=600,bbox_inches='tight');plt.close(fig)
fig,axes=plt.subplots(1,5,figsize=(7.2,3.4),sharey=True)
for ax,m,name,c in zip(axes,methods,display,colors):
    row=s.loc[m];sub=g.loc[g.method.eq(m)];vals=[132,int(row.method_scoreable_genes),int(sub.GE3_n.gt(0).sum()),int(row.genes_C4_ge70),int(row.balanced_genes)];ax.barh(np.arange(5),vals,color=c);ax.set_yticks(range(5),['U2','Scoreable','Any GE3 variant','C4 ≥70%','Balanced']);ax.set_xlim(0,145);ax.set_title(name);ax.set_xlabel('Genes')
    for y,n in enumerate(vals):ax.text(n+2,y,str(n),va='center',fontsize=7)
axes[0].invert_yaxis();fig.tight_layout();save(fig,'Figure_S1_attrition')
mat=g.loc[g.method.isin(methods)].pivot(index='gene',columns='method',values='C4_fraction').reindex(columns=methods).sort_index();fig,axes=plt.subplots(1,2,figsize=(7,17),gridspec_kw={'width_ratios':[5,2]},sharey=True)
im=axes[0].imshow(mat,vmin=0,vmax=1,aspect='auto',cmap='viridis');axes[0].set_xticks(range(5),display,rotation=45,ha='right',rotation_mode='anchor');axes[0].set_yticks(range(len(mat)),mat.index,fontsize=6);fig.colorbar(im,ax=axes[0],fraction=.045,pad=.03,label='C4 fraction')
cls={'SEVERE_IMBALANCE':0,'MATERIAL_IMBALANCE':1,'MINOR_IMBALANCE':2};b=g.loc[g.method.isin(methods[-2:])].pivot(index='gene',columns='method',values='balance_class').reindex(index=mat.index,columns=methods[-2:]);a=np.array([[cls.get(x,np.nan) for x in r] for r in b.values]);im=axes[1].imshow(a,vmin=-.5,vmax=2.5,aspect='auto',cmap=ListedColormap(['#D49A9A','#D9C58D','#8FC1AE']));axes[1].set_xticks([0,1],['SPPIDER','X-PAIR'],rotation=45,ha='right',rotation_mode='anchor');cb=fig.colorbar(im,ax=axes[1],ticks=[0,1,2],fraction=.045,pad=.03);cb.ax.set_yticklabels(['Severe','Material','Minor']);fig.tight_layout();save(fig,'Figure_S2_coverage_heatmap')
paired=g.loc[g.U1&g.method.isin(methods)].pivot(index='gene',columns='method',values='C4_fraction').reindex(columns=methods).sort_index();fig,ax=plt.subplots(figsize=(7,8));
for i,(gene,r) in enumerate(paired.iterrows()):
    xs=r.to_numpy();ys=i+np.linspace(-.23,.23,5);ax.plot(xs,ys,color='#BBBBBB',lw=.5)
    for j,c in enumerate(colors):ax.scatter(xs[j],ys[j],color=c,s=13,label=display[j] if i==0 else None)
ax.set_yticks(range(len(paired)),paired.index,fontsize=7);ax.set_xlim(-.03,1.03);ax.set_xlabel('C4 fraction');ax.axvline(.7,color='#555555',ls='--',lw=.7);ax.legend(loc='upper center',bbox_to_anchor=(.5,1.08),ncol=3);fig.tight_layout();save(fig,'Figure_S3_U1_paired')
labels=['Historical\nMutPred']+display+['Both available','Consensus high'];vals=[3]+[int(s.loc[m,'balanced_genes']) for m in methods]+[int(s.loc[m,'balanced_genes']) for m in ['SEQUENCE_CONSENSUS_AVAILABILITY','CONSENSUS_HIGH_SECONDARY']];fig,ax=plt.subplots(figsize=(7.2,3.5));ax.bar(range(len(vals)),vals,color=['#BBBBBB']+colors+['#96BBA9','#608D7B']);ax.set_xticks(range(len(vals)),labels,rotation=35,ha='right',rotation_mode='anchor');ax.set_ylabel('Balanced genes in U2');ax.set_ylim(0,max(vals)*1.18+1)
for i,n in enumerate(vals):ax.text(i,n+1,str(n),ha='center',fontsize=7)
fig.tight_layout();save(fig,'Figure_S4_balanced_counts')
fig,axes=plt.subplots(1,2,figsize=(6.5,3),sharey=True)
for ax,metric in zip(axes,['AUROC','AUPRC']):
    for i,m in enumerate(['sppider','xpair']):
        rr=v.loc[v.method.eq(m)&v.aggregation.eq('macro')&v.metric.eq(metric)]
        if len(rr):r=rr.iloc[0];ax.errorbar(i,r.estimate,yerr=[[max(0,r.estimate-r.CI95_low)],[max(0,r.CI95_high-r.estimate)]],fmt='o',color=colors[i+3],capsize=3)
    ax.set_xticks([0,1],['SPPIDER','X-PAIR']);ax.set_ylim(0,1);ax.set_title('Macro '+metric);ax.set_ylabel('Estimate (95% protein bootstrap CI)')
fig.tight_layout();save(fig,'Figure_S5_structural_validation')
fig,axes=plt.subplots(2,3,figsize=(7.2,5.5));names=['query_length','credible_partner_n','PDB_variant_mapping_fraction','canonical_mapping_completeness','structurally_mappable_variant_fraction','query_XPAIR_length_supported']
for ax,name in zip(axes.flat,names):
    for m,label,c in zip(methods,display,colors):
        d=g.loc[g.method.eq(m)].merge(features,on='gene');ax.scatter(d[name],d.C4_fraction,s=9,alpha=.6,color=c,label=label)
    ax.set_xlabel({'query_length':'Query length (aa)','credible_partner_n':'Credible partner count','PDB_variant_mapping_fraction':'PDB variant mapping fraction','canonical_mapping_completeness':'Canonical mapping fraction','structurally_mappable_variant_fraction':'Structurally mappable fraction','query_XPAIR_length_supported':'X-PAIR query length support'}[name]);ax.set_ylabel('C4 fraction');ax.axhline(.7,color='#AAAAAA',lw=.5)
handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,1.01),fontsize=7,ncol=5);fig.tight_layout(rect=[0,0,1,.95]);save(fig,'Figure_S6_failure_architecture');print('FIGURES_GENERATED_VISUAL_REVIEW_PENDING',flush=True)
