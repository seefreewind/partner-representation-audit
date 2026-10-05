from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'results/phase2P1';F=O/'figures';F.mkdir(exist_ok=True)
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'svg.fonttype':'none','pdf.fonttype':42,'font.size':8,'axes.titlesize':9,'axes.spines.right':False,'axes.spines.top':False,'axes.linewidth':.8,'legend.frameon':False,'xtick.labelsize':7,'ytick.labelsize':8})
blue='#27698E';gray='#8B9298';orange='#B87836';genes=pd.read_csv(O/'xpair_primary_per_gene.tsv',sep='\t').gene.tolist();y=np.arange(10)
x=pd.read_csv(O/'xpair_primary_per_gene.tsv',sep='\t').set_index('gene').loc[genes];s=pd.read_csv(O/'sppider_confirmatory_per_gene.tsv',sep='\t').set_index('gene').loc[genes];b=pd.read_csv(O/'baseline_comparisons.tsv',sep='\t');ctx=pd.read_csv(O/'residue_ceiling_context.tsv',sep='\t').set_index('gene').loc[genes]
contract={'backend':'python','archetype':'quantitative grid','export':'SVG editable text + PNG 600 dpi; no PDF requested','width_mm':183,'glyph_min_pt':6.5,'selection':'All 10 genes in frozen order; permutation example alphabetically first eligible COL3A1','panels':{'P1':'All frozen partner-identity effects with zero reference; no per-gene uncertainty was preregistered','P2':'Within-gene paired labeled and sorted R²','P3':'Cross-method effect concordance without revising primary claim','P4':'Frozen generic baselines, missingness and row-matched X-PAIR contrasts','P5':'Frozen same-residue discordance and ceiling alongside primary effect; separate axes; descriptive','P6':'Neutral COL3A1 paired permutation null with observed primary ΔR²'},'data_source':'results/phase2P1/*.tsv and xpair_permutations/COL3A1.npz','replicate_unit':'One gene-residue observation; global summaries weight genes equally','uncertainty':'No per-gene CI; global gene-bootstrap CI is reported in text'}
(F/'figure_contract.json').write_text(json.dumps(contract,indent=2)+'\n')
def save(fig,name):
 fig.savefig(F/(name+'.svg'));fig.savefig(F/(name+'.png'),dpi=600);plt.close(fig)
def baseaxis(ax,title,xlabel):
 ax.set_yticks(y,genes);ax.invert_yaxis();ax.set_title(title,loc='left',pad=12);ax.set_xlabel(xlabel);ax.grid(axis='x',color='#E9ECEF',linewidth=.5);ax.set_axisbelow(True)
fig,ax=plt.subplots(figsize=(183/25.4,110/25.4));fig.subplots_adjust(left=.15,right=.94,bottom=.16,top=.88);d=x.DELTA_R2_PARTNER_ID.to_numpy();ax.axvline(0,color=gray,lw=.9);ax.scatter(d,y,c=blue,s=25);baseaxis(ax,'P1  Partner-identity effect across the frozen cohort','ΔR² partner identity (X-PAIR labeled − sorted)');ax.set_xlim(min(d)-.006,max(d)+.009)
for a,j in zip(d,y):ax.annotate(f'{a:.5f}',(a,j),xytext=(7,0),textcoords='offset points',va='center',fontsize=7)
save(fig,'FigureP1')
fig,ax=plt.subplots(figsize=(183/25.4,110/25.4));fig.subplots_adjust(left=.15,right=.94,bottom=.16,top=.85)
for j,gene in enumerate(genes):ax.plot([x.loc[gene,'R2_XPAIR_SORTED'],x.loc[gene,'R2_XPAIR_LABELED']],[j,j],color='#C7CDD1',lw=1)
ax.scatter(x.R2_XPAIR_LABELED,y,color=blue,s=26,label='X-PAIR labeled');ax.scatter(x.R2_XPAIR_SORTED,y,color=orange,s=24,marker='s',label='X-PAIR sorted');baseaxis(ax,'P2  Paired phenotype information','PERMANOVA R²');ax.legend(loc='lower right',fontsize=8);save(fig,'FigureP2')
fig,ax=plt.subplots(figsize=(183/25.4,112/25.4));fig.subplots_adjust(left=.13,right=.94,bottom=.16,top=.84);a=x.DELTA_R2_PARTNER_ID.to_numpy();bb=s.DELTA_R2_SPPIDER.to_numpy();lo=min(a.min(),bb.min())-.003;hi=max(a.max(),bb.max())+.004;ax.axvline(0,color=gray,lw=.7);ax.axhline(0,color=gray,lw=.7);ax.plot([lo,hi],[lo,hi],color='#C7CDD1',ls='--',lw=.8);ax.scatter(a,bb,color=blue,s=27)
# All ten labels repositioned only for layout; no observation removed.
label_positions=[(-.008,.006),(.009,.002),(.002,.007),(-.023,.003),(-.004,.003),(-.002,-.004),(.005,.003),(.001,-.006),(.011,-.020),(.006,-.005)]
for name,aa,bbb,pos in zip(genes,a,bb,label_positions):ax.annotate(name,(aa,bbb),xytext=pos,textcoords='data',fontsize=7,arrowprops={'arrowstyle':'-','color':'#ADB5BD','lw':.6})
c=pd.read_csv(O/'cross_method_concordance.tsv',sep='\t').iloc[0];ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel='ΔR² X-PAIR',ylabel='ΔR² SPPIDER');ax.set_title(f'P3  Cross-method concordance (Spearman ρ = {c.Spearman:.3f}; n = 10 genes)',loc='left',pad=12);save(fig,'FigureP3')
fig,axes=plt.subplots(1,2,figsize=(183/25.4,135/25.4),gridspec_kw={'width_ratios':[1.15,1]});fig.subplots_adjust(left=.14,right=.97,bottom=.20,top=.84,wspace=.46)
cols=['X-PAIR labeled','X-PAIR sorted','Position','DOMAIN_PROXY','Generic interface','Generic 3D'];M=np.full((10,6),np.nan);M[:,0]=x.R2_XPAIR_LABELED;M[:,1]=x.R2_XPAIR_SORTED;N=np.full((10,4),np.nan);counts={}
for i,gene in enumerate(genes):
 for j,name in enumerate(['POSITION','DOMAIN_PROXY','GENERIC_INTERFACE','GENERIC3D']):
  row=b.loc[b.gene.eq(gene)&b.baseline.eq(name)].iloc[0];counts[i,j]=(row.n_complete_case,row.n_primary,row.status)
  if row.status=='EVALUABLE':M[i,j+2]=row.R2_BASELINE;N[i,j]=row.DELTA_R2_BASELINE
cm=plt.colormaps['Blues'].copy();cm.set_bad('#F0F0F0');im=axes[0].imshow(np.ma.masked_invalid(M),aspect='auto',cmap=cm,vmin=0,vmax=np.nanmax(M));axes[0].set_yticks(y,genes);axes[0].set_xticks(np.arange(6),cols,rotation=65,ha='right',rotation_mode='anchor',fontsize=7);axes[0].set_title('a  Marginal R² with availability',loc='left',pad=12)
for i in range(10):
 for j in range(6):
  text='NA' if not np.isfinite(M[i,j]) else f'{M[i,j]:.3f}';axes[0].text(j,i,text,ha='center',va='center',fontsize=6.5,color='white' if np.isfinite(M[i,j]) and M[i,j]>.07 else '#222222')
axes[1].axvline(0,color=gray,lw=.8);symbols=['o','s','^','D'];colors=['#63737E','#B87836','#5999A4','#A596B1']
for j,name in enumerate(cols[2:]):axes[1].scatter(N[:,j],y+(j-1.5)*.13,s=18,marker=symbols[j],color=colors[j],label=name)
axes[1].set_yticks(y,genes);axes[1].invert_yaxis();axes[1].set_xlabel('Row-matched X-PAIR − baseline R²');axes[1].set_title('b  Complete-case paired contrasts',loc='left',pad=12);axes[1].legend(loc='upper center',bbox_to_anchor=(.5,-.18),ncol=2,fontsize=7);fig.colorbar(im,ax=axes[0],fraction=.025,pad=.03,label='R²');save(fig,'FigureP4');b.to_csv(F/'FigureP4_source_data.tsv',sep='\t',index=False)
fig,axes=plt.subplots(1,3,figsize=(183/25.4,110/25.4));fig.subplots_adjust(left=.13,right=.97,bottom=.18,top=.82,wspace=.40)
for j,(col,title,xlabel) in enumerate([('DiscordantResidueRate','a  Same-residue discordance','Discordant residue fraction'),('RESIDUE_INFORMATION_CEILING','b  Frozen information ceiling','Variant-mass ceiling'),('DELTA_R2_PARTNER_ID','c  Partner-identity ΔR²','ΔR² partner identity')]):
 ax=axes[j];ax.scatter(ctx[col],y,s=24,color=blue if j==2 else gray);ax.set_yticks(y,genes if j==0 else ['']*10);ax.invert_yaxis();ax.set_title(title,loc='left',fontsize=8,pad=12);ax.set_xlabel(xlabel,fontsize=8)
 if j<2:ax.set_xlim(-.03,1.05)
 else:ax.axvline(0,color=gray,lw=.7)
fig.suptitle('P5  Descriptive residue-information limits',x=.13,ha='left',y=.96,fontsize=9);save(fig,'FigureP5')
# Alphabetical first eligible gene is a neutral rule fixed before figure rendering.
gene=genes[0];assert gene=='COL3A1';z=np.load(O/'xpair_permutations'/(gene+'.npz'))['delta_R2_perm'];obs=x.loc[gene,'DELTA_R2_PARTNER_ID'];fig,ax=plt.subplots(figsize=(183/25.4,95/25.4));fig.subplots_adjust(left=.12,right=.96,bottom=.19,top=.84);ax.hist(z,bins=45,color='#CAD8E1',edgecolor='white',lw=.3);ax.axvline(0,color=gray,lw=.8,ls='--',label='Zero');ax.axvline(obs,color=blue,lw=1.5,label=f'Observed ΔR² = {obs:.5f}');ax.set(xlabel='Permuted ΔR² partner identity',ylabel='Number of assignments');ax.set_title('P6  COL3A1: neutral alphabetical exemplar, 10,000 paired assignments',loc='left',pad=12,fontsize=9);ax.legend(fontsize=8);save(fig,'FigureP6')
print('ALL_SIX_FIGURES_EXPORTED')
