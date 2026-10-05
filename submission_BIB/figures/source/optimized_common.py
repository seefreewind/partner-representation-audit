from pathlib import Path
import json, re, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'submission_BIB/figures'
QA_PREVIEW_DPI=600/4 # QA only; all deliverable PNGs use 600 dpi
S=ROOT/'results/phase2S_full'; P=ROOT/'results/phase2P1'
plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.titlesize':9.5,'axes.labelsize':8.5,'xtick.labelsize':8.5,'ytick.labelsize':8.5,'legend.fontsize':8.5,'legend.frameon':False,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'figure.facecolor':'white','axes.facecolor':'white'})
COL=['#444444','#777777','#AAAAAA','#B57900','#0072B2']; METHODS=['MUTPRED_PPI_DISRUPTION','STRUCTURAL_BINARY_INTERFACE','STRUCTURAL_INTERFACE_DISTANCE','SPPIDER_SEQ_PARTNER_INTERFACE','XPAIR_PARTNER_INTERFACE']; NAMES=['MutPred-PPI','Structural binary','Structural distance','SPPIDER-seq','X-PAIR']
summary=json.loads((P/'reporting_summary.json').read_text()); cohort=json.loads((ROOT/'results/phase2P0/cohort_summary.json').read_text())
coverage=pd.read_csv(S/'method_comparison_summary.tsv',sep='\t').query("universe == 'U2'").set_index('method').loc[METHODS]
xp=pd.read_csv(P/'xpair_primary_per_gene.tsv',sep='\t').set_index('gene').sort_index(); sp=pd.read_csv(P/'sppider_confirmatory_per_gene.tsv',sep='\t').set_index('gene').reindex(xp.index)
# Historical numbers occur in frozen narrative reports and figure legends, not a structured summary.
# Extract them from the existing, audited caption instead of repeating scientific literals.
oldleg=json.loads((ROOT/'submission_BIB/source_data/figure_legends.json').read_text())
hist=[int(x) for x in re.search(r'attrition is ([\d→]+) genes',oldleg['1'][1]).group(1).split('→')]
scored,total,scorepct=re.search(r'scored ([\d,]+)/([\d,]+) submitted rows \(([\d.]+)%\)',oldleg['1'][1]).groups(); scored=int(scored.replace(',',''));total=int(total.replace(',',''));scorepct=float(scorepct)
tracks=[float(x.rstrip('.')) for x in re.search(r'summaries, ([\d.]+) and ([\d.]+)',oldleg['3'][1]).groups()]

def save(fig,number,supp=False):
    dest=OUT/('supplementary' if supp else 'main'); stem=f'Figure{number}'
    fig.savefig(dest/f'{stem}.pdf',format='pdf')
    fig.savefig(dest/f'{stem}.svg',format='svg')
    fig.savefig(dest/f'{stem}.png',dpi=600)
    fig.savefig(ROOT/'submission_BIB/qa_figures_optimized'/f'{stem}_preview.png',dpi=QA_PREVIEW_DPI)
    plt.close(fig)
def panel(ax,label,title):
    ax.set_title(title,loc='left',pad=13)
    ax.text(-.14,1.08,label,transform=ax.transAxes,fontweight='bold',fontsize=11)
def zero(ax):ax.axvline(0,color='#777777',lw=.8,ls='--')
def text(ax,x,y,s,size=9,**kw):ax.text(x,y,s,fontsize=size,**kw)
def arrow(ax,a,b):ax.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'->','color':'#777777','lw':1},annotation_clip=False)

def figure1():
    fig,axes=plt.subplots(3,1,figsize=(170/25.4,155/25.4));fig.subplots_adjust(left=.05,right=.97,top=.94,bottom=.12,hspace=.37)
    for ax in axes:ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
    a,b,c=axes
    text(a,0,1.04,'a   STAGE 1  Technical scoreability',11,fontweight='bold')
    for x,n,l in zip([.02,.36,.72],hist[:3],['strict genes','not-known-incompatible','structure-ready']):text(a,x,.62,str(n),17,fontweight='bold');text(a,x,.39,l)
    arrow(a,(.18,.64),(.32,.64));arrow(a,(.6,.64),(.68,.64))
    text(a,.02,.02,f'{scorepct:.2f}% row scoreability',13,fontweight='bold');text(a,.57,.02,f'≠  {hist[-1]} balanced genes',13,fontweight='bold')
    text(a,.02,-.2,f'MutPred-PPI: {scored:,}/{total:,} submitted rows scored',8.5)
    text(b,0,1.04,'b   STAGE 2  Phenotype-balanced representation',11,fontweight='bold')
    for y,i in zip([.68,.36,.04],[0,3,4]):
        r=coverage.iloc[i];text(b,.02,y,NAMES[i],10,color=COL[i]);text(b,.38,y,f'{int(r.genes_C4_ge70)} with C4 ≥70%',10);arrow(b,(.66,y+.035),(.74,y+.035));text(b,.78,y,f'{int(r.balanced_genes)} balanced',10,fontweight='bold')
    text(c,0,1.04,'c   STAGE 3  Incremental phenotype information',11,fontweight='bold')
    text(c,.02,.7,f'{cohort["shared_genes"]} shared → {summary["primary_genes"]} phenotype-ready → {summary["primary_residues"]} residues',11)
    text(c,.02,.35,f'Median ΔR² = {summary["median_delta_R2"]:.4f}\n95% CI includes zero',11,va='center')
    text(c,.66,.23,'GLOBAL FAIL',17,fontweight='bold')
    fig.text(.05,.045,'Evaluation order, not causal flow. Failed gates remain frozen.',fontsize=9)
    save(fig,1)

def figure2():
    fig,ax=plt.subplots(figsize=(170/25.4,94/25.4));fig.subplots_adjust(left=.3,right=.93,top=.9,bottom=.25)
    ys=np.array([0,1,2,3.7,4.7]);h=.28
    for y,(i,r) in zip(ys,enumerate(coverage.itertuples())):
        for yy,n,hatch in [(y-h/2,r.genes_C4_ge70,None),(y+h/2,r.balanced_genes,'////')]:
            ax.barh(yy,n,height=h,color=COL[i],hatch=hatch,edgecolor='white' if hatch is None else '#222222',linewidth=.4)
            ax.text(n+1.7,yy,str(int(n)),va='center')
    ax.set_yticks(ys,NAMES);ax.invert_yaxis();ax.set_xlim(0,140);ax.set_xlabel('Genes among the common 132-gene denominator')
    ax.axhline(2.8,color='#DDDDDD',lw=.7);ax.text(.02,1.07,'Structure-dependent / structural',transform=ax.transAxes,fontweight='bold');ax.text(.02,.36,'Sequence-based',transform=ax.transAxes,fontweight='bold')
    fig.legend(handles=[Patch(facecolor='#888888',label='C4 ≥70%'),Patch(facecolor='#888888',edgecolor='#222222',hatch='////',label='Phenotype-balanced')],loc='lower center',bbox_to_anchor=(.53,.08),ncol=2)
    fig.text(.04,.025,'Representation availability; score semantics differ across methods.',fontsize=8.5)
    save(fig,2)

def figure3():
    fig=plt.figure(figsize=(170/25.4,133/25.4));gs=fig.add_gridspec(2,2,height_ratios=[1,1],left=.18,right=.95,top=.86,bottom=.14,hspace=.8,wspace=.55)
    a=fig.add_subplot(gs[0,0]);b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[1,:]);panel(a,'a','Internal structural validity');panel(b,'b','Partner-track similarity');panel(c,'c','Input-contract availability')
    v=pd.read_csv(S/'Table_S7_structural_validation.tsv',sep='\t');v=v[v.aggregation.eq('macro')]
    for i,m in enumerate(['sppider','xpair']):
        for j,metric in enumerate(['AUROC','AUPRC']):
            r=v[v.method.eq(m)&v.metric.eq(metric)].iloc[0];y=j+(.1 if i==0 else -.1)
            a.errorbar(r.estimate,y,xerr=[[r.estimate-r.CI95_low],[r.CI95_high-r.estimate]],fmt=['o','s'][i],color=COL[i+3],capsize=3,label=NAMES[i+3] if j==0 else None)
    a.set_yticks([0,1],['AUROC','AUPRC']);a.set_xlim(0,1);a.set_ylim(1.5,-.5);a.set_xlabel('Estimate and 95% bootstrap CI');fig.legend(*a.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.55,.995),ncol=2)
    for i,n in enumerate(tracks):b.scatter(n,i,s=45,marker=['o','s'][i],color=COL[i+3]);b.text(n-.02,i+.26,f'{n:.3f}',ha='center')
    b.set_yticks([0,1],NAMES[3:]);b.set_ylim(1.6,-.6);b.set_xlim(0,1);b.set_xlabel('Median within-gene Spearman')
    res=pd.read_csv(S/'Table_S9_computational_resources.tsv',sep='\t').set_index('method')
    length=pd.read_csv(S/'xpair_length_constraint_audit.tsv',sep='\t')
    # Pair status provides deterministic counts; no biological or statistical calculation.
    statuses=pd.read_csv(S/'pair_status_xpair.tsv',sep='\t');lc=int(statuses.status.eq('MODEL_LENGTH_UNSUPPORTED').sum())
    for i,m in enumerate(['sppider','xpair']):
        r=res.loc[m];n=int(r.pairs_scored);tot=int(r.frozen_pairs_total);other=tot-n-(lc if i else 0)
        c.barh(i,n,color=COL[i+3],height=.48);c.text(n/2,i,f'{n:,} scored',ha='center',va='center',color='white')
        start=n
        if i:c.barh(i,lc,left=start,color='#CCCCCC',edgecolor='#444444',hatch='////',height=.48);c.text(start+lc/2,i,f'{lc:,}',ha='center',va='center');start+=lc
        c.barh(i,other,left=start,color='#F0F0F0',edgecolor='#777777',height=.48);c.text(tot+180,i,f'{other:,} other',va='center')
    c.set_yticks([0,1],NAMES[3:]);c.set_ylim(1.65,-.7);c.set_xlim(0,18400);c.set_xticks([0,8000,15942],['0','8,000','15,942']);c.set_xlabel('Pairs in the frozen plan (15,942 per model)')
    fig.text(.12,.018,'Hatched segment: X-PAIR length unsupported. Higher track correlation means less variation.',fontsize=8.5)
    save(fig,3)

def figure4():
    fig=plt.figure(figsize=(170/25.4,163/25.4));gs=fig.add_gridspec(2,2,height_ratios=[.8,1.4],left=.15,right=.95,top=.92,bottom=.16,hspace=.46,wspace=.5)
    a=fig.add_subplot(gs[0,:]);a.axis('off');a.set_xlim(0,1);a.set_ylim(0,1);panel(a,'a','Phenotype cohort attrition')
    text(a,.5,.88,f'{cohort["shared_genes"]} shared balanced genes → {cohort["total_unique_residues"]:,} unique pathogenic residues',10,ha='center')
    branches=[('single-disease_residues','single disease'),('cross_disease_residues','cross disease'),('multilabel_ambiguous_residues','multilabel / ambiguous')]
    for x,(k,l) in zip([.18,.5,.82],branches):
        k=k.replace('single-disease','single_disease');arrow(a,(.5,.75),(x,.59));text(a,x,.43,f'{cohort[k]}\n{l}',10,ha='center',va='center')
    arrow(a,(.18,.29),(.18,.1));text(a,.18,-.05,'Unique-residue group gate',9,ha='center');text(a,.62,-.05,f'{summary["primary_genes"]} genes · {summary["primary_residues"]} primary residues\nMedian common partners: {summary["median_common_partners"]:g}',10,ha='center');arrow(a,(.36,-.03),(.4,-.03))
    b=fig.add_subplot(gs[1,0]);c=fig.add_subplot(gs[1,1]);panel(b,'b','Primary single-disease cohort');panel(c,'c','Broader annotated registry:\ncross-disease-linked residues')
    y=np.arange(len(xp));b.barh(y,xp.n_group1,color='#555555',label='Group 1');b.barh(y,xp.n_group2,left=xp.n_group1,color='#D6D6D6',edgecolor='#777777',label='Group 2')
    for i,r in enumerate(xp.itertuples()):
        b.text(r.n_group1/2,i,str(r.n_group1),ha='center',va='center',color='white');b.text(r.n_group1+r.n_group2/2,i,str(r.n_group2),ha='center',va='center')
    b.set_yticks(y,xp.index);b.invert_yaxis();b.set_xlim(0,135);b.set_xlabel('Independent residues');b.legend(loc='upper center',bbox_to_anchor=(.5,-.15),ncol=2)
    ctx=pd.read_csv(P/'residue_ceiling_context.tsv',sep='\t').set_index('gene').reindex(xp.index);c.scatter(ctx.DiscordantResiduePct,y,color='#444444',s=24)
    for i,r in enumerate(ctx.itertuples()):c.text(r.DiscordantResiduePct+2,i,f'{r.DiscordantResiduePct:.1f}%',va='center')
    c.set_yticks(y,xp.index);c.invert_yaxis();c.set_xlim(-2,110);c.set_xlabel('Annotated residues linked to ≥2 groups (%)')
    fig.text(.05,.025,'WT residue × partner vectors tie substitutions at the same position: a representation ceiling.',fontsize=8.5)
    save(fig,4)

def figure5():
    fig=plt.figure(figsize=(170/25.4,175/25.4));gs=fig.add_gridspec(2,2,left=.15,right=.95,top=.93,bottom=.11,hspace=.51,wspace=.55,height_ratios=[1.1,1])
    a=fig.add_subplot(gs[0,0]);b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[1,0]);d=fig.add_subplot(gs[1,1])
    panel(a,'a','All primary genes');panel(b,'b','Labeled versus sorted');panel(c,'c','Preregistered global result');panel(d,'d','Cross-model concordance')
    y=np.arange(len(xp));a.scatter(xp.DELTA_R2_PARTNER_ID,y,s=26,color=COL[4]);zero(a);a.set_yticks(y,xp.index);a.invert_yaxis();a.set_xlim(min(-.005,float(xp.DELTA_R2_PARTNER_ID.min())*1.25),max(.005,float(xp.DELTA_R2_PARTNER_ID.max())*1.25));a.set_xticks([-.02,0,.01]);a.set_xlabel('ΔR² (labeled − sorted)')
    for i,r in enumerate(xp.itertuples()):b.plot([r.R2_XPAIR_SORTED,r.R2_XPAIR_LABELED],[i,i],color='#BBBBBB',lw=.8)
    b.scatter(xp.R2_XPAIR_LABELED,y,color=COL[4],s=32,marker='o',label='Labeled');b.scatter(xp.R2_XPAIR_SORTED,y,color='#555555',s=32,marker='s',label='Sorted');b.set_yticks(y,xp.index);b.invert_yaxis();b.set_xlim(-.003,max(float(xp.R2_XPAIR_LABELED.max()),float(xp.R2_XPAIR_SORTED.max()))*1.15);b.set_xticks([0,.04,.08]);b.set_xlabel('R²');b.legend(loc='upper center',bbox_to_anchor=(.5,-.2),ncol=2)
    c.set_xlim(-.007,.007);c.set_ylim(0,1);c.set_yticks([]);c.axvline(0,color='#777777',ls='--',lw=.8,ymin=.5)
    med=summary['median_delta_R2'];lo=summary['bootstrap_CI_low'];hi=summary['bootstrap_CI_high'];c.errorbar(med,.81,xerr=[[med-lo],[hi-med]],fmt='o',color=COL[4],capsize=4)
    c.text(.5,.63,f'Median ΔR² = {med:.4f}\n95% CI {lo:.4f} to {hi:.4f}\n{summary["positive_gene_n"]}/{summary["total_gene_n"]} positive',transform=c.transAxes,ha='center',va='top',fontsize=10)
    c.text(.5,.25,'Condition A: FAIL\nCondition B: PASS',transform=c.transAxes,ha='center',fontsize=10)
    c.text(.5,.03,'GLOBAL FAIL',transform=c.transAxes,ha='center',fontsize=17,fontweight='bold');c.set_xlabel('Global median ΔR²');c.set_xticks([-.005,0,.005])
    d.scatter(xp.DELTA_R2_PARTNER_ID,sp.DELTA_R2_SPPIDER,s=30,color='#555555');zero(d);d.axhline(0,color='#777777',ls='--',lw=.8);d.set_xlabel('X-PAIR ΔR²');d.set_ylabel('SPPIDER ΔR²');d.yaxis.set_major_formatter(__import__('matplotlib').ticker.FormatStrFormatter('%.3f'));d.set_xlim(float(xp.DELTA_R2_PARTNER_ID.min())*1.3,float(xp.DELTA_R2_PARTNER_ID.max())*1.3);d.set_ylim(float(sp.DELTA_R2_SPPIDER.min())*1.4,float(sp.DELTA_R2_SPPIDER.max())*1.6)
    for g in ['RS1','COL3A1']:d.annotate(g,(xp.loc[g,'DELTA_R2_PARTNER_ID'],sp.loc[g,'DELTA_R2_SPPIDER']),xytext=(-25,8) if g=='RS1' else (3,5),textcoords='offset points',fontsize=8.5)
    d.text(.03,-.31,f'Spearman = {summary["cross_method_Spearman"]:.3f}\n{summary["positive_both"]}/10 positive in both',transform=d.transAxes,fontsize=8.5)
    save(fig,5)

def supplementary(n):
    if n==1:
        fig,ax=plt.subplots(figsize=(170/25.4,85/25.4));fig.subplots_adjust(left=.14,bottom=.22,top=.9,right=.96);arr=np.load(P/'xpair_permutations/COL3A1.npz')['delta_R2_perm'];ax.hist(arr,bins=45,color='#CCCCCC',edgecolor='white');ax.axvline(xp.loc['COL3A1','DELTA_R2_PARTNER_ID'],color=COL[4],label='Observed',lw=1.5);zero(ax);ax.set_xlabel('Paired permutation ΔR² — COL3A1');ax.set_ylabel('Alternative assignments');ax.legend()
    elif n==2:
        r=pd.read_csv(P/'restricted_permutation_results.tsv',sep='\t').set_index('gene').reindex(xp.index);fig,ax=plt.subplots(figsize=(170/25.4,104/25.4));fig.subplots_adjust(left=.17,right=.95,bottom=.2,top=.92)
        for i,row in enumerate(r.itertuples()):ax.plot([row.unrestricted_p_DELTA,row.restricted_p_DELTA],[i,i],color='#BBBBBB',lw=.8)
        ax.scatter(r.unrestricted_p_DELTA,range(10),marker='o',color=COL[4],label='Unrestricted');ax.scatter(r.restricted_p_DELTA,range(10),marker='s',color='#555555',label='Restricted');ax.set_yticks(range(10),xp.index);ax.invert_yaxis();ax.axvline(.05,color='#777777',ls='--');ax.set_xlim(0,1);ax.set_xlabel('Paired-contrast permutation P');ax.legend(loc='upper center',bbox_to_anchor=(.5,-.15),ncol=2)
    elif n==3:
        r=pd.read_csv(S/'Table_S8_length_constraints.tsv',sep='\t');t=pd.read_csv(S/'sppider_long_sequence_audit.tsv',sep='\t');fig,axs=plt.subplots(1,2,figsize=(170/25.4,85/25.4));fig.subplots_adjust(left=.1,right=.96,bottom=.24,top=.83,wspace=.42)
        for ax,data,key,col,title,l in [(axs[0],r,'length_supported_pairs',COL[4],'X-PAIR contract','a'),(axs[1],t,'sppider_scored_partners',COL[3],'SPPIDER windowing','b')]:
            ax.scatter(data.query_length,data[key],s=13,color=col,marker='s' if l=='a' else 'o');ax.axvline(2000,color='#777777',ls='--',lw=.8);ax.set_xlabel('Query length (aa)');ax.set_ylabel('Partner pairs');panel(ax,l,title)
    elif n==4:
        r=pd.read_csv(S/'coverage_by_gene_all_methods.tsv',sep='\t').pivot(index='gene',columns='method',values='C4_fraction').reindex(columns=METHODS).sort_index();fig,axs=plt.subplots(2,2,figsize=(170/25.4,250/25.4));fig.subplots_adjust(left=.12,right=.97,top=.86,bottom=.06,hspace=.24,wspace=.35)
        from matplotlib.colors import ListedColormap
        cmap=ListedColormap(['#F0F0F0','#DADADA','#BBBBBB','#888888','#555555','#222222']);cmap.set_bad('white')
        for i,ax in enumerate(axs.flat):
            sub=r.iloc[i*33:(i+1)*33];im=ax.imshow(sub,vmin=0,vmax=1,cmap=cmap,aspect='auto');ax.set_yticks(range(len(sub)),sub.index);ax.set_xticks(range(5),['MP','Binary','Dist.','SPP','XP'],rotation=45,ha='right',rotation_mode='anchor');panel(ax,chr(97+i),f'Genes {i*33+1}–{(i+1)*33}')
        ca=fig.add_axes([.3,.96,.5,.013]);cb=fig.colorbar(im,cax=ca,orientation='horizontal',label='C4 fraction; white also denotes unavailable');cb.ax.xaxis.set_label_position('top')
    elif n==5:
        r=pd.read_csv(S/'shuffled_partner_control.tsv',sep='\t');fig,axs=plt.subplots(1,2,figsize=(170/25.4,87/25.4));fig.subplots_adjust(left=.11,right=.96,bottom=.24,top=.83,wspace=.42)
        for ax,k,l,tit in zip(axs,['Spearman','mean_absolute_track_change'],['a','b'],['True versus replacement\nsimilarity','Absolute track change']):
            for i,m in enumerate(['sppider','xpair']):
                v=r[r.method.eq(m)];ax.scatter(range(1,len(v)+1),v[k],color=COL[3+i],marker=['o','s'][i],s=13,label=NAMES[3+i])
            panel(ax,l,tit);ax.set_xlabel('Frozen control order');ax.set_ylabel('Spearman' if k=='Spearman' else 'Mean absolute difference')
        axs[0].legend(loc='upper center',bbox_to_anchor=(1.1,-.2),ncol=2)
    elif n==6:
        r=pd.read_csv(S/'Table_S9_computational_resources.tsv',sep='\t').set_index('method');fig,axs=plt.subplots(1,3,figsize=(170/25.4,86/25.4));fig.subplots_adjust(left=.08,right=.98,bottom=.23,top=.83,wspace=.5)
        for ax,k,l,tit,scale in zip(axs,['CPU_hours','peak_RAM_bytes','embedding_cache_bytes'],['a','b','c'],['Recorded CPU hours','Peak RAM (GB)','Embedding cache (GB)'],[1,1e9,1e9]):
            for i,m in enumerate(['sppider','xpair']):
                v=r.loc[m,k]/scale;ax.scatter(i,v,color=COL[i+3],marker=['o','s'][i],s=35);ax.text(i,v+.035*max(r[k]/scale),f'{v:.2f}',ha='center')
            ax.set_xticks([0,1],['SPPIDER','X-PAIR'],rotation=25,ha='right',rotation_mode='anchor');ax.set_ylim(0,max(r[k]/scale)*1.25);ax.set_xlim(-.5,1.5);panel(ax,l,tit)
    elif n==7:
        r=pd.read_csv(P/'baseline_summary.tsv',sep='\t').set_index('baseline').loc[['POSITION','GENERIC_INTERFACE','GENERIC3D','DOMAIN_PROXY']];fig,ax=plt.subplots(figsize=(170/25.4,85/25.4));fig.subplots_adjust(left=.35,right=.96,bottom=.25,top=.9);ax.scatter(r.median_XPAIR_increment,range(4),s=42,color='#555555');zero(ax);ax.set_yticks(range(4),[f'{i}  n={int(row.n_genes)}' for i,row in r.iterrows()]);ax.invert_yaxis();ax.set_xlabel('Paired X-PAIR − baseline R² (existing median)');fig.text(.05,.04,'Different complete-case gene sets prevent global ranking.',fontsize=10,fontweight='bold')
    save(fig,'S'+str(n),True)
