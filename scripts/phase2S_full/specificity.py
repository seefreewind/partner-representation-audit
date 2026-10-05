from common import *
import numpy as np,itertools
from scipy.stats import rankdata
rows=[];pairrows=[]
for method in ['sppider','xpair']:
    d=rd('results/phase2S_full/pair_status_'+method+'.tsv');d=d.loc[d.status.eq('SCORED')]
    for gene,group in d.groupby('gene'):
        if group.partner.nunique()<3:continue
        group=group.sort_values('partner');tracks=[]
        for r in group.itertuples(index=False):
            z=np.load(R/r.score_array_path);tracks.append(z['receptor'] if method=='sppider' else z['query'])
        a=np.stack(tracks).astype(float);vary=np.var(a,axis=0);pear=np.corrcoef(a);spearman=np.corrcoef(rankdata(a,axis=1));indices=np.triu_indices(len(a),1)
        for i,j in zip(*indices):pairrows.append(dict(method=method,gene=gene,partner1=group.iloc[i].partner,partner2=group.iloc[j].partner,Pearson=pear[i,j],Spearman=spearman[i,j],constant_track=bool(np.ptp(a[i])==0 or np.ptp(a[j])==0)))
        rows.append(dict(method=method,gene=gene,scored_partners=len(a),residue_count=a.shape[1],pairwise_Pearson_median=np.nanmedian(pear[indices]),pairwise_Spearman_median=np.nanmedian(spearman[indices]),pairwise_Pearson_mean=np.nanmean(pear[indices]),pairwise_Spearman_mean=np.nanmean(spearman[indices]),mean_across_partner_residue_variance=vary.mean(),meaningful_dispersion_fraction=(np.sqrt(vary)>=.05).mean(),dispersion_rule='population_residue_SD_ge0.05'))
        wr(pd.DataFrame(rows),'partner_specificity_metrics.tsv')
wr(pd.DataFrame(pairrows),'partner_specificity_pairwise.tsv');print('PARTNER_SPECIFICITY_COMPLETE',flush=True)
