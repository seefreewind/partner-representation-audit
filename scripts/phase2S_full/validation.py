from common import *
import numpy as np
from sklearn.metrics import roc_auc_score,average_precision_score,matthews_corrcoef,confusion_matrix
from collections import defaultdict
labels=rd('results/phase2S_full/experimental_validation_label_manifest.tsv');rng=np.random.default_rng(20261003)
def metrics(y,p,weights=None):
    pred=p>=.5;one=len(np.unique(y))<2;tn,fp,fn,tp=confusion_matrix(y,pred,labels=[0,1],sample_weight=weights).ravel()
    den=np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return {'AUROC':np.nan if one else roc_auc_score(y,p,sample_weight=weights),'AUPRC':np.nan if one else average_precision_score(y,p,sample_weight=weights),'MCC':(tp*tn-fp*fn)/den if den else np.nan,'sensitivity':tp/(tp+fn) if tp+fn else np.nan,'specificity':tn/(tn+fp) if tn+fp else np.nan}
summary=[];boot=[];groupmetrics=[]
for method in ['sppider','xpair']:
    pairs=rd('results/phase2S_full/pair_status_'+method+'.tsv').set_index(['gene','partner']);groups=defaultdict(list);paired=[]
    # Multiple experimental complexes for one residue/pair are merged by the
    # union contact label, preventing replicate observations from duplicating units.
    bypair=defaultdict(dict)
    for r in labels.loc[labels.status.eq('EXPERIMENTAL_LABELS')].itertuples(index=False):
        z=np.load(R/r.label_path)
        for pos,y in zip(z['positions'],z['label']):bypair[(r.gene,r.partner)][int(pos)]=max(bypair[(r.gene,r.partner)].get(int(pos),0),int(y))
    for key,label in sorted(bypair.items()):
        if key not in pairs.index or pairs.loc[key,'status']!='SCORED':continue
        r=pairs.loc[key];z=np.load(R/r.score_array_path);arr=z['receptor'] if method=='sppider' else z['query'];positions=sorted(label);y=np.array([label[x] for x in positions]);scores=arr[np.array(positions)-1]
        groups[key[0]].append((y,scores));paired.append({'method':method,'gene':key[0],'partner':key[1],'residues':len(y),'positive_residues':int(y.sum()),**metrics(y,scores),'training_overlap':'UNKNOWN_VALIDATION_NOT_HELD_OUT'})
    data={g:(np.concatenate([a for a,b in x]),np.concatenate([b for a,b in x])) for g,x in groups.items()};gs=sorted(data)
    if not gs:summary.append(dict(method=method,status='NO_EXPERIMENTAL_EVALUABLE_RESIDUES'));continue
    gm={g:metrics(*data[g]) for g in gs};groupmetrics.extend(dict(method=method,gene=g,residue_n=len(data[g][0]),**gm[g],status='TWO_CLASS' if len(np.unique(data[g][0]))==2 else 'ONE_CLASS_GROUP') for g in gs)
    y=np.concatenate([data[g][0] for g in gs]);scores=np.concatenate([data[g][1] for g in gs]);idx=np.concatenate([np.full(len(data[g][0]),i) for i,g in enumerate(gs)])
    pointmicro=metrics(y,scores);pointmacro={k:np.nanmean([gm[g][k] for g in gs]) for k in pointmicro}
    for iteration in range(1000):
        sample=rng.integers(0,len(gs),len(gs));weights=np.bincount(sample,minlength=len(gs));keep=weights[idx]>0;micro=metrics(y[keep],scores[keep],weights[idx][keep]);macro={k:np.nanmean([gm[gs[i]][k] for i in sample]) for k in pointmicro}
        for mode,values in [('micro',micro),('macro',macro)]:
            for metric,value in values.items():boot.append(dict(method=method,bootstrap=iteration,aggregation=mode,metric=metric,value=value,group_unit='query_protein',group_n=len(gs)))
    for mode,values in [('micro',pointmicro),('macro',pointmacro)]:
        for k,value in values.items():
            vals=np.array([r['value'] for r in boot if r['method']==method and r['aggregation']==mode and r['metric']==k]);vals=vals[np.isfinite(vals)]
            lo,hi=np.quantile(vals,[.025,.975]) if len(vals) else (np.nan,np.nan)
            summary.append(dict(method=method,aggregation=mode,metric=k,estimate=value,CI95_low=lo,CI95_high=hi,bootstrap_n=1000,valid_bootstrap_n=len(vals),protein_groups=len(gs),residue_pair_n=len(y),status='PHENOTYPE_BLIND_INTERNAL_VALIDATION_TRAINING_OVERLAP_UNKNOWN'))
wr(pd.DataFrame(summary),'structural_validation_metrics.tsv');wr(pd.DataFrame(boot),'structural_validation_bootstrap.tsv');wr(pd.DataFrame(groupmetrics),'structural_validation_by_protein.tsv');wr(pd.DataFrame(paired),'structural_validation_by_pair.tsv');print('STRUCTURAL_VALIDATION_COMPLETE',flush=True)
