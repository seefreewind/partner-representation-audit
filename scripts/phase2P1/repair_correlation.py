from pathlib import Path
import numpy as np,pandas as pd,shutil,json
from stats import *
R=Path(__file__).resolve().parents[2];O=R/'results/phase2P1';P=R/'results/phase2P0';A=O/'diagnostic_archive';A.mkdir(exist_ok=True)
old=O/'correlation_distance_sensitivity.tsv'
if not (A/old.name).exists():shutil.copyfile(old,A/old.name)
rows=[]
for gene in pd.read_csv(O/'xpair_primary_per_gene.tsv',sep='\t').gene:
 y=np.load(O/'paired_assignments'/(gene+'.npz'))['observed_labels'];perm=np.load(O/'paired_assignments'/(gene+'.npz'))['label_assignments'];a=np.load(P/'matrices'/(gene+'_xpair_primary.npz'))['scores'];b=np.load(P/'matrices'/(gene+'_xpair_unlabeled_control.npz'))['scores'];G,info=correlation_gram(a);H,other=correlation_gram(b);row=dict(gene=gene,**info,sorted_status=other['status'],sorted_Lingoes_constant=other.get('Lingoes_constant',np.nan));path=O/'sensitivity_permutations'/(gene+'_correlation.npz')
 if path.exists() and not (A/path.name).exists():shutil.copyfile(path,A/path.name)
 if G is not None and H is not None and info['status']=='EVALUABLE_LINGOES' and other['status']=='EVALUABLE_LINGOES':
  r=float(r2(G,y)[0]);rr=float(r2(H,y)[0]);nn=r2(G,perm);ns=r2(H,perm);row.update(R2_labeled=r,R2_sorted=rr,DELTA_R2=r-rr,p_labeled=empirical(F(nn,len(y),2),F(r,len(y),2)),p_DELTA=empirical(nn-ns,r-rr));np.savez_compressed(path,R2_labeled_perm=nn,R2_sorted_perm=ns,delta_R2_perm=nn-ns)
 else:
  if G is not None and info['status']=='EVALUABLE_LINGOES':row['R2_labeled']=float(r2(G,y)[0])
  row['status']='UNEVALUABLE_'+('ZERO_TOTAL_SORTED_DISTANCE' if other['status']=='DEGENERATE_ZERO_TOTAL_DISTANCE' else 'DIMENSION_OR_ROW_GATE')
  # Preserve old numerical artifacts in archive; explicitly mark present active arrays unusable.
  if path.exists():np.savez_compressed(path,R2_labeled_perm=np.full(len(perm),np.nan),R2_sorted_perm=np.full(len(perm),np.nan),delta_R2_perm=np.full(len(perm),np.nan),valid=np.array(False))
 rows.append(row)
pd.DataFrame(rows).to_csv(old,sep='\t',index=False,na_rep='NA')
(O/'correlation_implementation_correction.json').write_text(json.dumps(dict(reason='Mathematically identical rank vectors have exact zero distance. Floating point corrcoef introduced spurious tiny distances and a NaN delta was erroneously counted as an extreme tail. Exact identical rank pairs now have D=0, and preregistered zero-total-distance rule makes such controls unevaluable. No primary statistic, cohort, scaling, permutation or endpoint changed.',original_outputs_preserved='diagnostic_archive',new_analysis_choice=False),indent=2)+'\n')
print('CORRELATION_NUMERICAL_DEGENERACY_HANDLED')
