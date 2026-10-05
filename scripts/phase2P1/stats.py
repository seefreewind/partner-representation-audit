import numpy as np, math,itertools,hashlib
from scipy.stats import rankdata,binomtest

def standardize(x):
 x=np.asarray(x,dtype=float);mu=x.mean(axis=0);sd=x.std(axis=0,ddof=0);z=np.zeros_like(x);ok=sd>0;z[:,ok]=(x[:,ok]-mu[ok])/sd[ok];return z,mu,sd

def gram(x):
 x=np.asarray(x,dtype=float);c=x-x.mean(axis=0);return c@c.T

def r2(G,labels):
 Y=np.atleast_2d(labels);T=float(np.trace(G));out=np.zeros(len(Y))
 if T<=0:return np.full(len(Y),np.nan)
 for k in np.unique(Y[0]):
  W=(Y==k).astype(float);sizes=W.sum(axis=1);out+=np.einsum('ij,ij->i',W@G,W)/sizes
 assert np.isfinite(out).all(),'NONFINITE_PERMANOVA_NULL';return out/T

def F(r,n,k):return (r/(k-1))/((1-r)/(n-k))
def empirical(null,obs):return (1+int(np.count_nonzero(null>=obs)))/(len(null)+1)
def bh(p):
 p=np.asarray(p,float);out=np.full(len(p),np.nan);valid=np.isfinite(p);a=p[valid];order=np.argsort(a);s=a[order];v=np.minimum(1,np.minimum.accumulate((s*len(a)/np.arange(1,len(a)+1))[::-1])[::-1]);vals=np.empty(len(a));vals[order]=v;out[valid]=vals;return out

def multiset(y):
 y=np.asarray(y,dtype=np.int16);classes,counts=np.unique(y,return_counts=True);indices=list(range(len(y)))
 def recurse(level,remaining,a):
  if level==len(classes)-1:
   b=a.copy();b[remaining]=classes[level];yield b;return
  for chosen in itertools.combinations(remaining,int(counts[level])):
   b=a.copy();b[list(chosen)]=classes[level];chosen=set(chosen);yield from recurse(level+1,[i for i in remaining if i not in chosen],b)
 yield from recurse(0,indices,np.zeros(len(y),dtype=np.int16))

def permutations(y,seed,B=10000,strata=None):
 y=np.asarray(y,dtype=np.int16);sets=[np.arange(len(y))] if strata is None else [np.flatnonzero(np.asarray(strata)==s) for s in sorted(set(strata))];space=1
 for ix in sets:
  sizes=np.unique(y[ix],return_counts=True)[1];space*=math.factorial(len(ix))//math.prod(math.factorial(int(n)) for n in sizes)
 target=min(B,space-1);obs=y.tobytes();seen={obs};result=[]
 if space-1<=B:
  lists=[list(multiset(y[ix])) for ix in sets]
  for blocks in itertools.product(*lists):
   row=y.copy()
   for ix,block in zip(sets,blocks):row[ix]=block
   if row.tobytes()!=obs:result.append(row)
 else:
  rng=np.random.Generator(np.random.PCG64(seed))
  while len(result)<target:
   row=y.copy()
   for ix in sets:row[ix]=rng.permutation(y[ix])
   key=row.tobytes()
   if key not in seen:seen.add(key);result.append(row)
 out=np.array(result,dtype=np.int16);assert len(out)==target and len({row.tobytes() for row in out})==target;assert all(np.array_equal(np.unique(row,return_counts=True)[1],np.unique(y,return_counts=True)[1]) for row in out)
 return out,space

def analyze(a,b,y,perm):
 G=gram(a);H=gram(b);obs=float(r2(G,y)[0]);sort=float(r2(H,y)[0]);null=r2(G,perm);ns=r2(H,perm);delta=obs-sort;dn=null-ns
 return dict(labeled=obs,sorted=sort,delta=delta,p_labeled=empirical(F(null,len(y),len(np.unique(y))),F(obs,len(y),len(np.unique(y)))),p_delta=empirical(dn,delta),pseudoF=F(obs,len(y),len(np.unique(y)))),dict(R2_labeled_perm=null,R2_sorted_perm=ns,delta_R2_perm=dn)

def global_summary(delta):
 d=np.asarray(delta);rng=np.random.Generator(np.random.PCG64(20261005));meds=np.median(d[rng.integers(0,len(d),size=(10000,len(d)))],axis=1);low,high=np.quantile(meds,[.025,.975]);npos=int((d>0).sum());nties=int((d==0).sum());m=float(np.median(d));A=m>0 and low>0;B=npos/len(d)>=.7;return dict(median_delta_R2=m,Q1=np.quantile(d,.25),Q3=np.quantile(d,.75),minimum=d.min(),maximum=d.max(),bootstrap_CI_low=low,bootstrap_CI_high=high,positive_gene_n=npos,total_gene_n=len(d),positive_gene_fraction=npos/len(d),zero_gene_n=nties,sign_test_P=binomtest(npos,len(d)-nties,.5,alternative='greater').pvalue if len(d)>nties else 1.,condition_A_pass=bool(A),condition_B_pass=bool(B),global_support_rule_pass=bool(A and B and len(d)==10))

def correlation_gram(x):
 sd=x.std(axis=0);cols=sd>0
 if cols.sum()<10 or (np.ptp(x[:,cols],axis=1)==0).any():return None,dict(status='NONCONSTANT_DIMENSION_OR_ROW_GATE_FAIL',usable_dimensions=int(cols.sum()))
 ranks=rankdata(x[:,cols],axis=1);D=1-np.corrcoef(ranks);D=np.maximum(D,0);D[np.all(ranks[:,None,:]==ranks[None,:,:],axis=2)]=0;np.fill_diagonal(D,0);n=len(D);J=np.eye(n)-np.ones((n,n))/n;G=-.5*J@(D**2)@J;eig=np.linalg.eigvalsh(G);c=max(0,-float(eig.min()));corrected=D**2+2*c*(np.ones((n,n))-np.eye(n));G=-.5*J@corrected@J
 return G,dict(status='DEGENERATE_ZERO_TOTAL_DISTANCE' if np.trace(G)<=0 else 'EVALUABLE_LINGOES',usable_dimensions=int(cols.sum()),uncorrected_negative_eigenvalues=int((eig<-1e-10).sum()),Lingoes_constant=c)
