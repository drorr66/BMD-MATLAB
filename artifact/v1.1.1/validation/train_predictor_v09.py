import random, math, importlib.util
from pathlib import Path
import numpy as np
# helpers from explore file before baseline rows
s=open('/mnt/data/explore_v09d.py').read().split('# load actual baseline rows')[0]
ns={};exec(s,ns)
BMD=ns['BMD'];indicator=ns['indicator'];support_bits=ns['support_bits'];feats=ns['feats'];actual_new=ns['actual_new']

def make_case(A,B):
 m=BMD();a=indicator(m,support_bits(A));b=indicator(m,support_bits(B));ft=feats(m,a,b,1<<len(A),1<<len(B));new,res,cache=actual_new(m,a,b);return ft,new
rng=random.Random(20260820)
rows=[]
# include systematic fixed P low band with diverse Q patterns
P=tuple(range(8))
patterns=[]
for _ in range(500):
 B=tuple(sorted(rng.sample(range(0,22),8)));patterns.append((P,B))
# random pairs, varying 5-9 bits each to broaden
for _ in range(500):
 ka=rng.choice([5,6,7,8,9]);kb=rng.choice([5,6,7,8,9])
 A=tuple(sorted(rng.sample(range(0,20),ka)));B=tuple(sorted(rng.sample(range(0,20),kb)));patterns.append((A,B))
for i,(A,B) in enumerate(patterns):
 ft,new=make_case(A,B);rows.append((A,B,ft,new))
 if (i+1)%100==0: print('built',i+1)
# feature transforms
featnames=['p_nodes','q_nodes','common_levels','carry_runsum','carry_run2','carry_max','carry_massrun','adjacent_mass','invdist','suffix_product_sum','suffix_carry_sum','level_span']
def xvec(ft):return [math.log1p(float(ft[n])) for n in featnames]
X=np.array([xvec(ft) for _,_,ft,_ in rows]); y=np.log1p([new for *_,new in rows])
idx=np.arange(len(rows));rng.shuffle(idx.tolist())
# deterministic split by index modulo 5
tr=np.array([i%5!=0 for i in range(len(rows))]);te=~tr
for alpha in [0.001,0.01,0.1,1,10]:
 xx=np.c_[np.ones(tr.sum()),X[tr]];xt=np.c_[np.ones(te.sum()),X[te]];I=np.eye(xx.shape[1]);I[0,0]=0
 beta=np.linalg.solve(xx.T@xx+alpha*I,xx.T@y[tr]);pr=np.expm1(xt@beta);ac=np.expm1(y[te]);
 rel=np.abs(pr-ac)/ac;fac=np.maximum(pr/ac,ac/np.maximum(pr,1e-9))
 print('alpha',alpha,'medAPE',np.median(rel),'p90',np.quantile(rel,.9),'maxfac',max(fac),'beta',beta.tolist())
# choose .01 or .1 based p90
alpha=.01
xx=np.c_[np.ones(tr.sum()),X[tr]];I=np.eye(xx.shape[1]);I[0,0]=0
beta=np.linalg.solve(xx.T@xx+alpha*I,xx.T@y[tr])
print('FEATURES',featnames); print('BETA',','.join(f'{v:.17g}' for v in beta))
# evaluate named holdouts
hold=[
('H1',range(8),[8,10,11,12,13,14,15,16]),
('H2',range(8),[7,10,11,12,13,14,15,16]),
('H3',range(8),[5,8,10,12,14,16,18,20]),
('H4',range(8),[6,7,9,11,13,15,17,19]),
('H5',range(8),[1,4,9,11,13,15,17,19]),
]
for name,A,B in hold:
 ft,new=make_case(tuple(A),tuple(B));pv=math.expm1(np.r_[1,xvec(ft)]@beta)
 print(name,'pred',pv,'new',new,'ratio',pv/new)
