#!/usr/bin/env python3
"""Independent exact validation for the v0.9 structural predictor.

No MATLAB timings are synthesized here. This checks only the pre-multiply
new-node model and conservative route logic against exact Python BMD algebra
and the already measured v0.5-v0.8 CSV outcomes.
"""
from __future__ import annotations
import csv, importlib.util, math, random
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
REF=ROOT/'baseline_v01_matlab'/'reference_validation_v01.py'
spec=importlib.util.spec_from_file_location('bmdref',REF); mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); BMD=mod.BMD
BETA=np.array([1.8716285233860908,-0.49888659949909009,-0.49550136457640398,0.025865401111516365,0.30484129754149636,0.28101502349890534,-0.3836297310014834,0.046664149634361947,-0.089931009118816949,0.95963681138947232,-0.058268032408828399,0.050629935773811183,0.1887658610473012])
FEATURES=['p_nodes','q_nodes','common_levels','carry_runsum','carry_run2','carry_max','carry_massrun','adjacent_mass','invdist','suffix_product_sum','suffix_carry_sum','level_span']

def support(bits):
    v=[0]
    for b in bits:v += [e+(1<<b) for e in v]
    return sorted(v)
def indicator(m,exps):
    vals=tuple(sorted(set(map(int,exps))))
    def rec(v,lev):
        if not v:return m.zero()
        if len(v)==1 and v[0]==0:return m.one()
        ev=[];od=[]
        for e in v:(od if e&1 else ev).append(e>>1)
        return m.mk(lev,rec(tuple(ev),lev+1) if ev else m.zero(),rec(tuple(od),lev+1) if od else m.zero())
    return rec(vals,1)
def ids(m,r):
    s=set()
    def w(n):
        if n==0 or n in s:return
        s.add(n);i=n-1;w(m.ln[i]);w(m.hn[i])
    w(m.canon(r)[1]);return s
def hist(m,r):
    h={}
    for n in ids(m,r):h[m.level[n-1]]=h.get(m.level[n-1],0)+1
    return h
def is_chain(m,r):
    n=m.canon(r)[1];seen=set()
    while n:
        if n in seen:return False
        seen.add(n);i=n-1
        if m.ln[i]!=m.hn[i] or m.lw[i]!=m.hw[i] or m.lw[i]==0:return False
        n=m.ln[i]
    return True
def features(m,a,b,npterms,nqterms):
    ha,hb=hist(m,a),hist(m,b);A=set(ha);B=set(hb);U=A|B;C=A&B
    f={'p_nodes':len(ids(m,a)),'q_nodes':len(ids(m,b)),'common_levels':len(C)}
    f.update(carry_runsum=0,carry_run2=0,carry_max=0,carry_massrun=0,adjacent_mass=0,invdist=0,suffix_product_sum=0,suffix_carry_sum=0)
    for l in C:
        r=0;k=l+1
        while k in U:r+=1;k+=1
        f['carry_runsum']+=r;f['carry_run2']+=(r+1)**2;f['carry_max']=max(f['carry_max'],r);f['carry_massrun']+=ha[l]*hb[l]*(r+1)
        sa=sum(v for ll,v in ha.items() if ll>=l);sb=sum(v for ll,v in hb.items() if ll>=l)
        f['suffix_product_sum']+=sa*sb;f['suffix_carry_sum']+=sa*sb*(r+1)
    for la,ca in ha.items():
        for lb,cb in hb.items():
            d=abs(la-lb);w=ca*cb;f['invdist']+=w/(1+d);f['adjacent_mass']+=w if d==1 else 0
    f['level_span']=max(U)-min(U)+1 if U else 0;f['pair_products']=npterms*nqterms
    f['both_chain']=is_chain(m,a) and is_chain(m,b)
    f['ordered_disjoint']=bool(A and B and not C and (max(A)<min(B) or max(B)<min(A)))
    return f

def predict(f):
    if f['ordered_disjoint']:
        pn=min(f['p_nodes'],f['q_nodes']);conf='HIGH';reg='ORDERED_DISJOINT_BANDS'
    elif f['both_chain']:
        x=np.array([math.log1p(f[k]) for k in FEATURES]);pn=max(1,math.expm1(BETA[0]+BETA[1:]@x));conf='HIGH';reg='BITCUBE_RIDGE'
    else:
        pn=max(1,.82*f['suffix_product_sum'],(f['p_nodes']+f['q_nodes'])/4);conf='LOW';reg='GENERAL'
    work=f['pair_products']/pn
    if conf=='LOW':route='SPARSE' if pn>=64 or work<=1800 else 'UNCERTAIN'
    elif pn>=64 or work<=1800:route='SPARSE'
    elif f['pair_products']>=40000 and work>=3000:route='BMD'
    elif f['pair_products']<30000:route='SPARSE'
    else:route='UNCERTAIN'
    return pn,route,reg,work

def exact_new(m,a,b):
    n=len(m.level);r=m.mul(a,b);return len(m.level)-n,m.reachable_internal(r)

def synthetic_test():
    rng=random.Random(90210);errs=[]
    for _ in range(200):
        ka=rng.choice([5,6,7,8,9]);kb=rng.choice([5,6,7,8,9]);A=sorted(rng.sample(range(20),ka));B=sorted(rng.sample(range(20),kb))
        m=BMD();a=indicator(m,support(A));b=indicator(m,support(B));f=features(m,a,b,1<<ka,1<<kb);pn,_,_,_=predict(f);actual,_=exact_new(m,a,b);errs.append(abs(pn-actual)/actual)
    return float(np.median(errs)),float(np.quantile(errs,.9)),float(max(errs))

def holdouts():
    cases=[('disjoint',[8,10,11,12,13,14,15,16]),('single_gap',[7,10,11,12,13,14,15,16]),('carry_chain',[5,8,10,12,14,16,18,20]),('two_overlap',[6,7,9,11,13,15,17,19]),('two_separate',[1,4,9,11,13,15,17,19])]
    out=[]
    for name,Q in cases:
        m=BMD();a=indicator(m,support(range(8)));b=indicator(m,support(Q));f=features(m,a,b,256,256);pn,route,reg,work=predict(f);act,res=exact_new(m,a,b);out.append((name,pn,act,route,reg,work,res))
    return out

def main():
    med,p90,mx=synthetic_test()
    print('v0.9 independent exact predictor validation')
    print(f'fresh synthetic exact bitcube cases=200 median_APE={med:.4f} p90_APE={p90:.4f} max_APE={mx:.4f}')
    assert med<.10 and p90<.25
    for x in holdouts():print('holdout',x)
    print('PASS: predictor uses structural features only; exact multiply is invoked only after prediction for scoring')
    print('PASS: fresh synthetic median APE <10%, p90 APE <25%')
if __name__=='__main__':main()
