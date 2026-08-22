#!/usr/bin/env python3
"""Independent exact validation for the frozen v1.0 cross-family case set.

This does not synthesize MATLAB timings.  It reconstructs the 15 deterministic
sparse polynomial pairs independently, converts them to the Python exact BMD
reference, multiplies them, and checks exact coefficient dictionaries.
"""
from __future__ import annotations
import importlib.util, math
from pathlib import Path

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
REF=ROOT/'baseline_v01_matlab'/'reference_validation_v01.py'
spec=importlib.util.spec_from_file_location('bmdref',REF); mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); BMD=mod.BMD

MASK=[0,1,4,5,16,17,20,21]

def norm(d): return {int(k):int(v) for k,v in d.items() if v}
def geom(count,step=1): return {i*step:1 for i in range(count)}
def repeat_mask(blocks,stride=256): return {(j*stride+m):1 for j in range(blocks) for m in MASK}
def lcg_support(count,modulus,seed):
    used=set(); out=[]; x=seed
    while len(out)<count:
        x=(1664525*x+1013904223)%(2**32); v=x%modulus
        if v not in used: used.add(v); out.append(v)
    return sorted(out)
def irregular(count,modulus,seed): return {e:1 for e in lcg_support(count,modulus,seed)}
def binomial(n,sign=1):
    d={}; c=1
    for k in range(n+1):
        if k>0: c=c*(n-k+1)//k
        d[k]=c*(sign**k)
    return norm(d)
def dense_pm1(count,seed):
    d={}
    for k in range(count):
        c=1
        if (k+seed)%3==0:c=-1
        if (k+2*seed)%7==0:c=2
        d[k]=c
    return norm(d)
def weighted(count,pat): return norm({k:pat[k%len(pat)] for k in range(count)})

def cases():
    return [
      ('box_overlap_128','BOXCAR',geom(128),geom(128)),
      ('box_overlap_256','BOXCAR',geom(256),geom(256)),
      ('box_disjoint_128x320','BLOCK_GENERATING',geom(128),geom(320,256)),
      ('box_disjoint_256x224','BLOCK_GENERATING',geom(256),geom(224,512)),
      ('comb_3x5_128','COMB',geom(128,3),geom(128,5)),
      ('comb_5x9_192','COMB',geom(192,5),geom(192,9)),
      ('periodic_mask_overlap','PERIODIC_MASK',repeat_mask(64),geom(128)),
      ('periodic_mask_disjoint','PERIODIC_MASK',repeat_mask(64),geom(128,32768)),
      ('irregular_sparse_64','IRREGULAR_CONTROL',irregular(64,2048,17),irregular(64,2048,91)),
      ('irregular_sparse_96','IRREGULAR_CONTROL',irregular(96,4096,123),irregular(96,4096,991)),
      ('binomial_12x10','BINOMIAL',binomial(12,1),binomial(10,1)),
      ('binomial_18x14','BINOMIAL',binomial(18,1),binomial(14,1)),
      ('difference_12x10','FINITE_DIFFERENCE',binomial(12,-1),binomial(10,-1)),
      ('dense_pm1_64','DENSE_CONTROL',dense_pm1(64,7),dense_pm1(64,19)),
      ('weighted_periodic_96','WEIGHTED_PERIODIC',weighted(96,[1,-1,2,-2,3,-3]),weighted(96,[2,1,-2,-1])),
    ]

def conv(a,b):
    out={}
    for ea,ca in a.items():
        for eb,cb in b.items(): out[ea+eb]=out.get(ea+eb,0)+ca*cb
    return norm(out)

def indicator(m,exps):
    vals=tuple(sorted(set(map(int,exps))))
    def rec(v,lev):
        if not v:return m.zero()
        if len(v)==1 and v[0]==0:return m.one()
        ev=[]; od=[]
        for e in v:(od if e&1 else ev).append(e>>1)
        lo=rec(tuple(ev),lev+1) if ev else m.zero(); hi=rec(tuple(od),lev+1) if od else m.zero()
        return m.mk(lev,lo,hi)
    return rec(vals,1)

def to_bmd(m,p):
    vals=set(p.values())
    if vals=={1}: return indicator(m,p.keys())
    if len(vals)==1:
        c=next(iter(vals)); return m.scale(indicator(m,p.keys()),c)
    r=m.zero()
    for e in sorted(p): r=m.add(r,m.monomial(e,p[e]))
    return r

def main():
    print('v1.0 independent exact cross-family validation')
    max_nodes=0
    for name,fam,p,q in cases():
        expected=conv(p,q); m=BMD(); a=to_bmd(m,p); b=to_bmd(m,q); before=len(m.level); r=m.mul(a,b); got=m.to_terms_univariate(r)
        assert got==expected,(name,len(got),len(expected))
        new=len(m.level)-before; rn=m.reachable_internal(r); max_nodes=max(max_nodes,len(m.level))
        max_coeff=max(map(abs,expected.values()),default=0)
        assert max_coeff <= 9007199254740991,(name,max_coeff)
        print(f'{name:28s} family={fam:18s} terms={len(p):4d}x{len(q):4d} new={new:6d} result_nodes={rn:6d} result_terms={len(expected):6d}')
    print(f'PASS: 15/15 exact products; max total internal nodes observed={max_nodes}')
    print('PASS: all result coefficients remain within MATLAB v1.0 exact-integer guard')
if __name__=='__main__': main()
