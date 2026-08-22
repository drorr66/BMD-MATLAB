#!/usr/bin/env python3
"""Independent exact validation for the frozen v1.1 publication workload set.

No MATLAB timings are synthesized. This independently reconstructs all 60
polynomial pairs, computes exact Python-integer convolution, converts operands
to the exact BMD reference, and verifies coefficient dictionaries and structural
counts. It also writes expected_structure_v11.csv for audit only.
"""
from __future__ import annotations
import csv, importlib.util
from math import comb
from pathlib import Path
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
REF=ROOT/'baseline_v01_matlab'/'reference_validation_v01.py'
spec=importlib.util.spec_from_file_location('bmdref',REF); mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); BMD=mod.BMD
MAX_EXACT=9007199254740991
M1=[0,1,4,5,16,17,20,21]; M2=[0,2,3,10,11,18,24,25]

def norm(d): return {int(k):int(v) for k,v in d.items() if v}
def geom(n,step=1): return {i*step:1 for i in range(n)}
def conv(a,b):
    out={}
    for ea,ca in a.items():
        for eb,cb in b.items(): out[ea+eb]=out.get(ea+eb,0)+ca*cb
    return norm(out)
def repeat_mask(mask,blocks,stride=256): return {j*stride+m:1 for j in range(blocks) for m in mask}
def lcg_support(count,modulus,seed):
    used=set(); out=[]; x=seed
    while len(out)<count:
        x=(1664525*x+1013904223)%(2**32); v=x%modulus
        if v not in used: used.add(v); out.append(v)
    return sorted(out)
def irregular(n,modulus,s1): return {e:1 for e in lcg_support(n,modulus,s1)}
def subset(weights):
    p={0:1}
    for w in weights: p=conv(p,{0:1,w:1})
    return p
def bounded(factors):
    p={0:1}
    for d,m in factors: p=conv(p,{j*d:1 for j in range(m+1)})
    return p
def binomial(n,sign): return {k:comb(n,k)*(sign**k) for k in range(n+1)}
def weighted(n,pat,step=1): return norm({k*step:pat[k%len(pat)] for k in range(n)})

def cases():
    out=[]
    for a,b in [(32,32),(64,64),(128,64),(128,128),(256,128),(256,256)]: out.append((f'fir_box_{a}x{b}','BOXCAR_FIR',geom(a),geom(b)))
    for a,b,d1,d2 in [(32,32,2,3),(48,48,2,4),(64,64,3,5),(64,96,4,8),(96,96,5,9),(128,64,8,16),(128,128,7,11),(192,96,16,32)]: out.append((f'cic_comb_{a}x{b}_d{d1}_d{d2}','CIC_COMB',geom(a,d1),geom(b,d2)))
    for a,b,s in [(64,64,128),(64,128,128),(128,64,256),(128,128,256),(128,256,256),(256,64,512),(256,128,512),(256,256,512)]: out.append((f'polyphase_sep_{a}x{b}_s{s}','POLYPHASE_SEPARATED',geom(a),geom(b,s)))
    P=[(M1,16,64,1),(M1,32,64,1),(M1,32,128,1),(M2,32,96,1),(M1,16,64,8192),(M1,32,64,16384),(M1,32,128,32768),(M2,32,96,32768)]
    for i,(mask,blocks,b,s) in enumerate(P,1): out.append((f'periodic_mask_{i:02d}','PERIODIC_MASK',repeat_mask(mask,blocks),geom(b,s)))
    W=[([1,2,4,8,16,32],[64,128,256,512,1024,2048]),([1,3,7,15,31,63],[127,255,511,1023,2047]),([2,5,11,23,47,95],[191,383,767,1535]),([3,7,13,29,59],[5,11,17,37,73]),([4,9,19,39,79],[6,13,27,55,111]),([5,10,20,40,80,160],[7,14,28,56,112,224]),([1,5,10,25,50],[2,10,20,50,100]),([6,13,27,54,109],[8,17,35,71,143])]
    for i,(w1,w2) in enumerate(W,1): out.append((f'subset_sum_{i:02d}','SUBSET_SUM',subset(w1),subset(w2)))
    F=[([(1,3),(4,3),(16,3)],[(64,3),(256,3)]),([(1,4),(5,2),(25,2)],[(125,2),(625,2)]),([(2,3),(7,3),(29,2)],[(5,3),(19,2),(83,2)]),([(3,4),(12,3),(48,2)],[(192,2),(768,2)]),([(1,5),(10,3),(100,2)],[(1000,2),(10000,1)]),([(4,3),(20,3),(100,2)],[(6,3),(30,3),(150,2)])]
    for i,(f1,f2) in enumerate(F,1): out.append((f'bounded_resource_{i:02d}','BOUNDED_RESOURCE',bounded(f1),bounded(f2)))
    B=[(8,8,1,1),(12,10,1,1),(16,12,1,-1),(18,14,-1,-1),(20,8,1,-1),(22,10,-1,1)]
    for i,(n,m,s1,s2) in enumerate(B,1): out.append((f'binomial_diff_{i:02d}','BINOMIAL_DIFFERENCE',binomial(n,s1),binomial(m,s2)))
    I=[(24,1024,17,91),(32,2048,123,991),(40,4096,7,19),(48,4096,71,271),(56,8192,313,911),(64,8192,17,991)]
    for i,(n,md,s1,s2) in enumerate(I,1): out.append((f'irregular_control_{i:02d}','IRREGULAR_CONTROL',irregular(n,md,s1),irregular(n,md,s2)))
    WP=[(32,[1,-1,2,-2],[2,1,-2,-1],1,1),(48,[1,0,-1,2,0,-2],[1,2,-1,-2],1,1),(64,[1,-1,1,-1],[2,-2,3,-3],2,3),(96,[1,-1,2,-2,3,-3],[2,1,-2,-1],1,1)]
    for i,(n,p1,p2,s1,s2) in enumerate(WP,1): out.append((f'weighted_periodic_{i:02d}','WEIGHTED_PERIODIC',weighted(n,p1,s1),weighted(n,p2,s2)))
    assert len(out)==60
    return out

def indicator(m,exps):
    vals=tuple(sorted(set(map(int,exps)))); memo={}
    def rec(v,lev):
        key=(v,lev)
        if key in memo:return memo[key]
        if not v:return m.zero()
        if len(v)==1 and v[0]==0:return m.one()
        ev=[];od=[]
        for e in v:(od if e&1 else ev).append(e>>1)
        lo=rec(tuple(ev),lev+1) if ev else m.zero(); hi=rec(tuple(od),lev+1) if od else m.zero()
        r=m.mk(lev,lo,hi); memo[key]=r; return r
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
    rows=[]; maxnew=0; maxcoeff=0
    print('v1.1 independent exact publication validation')
    for idx,(name,fam,p,q) in enumerate(cases(),1):
        expected=conv(p,q); mc=max(map(abs,expected.values()),default=0); maxcoeff=max(maxcoeff,mc); assert mc<=MAX_EXACT,(name,mc)
        m=BMD(); a=to_bmd(m,p); b=to_bmd(m,q); before=len(m.level); r=m.mul(a,b); got=m.to_terms_univariate(r); assert got==expected,name
        new=len(m.level)-before; rn=m.reachable_internal(r); maxnew=max(maxnew,new)
        rows.append({'case_index':idx,'case_name':name,'family':fam,'p_terms':len(p),'q_terms':len(q),'pair_products':len(p)*len(q),'actual_new_nodes':new,'actual_result_nodes':rn,'result_terms':len(expected),'max_abs_coeff':mc})
        print(f'{idx:02d} {name:30s} {fam:22s} {len(p):4d}x{len(q):4d} new={new:6d} result={rn:5d} terms={len(expected):5d}')
    out=HERE/'expected_structure_v11.csv'
    with out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    print(f'PASS: 60/60 exact products; max new nodes={maxnew}; max abs coefficient={maxcoeff}')
    print(f'WROTE: {out}')
if __name__=='__main__': main()
