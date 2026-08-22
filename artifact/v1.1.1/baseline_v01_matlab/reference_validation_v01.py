#!/usr/bin/env python3
"""Independent exact-integer reference validator for BMD-MATLAB v0.1.

This is intentionally not the production MATLAB implementation. It mirrors the
algebra and reduction rules using Python arbitrary-precision integers so we can
validate correctness and structural scaling in an environment without MATLAB.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import gcd
import csv, json, random, time
from pathlib import Path

VAR_STRIDE = 128
TERMINAL = 0

class BMD:
    def __init__(self):
        self.level=[]; self.ln=[]; self.hn=[]; self.lw=[]; self.hw=[]
        self.unique={}; self.add_cache={}; self.mul_cache={}
        self.add_hits=self.add_misses=self.mul_hits=self.mul_misses=0

    @staticmethod
    def zero(): return (0, TERMINAL)
    @staticmethod
    def one(): return (1, TERMINAL)
    @staticmethod
    def const(c): return (int(c), TERMINAL) if c else (0, TERMINAL)

    def canon(self,e):
        w,n=e
        return (0,TERMINAL) if w==0 else (int(w),int(n))

    def decode_level(self, level):
        return (level-1)//VAR_STRIDE + 1, (level-1)%VAR_STRIDE + 1

    def top(self,e):
        w,n=self.canon(e)
        return 10**30 if n==TERMINAL else self.level[n-1]

    def mk(self, level, low, high):
        low=self.canon(low); high=self.canon(high)
        if high[0]==0:
            return low
        for e in (low,high):
            if e[0] and e[1]!=TERMINAL:
                assert self.level[e[1]-1] > level
        lw,ln=low; hw,hn=high
        g=gcd(abs(lw),abs(hw))
        if g==0: return self.zero()
        first=lw if lw else hw
        factor=-g if first<0 else g
        nlw=lw//factor; nhw=hw//factor
        if nlw==0: ln=TERMINAL
        if nhw==0: hn=TERMINAL
        key=(level,ln,nlw,hn,nhw)
        nid=self.unique.get(key)
        if nid is None:
            self.level.append(level); self.ln.append(ln); self.hn.append(hn)
            self.lw.append(nlw); self.hw.append(nhw)
            nid=len(self.level); self.unique[key]=nid
        return (factor,nid)

    def split(self,e,top):
        e=self.canon(e); w,n=e
        if n==TERMINAL or self.top(e)>top:
            return e,self.zero()
        assert self.top(e)==top
        i=n-1
        return self.canon((w*self.lw[i],self.ln[i])), self.canon((w*self.hw[i],self.hn[i]))

    @staticmethod
    def ordered_pair(a,b): return (a,b) if a<=b else (b,a)

    def add(self,a,b):
        a=self.canon(a); b=self.canon(b)
        if a[0]==0: return b
        if b[0]==0: return a
        if a[1]==TERMINAL and b[1]==TERMINAL: return self.const(a[0]+b[0])
        aa,bb=self.ordered_pair(a,b); key=(aa,bb)
        if key in self.add_cache:
            self.add_hits+=1; return self.add_cache[key]
        self.add_misses+=1
        top=min(self.top(a),self.top(b))
        a0,a1=self.split(a,top); b0,b1=self.split(b,top)
        r=self.mk(top,self.add(a0,b0),self.add(a1,b1))
        self.add_cache[key]=r; return r

    def neg(self,a):
        a=self.canon(a); return self.zero() if a[0]==0 else (-a[0],a[1])
    def sub(self,a,b): return self.add(a,self.neg(b))
    def scale(self,a,c):
        a=self.canon(a); c=int(c)
        return self.zero() if a[0]==0 or c==0 else (a[0]*c,a[1])

    def mul(self,a,b):
        a=self.canon(a); b=self.canon(b)
        if a[0]==0 or b[0]==0: return self.zero()
        if a[1]==TERMINAL: return self.scale(b,a[0])
        if b[1]==TERMINAL: return self.scale(a,b[0])
        aa,bb=self.ordered_pair(a,b); key=(aa,bb)
        if key in self.mul_cache:
            self.mul_hits+=1; return self.mul_cache[key]
        self.mul_misses+=1
        ta,tb=self.top(a),self.top(b); top=min(ta,tb)
        a0,a1=self.split(a,top); b0,b1=self.split(b,top)
        p00=self.mul(a0,b0); p01=self.mul(a0,b1); p10=self.mul(a1,b0)
        if ta==tb:
            p11=self.mul(a1,b1)
            cross=self.add(p01,p10)
            base=self.mk(top,p00,cross)
            var_idx,bit_idx=self.decode_level(top)
            assert bit_idx < VAR_STRIDE
            carry_var=self.mk(top+1,self.zero(),self.one())
            r=self.add(base,self.mul(carry_var,p11))
        else:
            p11=self.mul(a1,b1)
            high=self.add(self.add(p01,p10),p11)
            r=self.mk(top,p00,high)
        self.mul_cache[key]=r; return r

    def monomial_var(self,var_idx,exponent,coeff=1):
        exponent=int(exponent); coeff=int(coeff)
        assert var_idx>=1 and exponent>=0
        r=self.const(coeff)
        if exponent == 0:
            return r
        max_bit=exponent.bit_length()
        for bit in range(max_bit,0,-1):
            if (exponent >> (bit-1)) & 1:
                level=(var_idx-1)*VAR_STRIDE+bit
                r=self.mk(level,self.zero(),r)
        return r
    def monomial(self,e,c=1): return self.monomial_var(1,e,c)
    def variable(self,v): return self.monomial_var(v,1,1)

    def power(self,a,n):
        n=int(n); assert n>=0
        result=self.one(); base=a
        while n:
            if n&1: result=self.mul(result,base)
            n >>= 1
            if n: base=self.mul(base,base)
        return result

    def from_dense(self,p_desc):
        r=self.zero(); deg=len(p_desc)-1
        for i,c in enumerate(p_desc):
            if c: r=self.add(r,self.monomial(deg-i,int(c)))
        return r

    def to_terms_univariate(self,e):
        memo={TERMINAL:{0:1}}
        def node(n):
            if n in memo: return memo[n]
            i=n-1; var,bit=self.decode_level(self.level[i]); assert var==1
            lo=node(self.ln[i]); hi=node(self.hn[i]); shift=1<<(bit-1)
            out={}
            for d,c in lo.items(): out[d]=out.get(d,0)+self.lw[i]*c
            for d,c in hi.items(): out[d+shift]=out.get(d+shift,0)+self.hw[i]*c
            out={d:c for d,c in out.items() if c}; memo[n]=out; return out
        w,n=self.canon(e)
        if w==0:return {}
        return {d:w*c for d,c in node(n).items() if w*c}

    def reachable_internal(self,e):
        seen=set()
        def walk(n):
            if n==TERMINAL or n in seen:return
            seen.add(n); i=n-1; walk(self.ln[i]); walk(self.hn[i])
        walk(self.canon(e)[1]); return len(seen)

    def stats(self,e):
        return dict(reachable_internal=self.reachable_internal(e),
                    final_nodes_including_terminal=self.reachable_internal(e)+1,
                    total_internal=len(self.level),unique_entries=len(self.unique),
                    add_cache=len(self.add_cache),mul_cache=len(self.mul_cache),
                    add_hits=self.add_hits,add_misses=self.add_misses,
                    mul_hits=self.mul_hits,mul_misses=self.mul_misses)

def dense_conv(a,b):
    # descending coefficient vectors, exact ints
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return trim(out)

def trim(p):
    i=0
    while i<len(p)-1 and p[i]==0:i+=1
    return p[i:]

def test_random_mult(count=300,seed=17):
    rnd=random.Random(seed)
    for t in range(count):
        da=rnd.randint(0,8); db=rnd.randint(0,8)
        pa=[rnd.randint(-3,3) for _ in range(da+1)]
        pb=[rnd.randint(-3,3) for _ in range(db+1)]
        if not any(pa): pa[-1]=1
        if not any(pb): pb[-1]=1
        pa=trim(pa); pb=trim(pb)
        m=BMD(); C=m.mul(m.from_dense(pa),m.from_dense(pb))
        terms=m.to_terms_univariate(C)
        deg=max(terms,default=0); got=[terms.get(d,0) for d in range(deg,-1,-1)] if terms else [0]
        exp=dense_conv(pa,pb)
        if got!=exp: raise AssertionError((t,pa,pb,got,exp))

def validate_paper():
    out={}
    m=BMD(); F=m.monomial(5); G=m.add(m.monomial(3),m.monomial(1,2)); H=m.mul(F,G)
    assert m.to_terms_univariate(H)=={8:1,6:2}; out['paper_example']='x^8 + 2*x^6'
    m=BMD(); x=m.monomial(1); lhs=m.mul(m.add(x,m.one()),m.sub(x,m.one())); rhs=m.sub(m.monomial(2),m.one())
    assert lhs==rhs; out['canonical_identity']=True
    mon={}; geo={}
    for n in [50,100,150,200,250,300,350,400,450]:
        m=BMD(); r=m.monomial(n); mon[n]=m.stats(r)['final_nodes_including_terminal']
        m=BMD(); r=m.zero()
        for k in range(n+1): r=m.add(r,m.monomial(k))
        geo[n]=m.stats(r)['final_nodes_including_terminal']
    out['table1_monomial_counts']=mon; out['table1_geometric_counts']=geo
    small=[]
    for rp in [1,4,8]:
        for n in range(1,6):
            m=BMD(); r=m.one()
            for k in range(1,n+1):
                r=m.mul(r,m.power(m.add(m.variable(k),m.one()),rp))
            nodes=m.stats(r)['final_nodes_including_terminal']
            assert nodes==rp*n+1,(rp,n,nodes)
            small.append({'r':rp,'n':n,'nodes':nodes})
    out['table2_small']=small
    return out

def benchmark_table2():
    plans={1:[10,25,50,100,150,200,300],4:[2,4,6,8,10,12,16,20],8:[1,2,3,4,5,6,8,10]}
    rows=[]
    for rp,ns in plans.items():
        for n in ns:
            samples=[]; best=None
            for _ in range(3):
                m=BMD(); r=m.one(); t0=time.perf_counter()
                for k in range(1,n+1):
                    factor=m.power(m.add(m.variable(k),m.one()),rp)
                    r=m.mul(r,factor)
                samples.append(time.perf_counter()-t0); best=(m,r)
            samples.sort(); dt=samples[len(samples)//2]
            m,r=best; s=m.stats(r)
            rows.append({'r':rp,'n':n,'time_s_median3':dt,'final_nodes':s['final_nodes_including_terminal'],
                         'expected_nodes':rp*n+1,'total_internal_nodes':s['total_internal'],
                         'mul_cache':s['mul_cache'],'mul_hits':s['mul_hits'],'mul_misses':s['mul_misses']})
            assert s['final_nodes_including_terminal']==rp*n+1
    return rows

def benchmark_univariate():
    rows=[]
    # structural families; exact Python BMD only, not a MATLAB performance claim
    for family,ns in [('monomial',[10,100,1000,10000,1000000]),('geometric',[50,100,200,400,800]),('binomial',[5,10,15,20,25])]:
        for n in ns:
            m=BMD(); t0=time.perf_counter()
            if family=='monomial': r=m.monomial(n)
            elif family=='geometric':
                r=m.zero()
                for k in range(n+1): r=m.add(r,m.monomial(k))
            else:
                r=m.one(); f=m.add(m.monomial(1),m.one())
                for _ in range(n): r=m.mul(r,f)
            dt=time.perf_counter()-t0; s=m.stats(r)
            rows.append({'family':family,'n':n,'time_s':dt,'final_nodes':s['final_nodes_including_terminal'],
                         'total_internal_nodes':s['total_internal'],'add_cache':s['add_cache'],'mul_cache':s['mul_cache']})
    return rows

def main():
    root=Path(__file__).resolve().parent; outdir=root/'results'; outdir.mkdir(exist_ok=True)
    test_random_mult(300)
    proof=validate_paper(); proof['random_multiplications']=300
    (outdir/'reference_validation.json').write_text(json.dumps(proof,indent=2,sort_keys=True),encoding='utf-8')
    t2=benchmark_table2(); uni=benchmark_univariate()
    for name,rows in [('reference_table2_scaling.csv',t2),('reference_univariate_scaling.csv',uni)]:
        with (outdir/name).open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print('PASS: 300 random exact multiplications')
    print('PASS: paper example and canonical identity')
    print('PASS: Table 2 node identity r*n+1 for all reference cases')
    print('Wrote',outdir/'reference_validation.json')
    print('Wrote',outdir/'reference_table2_scaling.csv')
    print('Wrote',outdir/'reference_univariate_scaling.csv')
    print('\nTable 2 reference scaling:')
    for row in t2: print(row)

if __name__=='__main__': main()
