#!/usr/bin/env python3
"""Independent exact-integer structural validation for BMD-MATLAB v0.7.

Validates the fixed-input operation-closure family without MATLAB timing.
Both operands always have 256 explicit terms and 8 BMD nodes; only the bit
band position of Q changes. Exact BMD products are compared with an
independent sparse convolution.
"""
from __future__ import annotations
import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
REF = HERE.parent / 'baseline_v01_matlab' / 'reference_validation_v01.py'
spec = importlib.util.spec_from_file_location('bmd_ref_v01', REF)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
BMD = mod.BMD


def indicator_exponents(m: BMD, exponents):
    vals = tuple(sorted(set(int(e) for e in exponents)))
    def rec(v, level):
        if not v:
            return m.zero()
        if len(v) == 1 and v[0] == 0:
            return m.one()
        even, odd = [], []
        for e in v:
            (odd if e & 1 else even).append(e >> 1)
        lo = rec(tuple(even), level + 1) if even else m.zero()
        hi = rec(tuple(odd), level + 1) if odd else m.zero()
        return m.mk(level, lo, hi)
    return rec(vals, 1)


def reachable_union(m: BMD, refs):
    seen = set()
    def walk_node(n):
        if n == 0 or n in seen:
            return
        seen.add(n)
        i = n - 1
        walk_node(m.ln[i]); walk_node(m.hn[i])
    for r in refs:
        walk_node(m.canon(r)[1])
    return len(seen)


def sparse_product(a, b):
    out = {}
    for ea in a:
        for eb in b:
            e = ea + eb
            out[e] = out.get(e, 0) + 1
    return {e: c for e, c in out.items() if c}


def main():
    n = 256
    shifts = [0,1,2,3,4,5,6,7,8,9,10,12]
    p_exps = list(range(n))
    rows = []
    print('v0.7 independent exact operation-closure validation')
    print('s overlap Pnodes Qnodes union result new cache sparse_out max_coeff exact')
    for s in shifts:
        q_exps = [j * (1 << s) for j in range(n)]
        m = BMD()
        p = indicator_exponents(m, p_exps)
        q = indicator_exponents(m, q_exps)
        pnodes = m.reachable_internal(p)
        qnodes = m.reachable_internal(q)
        union = reachable_union(m, [p,q])
        total_before = len(m.level)
        prod = m.mul(p,q)
        new_nodes = len(m.level) - total_before
        result_nodes = m.reachable_internal(prod)
        got = m.to_terms_univariate(prod)
        expected = sparse_product(p_exps, q_exps)
        assert got == expected, s
        assert pnodes == 8 and qnodes == 8, (s,pnodes,qnodes)
        if s >= 8:
            assert result_nodes == 16, (s,result_nodes)
            assert new_nodes == 8, (s,new_nodes)
            assert len(expected) == 65536 and max(expected.values()) == 1
        rows.append((s,max(0,8-s),pnodes,qnodes,union,result_nodes,new_nodes,len(m.mul_cache),len(expected),max(expected.values())))
        print(f'{s:2d} {max(0,8-s):7d} {pnodes:6d} {qnodes:6d} {union:5d} {result_nodes:6d} {new_nodes:4d} {len(m.mul_cache):5d} {len(expected):10d} {max(expected.values()):9d} PASS')

    # Strong endpoint separation: same input node counts, radically different
    # operation closure after changing only bit-band alignment.
    assert rows[0][5] > 20 * rows[8][5], (rows[0],rows[8])
    assert rows[0][6] > 50 * rows[8][6], (rows[0],rows[8])
    print('PASS fixed explicit input cardinality: 256 x 256 terms at every shift')
    print('PASS fixed individual BMD input size: 8 nodes each at every shift')
    print(f'PASS result-node endpoint separation: {rows[0][5]} -> {rows[8][5]}')
    print(f'PASS new-node endpoint separation: {rows[0][6]} -> {rows[8][6]}')
    print('PASS exact BMD product equals independent sparse convolution at every shift')


if __name__ == '__main__':
    main()
