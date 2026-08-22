#!/usr/bin/env python3
"""Independent exact-integer validation for BMD-MATLAB v0.8.

The experiment fixes both operands at 256 terms and 8 BMD nodes and chooses
Q bit signatures so cold multiplication creates exactly 8,12,...,40 new BMD
workspace nodes. BMD products are compared exactly with independent sparse
convolution.
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


def support_from_bits(bits):
    vals = [0]
    for b in bits:
        vals += [e + (1 << b) for e in vals]
    vals = sorted(vals)
    assert len(vals) == len(set(vals)) == (1 << len(bits))
    return vals


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
    def walk(n):
        if n == 0 or n in seen:
            return
        seen.add(n)
        i = n - 1
        walk(m.ln[i]); walk(m.hn[i])
    for r in refs:
        walk(m.canon(r)[1])
    return len(seen)


def sparse_product(a, b):
    out = {}
    for ea in a:
        for eb in b:
            e = ea + eb
            out[e] = out.get(e, 0) + 1
    return {e: c for e, c in out.items() if c}


def main():
    p_bits = tuple(range(8))
    high7 = tuple(range(9, 16))
    overlap_bits = [None, 7, 6, 5, 4, 3, 2, 1, 0]
    expected_new = [8, 12, 16, 20, 24, 28, 32, 36, 40]
    expected_result = [16, 16, 17, 18, 19, 20, 21, 22, 23]
    p_exps = support_from_bits(p_bits)

    print('v0.8 independent exact threshold-calibration validation')
    print('target overlap qbits Pnodes Qnodes union result new cache sparse_out max_coeff exact')
    for ov, enew, eres in zip(overlap_bits, expected_new, expected_result):
        q_bits = tuple(range(9,17)) if ov is None else (ov,) + high7
        q_exps = support_from_bits(q_bits)
        assert len(p_exps) == len(q_exps) == 256

        m = BMD()
        p = indicator_exponents(m, p_exps)
        q = indicator_exponents(m, q_exps)
        pn = m.reachable_internal(p)
        qn = m.reachable_internal(q)
        union = reachable_union(m, [p, q])
        before = len(m.level)
        prod = m.mul(p, q)
        new_nodes = len(m.level) - before
        result_nodes = m.reachable_internal(prod)
        got = m.to_terms_univariate(prod)
        expected = sparse_product(p_exps, q_exps)

        assert pn == 8 and qn == 8 and union == 16, (ov,pn,qn,union)
        assert new_nodes == enew, (ov,new_nodes,enew)
        assert result_nodes == eres, (ov,result_nodes,eres)
        assert got == expected, ov
        assert len(p_exps)*len(q_exps) == 65536
        label = -1 if ov is None else ov
        bits_text = ','.join(map(str,q_bits))
        print(f'{enew:6d} {label:7d} {bits_text:23s} {pn:6d} {qn:6d} {union:5d} {result_nodes:6d} {new_nodes:4d} {len(m.mul_cache):5d} {len(expected):10d} {max(expected.values()):9d} PASS')

    print('PASS fixed explicit input cardinality: 256 x 256 terms in all 9 cases')
    print('PASS fixed individual BMD input size: 8 nodes each; union 16 nodes')
    print('PASS exact new-node calibration grid: 8,12,16,20,24,28,32,36,40')
    print('PASS exact BMD product equals independent sparse convolution in every case')


if __name__ == '__main__':
    main()
