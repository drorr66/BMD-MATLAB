#!/usr/bin/env python3
"""Independent exact-integer reference validation for BMD-MATLAB v0.6.

Validates the controlled structural-sharing family and the arbitrary exponent
set recursion without requiring MATLAB. It does NOT make MATLAB timing claims.
"""
from __future__ import annotations
import importlib.util
import math
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


def template_bank(blocks=64, universe=256, k=128):
    bank = [tuple(range(k))]
    for t in range(2, blocks + 1):
        scored = []
        for i in range(1, universe + 1):
            x = math.sin(i * 12.9898 + t * 78.233) * 43758.5453123
            score = x - math.floor(x)
            scored.append((score, i - 1))
        scored.sort()
        bank.append(tuple(sorted(i for _, i in scored[:k])))
    assert len(set(bank)) == blocks
    assert all(len(x) == k for x in bank)
    return bank


def sharing_case(bank, unique_templates, blocks=64, stride_power=9):
    stride = 1 << stride_power
    exps = []
    ids = []
    for j in range(blocks):
        tid = j % unique_templates
        ids.append(tid)
        exps.extend(j * stride + i for i in bank[tid])
    assert len(exps) == len(set(exps))
    return exps, ids


def sparse_multiply_terms(a, b):
    out = {}
    for ea in a:
        for eb in b:
            e = ea + eb
            out[e] = out.get(e, 0) + 1
    return {e: c for e, c in out.items() if c}


def main():
    bank = template_bank()
    counts = [1, 2, 4, 8, 16, 24, 32, 48, 64]
    q_exp = list(range(32))
    prev_nodes = None
    print('v0.6 independent exact validation')
    print('T  terms  P_nodes  result_nodes  exact_product')
    for t in counts:
        exps, ids = sharing_case(bank, t)
        assert len(exps) == 8192
        assert len(set(ids)) == t

        m = BMD()
        p = indicator_exponents(m, exps)
        q = indicator_exponents(m, q_exp)
        p_nodes = m.reachable_internal(p)
        prod = m.mul(p, q)
        result_nodes = m.reachable_internal(prod)
        got = m.to_terms_univariate(prod)
        expected = sparse_multiply_terms(exps, q_exp)
        assert got == expected, t

        # The chosen deterministic family is designed so loss of deliberate
        # block-template reuse increases final P complexity overall. A few
        # local non-monotonicities would not invalidate the experiment, but
        # endpoints must be strongly separated.
        if t == 1:
            first_nodes = p_nodes
        if t == 64:
            last_nodes = p_nodes
        print(f'{t:2d} {len(exps):6d} {p_nodes:7d} {result_nodes:12d} PASS')
        prev_nodes = p_nodes

    assert last_nodes > 50 * first_nodes, (first_nodes, last_nodes)
    print(f'PASS fixed sparse cardinality: 8192 terms at all {len(counts)} levels')
    print(f'PASS sharing endpoint separation: {first_nodes} -> {last_nodes} BMD nodes')
    print('PASS exact BMD products equal independent sparse convolution at every level')


if __name__ == '__main__':
    main()
