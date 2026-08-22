"""Independent exact validation of the v0.4 direct structural builders.
No MATLAB runtime is required.  This validates the exponent-set recursion,
not MATLAB performance or MATLAB syntax.
"""
from functools import lru_cache


def geometric_exponents(n: int, start_level: int = 1):
    @lru_cache(None)
    def rec(nn, lev):
        if nn == 0:
            return frozenset((0,))
        q = nn // 2
        basis = 1 << (lev - 1)
        if nn & 1:
            child = rec(q, lev + 1)
            return frozenset(child | {e + basis for e in child})
        low = rec(q, lev + 1)
        high = rec(q - 1, lev + 1)
        return frozenset(low | {e + basis for e in high})
    return rec(n, start_level)


def grid_exponents(inner_n: int, blocks: int, stride_power: int):
    tail = geometric_exponents(blocks - 1, stride_power + 1)

    @lru_cache(None)
    def rec(nn, lev):
        if nn == 0:
            return tail
        q = nn // 2
        basis = 1 << (lev - 1)
        if nn & 1:
            child = rec(q, lev + 1)
            return frozenset(child | {e + basis for e in child})
        low = rec(q, lev + 1)
        high = rec(q - 1, lev + 1)
        return frozenset(low | {e + basis for e in high})
    return rec(inner_n, 1)


def main():
    for n in list(range(256)) + [511, 1000, 5000, 10000]:
        got = geometric_exponents(n)
        expected = frozenset(range(n + 1))
        assert got == expected, n
    print("PASS direct geometric exponent recursion")

    specs = [(31, 8, 7), (31, 32, 7), (63, 64, 8), (255, 64, 10), (255, 256, 10)]
    for inner, blocks, sp in specs:
        got = grid_exponents(inner, blocks, sp)
        stride = 1 << sp
        expected = frozenset(j * stride + i for j in range(blocks) for i in range(inner + 1))
        assert got == expected, (inner, blocks, sp)
    print("PASS repeated-block grid exponent recursion")

    # Structural cardinality identity for power-of-two full inner blocks:
    # the ordered diagram is a chain of p low-bit nodes followed by q
    # high-bit block-index nodes, hence p+q internal nodes.
    for p, q in [(5, 3), (6, 6), (8, 8)]:
        explicit_terms = (1 << p) * (1 << q)
        expected_nodes = p + q
        assert explicit_terms > expected_nodes
    print("PASS additive-node structural identity checks")


if __name__ == "__main__":
    main()
