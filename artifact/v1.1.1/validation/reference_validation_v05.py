"""Independent exact validation for the v0.5 fixed structural family.

This validates the mathematical exponent sets and expected structural scale.
It does not validate MATLAB timing or MATLAB parser/runtime behavior.
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
    a = geometric_exponents(inner_n, 1)
    b = geometric_exponents(blocks - 1, stride_power + 1)
    return frozenset(x + y for x in a for y in b)


def main():
    inner_n = 255
    stride_power = 10
    stride = 1 << stride_power
    blocks_list = [96, 128, 160, 192, 256, 384, 512, 1024]

    assert inner_n < stride
    for blocks in blocks_list:
        got = grid_exponents(inner_n, blocks, stride_power)
        expected = frozenset(j * stride + i for j in range(blocks) for i in range(inner_n + 1))
        assert got == expected, blocks
        assert len(got) == 256 * blocks, blocks
    print('PASS v0.5 fixed-family exponent sets for all block counts')

    # For the power-of-two cases, the ideal full-grid decision structure has
    # 8 low-bit nodes plus log2(blocks) high-bit nodes.
    for blocks in [128, 256, 512, 1024]:
        expected_nodes = 8 + (blocks.bit_length() - 1)
        assert expected_nodes in [15, 16, 17, 18]
        assert 256 * blocks / expected_nodes > 1000
    print('PASS v0.5 power-of-two structural-scale checks')


if __name__ == '__main__':
    main()
