# v0.6 Experiment Notes - Structural-Sharing Map

## Research question

v0.5 established a stable cold-multiplication crossover against the explicit sparse baseline for a maximally repeated-block family. v0.6 does **not** search for another size crossover. It asks what happens when the problem size is held fixed and deliberate structural sharing is progressively removed.

## Controlled family

For block stride `S = 512`, define 64 blocks. Each block contains exactly 128 exponents selected from a low universe of 256 positions. Thus every `P` contains exactly 8,192 monomials.

A deterministic bank of 64 equal-cardinality templates is used. The first template is deliberately structured (`0:127`); the remaining templates are deterministic pseudo-random-looking 128-of-256 masks.

For a requested template count `T`, block `j` uses template `1 + mod(j-1,T)`. Hence:

- `T=1`: one block pattern repeated 64 times.
- `T=2`: two patterns, each repeated about 32 times.
- ...
- `T=64`: one distinct pattern per block.

The explicit sparse term count is unchanged at every T.

## Fixed multiplier

`Q = sum_{k=0}^{31} x^k`.

This intentionally overlaps the low exponent-bit levels used by P. Unlike the disjoint-bit multiplication used to locate the v0.5 crossover, multiplication must now interact with P's internal BMD structure. The expectation is therefore that BMD multiplication cost will rise as the BMD representation loses sharing, while sparse outer-product cardinality remains approximately fixed.

## Measurements

For every T:

- BMD P reachable node count
- BMD P lower-bound numeric payload
- explicit sparse P bytes
- sparse/BMD compression ratio
- cold BMD multiplication median, Q25, Q75, CV
- sparse multiplication median, Q25, Q75, CV
- median speed ratio
- conservative robust speed interval
- BMD result nodes
- BMD workspace nodes created by multiplication
- multiplication cache entries
- sparse result term count
- numerical agreement at two evaluation points

Seven independent trials are retained in a separate raw-trial CSV.

## Pre-registered interpretation

A useful outcome is not merely a single faster/slower point. We are looking for a **transition region** in which:

1. the sparse workload remains fixed,
2. BMD node count rises strongly as template reuse is removed,
3. BMD cold-multiplication advantage shrinks and eventually disappears,
4. the transition can be associated with both a template-count range and an approximate BMD node-count range.

If BMD remains faster even at 64 unique templates, the experiment will report that no threshold was reached in this family. If sparse is faster even at one template, the prior v0.5 advantage does not transfer to this overlapping-bit workload. Either result is informative.

## Independent validation before MATLAB

The exact-integer Python reference validates all nine sharing levels:

- exactly 8,192 input terms at every level,
- exact BMD product equals independent sparse convolution,
- P BMD nodes increase from 13 at T=1 to 2,288 at T=64 in the reference implementation.

These are structural/correctness checks only, not MATLAB timing claims.
