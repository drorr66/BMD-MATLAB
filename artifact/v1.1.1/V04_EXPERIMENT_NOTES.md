# BMD-MATLAB v0.4 experiment notes

## Decision entering v0.4

v0.3 produced a genuine but qualified signal:

- conventional sparse storage remains the right baseline for isolated monomials and very low term counts;
- BMD shows strong compression when many explicit terms share ordered binary structure;
- canonical equality beat sparse at sufficiently large geometric sums;
- exact memoized repeats can be much faster than recomputing sparse products;
- sequential MATLAB construction was wasting most of its workspace on unreachable nodes.

The next question is therefore structural, not merely “BMD vs MATLAB dense”.

## H1 - structural representation

As block repetition increases, sparse term count should grow roughly with the number of explicit terms while BMD node count should grow much more slowly if the diagram shares the high-bit block index.

Primary measures:

- `sparse_terms_per_bmd_node`
- `sparse_bytes_over_bmd`
- BMD node growth as `blocks` increases

## H2 - cold multiplication

For factors with disjoint bit ranges, generic BMD multiplication may exploit small operand DAGs. If BMD cold multiplication crosses sparse multiplication as structural sharing grows, that is a substantially stronger result than beating dense `conv`.

## H3 - canonical equality

The same polynomial is produced via generic multiplication and via the direct structural builder. Canonical BMD equality is an O(1) reference comparison in one manager, while dense/sparse equality scales with explicit representation size.

## H4 - related cache reuse

The benchmark first computes `A*B`, then changes `B` by extending the block-index geometric sum by one block and computes `A*Bplus`.

`related_cache_speedup > 1` indicates that computed-table entries from the old product accelerate a new, structurally related product. This is a stronger memoization test than timing the identical top-level multiplication twice.

## Builder remediation

`geometricSum(n)` and `geometricGrid(...)` use direct bit-recursive construction. For these families the expected invariant is:

`total_internal_nodes == reachable_nodes`

Any violation fails the v0.4 tests/benchmark row.

## What would count as a strong positive result

A compelling v0.4 signal would combine several of these:

1. BMD nodes grow logarithmically or near-logarithmically while sparse terms grow linearly with repeated blocks.
2. Direct builder eliminates the large garbage-node multiplier seen in v0.3.
3. Cold BMD multiplication becomes competitive with or faster than sparse multiplication for strongly shared cases.
4. BMD equality advantage grows with explicit term count.
5. Related (not identical) multiplication benefits measurably from prior computed-table state.

No single warm-cache number is sufficient by itself.
