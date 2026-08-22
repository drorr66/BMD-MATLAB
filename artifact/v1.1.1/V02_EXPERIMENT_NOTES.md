# BMD-MATLAB v0.2 — experiment notes

Date: 2026-08-20

## Baseline from the first real MATLAB Online run

The v0.1 run completed successfully in MATLAB Online. All supplied benchmark rows returned `OK`, and the Table 2 structural replication matched the published final-node formula in every tested row.

The baseline also exposed the key measurement problem addressed by v0.2: a final BMD may be very small while the construction manager retains many historical/unreachable nodes and computed-table entries. Comparing that entire construction process with native MATLAB vector creation or `conv` obscures the operation-level question.

The original CSVs are stored in `baseline_v01_matlab/`.

## v0.2 measurement model

### Representation track

Measure separately:

- BMD construction time;
- final reachable node count;
- total workspace nodes after construction;
- unreachable/garbage node count;
- workspace/reachable ratio;
- compacted node count and lower-bound payload;
- dense coefficient count and actual MATLAB bytes.

### Operation track

Operands are built before timing.

For `add` and `multiply`:

1. compact the operand roots into a fresh manager **before** the timer starts;
2. start with empty computed caches and no prior result nodes;
3. time only the operation;
4. repeat from another fresh compact manager for the cold median;
5. separately use `timeit` after prewarming the exact operation to measure computed-table reuse.

For `evaluate` and `equality`, prebuilt operands are timed directly with `timeit`.

### Cases selected to expose representation effects

- high-degree monomial addition;
- high-degree monomial multiplication;
- geometric-sum multiplication;
- evaluation of highly compressed monomials/geometric sums;
- canonical equality versus dense coefficient-vector equality.

This deliberately includes both favorable and unfavorable BMD structures. The goal is to find crossover behavior, not to construct a benchmark that always favors BMD.
