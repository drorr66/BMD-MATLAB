# BMD-MATLAB v0.3 - experiment notes

## v0.2 evidence motivating this run

MATLAB Online v0.2 produced correct results for all measured operations and exact Table 2 node-count replication. Selected measured ratios suggested BMD crossovers against dense MATLAB, including large-degree monomial addition/evaluation and canonical equality.

However, the strongest monomial results are also cases where a one-term sparse polynomial is the obvious conventional alternative. v0.3 adds that baseline before making any claim that the advantage is specifically due to BMD.

## Experimental controls added

### Materialization barrier

Dense/sparse timed builders and binary operations store the complete result in persistent storage. This intentionally creates a visible side effect and retains the full output array/term list. Absolute timings therefore include the assignment/retention boundary, but the benchmark no longer relies on a small scalar checksum as evidence that the entire result was created.

### Sparse term-list control

The sparse representation uses exact integer exponents (`uint64`) and floating coefficients (`double`). It is normalized by exponent and combines duplicate terms. It is not intended as the fastest possible sparse polynomial package; it is a transparent MATLAB baseline for the core representation question.

### Crossover sampling

Quick-run points are concentrated around v0.2 transition regions rather than logarithmically sparse endpoints. `summarize_crossovers_v03` reports a crossover only when measured adjacent points actually bracket ratio 1, and interpolates in log(n)-log(ratio) space.

## Expected interpretation

A result is especially interesting when BMD beats **both** dense and sparse baselines. A monomial case where BMD beats dense but loses to sparse should be interpreted as a dense-representation failure, not evidence of uniquely superior BMD compression.

The geometric-sum family is more diagnostic: it has `N+1` nonzero terms, so ordinary sparse storage is still O(N), while the canonical BMD can have O(log N)-scale reachable structure for this family.
