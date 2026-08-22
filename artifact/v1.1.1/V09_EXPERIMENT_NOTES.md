# BMD-MATLAB v0.9 — Static Pre-Multiply Predictor

## Research question

Can MATLAB decide **before multiplication** whether a polynomial product should be routed to BMD or to the explicit sparse baseline, using only the existing operand DAGs plus cheap metadata?

v0.9 is deliberately a **predictor experiment**, not another search for a crossover point. The prediction function never calls `BMDManager.multiply`.

## What v0.8 established

For the calibrated 256×256-term family, BMD multiplication cost tracked newly-created BMD nodes closely. The measured median crossover was near 30 new nodes and the conservative quartile crossover was near 22–23 nodes. Earlier v0.5 data also showed that the same 8-new-node BMD operation only became worthwhile once the explicit sparse workload was large enough.

Therefore a practical router needs both:

1. an estimate of **operation closure** (new BMD work), and
2. an estimate of **explicit sparse work**.

## v0.9 predictor contract

`predict_multiply_route_v09(mgr,p,q,pTerms,qTerms)` receives:

- the two canonical BMD operands,
- explicit term counts `pTerms`, `qTerms` supplied as metadata.

It returns:

- structural features,
- predicted new BMD nodes,
- prediction regime/confidence,
- `pairs_per_predicted_new_node`,
- route: `BMD`, `SPARSE`, or `UNCERTAIN`.

`UNCERTAIN` is intentional. v0.9 prioritizes avoiding a false BMD route over forcing 100% coverage.

### Important limitation

v0.9 does **not** yet infer sparse term counts cheaply from an arbitrary BMD. A production adaptive engine will need term-count metadata, a safe estimator, or a representation object that keeps this metadata as the polynomial evolves.

## Structural regimes

### 1. Ordered disjoint bands — high confidence

If the reachable BMD levels of the two operands are disjoint and one entire band precedes the other, there is no same-level carry collision. This includes the successful v0.5 grid regime.

### 2. Bit-cube chains — high confidence

For operands of the form

`prod_{b in B} (1 + x^(2^b))`

the canonical operand BMD is a simple equal-branch chain. v0.9 detects this topology and uses a ridge predictor trained only on **exact structural new-node counts**, not MATLAB timings.

Training corpus:

- 1,000 deterministic exact-reference bit-cube products,
- 800 model-training cases,
- 200 model-test cases.

Model-test result from the build process:

- median absolute percentage error ≈ 5.35%,
- 90th percentile ≈ 14.34%,
- worst multiplicative error ≈ 1.48×.

A second, separately-seeded 200-case exact reference validation produced:

- median APE ≈ 5.08%,
- 90th percentile ≈ 13.64%.

### 3. General heavy overlap — low confidence

For large, overlapping DAGs outside the calibrated bit-cube topology, v0.9 uses a conservative structural estimate and **will not automatically route to BMD**. It returns SPARSE only when the predicted DAG work is clearly high; otherwise it may abstain.

## Provisional routing policy

The router combines predicted new nodes with explicit pair products `pTerms*qTerms`.

The policy was distilled from the measured MATLAB v0.5–v0.8 data and is **not a universal theorem**:

- clear high DAG work / low work-per-node → SPARSE,
- sufficiently large explicit workload and high work-per-predicted-node → BMD,
- near the crossover or outside calibrated structures → UNCERTAIN.

The constants are intentionally conservative and should be re-calibrated if the implementation moves to MEX/C++ or if the sparse baseline changes.

## Two validation layers

### Calibration replay

`run_predictor_replay_v09` rebuilds the operands from v0.5–v0.8 and applies the predictor without multiplication, then compares predictions to the already-measured MATLAB outcomes. This is **in-sample replay** and is not claimed as independent validation.

### Fresh MATLAB holdout

`run_predictor_holdout_v09` defines 11 cases that were not in the v0.5–v0.8 calibration tables:

- five new bit-cube signatures,
- three new repeated-grid sizes,
- three new sharing-map template counts.

For each case:

1. build and compact operands,
2. run and freeze the predictor,
3. only then warm up multiplication,
4. perform 7 alternating cold BMD/sparse timing trials,
5. compare predicted route with actual winner,
6. compare predicted vs actual new-node count.

This holdout is the main scientific result of v0.9.

## Success criteria

The most important failure is **false BMD routing**: automatically choosing BMD when sparse is actually faster. Desired outcome:

- zero or very few false-BMD holdout decisions,
- useful nonzero routing coverage,
- high accuracy on cases where the predictor does not abstain,
- low predictor overhead relative to the operation,
- acceptable new-node prediction error in high-confidence regimes.

## Next step if v0.9 succeeds

v1.0 should validate the router on broader polynomial families and add a third candidate, dense MATLAB `conv`. Only after the routing criterion is stable should a MEX/C++ BMD backend be treated as the primary performance optimization.
