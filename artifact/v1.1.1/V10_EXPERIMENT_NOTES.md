# BMD-MATLAB v1.0 — Cross-Family Validation Contract

## Question

Can the static v0.9 pre-multiply predictor generalize beyond its calibration families, and is BMD still economically useful after conversion cost is included?

## Frozen decision rule

`predict_multiply_route_v09.m` is reused unchanged. No v1.0 result is used to alter its coefficients or routing thresholds before scoring.

This is essential: v1.0 is validation, not another calibration pulse.

## New case families

The frozen 15-case suite contains:

- **BOXCAR:** overlapping moving-average/FIR polynomial products.
- **BLOCK_GENERATING:** digit-separated/blocked generating functions.
- **COMB:** arithmetic-stride polynomials representative of comb and multirate kernels.
- **PERIODIC_MASK:** repeated low-bit masks, with both overlapping and separated scales.
- **IRREGULAR_CONTROL:** deterministic pseudo-random sparse supports as a low-sharing control.
- **BINOMIAL:** repeated first-order factors with non-unit coefficients.
- **FINITE_DIFFERENCE:** signed binomial finite-difference stencil.
- **DENSE_CONTROL:** small dense signed-coefficient control.
- **WEIGHTED_PERIODIC:** periodic integer-coefficient kernel.

The cases are deterministic and are created without using MATLAB timing information.

## What is timed

For each case:

1. Existing sparse term lists are the input baseline.
2. Sparse -> BMD conversion is timed three times; conversion includes final `compact()`.
3. The frozen v0.9 predictor runs before any multiply/warm-up.
4. Cold BMD multiply, sparse term-list multiply, and (when bounded) dense `conv` are timed five times with rotating order.
5. Exact/numerical consistency is checked after timing.

Dense materialization is attempted only when the result coefficient-vector length is <= 250,000.

## Key outputs

### Operation-only / BMD-resident

- BMD versus sparse median time.
- v0.9 predicted route and routed accuracy.
- false-BMD and false-sparse counts.
- actual new BMD nodes versus predicted new nodes.

### Sparse-input economics

`bmd_single_shot_from_sparse_s = conversion + predictor + BMD multiply`

This answers whether converting only for one multiplication makes sense.

If BMD multiply is faster than sparse multiply, v1.0 also reports:

`break_even_reuses_route_once = ceil((conversion + predictor) / (sparse_time - BMD_time))`

This estimates how many repeated operations are needed to amortize the one-time conversion/prediction cost.

### Dense diagnostic

For bounded degrees, v1.0 also measures coefficient-vector construction and MATLAB `conv`. Dense is not yet part of the v0.9 route; it is included to prevent a false BMD-vs-sparse victory when standard dense MATLAB is actually faster.

## Caveats

- Pure MATLAB BMD node/object/hash overhead remains implementation-specific.
- Sparse baseline is a simple `(exponent, coefficient)` term-list implementation, not Symbolic Math Toolbox.
- The v0.9 predictor still needs an existing BMD DAG plus caller-supplied term counts. v1.0 therefore does not solve cheap routing directly from arbitrary sparse input; it measures the consequence of that limitation.
- BMD payload-byte metrics are lower bounds and are not MATLAB heap usage.
- A successful v1.0 would justify building a three-way adaptive engine; it would not establish a universal mathematical threshold.

## Representation-availability caveat

The operation-only predictor score assumes both BMD and sparse representations are available for the timed alternative.  A system that stores only BMD would also need to account for BMD -> sparse conversion when the route says `SPARSE`; v1.0 does not implement that reverse conversion.  Conversely, for sparse-only input, the current v0.9 predictor cannot run until a BMD has already been built.  `current_router_from_sparse_s` therefore exposes that architectural cost explicitly rather than hiding it.
