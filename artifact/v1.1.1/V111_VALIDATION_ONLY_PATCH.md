# BMD-MATLAB v1.1.1 - Validation-Only Patch

Date: 2026-08-21

## Scope

v1.1.1 is not a new performance experiment and does not retune the router.
It exists only to replace the single ill-conditioned validation gate that
caused `binomial_diff_04` to be reported as `ERROR:BMD:V11Numeric` in the
first v1.1 MATLAB run.

Frozen from v1.1 without changes:

- `BMDManager.m`
- v0.9 predictor and features
- all 60 workload definitions in `build_publication_case_v11.m`
- both sparse multiplication baselines
- dense eligibility cap
- operation/build/predict trial counts
- operation ordering and all timed regions
- routing thresholds and `UNCERTAIN` policy
- memory metric
- preregistered claims/questions

The first v1.1 result set is archived verbatim under
`baseline_v11_first_environment/`, including the one validation error and
all four false-sparse routes.

## Why the v1.1 gate failed

`binomial_diff_04` is

`(1-x)^18 * (1-x)^14 = (1-x)^32`.

v1.1 compared BMD and explicit sparse evaluations at `x=0.999`.  The true
value is approximately `1e-96`, while direct double summation of the 33
alternating coefficients is approximately `-1.046e-7` in an independent
reference calculation.  This is catastrophic cancellation, not evidence of
an incorrect polynomial product.

The independent exact v1.1 reference validator still passes all 60/60
products coefficient-for-coefficient.

## New validation gate

Validation remains outside every timed region.

For dense-eligible results (`result coefficient count <= 250000`):

1. exact equality between the two explicit sparse result term lists;
2. coefficient-for-coefficient equality between `BMDManager.toDense` and
   the explicit sparse result;
3. the existing tolerance-based `conv` versus sparse-result check.

For dense-ineligible results:

1. exact equality between the two explicit sparse result term lists;
2. three deterministic modular fingerprints comparing BMD and sparse
   results.  The moduli are about 1e6, so all modular products are exactly
   representable in MATLAB double arithmetic.

A stable `x=0.5` evaluation is retained only as a diagnostic column and is
not a validation gate.

## New output columns

`publication_validation_results_v111.csv` reports:

- `validation_mode`
- `sparse_crosscheck`
- `bmd_coefficient_exact`
- `bmd_modular_fingerprint`
- `dense_crosscheck`
- `stable_eval_diagnostic`
- `validation_failure_component`
- `numeric_check`

This makes any future validation failure attributable to a specific check
rather than a single opaque Boolean.

## Required sequence

1. Run v1.1.1 in the same MATLAB environment used for v1.1.
2. Require 60/60 successful validation rows.  Do not retune the router.
3. Archive the first-environment v1.1.1 outputs.
4. Run the exact same v1.1.1 ZIP, unchanged, in a second MATLAB environment.
5. Only then freeze publication figures/tables.
