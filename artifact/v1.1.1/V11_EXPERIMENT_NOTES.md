# BMD-MATLAB v1.1 - Publication-Readiness Validation

## Purpose

v1.1 is a falsification-oriented validation package. It does not retune the v0.9 predictor and does not modify the BMD multiplication algorithm. It asks whether the conclusions from v0.5-v1.0 survive a larger and more application-oriented workload set and a stronger explicit-sparse baseline.

## Frozen components

The following files are byte-identical to v1.0/v0.9:

- `BMDManager.m`
- `predict_multiply_route_v09.m`
- `bmd_predictor_features_v09.m`

No routing threshold is changed in v1.1.

## Workload set

60 deterministic cases across 9 families:

- BOXCAR_FIR - moving-average/FIR convolution formulas
- CIC_COMB - multirate/CIC/comb polynomial supports
- POLYPHASE_SEPARATED - separated-scale/polyphase or blocked generating kernels
- PERIODIC_MASK - repeated periodic/spectral-mask supports
- SUBSET_SUM - 0-1 knapsack/subset-sum generating functions
- BOUNDED_RESOURCE - bounded knapsack/resource generating functions
- BINOMIAL_DIFFERENCE - binomial and finite-difference transforms
- WEIGHTED_PERIODIC - integer-weighted periodic DSP kernels
- IRREGULAR_CONTROL - deterministic negative controls outside the structured families

54 cases are application-formula cases and 6 are controls. These are canonical deterministic formulas, not claims of external production traces.

## Stronger sparse baseline

BMD is compared against the faster median of two independent explicit sparse implementations:

1. the existing vectorized `unique + accumarray` implementation;
2. a new explicit `sort + run-reduce` implementation.

The selected best sparse baseline is used for BMD-vs-sparse routing accuracy, speedup, regret and break-even calculations.

## Dense baseline

MATLAB `conv` is measured when the result coefficient count is at most 250,000. It remains diagnostic; v0.9 itself is only a BMD-vs-sparse router.

## Memory reporting

v1.1 reports actual MATLAB `whos` bytes for:

- explicit sparse term-list values;
- dense coefficient vectors when eligible;
- a concrete packed numeric BMD DAG containing only nodes reachable from the selected roots.

The packed BMD figure is more concrete than the old theoretical payload lower bound, but it is NOT the live `BMDManager` footprint: `containers.Map` unique-table and computed-cache overhead is deliberately excluded and must not be described as full manager memory.

## Conversion / reuse economics

Sparse-to-BMD conversion includes `compact()`. The package reports single-shot cost and the number of repeated multiplications required to amortize conversion plus one predictor call when BMD operation time beats the strongest sparse baseline.

## Reproducibility

Run the exact same ZIP unchanged on another MATLAB installation for the cross-machine check. `run_metadata_v11.txt` records MATLAB version, platform, trial counts and frozen benchmark settings without recording machine identity.

## Independent validation

`validation/reference_validation_v11.py` reconstructs all 60 cases independently using Python integers and the exact Python BMD reference. It verifies all products coefficient-for-coefficient and emits `validation/expected_structure_v11.csv`.

## Claims v1.1 can test

- whether the frozen v0.9 BMD/sparse router continues to avoid false-BMD choices;
- whether BMD operation wins survive a stronger sparse baseline;
- whether application-formula families exhibit the same structural-closure behavior;
- whether BMD representation compression survives when measured as a concrete packed MATLAB DAG;
- whether conversion amortization remains the limiting factor for single-shot use;
- whether conclusions remain stable on a second MATLAB environment.

v1.1 does not by itself establish a universal polynomial-multiplication theorem or full live-manager memory usage.
