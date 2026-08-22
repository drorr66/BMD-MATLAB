# v1.1 Pre-Registered Analysis Contract

This file is written before the MATLAB v1.1 timings are observed. It prevents post-hoc threshold tuning.

## Frozen algorithmic state

- BMD multiplication: unchanged from v1.0.
- Predictor: `v0.9-static-20260820`, unchanged.
- Predictor routing thresholds: unchanged.
- Dense eligibility cap: result coefficient count <= 250,000.
- Sparse comparison: faster median of `unique+accumarray` and `sort+run-reduce`.
- Operation trials: 5 per case.

## Primary questions

1. **Safety:** How many false-BMD and false-sparse routes occur against the stronger sparse baseline?
2. **Coverage:** What fraction of cases receive BMD/SPARSE rather than UNCERTAIN?
3. **Operation niche:** In which application-formula families does BMD beat the strongest sparse baseline on a cold BMD-resident multiplication?
4. **Conversion economics:** When BMD operation time wins, how many reuses are required to amortize sparse-to-BMD conversion plus one predictor call?
5. **Dense competition:** When dense representation is feasible under the frozen cap, which representation wins operation time and single-shot time?
6. **Representation size:** How often is the concrete packed BMD result smaller than the explicit sparse result, and by how much?
7. **Predictor generality:** What are median and p90 errors for predicted new BMD nodes on these new families?

## Interpretation rules

- `false_bmd > 0` is reported prominently; it is not tuned away in v1.1.
- BMD speed claims use `sparse_best_median_s`, not the older single sparse implementation.
- A BMD operation win is distinct from a sparse-input single-shot win.
- Packed BMD bytes are not called full live-manager memory.
- Application-formula cases are not described as production traces.
- Results from the current machine are not called cross-platform until the unchanged ZIP is run on a second MATLAB environment.

## Publication-readiness target (guideline, not a theorem)

A strong outcome would be: low/zero false-BMD routing, useful nonzero coverage of genuine BMD operation wins on application-formula cases, consistent conversion break-even behavior, and qualitatively similar conclusions on a second MATLAB environment. Failure of any of these is scientifically reportable and must remain in the results.
