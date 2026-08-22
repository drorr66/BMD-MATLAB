# Publication Readiness Checklist - v1.1.1

## Evidence already frozen

- [x] v0.9 predictor frozen before v1.0/v1.1 validation.
- [x] v1.1 preregistered analysis contract written before timings.
- [x] v1.1 first-environment run archived verbatim.
- [x] v1.1 result: 60 attempted, 59 passed the old numeric gate.
- [x] v1.1 result: zero false-BMD and four false-sparse routes retained.
- [x] Independent exact reference validation passes 60/60 products.
- [x] Root cause of the single v1.1 validation error isolated to the
      ill-conditioned `x=0.999` evaluation of `binomial_diff_04`.

## v1.1.1 first-environment gate

- [ ] Run exact v1.1.1 ZIP without edits.
- [ ] 60/60 rows have `status=OK` and `numeric_check=1`.
- [ ] `binomial_diff_04` passes `bmd_coefficient_exact`.
- [ ] No performance/routing threshold is changed after seeing results.
- [ ] Archive all five v111 result files and the package SHA-256.
- [ ] Compare qualitative v1.1 vs v1.1.1 performance conclusions; validation
      patch must not be presented as a performance improvement.

## Second-environment replication gate

- [ ] Run the exact same v1.1.1 ZIP unchanged on a different MATLAB
      environment/machine.
- [ ] Archive raw trials and environment metadata.
- [ ] Compare BMD/sparse/dense winners by family.
- [ ] Compare false-BMD / false-sparse / UNCERTAIN outcomes.
- [ ] Compare conversion break-even ordering and the separated-scale niche.
- [ ] Explain any machine-dependent crossover shifts without retuning the
      frozen first-environment results.

## Publication freeze

- [ ] Freeze final dataset and figures only after both environments are archived.
- [ ] Decide venue-specific external CAS/library baseline separately.
- [ ] Keep packed-BMD memory explicitly labeled as packed numeric DAG bytes,
      excluding live `containers.Map`/cache overhead.
- [ ] Do not claim BMD generally replaces dense convolution, sparse
      multiplication, FFT/NTT or CAS algorithms.
