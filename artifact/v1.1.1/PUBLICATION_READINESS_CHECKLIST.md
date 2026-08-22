# Publication Readiness Checklist

## Completed before v1.1 MATLAB run

- [x] Core BMD arithmetic independently reference-validated.
- [x] Dense, explicit sparse and BMD comparisons established in prior versions.
- [x] Cold-operation crossover isolated under fixed explicit workload.
- [x] Pre-multiply v0.9 router frozen before v1.0/v1.1 validation.
- [x] v1.0 cross-family validation completed on 15 cases / 9 families.
- [x] v1.1 enlarged to 60 frozen cases / 9 families.
- [x] Second independent explicit-sparse algorithm added.
- [x] Concrete packed MATLAB representation-byte metric added.
- [x] Conversion and reuse break-even included.
- [x] Analysis questions preregistered before v1.1 timings.
- [x] Independent exact Python validation passed 60/60 v1.1 products.

## Required after first v1.1 run

- [ ] Review all 60 rows; no failed numeric checks.
- [ ] Report false-BMD / false-sparse without retuning.
- [ ] Compare old sparse baseline with stronger sort-reduce baseline.
- [ ] Analyze results by application family, not only aggregate accuracy.
- [ ] Report conversion-aware single-shot and reuse break-even separately.
- [ ] Report packed-BMD memory honestly as packed representation, not live-manager footprint.

## Required before final submission

- [ ] Run the exact same ZIP unchanged on a second MATLAB environment/machine.
- [ ] Compare qualitative winners, router errors, speedup ordering and crossover/reuse conclusions across environments.
- [ ] Decide whether an optional external CAS/library baseline is needed for the target venue; do not claim state-of-the-art multiplication performance without it.
- [ ] Freeze final figures/tables only after both environments are archived.
- [ ] Make code, case definitions, raw trial CSVs and environment metadata available as supplementary material where venue rules allow.

## Scope discipline

The intended publishable claim is an adaptive-representation result for structured polynomial workloads, not a claim that BMD generally replaces MATLAB `conv`, sparse multiplication, FFT/NTT, or computer-algebra systems.
