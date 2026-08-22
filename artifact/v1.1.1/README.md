# BMD-MATLAB v1.1.1

Publication validation-only patch for the BMD adaptive polynomial-routing
experiments.

v1.1.1 does **not** change the BMD implementation, predictor, thresholds,
60 frozen workloads, sparse/dense baselines, timing regions, trial counts or
memory methodology from v1.1.  It replaces only the v1.1 numeric validation
gate that evaluated an alternating binomial polynomial near `x=1` and was
susceptible to catastrophic cancellation.

See:

- `V111_VALIDATION_ONLY_PATCH.md` - exact patch scope and rationale
- `V11_PREREGISTERED_CLAIMS.md` - original preregistered v1.1 analysis contract
- `PUBLICATION_READINESS_CHECKLIST_v111.md` - remaining publication gates
- `FROZEN_COMPONENT_HASHES_v111.txt` - byte-identity evidence for frozen components
- `baseline_v11_first_environment/` - verbatim first v1.1 MATLAB results

## Run

```matlab
unzip('BMD_MATLAB_v1_1_1.zip')
cd('BMD_MATLAB_v1_1_1')
run_all
```

Expected result files:

- `results/publication_validation_results_v111.csv`
- `results/publication_validation_trials_v111.csv`
- `results/publication_validation_summary_v111.csv`
- `results/publication_family_summary_v111.csv`
- `results/run_metadata_v111.txt`

After a clean 60/60 first-environment rerun, run this exact ZIP unchanged on
a second MATLAB environment before final publication freeze.
