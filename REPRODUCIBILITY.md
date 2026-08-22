# Reproducibility Notes

This document accompanies the research artifact for the manuscript:

**Manipulation of Large-Scale Polynomials Using BMDs Revisited: Structural Closure, Representation-State-Aware Evaluation, and Cross-Environment Replication**

## Artifact identity

The manuscript identifies the final validation package as **BMD-MATLAB v1.1.1**.

The preserved source ZIP was checked before repository publication and has SHA-256:

`15cdcadb00f9bceb0921e430befc65e6463d271873e0f20a277d1dc09d59844c`

The extracted package is stored under `artifact/v1.1.1/` and contains 145 files. The package was committed without renaming or restructuring files inside the frozen artifact.

The Git commit that first added the frozen artifact and both final result sets is:

`64f7104f46da1e8bfcd64a1cbb7e6906c4b0d3ad`

Git tag `v1.1.1` has been created and verified. It identifies the frozen reproducibility state for this research package prior to GitHub Release and archival publication.

## Final validation design

The final study uses:

- 60 deterministic workloads in nine families;
- five measured trials per multiplication method;
- three trials for sparse-to-BMD conversion;
- rotated execution order across multiplication methods;
- correctness checks outside timed regions;
- a frozen dense-size policy;
- a frozen pre-operation BMD-versus-sparse screening rule.

## Execution environments

The complete validation package was run in two MATLAB R2026a environments:

1. **Linux / MATLAB Online** — GLNXA64, MATLAB R2026a Update 5
2. **Windows / MATLAB Desktop** — PCWIN64, MATLAB R2026a Update 4

The frozen research artifact is shared between the two environments. Final outputs from each execution environment are preserved separately in:

- `results/linux-matlab-online/`
- `results/windows-matlab-desktop/`

Each directory contains:

- `publication_family_summary_v111.csv`
- `publication_validation_results_v111.csv`
- `publication_validation_summary_v111.csv`
- `publication_validation_trials_v111.csv`
- `run_metadata_v111.txt`

## Artifact/result separation

The frozen `artifact/v1.1.1/` directory contains the implementation, deterministic workload construction, validation machinery, historical experiment notes, integrity information, and v1.1.1 publication-validation scripts.

The final CSV and metadata outputs are stored outside that frozen package under `results/` because they are products of the final executions rather than components of the source ZIP itself. This separation is intentional.

## Timing boundaries

Correctness validation is not part of measured timing regions. Operation timing, sparse-to-BMD conversion timing, and diagnostic validation are kept separate in the reported methodology.

## Dense baseline scope

The dense comparison in the manuscript uses MATLAB `conv` only for cases admitted by the study's frozen dense-size policy. The artifact does not claim to benchmark all possible dense methods such as FFT-, NTT-, or multiprecision-based polynomial multiplication.

## Sparse baseline scope

The explicit sparse baseline is the faster median of two independently implemented pair-product sparse multipliers. These are controlled reproducible baselines, not claims to represent the best available sparse polynomial algorithms.

## Validation entry points

The frozen artifact should be treated as immutable. Reproduction and validation should begin with the package's own instructions and readiness material under:

`artifact/v1.1.1/`

Relevant v1.1.1 files include the publication validation runner, publication summary runner, result validator, metadata writer, readiness checklist, and integrity files contained in the package. Follow those package-native instructions rather than editing the frozen files in place.

## Release checklist

Before making this repository public:

- verify the repository checkout against the frozen package identity recorded above;
- verify both preserved final result sets;
- reconstruct the manuscript's aggregate tables and key figures from the preserved outputs;
- choose and add an explicit software/data license;
- verify that Git tag `v1.1.1` points to the intended final documentation state;
- create the GitHub Release for `v1.1.1`;
- archive the release with a DOI-providing service such as Zenodo if desired;
- update the manuscript's Data and Code Availability statement with the public repository URL and archival DOI.

## Change policy

Do not modify files inside `artifact/v1.1.1/`. Any subsequent code changes, corrected experiments, or new validation protocols should be introduced outside the frozen directory and released under a new version identifier.
