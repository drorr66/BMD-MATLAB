# BMD-MATLAB

Reproducibility repository for the study:

**Manipulation of Large-Scale Polynomials Using BMDs Revisited: Structural Closure, Representation-State-Aware Evaluation, and Cross-Environment Replication**

This repository preserves the MATLAB research artifact, deterministic validation material, and the two final result sets used by the manuscript.

## Status

The repository is currently private while the reproducibility package is being finalized. The manuscript identifies the frozen validation package as **BMD-MATLAB v1.1.1**.

The preserved source ZIP for that package was independently checked before repository publication and has SHA-256:

`15cdcadb00f9bceb0921e430befc65e6463d271873e0f20a277d1dc09d59844c`

The extracted frozen artifact is stored under `artifact/v1.1.1/`. It contains 145 files and was committed without renaming or restructuring the files inside the package.

## Repository structure

```text
BMD-MATLAB/
├── artifact/
│   └── v1.1.1/                  frozen BMD-MATLAB research artifact
├── results/
│   ├── linux-matlab-online/     final Linux / MATLAB Online result set
│   └── windows-matlab-desktop/  final Windows / MATLAB Desktop result set
├── figures/                     figure-related material, if added separately
├── docs/                        supplementary reproducibility documentation
├── CITATION.cff
├── REPRODUCIBILITY.md
└── README.md
```

The final result files are preserved separately from the frozen artifact because the artifact contains the code and validation machinery used to generate and check them, while the two `results/` directories contain the preserved outputs from the final execution environments.

## Experimental scope

The final validation set contains 60 deterministic workloads across nine families. The study compares three representation states for polynomial multiplication:

- dense coefficient vectors using MATLAB `conv` when admitted by the frozen dense-size policy;
- explicit sparse term lists, implemented independently in two pair-product forms;
- canonical multiplicative Binary Moment Diagrams (*BMDs).

The two final execution environments are:

- **Linux / MATLAB Online:** MATLAB R2026a Update 5, GLNXA64
- **Windows / MATLAB Desktop:** MATLAB R2026a Update 4, PCWIN64

## Preserved final result files

Each environment directory contains the final publication-oriented outputs:

- `publication_family_summary_v111.csv`
- `publication_validation_results_v111.csv`
- `publication_validation_summary_v111.csv`
- `publication_validation_trials_v111.csv`
- `run_metadata_v111.txt`

## Reproducibility principles

The repository preserves:

- the exact frozen research artifact used for the final validation stage;
- the workload definitions and implementation code contained in that artifact;
- the final raw/publication result sets from both execution environments;
- correctness and validation machinery kept outside the measured timing regions;
- explicit separation between operation timing, sparse-to-BMD conversion timing, and correctness validation;
- integrity metadata and historical validation material retained inside the frozen artifact.

## Versioning and integrity

The manuscript refers to the package as:

`BMD-MATLAB v1.1.1`

Source ZIP SHA-256:

`15cdcadb00f9bceb0921e430befc65e6463d271873e0f20a277d1dc09d59844c`

The Git commit that first added the frozen artifact and the two preserved final result sets is:

`64f7104f46da1e8bfcd64a1cbb7e6906c4b0d3ad`

A Git tag/release `v1.1.1` should be created after the remaining repository checks are complete and before public archival.

## Reproducing and validating the study

Start with the instructions and validation material inside:

`artifact/v1.1.1/`

In particular, consult the package's own run/readiness notes and v1.1.1 publication-validation scripts before modifying any file. The frozen artifact should be treated as immutable; any future experimental changes should be made outside `artifact/v1.1.1/` or in a later versioned artifact.

## Citation

Citation metadata is provided in `CITATION.cff`. The manuscript does not yet have a publication DOI.

## License

No software/data license has been selected yet. A license should be added before the repository is made public.
