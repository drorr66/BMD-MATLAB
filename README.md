# BMD-MATLAB

Reproducibility repository for the study:

**Manipulation of Large-Scale Polynomials Using BMDs Revisited: Structural Closure, Representation-State-Aware Evaluation, and Cross-Environment Replication**

This repository is intended to contain the MATLAB implementation, deterministic validation workloads, validation scripts, preserved raw result files, and supporting reproducibility material for the study.

## Status

The repository is currently being prepared as the reproducibility package for the final validation artifact. The manuscript reports the frozen package as **BMD-MATLAB v1.1.1**.

The repository will remain private while the artifact is assembled and checked. It can be made public before submission or when the manuscript is released.

## Planned repository structure

```text
BMD-MATLAB/
├── src/                         MATLAB implementation
├── workloads/                   deterministic workload definitions
├── validation/                  correctness and validation scripts
├── results/
│   ├── linux-matlab-online/     preserved Linux / MATLAB Online results
│   └── windows-matlab-desktop/  preserved Windows / MATLAB Desktop results
├── figures/                     scripts/data used to reproduce manuscript figures
├── docs/                        supplementary reproducibility documentation
├── CITATION.cff
├── REPRODUCIBILITY.md
└── README.md
```

## Experimental scope

The final validation set contains 60 deterministic workloads across nine families. The study compares three representation states for polynomial multiplication:

- dense coefficient vectors using MATLAB `conv` when admitted by the frozen dense-size policy;
- explicit sparse term lists, implemented independently in two pair-product forms;
- canonical multiplicative Binary Moment Diagrams (*BMDs).

The paper reports results from two MATLAB R2026a environments: Linux / MATLAB Online and Windows / MATLAB Desktop.

## Reproducibility principles

The published artifact is intended to preserve:

- the exact workload definitions used in the final validation runs;
- the implementation and validation code used for the reported experiments;
- raw outputs from both execution environments;
- enough documentation to reproduce the reported aggregate tables and figures;
- explicit separation between operation timing, representation-conversion timing, and correctness validation.

## Versioning

The manuscript currently identifies the final validation package as:

`BMD-MATLAB v1.1.1`

A release/tag corresponding to that package will be created only after the repository contents have been checked against the preserved research artifact.

## Citation

Citation metadata is provided in `CITATION.cff`. The manuscript is not yet assigned a publication DOI.

## License

No software license has been selected yet. A license will be added before public release.
