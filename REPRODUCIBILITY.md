# Reproducibility Notes

This document will accompany the research artifact for the manuscript:

**Manipulation of Large-Scale Polynomials Using BMDs Revisited: Structural Closure, Representation-State-Aware Evaluation, and Cross-Environment Replication**

## Artifact identity

The manuscript identifies the final validation package as **BMD-MATLAB v1.1.1**.

Before public release, the repository contents will be checked against the preserved final package and the manuscript's reported SHA-256 value. The GitHub release/tag will be created only after that consistency check.

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

1. Linux / MATLAB Online (GLNXA64, Update 5)
2. Windows / MATLAB Desktop (PCWIN64, Update 4)

The same workload definitions and research code are intended to be preserved here together with the raw outputs from both environments.

## Intended artifact contents

The public artifact should include, at minimum:

- MATLAB source code required for BMD construction, canonicalization, and multiplication;
- explicit sparse multiplication implementations used in the study;
- deterministic workload definitions;
- validation and correctness-check scripts;
- frozen screening-rule implementation and configuration;
- raw result files from both environments;
- scripts or documented procedures used to reconstruct reported tables and figures;
- version and integrity metadata.

## Timing boundaries

Correctness validation is not part of measured timing regions. Operation timing, sparse-to-BMD conversion timing, and diagnostic validation are kept separate in the reported methodology.

## Dense baseline scope

The dense comparison in the manuscript uses MATLAB `conv` only for cases admitted by the study's frozen dense-size policy. The artifact does not claim to benchmark all possible dense methods such as FFT-, NTT-, or multiprecision-based polynomial multiplication.

## Sparse baseline scope

The explicit sparse baseline is the faster median of two independently implemented pair-product sparse multipliers. These are controlled reproducible baselines, not claims to represent the best available sparse polynomial algorithms.

## Release checklist

Before making this repository public:

- verify that committed source and workload files match the final validation package;
- verify both preserved raw result sets;
- verify the manuscript SHA-256 statement against the intended frozen package;
- reconstruct the manuscript's aggregate tables and key figures;
- choose and add an explicit software/data license;
- create and verify Git tag/release `v1.1.1`;
- archive the release with a DOI-providing service such as Zenodo if desired;
- update the manuscript's Data and Code Availability statement with the public repository URL and archival DOI.
