# Reproducibility Notes

This document accompanies the research artifact for:

**Manipulation of Large-Scale Polynomials Using BMDs Revisited: Structural Closure, Representation-State-Aware Evaluation, and Cross-Environment Replication**

## Preserved archive

The complete reproducibility archive is stored at:

`reproducibility/BMD_v15_reproducibility_2026-08-23.zip`

SHA-256:

`b50775f3ac54b8d1eeb30e589e31c24572628a70981d6ce6de2f9ffac1114a69`

The archive is self-contained and includes:

- the MATLAB implementation and deterministic workload builders;
- 60 workloads in nine structural families;
- raw 50-trial operation measurements for Linux and Windows;
- BMD construction, sparse conversion, dense materialization, and FFT precomputation measurements;
- correctness-validation outputs;
- execution-environment metadata;
- per-case recovery data;
- a complete internal SHA-256 manifest.

## Execution environments

1. **Linux / MATLAB Online** — GLNXA64, MATLAB R2026a Update 5
2. **Windows / MATLAB Desktop** — PCWIN64, MATLAB R2026a Update 4

Both environments use the same workload definitions, representation implementations, validation checks, and measurement policy.

## Measurement design

The benchmark compares:

- *BMD multiplication;
- the faster of two independent explicit sparse multipliers;
- direct MATLAB `conv` when eligible;
- FFT_COLD, which recomputes the operand transforms;
- FFT_RESIDENT, which reuses precomputed operand spectra.

Each operation route is measured in 50 trials per workload. Correctness validation remains outside timed regions. Sparse-to-*BMD conversion, dense materialization, and FFT spectrum preparation are measured separately for single-shot accounting.

Each timed *BMD multiplication begins from the two compact operand DAGs in a fresh manager so that product or intermediate nodes from an earlier trial cannot be reused.

## Verification procedure

Verify the downloaded archive before extraction:

```bash
sha256sum BMD_v15_reproducibility_2026-08-23.zip
unzip -t BMD_v15_reproducibility_2026-08-23.zip
```

After extraction, verify the files against `SHA256_MANIFEST.txt` from the archive root.

## Interpretation boundaries

The direct dense route is limited by the study's 250,000-coefficient eligibility threshold. The explicit sparse implementations are controlled reproducible baselines and are not claimed to represent every optimized sparse-polynomial library. Results should be interpreted across both execution environments because several winner identities are platform-sensitive.

## Release checklist

Before public release:

- choose and add an explicit software and data license;
- create a GitHub release containing the preserved archive;
- archive the release with a DOI-providing service if desired;
- add the public repository URL and archival DOI to the manuscript.
