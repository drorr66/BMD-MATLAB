# BMD-MATLAB

Reproducibility repository for the study:

**Manipulation of Large-Scale Polynomials Using BMDs Revisited: Structural Closure, Representation-State-Aware Evaluation, and Cross-Environment Replication**

This repository preserves the MATLAB implementation, deterministic workload definitions, validation material, and complete Linux and Windows result sets used by the manuscript.

## Current reproducibility bundle

The complete self-contained archive is:

`reproducibility/BMD_v15_reproducibility_2026-08-23.zip`

SHA-256:

`b50775f3ac54b8d1eeb30e589e31c24572628a70981d6ce6de2f9ffac1114a69`

The archive contains the implementation, 60 deterministic workloads, raw 50-trial operation measurements, construction and conversion measurements, validation outputs, environment metadata, recovery data, and a complete SHA-256 manifest.

## Experimental scope

The study compares five operation routes across 60 workloads in nine families:

- canonical multiplicative Binary Moment Diagrams (*BMDs);
- two independent explicit sparse implementations;
- direct MATLAB `conv` when eligible;
- FFT convolution with fresh transforms;
- FFT convolution with resident operand spectra.

Each operation route was measured in 50 trials in both environments:

- **Linux / MATLAB Online:** MATLAB R2026a Update 5, GLNXA64
- **Windows / MATLAB Desktop:** MATLAB R2026a Update 4, PCWIN64

## Repository structure

```text
BMD-MATLAB/
├── reproducibility/
│   └── BMD_v15_reproducibility_2026-08-23.zip
├── artifact/                    earlier preserved artifact
├── results/                     earlier preserved result layout
├── docs/
├── figures/
├── CITATION.cff
├── REPRODUCIBILITY.md
└── README.md
```

## Verification

After downloading the archive:

```bash
sha256sum BMD_v15_reproducibility_2026-08-23.zip
unzip -t BMD_v15_reproducibility_2026-08-23.zip
```

The expected SHA-256 value is shown above. The internal `SHA256_MANIFEST.txt` verifies every preserved file after extraction.

## Citation and release status

Citation metadata is provided in `CITATION.cff`. The repository remains private while release metadata, licensing, and archival publication are finalized.
