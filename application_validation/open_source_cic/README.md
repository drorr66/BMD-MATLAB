# Open-source FPGA CIC DBB application validation

This experiment is external to the frozen 60-case benchmark and uses the
hardware-test configuration from the public repository:

`davemuscle/sigma_delta_converters`

Upstream commit inspected: `198bdecf66dd147b26ec9f4196e8bb03c9abfb53`.

The upstream MAX10 hardware-test top-level explicitly sets:

```systemverilog
localparam OSR = 128;
localparam CIC = 2;
...
sigma_delta_adc #(
    .OVERSAMPLE_RATE(OSR),
    .CIC_STAGES(CIC),
    ...
)
```

The upstream `sigma_delta_adc.sv` instantiates `CIC_STAGES` cascaded
integrator/comb stages and decimates by `OVERSAMPLE_RATE`.

For differential delay M=1, one equivalent FIR factor of this CIC decimator is

`S_128(x) = 1 + x + ... + x^127`.

We test the exact DBB factorization

`S_128(x) = S_8(x) * S_16(x^8)`

with `s=3`, stride `2^s=8`.  Thus the low factor has eight terms and degree 7,
while the high factor has sixteen terms at exponents 0,8,...,120.

The factorization is not claimed to be the internal RTL architecture.  It is
an exact factorization of one FIR-equivalent stage implied by the documented
open-source CIC configuration.

## Run

From the BMD-MATLAB repository root:

```matlab
run('application_validation/open_source_cic/run_open_source_cic_dbb.m')
```

The script checks:

- fixed upstream provenance and configuration metadata;
- DBB precondition before multiplication;
- exact sparse identity with `S_128`;
- exact BMD/sparse product equality;
- `N_new = |V(F)|`;
- `|V(FG)| = |V(F)| + |V(G)|`;
- nontriviality of the low factor.

Output files are written under `application_validation/open_source_cic/results/`.
