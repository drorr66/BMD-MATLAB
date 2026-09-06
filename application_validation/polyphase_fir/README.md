# Polyphase FIR DBB application validation

This experiment is intentionally **outside** the frozen `artifact/v1.1.1` benchmark.
It does not modify the 60 frozen workloads or any v1.1.1 implementation file.

## Purpose

Test whether the exact disjoint-bit-block (DBB) closure result from the revised
*BMD paper appears in a public, application-level multirate FIR design.

The source design is the MathWorks **FIR Decimation for FPGA** example:

```matlab
decimFactor = 8;
coeffs = firpm(30,[0 0.1 0.2 0.5]*2,[1 1 0 0]);
inputVecSize = 4;
```

Source page:
https://www.mathworks.com/help/dsphdl/ug/fir-decim-hdl.html

The 31 floating-point FIR coefficients are quantized to signed Q15 integers
(`round(coeffs*2^15)`) because the frozen BMD implementation requires exact
integer edge weights.

## What is tested

1. **Whole-filter nontrivial DBB search.**  For `s=1..4`, the script tests
   whether the quantized FIR polynomial can be factored exactly as

   `P(x) = F(x) * G(x^(2^s))`, with `deg(F) < 2^s`.

   The test is performed as an exact rank-one test on the coefficient matrix
   indexed by residue modulo `2^s`.  This is a strict test; failure is a valid
   negative result and is not repaired by changing the filter.

2. **Natural polyphase-branch DBB validation.**  For the eight natural
   polyphase branches

   `x^k * H_k(x^8),  k=0,...,7`,

   the script builds the two operands in a compact operand-only BMD manager,
   multiplies them, and checks the theorem predictions

   `N_new = |V(F)|`

   and

   `|V(FG)| = |V(F)| + |V(G)|`.

3. **Nontriviality gate.**  The script reports whether the low-bit factor has
   at least two terms and at least two reachable BMD nodes.  A monomial phase
   shift is structurally valid DBB but is deliberately labelled *trivial* for
   publication purposes.

4. **Exact product validation.**  Every BMD product is checked against the
   frozen explicit sparse multiplier.

## Run

From the repository root in MATLAB R2026a with Signal Processing Toolbox:

```matlab
run('application_validation/polyphase_fir/run_polyphase_fir_dbb.m')
```

Outputs are written to:

- `application_validation/polyphase_fir/results/polyphase_branch_results.csv`
- `application_validation/polyphase_fir/results/whole_filter_dbb_search.csv`
- `application_validation/polyphase_fir/results/summary.txt`

## Interpretation rule fixed before execution

- If a **nontrivial whole-filter DBB factorization** is found, it can serve as a
  strong external application validation candidate.
- If only the eight monomial polyphase branches satisfy DBB, report them as an
  application-level structural sanity check, **not** as strong evidence that
  DBB is common in real FIR filters.
- If any exact theorem prediction fails, do not use the application result in
  the paper until the discrepancy is explained.

This experiment is designed to permit a negative result without changing the
benchmark, filter, quantization rule, or DBB criterion after seeing outcomes.
