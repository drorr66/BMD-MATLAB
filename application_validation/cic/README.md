# CIC DBB application validation

This experiment is separate from the frozen `artifact/v1.1.1` 60-case benchmark.

## Application basis

For a CIC decimator with rate change `R`, differential delay `M`, and `N`
sections, the standard transfer function can be written

`H_CIC(z) = [(1-z^(-R*M))/(1-z^(-1))]^N`

or equivalently

`H_CIC(z) = [sum_{k=0}^{R*M-1} z^(-k)]^N`.

Thus the equivalent response is a cascade of `N` identical length-`R*M`
all-ones FIR stages.

The experiment fixes `R=256`, `M=1`, `N=3` and validates **one** of those
identical equivalent FIR stages. The order `N` does not change the polynomial
of an individual stage.

## Nontrivial DBB case

With `x=z^-1`, one stage is

`S_256(x)=1+x+...+x^255`.

At the natural radix-16 exponent boundary,

`S_256(x) = S_16(x) * S_16(x^16)`.

Therefore

- `F(x)=S_16(x)` has 16 terms and degree 15 `< 16=2^4`;
- `G(x)=S_16(x^16)` has 16 terms and every exponent is a multiple of 16;
- the exact product has 256 terms.

This is a nontrivial DBB multiplication, unlike the monomial low factors in
the basic polyphase FIR branch experiment.

## Scientific interpretation boundary

The factorization is an exact algebraic factorization of the standard CIC
equivalent FIR stage. We **do not** claim that a particular deployed CIC FPGA
core internally implements this radix-16 factorization. The application claim
is about the transfer polynomial arising from the CIC design, not its vendor
microarchitecture.

## Run

From the repository root:

```matlab
run('application_validation/cic/run_cic_dbb.m')
```

No additional MATLAB toolbox is required beyond what the frozen BMD artifact
already requires.

Outputs:

- `application_validation/cic/results/cic_dbb_result.csv`
- `application_validation/cic/results/summary.txt`

The publication gate is positive only if the DBB precondition, exact sparse
identity, exact BMD product, and both theorem predictions all pass.
