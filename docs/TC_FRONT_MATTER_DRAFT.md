# IEEE Transactions on Computers front-matter draft

Status: manuscript framing draft after proof hardening and application validation.

## Working title options

1. **Exact Structural Predictability in *BMD Polynomial Multiplication**
2. **Disjoint-Bit-Block Closure in *BMDs: Exact Node Allocation and Product Size**
3. **Predictable Polynomial Multiplication with *BMDs under Disjoint Bit Blocks**

Preferred working title: **Disjoint-Bit-Block Closure in *BMDs: Exact Node Allocation and Product Size**.

## Abstract

Multiplicative binary moment diagrams (*BMDs) provide a canonical weighted decision-diagram representation for integer polynomials, but their behavior under polynomial multiplication has remained difficult to predict from operand structure alone. This paper identifies a structurally recognizable class of products for which the multiplication closure can be characterized exactly. For univariate integer polynomials `F` and `G`, if `deg(F) < 2^s` and `G(x)=H(x^(2^s))` for some nonconstant `H`, then the operand supports occupy disjoint exponent-bit blocks. We prove an **Exact Disjoint-Bit-Block (DBB) Closure Theorem**: under canonical normalized weighted *BMDs and exact arithmetic, multiplication in an operand-only manager creates exactly `|V(F)|` new nonterminal nodes and the resulting product contains exactly `|V(F)|+|V(G)|` reachable nonterminal nodes. The theorem is independent of the empirical study and gives an a priori structural prediction before multiplication is executed. We evaluate the result on a frozen 60-case benchmark suite and find 14 cases satisfying the DBB condition; both exact predictions hold in all 14 cases. We further validate the theorem on a configuration derived from an independently developed open-source FPGA CIC decimator, where the equivalent FIR stage `S_128(x)` admits the exact factorization `S_8(x)S_16(x^8)`: the theorem predicts three new nodes and a seven-node product DAG, and both values are reproduced exactly in MATLAB R2026a on Windows and MATLAB Online/Linux. A public polyphase FIR design is also used as a negative/control case: its complete quantized filter does not satisfy a nontrivial DBB factorization under the preregistered search, while its natural monomial polyphase branches satisfy only the expected trivial structural cases. These results establish a precise sufficient condition under which *BMD multiplication is not merely compact empirically but structurally predictable.

## Introduction

Polynomial multiplication is a basic operation in symbolic algebra, arithmetic verification, digital-signal-processing analysis, and hardware-oriented computation. Decision-diagram representations can exploit repeated structure in such polynomials, but multiplication may also create intermediate and final graph structure that is hard to infer from the operand diagrams. For a representation intended to support large structured polynomials, this predictability question matters independently of raw running time: can one determine, from the operand structure itself, how much canonical graph structure multiplication must create?

This paper revisits multiplication in multiplicative binary moment diagrams (*BMDs) from that perspective. Rather than asking only whether a particular benchmark compresses well or whether one implementation is faster than another, we ask whether there exists a recognizable class of products for which the canonical closure of multiplication is exactly determined before the operation is performed.

The key observation is that exponent bits can separate the structural roles of the two operands. Let `s >= 1`. If a low-degree factor `F` satisfies

`deg(F) < 2^s`,

then all exponent bits used by `F` lie below bit `s`. If the second factor has the form

`G(x)=H(x^(2^s))`,

then every exponent in `G` is divisible by `2^s`, so the lowest `s` exponent bits are zero. We call this the **disjoint-bit-block (DBB) condition**. In a canonical *BMD, the nonterminal levels reachable from `F` and `G` therefore occupy disjoint level ranges.

This separation turns the recursive multiplication rule into a one-sided closure process over the low-bit block. At each node of `F`, the high-bit operand has no component at that level, so the equal-top-level multiplication branch and its carry term are never invoked in the low block. Canonical normalization then yields a one-to-one correspondence between the reachable nonterminal nodes of `F` and the new low-block product nodes. We formalize this as the **Exact DBB Closure Theorem**. For nonzero `F` and nonconstant `G` satisfying DBB, exact multiplication in an operand-only manager creates

`N_new = |V(F)|`

new nonterminal canonical nodes, while the final product has

`|V(FG)| = |V(F)| + |V(G)|`

reachable nonterminal nodes.

These are structural equalities, not asymptotic tendencies or fitted empirical models. The allocation prediction depends on an operand-only or equivalently fresh manager, because a warm unique table may already contain some product nodes; the final reachable-size identity is a canonical representation statement. The theorem does not claim a machine-independent wall-clock complexity bound. In particular, we deliberately separate exact structural closure from implementation-sensitive recursion and cache behavior.

The theorem arose from a reproducible empirical study but is proved independently of it. We evaluate a frozen 60-case benchmark suite designed to exercise multiplication across sparse, dense, periodic, and structured polynomial families. Fourteen cases satisfy the DBB condition. In all 14, the measured number of newly allocated nodes equals `|V(F)|`, and the measured reachable product size equals `|V(F)|+|V(G)|`. Thus the empirical suite provides exact validation of two theorem predictions rather than evidence from which the theorem is inferred.

To test whether the phenomenon extends beyond the controlled benchmark suite, we use two external application-derived cases. The positive case is derived from a public open-source FPGA sigma-delta converter containing a CIC decimator with hardware-test parameters `R=128`, differential delay `M=1`, and two CIC stages. One FIR-equivalent stage is

`S_128(x)=1+x+...+x^127`,

which admits the exact DBB factorization

`S_128(x)=S_8(x) S_16(x^8)`.

For this independently selected configuration, the theorem predicts `|V(F)|=3`, `N_new=3`, `|V(G)|=4`, and `|V(FG)|=7`. These values, together with exact sparse and *BMD product equality, are reproduced in both MATLAB R2026a on Windows and MATLAB Online/Linux. We do not claim that the upstream RTL internally implements this factorization; the validation concerns the FIR-equivalent polynomial implied by the documented CIC configuration.

The second external case is deliberately conservative. We evaluate a public MathWorks FIR decimation design using a fixed Q15 quantization and preregistered search over candidate DBB splits. No nontrivial whole-filter DBB factorization is found. Its natural polyphase branches satisfy only monomial low-factor cases, which we classify as structural sanity checks rather than strong application evidence. This negative/control result is important because it demonstrates that the DBB criterion is selective: the analysis does not force arbitrary application polynomials into the theorem's scope.

The resulting contribution is therefore not a claim that all polynomial multiplication is simple in *BMDs, nor that DBB is a necessary condition for compact products. It is a precise sufficient condition under which two quantities that are normally discovered only after multiplication—new canonical-node allocation and final reachable product size—are exactly predictable from the operand diagrams.

## Contributions

This work makes four contributions:

1. **A structurally testable sufficient condition for exact multiplication closure.** We define the disjoint-bit-block (DBB) condition, `deg(F)<2^s` and `G(x)=H(x^(2^s))`, which can be checked before multiplication from the operand support structure.

2. **An exact closure theorem for canonical weighted *BMDs.** We prove that, for nonzero `F` and nonconstant `G` satisfying DBB, multiplication in an operand-only manager creates exactly `|V(F)|` new nonterminal canonical nodes and the resulting product has exactly `|V(F)|+|V(G)|` reachable nonterminal nodes. The proof establishes level separation, an injective node correspondence, absence of collisions with either operand DAG, absence of additional low-block nodes, and reachability of the complete `G` subgraph.

3. **Exact empirical validation on a frozen reproducible benchmark suite.** Among 60 benchmark cases, 14 satisfy DBB, and both theorem predictions are exact in all 14 cases. The validation is kept logically separate from the proof and uses the frozen implementation and benchmark artifacts.

4. **External positive and negative/control validation.** We reproduce the exact predictions on a polynomial derived from an open-source FPGA CIC hardware-test configuration (`R=128`, `M=1`, `N=2`) in both Windows and Linux/MATLAB Online environments. We also evaluate a public polyphase FIR design that fails the nontrivial whole-filter DBB test under a fixed search protocol, thereby delimiting rather than broadening the theorem's scope after observing outcomes.

## Claims intentionally excluded from the front matter

The following statements should not appear in the Abstract or Introduction unless separately established later:

- DBB is necessary for compact or efficient *BMD multiplication.
- The number of recursive multiplication subproblems is unconditionally `O(|V(F)|)`.
- The theorem provides a machine-independent runtime bound.
- The open-source CIC RTL internally implements the DBB factorization used for validation.
- One external CIC case establishes that DBB is common across deployed DSP systems.
- The MathWorks whole-filter FIR case is a positive nontrivial application validation.

## Suggested transition into the technical sections

After the Introduction, the paper should proceed in the following order:

1. *BMD preliminaries and canonical normalization.
2. Multiplication rule and exact-arithmetic implementation model.
3. DBB definition and Exact DBB Closure Theorem.
4. Proof.
5. Frozen benchmark methodology and 14-case theorem validation.
6. External application-derived validation: open-source FPGA CIC.
7. Negative/control application case: polyphase FIR.
8. Performance/representation results that remain relevant after the theorem reframing.
9. Discussion, limitations, and threats to validity.
10. Conclusion.
