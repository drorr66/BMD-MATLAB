# Disjoint-Bit-Block Closure in *BMDs: Exact Node Allocation and Product Size

**Dror Rotter**  
I.J. Business Do Ltd., Israel  
dr​or@ij-do.com

## Abstract

Multiplicative binary moment diagrams (*BMDs) provide a canonical weighted decision-diagram representation for integer polynomials, but the graph structure produced by polynomial multiplication is generally difficult to predict from the operands alone. This paper identifies a structurally recognizable class for which multiplication closure can be characterized exactly. For univariate integer polynomials `F` and `G`, suppose that for some `s>=1`, `deg(F)<2^s` and `G(x)=H(x^(2^s))` for a nonconstant integer polynomial `H`. The operand supports then occupy disjoint exponent-bit blocks. We prove an Exact Disjoint-Bit-Block (DBB) Closure Theorem: under canonical normalized weighted *BMDs and exact arithmetic, multiplication in an operand-only manager creates exactly `|V(F)|` new nonterminal nodes, and the product contains exactly `|V(F)|+|V(G)|` reachable nonterminal nodes. The theorem is an a priori structural result and is independent of the empirical evaluation. In a frozen 60-case benchmark suite, 14 cases satisfy DBB and both predictions hold exactly in all 14. External validation uses an independently developed open-source FPGA CIC configuration with `R=128`, `M=1`, and two stages. Its FIR-equivalent stage admits `S_128(x)=S_8(x)S_16(x^8)`; the theorem predicts three new nodes and a seven-node product DAG, and both values are reproduced exactly on Windows and Linux/MATLAB Online. A separate public polyphase FIR case provides a negative/control example: its complete quantized filter does not admit a nontrivial DBB factorization under the fixed search used here. Broader timing experiments show that this exact structural predictability defines a representation niche rather than general *BMD wall-clock superiority.

**Index Terms—** decision diagrams, binary moment diagrams, *BMDs, polynomial multiplication, symbolic computation, canonical representations, multirate signal processing.

# I. Introduction

Polynomial multiplication is fundamental to symbolic computation and appears in arithmetic reasoning, signal-processing analysis, and hardware-oriented computation. Compact symbolic representations can exploit regularity that is invisible in a dense coefficient vector, but multiplication may create graph structure that is difficult to infer from the operand representations. For a canonical decision-diagram representation, this raises a question that is distinct from raw execution speed: can the canonical structure created by multiplication be predicted exactly before the operation is performed?

This paper studies that question for multiplicative binary moment diagrams (*BMDs). The emphasis differs from a conventional benchmark comparison. Rather than asking whether *BMD multiplication is generally faster or smaller than sparse, dense, or transform-based alternatives, we seek a recognizable structural class for which the result of canonical multiplication has an exact a priori graph-size characterization.

The key observation is a separation of exponent bits. Let `s>=1`. If

`deg(F)<2^s`,

then every exponent used by `F` lies entirely within the lowest `s` binary positions. If

`G(x)=H(x^(2^s))`,

then every exponent of `G` is divisible by `2^s`, and its lowest `s` exponent bits are zero. We call this the **disjoint-bit-block (DBB) condition**. Under the standard univariate *BMD level ordering, the reachable nonterminal levels of the two operands consequently occupy disjoint blocks.

This separation removes the common-level carry case while the multiplication recursion traverses the low-bit operand. Canonical normalization then induces a one-to-one correspondence between the reachable nonterminal nodes of `F` and newly created low-block nodes of the product. The resulting Exact DBB Closure Theorem states that, for nonzero `F` and nonconstant `G` satisfying DBB, multiplication in an operand-only manager creates exactly

`N_new = |V(F)|`,

while the canonical product has exactly

`|V(FG)| = |V(F)| + |V(G)|`

reachable nonterminal nodes.

These equalities are structural, not fitted performance models. The allocation identity requires an operand-only or equivalently fresh manager because a warm unique table may already contain transformed product nodes. The final reachable-size identity is a property of the canonical representation. We make no machine-independent runtime claim from the theorem.

The theoretical result is evaluated separately from its proof. A frozen suite contains 60 deterministic workloads across periodic, sparse, combinatorial, filter-derived, and irregular polynomial families. Fourteen satisfy DBB. In every qualifying case, both the predicted number of new nodes and the predicted final reachable product size are exact. The experiments therefore validate theorem consequences; they are not used to establish the theorem.

We also test the result outside the controlled workload suite. An independently developed open-source FPGA sigma-delta converter provides a documented CIC configuration with oversampling rate `R=128`, differential delay `M=1`, and two CIC stages. One FIR-equivalent CIC stage is

`S_128(x)=1+x+...+x^127`,

which admits the exact DBB factorization

`S_128(x)=S_8(x)S_16(x^8)`.

Here `s=3`, `|V(F)|=3`, and `|V(G)|=4`; the theorem therefore predicts `N_new=3` and `|V(FG)|=7`. Both predictions, as well as exact polynomial equality, are reproduced in MATLAB R2026a on Windows and MATLAB Online/Linux. This validation concerns the FIR-equivalent transfer polynomial; it does not assert that the upstream RTL internally realizes the same factorization.

A second external example is used deliberately as a control. A public 31-tap polyphase FIR decimation design is quantized to a fixed Q15 representation and tested for nontrivial whole-filter DBB structure at the preregistered splits `s=1,...,4`. No such factorization is found. The natural polyphase branches satisfy only monomial low-factor cases and are treated as structural sanity checks. Thus the external analysis includes both a positive nontrivial case and a case that remains outside the theorem's nontrivial whole-filter scope.

The contributions are fourfold:

1. We define the DBB condition, a support-level sufficient condition that is testable before multiplication.
2. We prove exact canonical closure under DBB: `N_new=|V(F)|` and `|V(FG)|=|V(F)|+|V(G)|`.
3. We validate both predictions exactly in all 14 DBB cases of a frozen 60-workload cross-environment study.
4. We provide external positive and negative/control validation using public DSP designs, including an independently selected FPGA CIC configuration.

The remainder of the paper develops the canonical *BMD model needed for the result, proves exact DBB closure, evaluates its predictions, and then places the structural result in the broader representation and performance tradeoffs observed in the frozen benchmark study.

# II. Canonical *BMD Model and Multiplication

## A. Univariate polynomial decomposition

For level `j>=1`, define the power variable

`T_j = x^(2^(j-1))`.

A polynomial whose exponents are represented in binary can be decomposed at level `j` as

`P = P_0 + T_j P_1`,

where `P_0` collects terms whose `(j-1)`st exponent bit is zero and `P_1` contains the corresponding residual terms whose bit is one. A reduced *BMD represents this decomposition as a directed acyclic graph with weighted references. Terminal references represent integer constants; nonterminal nodes identify a level and two weighted outgoing references.

The implementation studied here uses exact integer edge weights. Node construction applies redundant-node elimination and canonical normalization of outgoing integer weights by their greatest common divisor together with a fixed sign convention. The extracted common factor is carried on the incoming reference. A unique table then ensures that an equal normalized tuple `(level, low reference, high reference)` denotes the same internal node. Consequently scalar multiples of the same normalized subfunction can share an internal node while differing only in their incoming weights.

Let `V(P)` denote the set of nonterminal nodes reachable from the canonical root reference of `P`. This paper counts reachable canonical nodes rather than all nodes that may remain in a manager workspace after unrelated operations.

## B. Multiplication rule

Let

`A=A_0+T A_1`,

`B=B_0+T B_1`

at a common top level with power variable `T`. Their product is

`AB=A_0B_0 + T(A_0B_1+A_1B_0) + T^2 A_1B_1`.

Because `T^2` is the next binary power variable, the last term introduces a carry to the next level. This common-level case is one source of structural growth in general *BMD multiplication.

If only one operand has the current top level, the other operand has zero high cofactor there. For example, if

`A=A_0+T A_1`,  `B=B_0+T*0`,

then

`AB=A_0B_0+T(A_1B_0)`.

No `T^2` term occurs. The DBB result exploits a complete block of levels over which multiplication remains in this one-sided case.

## C. Operand-only manager

The allocation result below refers to multiplication beginning from an operand-only manager: the unique table contains exactly the nodes reachable from the two operand roots, with no product nodes left from previous operations. This isolates nodes created by the multiplication itself. In the frozen experimental protocol, the two operands are built and the manager is compacted to their reachable union before multiplication is timed and structural counters are recorded.

This freshness condition affects the number of newly inserted nodes, not the canonical identity of the final polynomial. A warm manager could reuse a transformed node that happens to have been created earlier and would therefore report a smaller allocation count even though the canonical product DAG is unchanged.

# III. Exact Disjoint-Bit-Block Closure

## A. Definition

**Definition 1 (DBB condition).** Let `F,G in Z[x]`. The ordered pair `(F,G)` satisfies the disjoint-bit-block condition at split `s>=1` if

`deg(F)<2^s`

and

`G(x)=H(x^(2^s))`

for some polynomial `H`. For the closure theorem, `F` is required to be nonzero and `H` (equivalently `G`) nonconstant.

For a sparse representation, the high-block condition can be checked from exponent support: all exponents of `G` must be divisible by `2^s`. Equivalently, for nonconstant `G`, the 2-adic valuation of the gcd of its support exponents is at least `s`.

## B. Level separation

**Lemma 1.** If `(F,G)` satisfies DBB at split `s`, every nonterminal level reachable from `F` lies in `{1,...,s}`, whereas every nonterminal level reachable from `G` is greater than `s`.

**Proof.** Since `deg(F)<2^s`, every exponent occurring in `F` has zero binary digits at positions `s` and above. Hence `F` uses only the first `s` power variables. Conversely, every exponent of `G(x)=H(x^(2^s))` is divisible by `2^s`, so its lowest `s` binary digits are zero. Thus every nonterminal level of `G` lies above `s`. ∎

## C. Canonical node correspondence

**Lemma 2.** For every reachable nonterminal node `u in V(F)`, multiplication by nonzero `G` creates one normalized product node `Phi(u)` at the same low-block level. The map `Phi` is injective, and no `Phi(u)` coincides with a node of either operand DAG.

**Proof.** Let the normalized polynomial represented by `u` be

`A_u=A_0+T A_1`,

where `T` is the power variable at the top level of `u`. Since the DAG is reduced, `A_1` is nonzero; otherwise the node would be redundant. By Lemma 1, `G` has zero high cofactor at this level. Therefore

`A_uG=A_0G+T(A_1G)`.

The common-level carry branch is not entered, and one normalized node at the current level represents the two recursively obtained child products.

Suppose `Phi(u)=Phi(v)` for normalized nodes `u` and `v`. Equality of normalized weighted nodes means that for some nonzero integers `alpha` and `beta`,

`alpha A_uG=beta A_vG`.

Since `Z[x]` is an integral domain and `G` is nonzero, cancellation yields `alpha A_u=beta A_v`. Canonical weighted normalization then implies `u=v`; hence `Phi` is injective. Different incoming scalar weights to the same normalized node do not create additional internal nodes because scalar factors are carried on the incoming reference.

No `Phi(u)` can coincide with a node of `F`: every subfunction represented in `F` has degree below `2^s`, whereas nonconstant `G` has degree at least `2^s`, and therefore `deg(A_uG)>=2^s`. Nor can `Phi(u)` coincide with a node of `G`: the high cofactor `A_1G` is nonzero, so `Phi(u)` retains a top level at most `s`, while every nonterminal node of `G` lies above `s`. ∎

## D. Exact closure theorem

**Theorem 1 (Exact DBB Closure).** Let `F,G in Z[x]`, with `F!=0` and `G` nonconstant. Suppose `(F,G)` satisfies DBB at split `s`. Under canonical normalized weighted *BMDs and exact arithmetic, multiplication performed in an operand-only manager creates exactly

`N_new=|V(F)|`

new nonterminal canonical nodes, and the resulting product has exactly

`|V(FG)|=|V(F)|+|V(G)|`

reachable nonterminal nodes.

**Proof.** By Lemma 2, each `u in V(F)` produces a distinct transformed node `Phi(u)`, and none is present in either operand DAG. Thus at least `|V(F)|` new nodes are inserted.

There are no other insertions. At every low-block level, `G` has zero high cofactor, so the multiplication uses only the child products `A_0G` and `A_1G`, operations involving zero, and one node construction at the current level. The equal-level carry branch is never entered while traversing the low-bit structure of `F`. At a constant leaf `c`, multiplication returns a scalar reference to the existing `G` DAG and creates no internal node. Hence every new nonterminal node is one of the `Phi(u)`, proving `N_new=|V(F)|`.

Every transformed node is reachable from the product root. Because `F!=0`, at least one reachable branch of `F` terminates at a nonzero constant; multiplication at that leaf returns a nonzero scalar reference to the root of `G`, so the entire reachable DAG `V(G)` is contained in the product DAG. By Lemma 1 and Lemma 2, the transformed low-block nodes and `V(G)` are disjoint. Therefore

`V(FG)=Phi(V(F)) dot-union V(G)`,

and `|V(FG)|=|V(F)|+|V(G)|`. ∎

## E. Scope

The theorem gives a sufficient condition, not a necessary characterization of compact *BMD multiplication. Products outside DBB may also have small canonical diagrams or favorable execution time. Conversely, DBB does not imply that *BMD multiplication is the fastest available numerical or symbolic method. The theorem predicts canonical allocation and reachable result size under its stated representation assumptions; it does not assert an unconditional `O(|V(F)|)` bound on distinct weighted recursive calls or on wall-clock time.

# IV. Experimental Methodology

The evaluation uses the frozen v1.1.1 artifact associated with repository commit `f21a1570e2587ce096ae021d6dbcd32b53706596`. The benchmark contains 60 deterministic polynomial-multiplication workloads across nine workload families: BOXCAR_FIR (6), CIC_COMB (8), POLYPHASE_SEPARATED (8), PERIODIC_MASK (8), SUBSET_SUM (8), BOUNDED_RESOURCE (6), BINOMIAL_DIFFERENCE (6), IRREGULAR_CONTROL (6), and WEIGHTED_PERIODIC (4).

All *BMD arithmetic uses exact integer weights subject to a conservative exact-range policy. For each workload, both operands are built in one manager and compacted to their reachable union before multiplication. This realizes the operand-only condition used by Theorem 1. Structural quantities are recorded from exact counters and reachable-node traversal.

Performance comparisons include two explicit sparse multipliers, direct coefficient convolution when the dense coefficient-position count does not exceed 250,000, a cold FFT route, and a resident-spectrum FFT route that excludes reusable forward-spectrum work. The latter is intentionally favorable to repeated-use FFT scenarios and is reported separately. Each timed route is warmed up and then measured over 50 trials.

Two MATLAB R2026a environments are used: MATLAB Online R2026a Update 5 on Linux GLNXA64, and MATLAB Desktop R2026a Update 4 on Windows with an Intel Core i5-1335U and 16 GB of memory. The application-derived structural experiments are deterministic exact validations rather than timing studies and therefore do not use the 50-trial timing protocol.

# V. Structural Validation

## A. DBB qualification and exact predictions

Fourteen of the 60 frozen workloads satisfy the DBB condition. They comprise all eight POLYPHASE_SEPARATED cases, PERIODIC_MASK cases 05--08, SUBSET_SUM case 01, and BOUNDED_RESOURCE case 01. In every one of these 14 cases, both predictions of Theorem 1 are exact: the observed number of multiplication-induced nonterminal nodes equals `|V(F)|`, and the reachable product DAG contains exactly `|V(F)|+|V(G)|` nonterminal nodes.

This 14/14 result is an exact structural validation, not a statistical fit. DBB qualification is determined from the operand structure before multiplication, and the predicted quantities are integer identities checked against the resulting canonical DAG.

## B. Anchor case: periodic_mask_07

The `periodic_mask_07` workload illustrates the scale separation that motivated the structural analysis. One operand contains an eight-position low-bit mask; the second contains 128 terms separated by stride `2^15`. Their product contains 32,768 nonzero terms and has degree 4,169,493, yet its canonical *BMD product contains only 15 reachable nonterminal nodes. The low-bit operand contains eight reachable nodes, and multiplication creates exactly eight new nodes, as predicted by Theorem 1.

On Linux, median *BMD multiplication time for this case is 0.430 ms, compared with 0.875 ms for the best sparse route and 68.398 ms for the resident-spectrum FFT route. On Windows the corresponding medians are 0.494 ms, 0.659 ms, and 95.830 ms. The packed reachable *BMD representation is approximately 431 times smaller than the sparse representation under the accounting used in the frozen study. This example is intentionally not generalized into a universal performance claim; it shows that exact structural closure can coincide with a large support expansion and a very small canonical result.

## C. Structural containment of observed *BMD-over-sparse wins

All observed *BMD-over-sparse multiplication wins in the frozen study occur inside the DBB subset: eight cases on Linux and two on Windows. This is empirical containment, not a converse theorem. DBB is not claimed to be necessary for a *BMD speed advantage, and the absence of observed wins outside DBB in this finite suite does not establish necessity.

Across the full suite, the number of newly created nodes is strongly associated with *BMD multiplication time. On a log-log scale, Pearson correlations are 0.987 on Linux and 0.916 on Windows, with Spearman correlations of 0.975 and 0.913. Partial correlations controlling for sparse pair-product count remain 0.985 and 0.930. These measurements motivate node creation as an implementation-relevant structural quantity, but they are not part of the theorem.

# VI. External Application-Derived Validation

## A. Open-source FPGA CIC configuration

To test the theorem on structure not selected from the benchmark generator, we use the public `davemuscle/sigma_delta_converters` FPGA project, pinned to upstream commit `198bdecf66dd147b26ec9f4196e8bb03c9abfb53`. Its hardware-test configuration documents a CIC decimator with oversampling rate `R=128`, differential delay `M=1`, and `N=2` stages.

For `M=1`, one FIR-equivalent CIC stage has transfer polynomial

`S_128(x)=1+x+...+x^127`.

It admits the exact factorization

`S_128(x)=S_8(x)S_16(x^8)`.

Taking `F=S_8(x)`, `G=S_16(x^8)`, and `s=3` gives `deg(F)=7<8`, while every exponent of `G` is divisible by eight. Both factors are nontrivial. The frozen *BMD implementation gives `|V(F)|=3` and `|V(G)|=4`; Theorem 1 therefore predicts

`N_new=3`,  `|V(FG)|=7`.

Both values are obtained exactly in MATLAB Desktop R2026a on Windows and MATLAB Online/Linux. Exact sparse multiplication and exact *BMD reconstruction also agree with `S_128(x)` in both environments. The configuration was independently selected from public FPGA code rather than constructed as a benchmark instance.

This validation is deliberately limited to the FIR-equivalent transfer polynomial associated with the documented CIC parameters. It does not imply that the upstream SystemVerilog internally implements the DBB factorization.

## B. Polyphase FIR negative/control case

A separate public MathWorks FIR-decimation example provides a control case. The design uses a 31-tap Parks--McClellan FIR filter with decimation factor eight. To satisfy the exact-integer arithmetic requirements of the frozen *BMD implementation, the coefficients are converted to a fixed Q15 integer representation. The frozen coefficient snapshot was independently checked against MATLAB's `firpm` result in MATLAB Online.

The complete quantized filter is tested for exact nontrivial DBB factorizations at splits `s=1,...,4`, corresponding to strides 2, 4, 8, and 16. None is found. Its eight natural polyphase branches do satisfy the formal support separation at stride eight, but their low factors are monomials and fail the preregistered nontriviality gate. All branch-level theorem predictions are nevertheless exact.

We therefore use this example as negative/control evidence rather than as a second positive application claim. It demonstrates that the external-validation procedure can reject a candidate instead of forcing a favorable DBB interpretation.

# VII. Representation and Performance Tradeoffs

The broader benchmark results delimit the practical meaning of exact structural closure. Numerical routes dominate single-shot execution across most of the frozen suite: the direct-convolution or FFT routes win 56 of 60 workloads in each environment. The full oracle winner counts are 27 direct-convolution, 29 resident-spectrum FFT, three *BMD, and one sparse case on Linux; on Windows they are 34, 22, one, and three, respectively.

Against the best explicit sparse multiplier, *BMD multiplication wins eight workloads on Linux and two on Windows. However, when the cost of constructing a *BMD from sparse-origin input is included, no workload produces a single-shot end-to-end *BMD win. This distinction is important: the structural theorem concerns multiplication of operands already represented canonically as *BMDs and should not be interpreted as eliminating representation-conversion cost.

Representation size tells a different story from execution time. Under the packed reachable-node accounting used in the study, the *BMD representation is smaller than the sparse representation in 31 of 60 workloads in both environments. Thus canonical structural compression occurs substantially more often than wall-clock superiority.

DBB itself does not guarantee the fastest route. Direct convolution beats *BMD multiplication in every DBB case for which direct convolution remains eligible under the study's coefficient-position cap. On Linux the DBB cases exhibit a clean observed crossover against the sparse baseline—cases with at most 8,192 sparse pair products lose while cases with at least 16,384 win—but Windows does not show a single equivalent threshold. This reinforces the distinction between the exact structural theorem and machine-dependent performance.

The resident-spectrum FFT baseline represents a favorable reuse scenario because reusable forward transforms are treated as already resident. Conversely, the sparse multipliers are transparent MATLAB implementations rather than claims about the best possible production sparse library. The performance section should therefore be read as a representation-state-aware comparison within the frozen experimental protocol, not as a universal ranking of polynomial multiplication algorithms.

# VIII. Discussion and Limitations

The Exact DBB Closure Theorem isolates a narrow but strong property: under an identifiable exponent-bit separation, two canonical graph quantities that are normally known only after multiplication can be predicted exactly beforehand. The result is useful even when another algorithm is faster, because it characterizes representation growth independently of timing noise and hardware effects.

The condition is sufficient rather than necessary. The benchmark contains products outside DBB that may still compress well, and nothing in the proof excludes other structural closure classes. A broader characterization of compact *BMD multiplication remains open.

The empirical study also has implementation limitations. The evaluated implementation is MATLAB rather than an optimized C++ or Rust library. Only two execution environments are used. Workloads are deterministic formulas rather than production traces. Arithmetic is restricted to exact integer coefficients within a conservative exact range. Packed reachable-node memory estimates exclude live unique tables, caches, and other manager workspace. The explicit sparse baselines are reproducible MATLAB algorithms, not optimized production libraries. The resident-spectrum FFT route assumes favorable reuse of forward spectra.

The external application evidence should likewise be interpreted narrowly. The FPGA CIC case establishes that a nontrivial DBB polynomial arises from parameters of an independently developed, hardware-tested public design; it does not establish prevalence of DBB across deployed DSP systems. The polyphase FIR control further shows that a plausible multirate application need not yield a nontrivial whole-filter DBB factorization.

Finally, the theorem is a structural allocation and canonical-size result. A weighted DAG may reach the same normalized node through different incoming scalar references, and multiplication caches operate on weighted references. We therefore do not infer an unconditional `O(|V(F)|)` bound on distinct recursive calls or a machine-independent wall-clock complexity theorem from the node-count identities alone.

# IX. Related Work

[To be rebuilt from verified primary references before submission. This section will cover the original BMD/*BMD representation and arithmetic-verification literature, related decision-diagram representations for arithmetic functions, polynomial multiplication baselines relevant to the evaluation, and the CIC/polyphase literature needed to establish application provenance. Historical prose from earlier drafts will be reused only after reference-level verification.]

# X. Conclusion

This paper identifies a class of polynomial products for which canonical *BMD multiplication is exactly predictable at the graph-structure level. If `deg(F)<2^s` and `G(x)=H(x^(2^s))` with nonzero `F` and nonconstant `G`, the operands occupy disjoint exponent-bit blocks. Under canonical normalized weighted *BMDs and exact arithmetic, multiplication in an operand-only manager creates exactly `|V(F)|` new nonterminal nodes, and the final product contains exactly `|V(F)|+|V(G)|` reachable nonterminal nodes.

Both predictions hold in all 14 DBB cases identified before multiplication within a frozen 60-workload suite. They also reproduce exactly for an FIR-equivalent polynomial derived from an independently developed open-source FPGA CIC configuration, while a separate public polyphase FIR case remains outside the nontrivial whole-filter DBB class. The broader timing study shows why the result should be interpreted as exact structural predictability rather than general performance superiority. The principal outcome is therefore a provable representation niche: a support condition that converts otherwise operation-dependent *BMD growth into an exact a priori quantity.

# Acknowledgments

[To be finalized. Any dedication or historical acknowledgment should be placed here rather than in the title/author block, consistent with the target publication format.]

# Reproducibility Note

The frozen experimental artifact, exact benchmark configuration, theorem-audit record, application-validation scripts, and pinned external-source provenance are maintained in the associated repository. Final manuscript references will identify the immutable artifact/commit used for each reported result.
