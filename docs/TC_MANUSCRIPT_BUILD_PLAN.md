# IEEE Transactions on Computers manuscript build plan

Status: approved framing and title lock for manuscript construction.

## Locked working title

**Disjoint-Bit-Block Closure in *BMDs: Exact Node Allocation and Product Size**

This title replaces the earlier `Manipulation of High-Degree Polynomials Using *BMDs Revisited` framing for the TC manuscript. The older title remains part of the historical artifact record and should not be retroactively changed in frozen materials.

## Manuscript objective

Build a 10--11 page IEEE two-column Transactions manuscript, leaving margin below the regular-page ceiling. The central contribution is exact structural predictability under DBB, not general *BMD performance superiority.

## Page and section budget

| Section | Target pages | Purpose |
|---|---:|---|
| Title, Abstract, Index Terms, I. Introduction | 1.3--1.5 | State exact structural problem and contributions |
| II. *BMD Model and Multiplication | 1.0--1.2 | Only preliminaries needed by theorem |
| III. Disjoint-Bit-Block Closure | 1.5--1.8 | DBB definition, theorem, publication proof, corollary |
| IV. Experimental Methodology | 0.8--1.0 | Frozen suite, environments, baselines, measurement protocol |
| V. Structural Validation | 1.0--1.2 | 14/60 exact predictions, anchor case, correlations where relevant |
| VI. External Application Validation | 0.8--1.0 | Open-source FPGA CIC positive case + FIR negative/control |
| VII. Representation and Performance Tradeoffs | 1.2--1.5 | Dense/FFT/sparse/*BMD results without overclaiming |
| VIII. Discussion and Limitations | 0.7--0.9 | Scope, sufficiency vs necessity, implementation limitations |
| IX. Related Work | 0.6--0.8 | *BMD history, symbolic arithmetic/DD context, multirate context as needed |
| X. Conclusion | 0.3--0.4 | Exact result and bounded interpretation |
| References | 0.8--1.0 | IEEE style |

Target total: approximately 10--11 formatted pages.

## Section construction rules

### I. Introduction

Use `docs/TC_FRONT_MATTER_DRAFT.md` as the source draft. Compress application detail so the Introduction motivates rather than duplicates Sections V--VII. Contributions remain four items: DBB condition, exact theorem, frozen 14/60 validation, external positive + negative/control validation.

### II. *BMD Model and Multiplication

Retain only notation necessary for the theorem:

- weighted reference and terminal semantics;
- level/power-variable interpretation `T_j=x^(2^(j-1))` for the univariate case;
- canonical gcd/sign normalization;
- redundant-node elimination;
- unique-table canonicity;
- multiplication decomposition, including the equal-level carry rule and the one-sided unequal-level rule;
- definition of `V(P)` and operand-only manager.

Do not turn this section into an implementation manual. MATLAB-specific details belong in methodology/reproducibility.

### III. Disjoint-Bit-Block Closure

Use `docs/DBB_THEOREM_PUBLICATION_PROOF.md` as authoritative proof source.

Recommended manuscript presentation:

1. Definition: DBB at split `s`.
2. Lemma: level separation.
3. Lemma: canonical node correspondence and noncollision.
4. Theorem: exact allocation and product-DAG size.
5. Corollary: unique support sums under the same bit-block separation.
6. Short scope paragraph: sufficient, not necessary; allocation identity assumes fresh/operand-only manager; no wall-clock complexity theorem.

The theorem equations are the visual center of the paper:

`N_new = |V(F)|`

and

`|V(FG)| = |V(F)| + |V(G)|`.

### IV. Experimental Methodology

Source of truth is the frozen v1.1.1 artifact and Research Record. Preserve:

- 60 deterministic workloads and family counts;
- exact-integer coefficient policy and conservative range checks;
- 50 timed trials after warmup;
- Linux/MATLAB Online R2026a Update 5 and Windows MATLAB R2026a Update 4 environment details;
- sparse baselines, direct convolution, cold FFT, resident-spectrum FFT;
- direct-convolution eligibility cap;
- operand-only compaction before multiplication timing;
- frozen repository commit provenance.

Methodology must distinguish structural exactness tests from timing experiments. CIC/polyphase structural validation does not require 50 timing repetitions.

### V. Structural Validation

Lead with theorem validation, not speed.

Required results:

- 14 of 60 workloads satisfy DBB;
- both theorem predictions exact in 14/14;
- family distribution of the 14 cases;
- `periodic_mask_07` as anchor example;
- every observed *BMD-over-sparse win lies inside DBB, explicitly labeled empirical containment rather than a necessity theorem;
- structural correlation with multiplication time may be reported as secondary evidence, not as part of the theorem.

Prefer one compact table for all DBB cases or a family-level table plus anchor row, depending on final page pressure.

### VI. External Application Validation

Positive case:

- independently selected `davemuscle/sigma_delta_converters` public FPGA configuration;
- upstream pinned commit and source file in reproducibility material;
- documented `R=128`, `M=1`, `N=2`;
- equivalent FIR stage `S_128(x)`;
- exact DBB factorization `S_128(x)=S_8(x)S_16(x^8)`, `s=3`;
- `|V(F)|=3`, `|V(G)|=4`, `N_new=3`, `|V(FG)|=7`;
- exact reproduction on Windows and MATLAB Online/Linux;
- mandatory caveat: the factorization describes the FIR-equivalent transfer polynomial and does not imply that the upstream RTL internally realizes the same decomposition.

Negative/control case:

- public MathWorks 31-tap FIR decimation design, M=8;
- fixed Q15 quantization;
- no nontrivial whole-filter DBB factorization for preregistered `s=1,...,4` search;
- natural polyphase branches formally satisfy DBB only through monomial low factors;
- classify as control/sanity evidence, not a positive nontrivial application result.

### VII. Representation and Performance Tradeoffs

Keep the broader benchmark because it establishes the practical boundary of the theorem.

Core results to preserve:

- numerical routes win 56/60 in each environment;
- oracle winners by method/environment;
- *BMD versus sparse win counts;
- packed *BMD smaller than sparse in 31/60 in each environment;
- no single-shot *BMD wins from sparse-origin inputs;
- DBB does not guarantee fastest execution;
- direct convolution beats *BMD in every DBB case where direct convolution is eligible;
- resident-spectrum FFT represents a favorable reuse scenario and must be described as such;
- sensitivity analysis can be retained if space permits, otherwise move to repository/supplementary record.

### VIII. Discussion and Limitations

Required limitations:

- MATLAB implementation rather than optimized production C++/Rust library;
- two execution environments, not a broad hardware survey;
- deterministic generated workloads rather than production traces;
- exact integer coefficients and bounded exact range;
- packed node-memory estimate excludes live manager tables/caches;
- transparent MATLAB sparse baselines are not claimed to represent best production sparse libraries;
- DBB is sufficient, not necessary;
- application evidence demonstrates external relevance but not prevalence across DSP systems.

### IX. Related Work

Rebuild from verified references. Do not copy the older narrative wholesale. Separate:

- original BMD/*BMD arithmetic-function representation and verification work;
- decision-diagram approaches to arithmetic/symbolic functions;
- sparse/dense/FFT polynomial multiplication background where directly relevant;
- CIC/polyphase references only to establish application provenance, not as a claim of DBB use in conventional DSP implementations.

### X. Conclusion

End on predictability rather than speed:

A recognizable support condition yields exact a priori prediction of canonical *BMD allocation and final product-DAG size. The result is exact on all qualifying frozen benchmark cases and on an independently derived FPGA-CIC transfer polynomial, while a negative FIR control demonstrates selectivity. Performance results delimit the result as a representation niche rather than general superiority.

## Figure/table budget

Aim for no more than 4--5 principal visual objects:

1. DBB level-separation schematic / theorem intuition.
2. Benchmark-family + DBB qualification summary table.
3. Exact theorem-validation table or compact plot.
4. CIC application factorization/result table or diagram.
5. Performance/representation tradeoff figure or table.

Additional sensitivity plots should be omitted from the main manuscript if they push the paper beyond the target page budget.

## Source hierarchy

When conflicts occur, use this priority:

1. Frozen code/data and exact execution outputs.
2. `BMD_TC_Research_Record_v1.0`.
3. `docs/DBB_THEOREM_ADVERSARIAL_AUDIT.md`.
4. `docs/DBB_THEOREM_PUBLICATION_PROOF.md`.
5. Application-validation READMEs/results on the validated branches.
6. Current six-page DBB manuscript for prose/results already verified.
7. Older 14-page manuscript only for recoverable background/reference leads, not for claims that conflict with newer evidence.

## Claim-control lock

Do not claim:

- general *BMD runtime superiority;
- DBB necessity;
- unconditional `O(|V(F)|)` recursive-call complexity;
- that the open-source CIC RTL implements the DBB decomposition internally;
- that the negative polyphase whole-filter case is a positive DBB application;
- memory totals beyond what the packed-node accounting actually measures.

## Immediate manuscript build sequence

1. Typeset title, abstract, index terms, and compressed Introduction.
2. Write Section II from frozen canonical semantics.
3. Insert Section III theorem/proof in final notation.
4. Build the structural-validation tables directly from frozen results.
5. Add CIC + negative/control application section.
6. Reframe performance results as tradeoffs after theorem validation.
7. Rebuild Related Work and bibliography from verified sources.
8. Add Discussion, limitations, conclusion, acknowledgments, and AI-use disclosure as required by target policy.
9. Format in IEEE Transactions two-column template and perform a page-budget pass.
10. Perform final claim-to-source audit before submission.
