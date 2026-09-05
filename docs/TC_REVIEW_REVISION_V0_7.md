# Review-driven manuscript revision v0.7

Date: 2026-09-05
Title retained: **Disjoint-Bit-Block Closure in *BMDs: Exact Node Allocation and Product Size**.

This revision follows the Gemini/Claude reviews and exact-integer checks. The complete edited Word manuscript and its Markdown source are delivered together with a companion verification package. The older `TC_MANUSCRIPT_DRAFT_V0_1.md`, front-matter draft, build plan, and earlier proof documents are historical drafting records; they should not be used as the current manuscript without the corrections below. This note does not claim that the complete v0.7 manuscript source has been committed here.

## Scientific changes

1. The main theorem is now the representation identity `|V(FG)| = |V(F)| + |V(G)|`, under fixed increasing exponent-bit order, gcd/sign normalization, exact integer arithmetic, nonzero F, and nonconstant G satisfying DBB.
2. `N_new = |V(F)|` is a separate corollary for the specified recursive procedure, operand-only initial manager, empty computed caches, and no node deletion during the operation. No unconditional linear recursion/runtime claim is restored.
3. Ordered decision-diagram grafting/composition is explicitly acknowledged as prior work, with Bryant (1986), Section 4.5, including the correction concerning reduction. The paper presents an arithmetic specialization rather than claiming to invent disjoint-support composition.
4. New propositions with proofs: for R=2^n, `|V(S_R)|=n`, while `|V(S_R^2)|=2R-2` in the same fixed order. Counts for R=8,32,128,4096 are tabulated.
5. Interleaved disjoint support is not enough: `F=1+x+x^4+2x^5`, `G=1+3x^2` have sizes 3,1,5. However, block separation is not necessary on individual pairs: `(1+x)(1+x^4)` times `(1+x^2)` has sizes 2,1,3. Both examples are included, preventing an incorrect necessity claim.
6. The external CIC factor is explicitly distinguished from its two-stage ideal square. The former has seven nodes, the latter 254. Neither is a verification of the finite-word-length converter or of the upstream hardware tests.
7. All 14 qualifying cases are presented as implementation-conformance checks, not a separate theorem contribution or a prevalence estimate.
8. The other 46 cases are analyzed using individual operand sizes, not operand-union size. For rho = product_size/(F_size+G_size), six equal one and 40 exceed one; min=1, Q1=2.4077574967405475, median=3.416883793642522, Q3=6.156377551020409, max=31.875. All six equalities are BINOMIAL_DIFFERENCE cases.
9. Four data tables replace long unstructured numerical prose: dyadic family counts; 14 DBB predicted/observed counts; non-DBB growth by family; retained performance summary.
10. Correlation claims without an adequately specified uncertainty model are removed from the main text. Fisher calculations from the previously reported win counts remain explicitly exploratory and conditional, not population or causal evidence.

## Evidence provenance

The new Python reanalysis independently reconstructs all 60 input pairs from `artifact/v1.1.1/build_publication_case_v11.m` at `f21a1570e2587ce096ae021d6dbcd32b53706596`. All 60 operand-union and product-node counts match the archived Linux v1.1.1 structural CSV, blob `4b339591646d0c7875650717db383a1cebaf6449`. The companion package identifies the transcribed archive columns, independent individual-root counts, exact arithmetic code, and resulting 60-row CSV.

The 50-trial timing summary is retained from the existing Research Record/manuscript. It is not recomputed from the older v1.1.1 timing columns. The later archive is `reproducibility/BMD_v15_reproducibility_2026-08-23.zip`; this revision did not newly execute MATLAB or obtain new timing results. Rechecking the complete later timing archive remains an explicit release-checklist task, not a silently asserted completion.

## Figures, wording, and attribution

The erroneous separate color infographic is withdrawn. It was not an embedded figure of v0.6. The valid exponent-separation and canonical-DAG figures are retained; the repetitive provenance/scope flowchart is omitted. The figure split is between x^4 and x^8. The external configuration remains R=128, M=1, N=2. Same-code Windows/MATLAB Online runs are described as cross-environment consistency checks. Dave Muscle retains authorship of the converter. The Professor Shuzo Yajima dedication is retained. Full source revisions remain in the availability section. Formulas are native editable Word equations; references are numbered in first-citation order.

## Diagnostic correction and archival status

The non-frozen application script's operand-union diagnostic was corrected in commit `b0dfcbfcba3331ff415836603f6331e840cd9503`. See `OPEN_SOURCE_CIC_DIAGNOSTIC_ERRATUM.md`. This code correction has not been represented as a new MATLAB execution. No frozen implementation or case-generator file was modified. No pull request was merged, no release was published, and no DOI was minted by this revision.

Remaining submission tasks include final target-format layout, collection of native application-run output files into the archive, explicit software-license approval, and DOI deposition of the approved artifact. The revision is not labeled submission-ready solely because the editorial corrections are implemented.
