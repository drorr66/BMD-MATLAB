# Publication-grade proof: Exact DBB Closure Theorem

This document is intended as the manuscript-ready proof source for the IEEE Transactions on Computers version. It incorporates the adversarial audit in `docs/DBB_THEOREM_ADVERSARIAL_AUDIT.md` and deliberately excludes the previously considered unconditional recursion-complexity claim.

## Definition (DBB condition)

For an integer `s >= 1`, a pair of nonzero univariate integer polynomials `(F,G)` satisfies the disjoint-bit-block (DBB) condition at split `s` if

`deg(F) < 2^s`

and

`G(x) = H(x^(2^s))`

for some integer polynomial `H`. In the theorem below, `G` is additionally required to be nonconstant.

Let `V(P)` denote the set of reachable nonterminal nodes in the canonical normalized weighted *BMD of polynomial `P`.

## Theorem (Exact DBB Closure)

Let `F,G in Z[x]`, with `F != 0` and `G` nonconstant. Suppose the pair `(F,G)` satisfies DBB at split `s`, i.e.

`deg(F) < 2^s`

and

`G(x) = H(x^(2^s))`

for some nonconstant `H in Z[x]`.

Assume canonical normalized weighted *BMDs, exact arithmetic, and multiplication performed in an operand-only manager containing exactly the nodes reachable from the two operands, with no pre-existing product nodes relevant to this multiplication.

Then multiplication creates exactly

`N_new = |V(F)|`

new nonterminal canonical nodes, and the product has exactly

`|V(FG)| = |V(F)| + |V(G)|`

reachable nonterminal nodes.

## Lemma 1 (Level separation)

All nonterminal levels reachable from `F` lie in the low block `{1,...,s}`, while all nonterminal levels reachable from `G` lie strictly above it.

### Proof

The *BMD power variables correspond to binary exponent bits. Since `deg(F) < 2^s`, every exponent occurring in `F` has zero binary digits at positions `s` and above. Hence every nonterminal node reachable from `F` is associated with one of the first `s` power variables.

Every exponent occurring in `G(x)=H(x^(2^s))` is divisible by `2^s`; therefore its lowest `s` binary digits are zero. Hence every nonterminal node reachable from `G` is associated with a level strictly above `s`.

Thus the nonterminal level sets of the two operand DAGs are disjoint. QED

## Lemma 2 (One transformed node per low-block node)

For every reachable nonterminal node `u in V(F)`, multiplication by `G` creates exactly one new normalized canonical node `Phi(u)` at the same low-block level as `u`, and distinct nodes of `F` yield distinct transformed nodes.

### Proof

Let the normalized polynomial represented by node `u` be

`A_u = A_0 + T A_1`,

where `T` is the power variable associated with the top level of `u`. Because `u` is a reachable nonterminal node in a reduced canonical DAG, its high cofactor is nonzero; otherwise redundant-node elimination would remove the node.

By Lemma 1, `G` has no occurrence of level `T`, so at this level its decomposition is simply

`G = G + T*0`.

Therefore the multiplication rule reduces to the one-sided identity

`A_u G = A_0 G + T(A_1 G)`.

No equal-top-level branch is entered and no `T^2` carry term is generated in the low block. Once the two child products have been obtained, exactly one normalized node at the current level is needed to represent the result.

Now suppose two distinct normalized nodes `u,v in V(F)` produced the same normalized product node. Equality of normalized weighted nodes means that the represented product polynomials differ by at most a nonzero scalar carried on the incoming reference. Hence for some nonzero integers `alpha,beta`,

`alpha A_u G = beta A_v G`.

Since `Z[x]` is an integral domain and `G != 0`, cancellation gives

`alpha A_u = beta A_v`.

Canonical weighted normalization then implies that `u` and `v` are the same normalized internal node, contradicting `u != v`. Thus the map `Phi:u -> normalized(A_u G)` is injective.

Finally, different incoming scalar weights to the same normalized node of `F` do not create additional normalized internal nodes: scalar multiplication is absorbed by the incoming edge weight under the canonical gcd/sign normalization. Therefore each normalized node of `F` contributes exactly one transformed normalized node. QED

## Lemma 3 (No collision with the operand DAGs)

For every `u in V(F)`, the transformed node `Phi(u)` is distinct from every node already present in `V(F) union V(G)`.

### Proof

First consider collision with `V(F)`. Every polynomial represented by a node of `F` has degree less than `2^s`. Since `G` is nonconstant and all its nonconstant exponents are multiples of `2^s`,

`deg(G) >= 2^s`.

For every nonzero `A_u`,

`deg(A_u G) = deg(A_u)+deg(G) >= 2^s`.

Hence `A_u G` cannot be scalar-equivalent to a polynomial represented by a node of `F`.

Now consider collision with `V(G)`. In the one-sided decomposition

`A_u G = A_0 G + T(A_1 G)`,

the high cofactor `A_1 G` is nonzero because `A_1 != 0`, `G != 0`, and `Z[x]` has no zero divisors. Therefore `Phi(u)` retains the same low-block top level as `u`, which is at most `s`. By Lemma 1, every nonterminal node of `G` has top level greater than `s`. Thus `Phi(u)` cannot coincide with a node of `G`.

Therefore every first construction of `Phi(u)` is a genuinely new unique-table insertion in an operand-only manager. QED

## Proof of the theorem

By Lemma 2, every node `u in V(F)` yields one transformed normalized node `Phi(u)`, and distinct nodes of `F` yield distinct transformed nodes. By Lemma 3, none of these transformed nodes is already present in either operand DAG. Thus multiplication creates at least `|V(F)|` new nonterminal nodes.

It remains to show that no additional nonterminal nodes are created.

At any low-block recursive step, the opposite operand `G` has zero high cofactor at that level. Consequently the multiplication routine performs only the two relevant child products `A_0 G` and `A_1 G`, additions involving zero, and one node construction at the current low-block level. The equal-top-level carry branch is never used while traversing the low-bit structure of `F`.

At a constant leaf `c`, multiplication returns a scalar multiple of the existing reference to `G`; canonical weighted representation absorbs `c` into the incoming edge weight and does not clone the internal DAG of `G`. Multiplication by zero similarly returns immediately. Therefore every newly inserted nonterminal node is one of the transformed nodes `Phi(u)`. Hence

`N_new = |V(F)|`.

Now count the nodes reachable from the product root. Because every node of `F` is reachable from the root of `F`, and the one-sided recursion preserves the low-block DAG structure, every transformed node `Phi(u)` is reachable from the root of `FG`. Thus the product contains exactly `|V(F)|` transformed low-block nodes.

Since `F != 0`, at least one reachable path through its DAG terminates at a nonzero constant `c`. At that leaf, multiplication returns a nonzero scalar reference to the existing root of `G`. Therefore the full reachable DAG `V(G)` occurs below the transformed low-block structure.

The two node sets are disjoint: each `Phi(u)` has top level at most `s`, while every nonterminal node of `G` has top level greater than `s`. Hence

`V(FG) = Phi(V(F)) dot-union V(G)`,

which gives

`|V(FG)| = |V(F)| + |V(G)|`.

QED

## Remarks for manuscript placement

1. The theorem is structural, not a wall-clock complexity theorem.
2. `F != 0` is necessary: if `F=0`, the product is zero and the `V(G)` term disappears.
3. `G` must be nonconstant: if `G` is constant, multiplication merely scales `F`, so the allocation identity does not hold for nonconstant `F`.
4. The operand-only-manager assumption is required for the allocation count `N_new`; a warm manager may already contain some transformed nodes. The final canonical reachable-size identity is representation-level and does not depend on allocation history.
5. Negative integer coefficients are allowed under the same exact canonical normalization and sign convention.
6. The theorem also covers nonconstant monomial `G`, provided its exponent support lies in the high block.

## Suggested manuscript compression

For the main paper, Lemma 1 can be stated inline, while Lemmas 2 and 3 can be combined into a single `node correspondence` lemma. The full expanded form above should be retained in the research record or supplementary material so every proof obligation remains explicit.
