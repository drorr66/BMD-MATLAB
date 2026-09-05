# Adversarial audit of the exact DBB closure theorem

Status: proof hardening for the IEEE Transactions on Computers manuscript.

## Candidate theorem under audit

Let F,G in Z[x], F != 0 and G nonconstant. Suppose for some s >= 1,

- deg(F) < 2^s, and
- G(x) = H(x^(2^s)) for a nonconstant H in Z[x].

Let V(P) be the reachable nonterminal-node set of the canonical normalized weighted *BMD for P. Multiplication starts in an operand-only manager containing exactly the nodes reachable from F and G, with exact arithmetic and empty computed caches.

The structural claims are:

1. N_new = |V(F)|, where N_new is the number of unique-table node insertions caused by multiplication; and
2. |V(FG)| = |V(F)| + |V(G)|.

A separate complexity claim was also considered and is audited below.

## Audit result

The two exact node-count claims survive the adversarial audit, subject to the assumptions stated above. The previously proposed O(|V(F)|) bound on the number of distinct recursive multiplication subproblems does **not** follow in full generality from node count alone for a weighted DAG and should not appear in the theorem without an additional hypothesis or a different complexity parameter.

## Proof-hardening lemmas

### Lemma 1: bit-block separation

Because deg(F) < 2^s, every exponent in F has zero binary digits at positions s and above. Hence every nonterminal node reachable from F has a BMD level in {1,...,s}.

Because every exponent of G is divisible by 2^s, every nonterminal node reachable from G has level at least s+1.

Thus the two operand DAGs have disjoint nonterminal level sets.

### Lemma 2: scalar invariance of a normalized weighted node

For any nonzero integer c and nonzero polynomial P, multiplying P by c changes only the incoming edge weight of its canonical normalized weighted *BMD root; it does not require a distinct normalized internal node for cP.

This is the property needed to handle the fact that one stored DAG node may be reached through different incoming weights. In the frozen implementation, makeNode removes the gcd (and fixes a sign convention) before unique-table lookup, so scalar multiples normalize to the same internal node whenever the underlying normalized polynomial is the same.

### Lemma 3: one-sided recursion at every low-bit node

For any reachable normalized node u of F at level j <= s, write its normalized polynomial as

A_u = A_0 + T_j A_1.

A_1 is nonzero; otherwise the node would be eliminated by redundant-node reduction.

At the same level G has decomposition G + T_j*0, because G has no low-bit levels. Therefore the multiplication rule reduces exactly to

A_u G = A_0 G + T_j(A_1 G).

The common-level carry branch containing T_j^2 is never entered in the low-bit block.

In the frozen multiply implementation, the opposite operand splits as b0=G, b1=0, so p01 and p11 are zero and the high branch reduces to p10. No additive or carry node is introduced beyond the transformed low-bit node itself.

### Lemma 4: one normalized product node per normalized F node

Define Phi(u) to be the normalized internal node representing the scalar-equivalence class of A_u G.

If Phi(u)=Phi(v), then canonical weighted normalization implies alpha A_u G = beta A_v G for some nonzero integers alpha,beta. Since Z[x] is an integral domain and G != 0, cancellation gives alpha A_u = beta A_v. Canonicity of the normalized weighted representation then gives u=v. Hence Phi is injective.

Different incoming scalar references to the same normalized node u do not generate additional distinct normalized product nodes by Lemma 2.

### Lemma 5: no collision with either operand DAG

Phi(u) cannot equal a node of F: every node of F represents a polynomial of degree < 2^s, while deg(A_u G) >= deg(G) >= 2^s because G is nonconstant and all nonconstant exponents of G are multiples of 2^s.

Phi(u) cannot equal a node of G: Phi(u) retains top level j <= s because its high cofactor A_1 G is nonzero, whereas every nonterminal node of G has top level >= s+1.

Therefore every first construction of Phi(u) is a genuinely new unique-table insertion in an operand-only manager.

### Lemma 6: no extra node insertions

At a low-bit recursive call the implementation performs only:

- recursive products A_0 G and A_1 G;
- multiplication by zero for the absent low-bit component of G, which terminates immediately;
- additions with zero, which return the nonzero operand directly; and
- one makeNode at the current low-bit level.

At a constant leaf c, multiply(c,G) scales the existing G reference and creates no internal nodes.

No carry variable is created because the equal-top-level branch is never entered while traversing F. Therefore every inserted internal node is one of the Phi(u), and there are exactly |V(F)| of them.

This establishes

N_new = |V(F)|.

### Lemma 7: exact reachable product size

Every Phi(u) is reachable from the product root because every u is reachable in F and the one-sided decomposition preserves the low-bit recursion structure.

Since F != 0, some reachable branch ends at a nonzero constant c. At that leaf multiply(c,G) returns a nonzero scalar reference to the existing root of G, so the whole reachable DAG V(G) occurs beneath the transformed low-bit structure.

The two sets are disjoint by top-level separation. Hence

V(FG) = Phi(V(F)) dot-union V(G),

so

|V(FG)| = |V(F)| + |V(G)|.

## Edge cases tested conceptually

### F is a nonzero constant

Then |V(F)|=0. Multiplication scales G, creates no internal nodes, and the product reachable size is |V(G)|. Both formulas remain valid.

### F = 0

The product is zero and V(G) is not reachable from the result, so the size formula fails. Therefore F != 0 is necessary.

### G is constant

For nonconstant F, multiplication merely scales F and creates no new nodes, so N_new=0 rather than |V(F)|. Therefore G must be nonconstant.

### G is a nonconstant monomial

The theorem still holds. Its nonterminal levels remain >= s+1, the product retains each low-bit level of F, and no carry is introduced.

### Negative coefficients

The structural argument is unchanged. The theorem must assume the canonical sign convention used by normalized weighted nodes, together with exact arithmetic.

### Zero low cofactor

Allowed. A node may have low cofactor zero. The transformed node still survives because its high cofactor is nonzero. Only a zero high cofactor would remove the node, and such a node is already absent from a canonical reduced operand DAG.

### Shared DAG node reached with multiple incoming weights

This does not change the exact unique-node count because scalar multiples normalize to the same internal product node. It does, however, matter for recursion-count complexity; see below.

### Pre-existing non-operand nodes or warm caches

The exact N_new statement is an allocation statement and therefore requires an operand-only manager (or an equivalent freshness condition). If the manager already contains some Phi(u), fewer insertions may occur even though the final reachable product size remains the same.

## Important correction: drop the unconditional O(|V(F)|) recursion bound

A previous proof sketch argued that a binary DAG with n nodes has at most 2n child references and therefore only O(n) recursive multiplication subproblems. That argument is not sufficient for a weighted DAG.

A normalized internal node can be reached with different incoming scalar weights along different paths. splitAt propagates the incoming weight into child references, and the multiplication cache key includes the full weighted reference. Thus the same node ID may participate in several distinct weighted recursive calls. In the worst case, the number of distinct weighted references generated from the root is not proved to be O(|V(F)|) by DAG node count alone.

Safe alternatives for the TC paper:

1. Remove the recursion-complexity claim from the theorem and keep only the exact structural allocation and reachable-size claims; or
2. State complexity in terms of R(F), the number of distinct weighted subreferences encountered by the recursive decomposition of the root reference, giving O(R(F)); or
3. Add a separate hypothesis that each reachable normalized node of F is encountered with only O(1) distinct incoming weights, in which case O(|V(F)|) follows.

Recommendation: use option 1 in the main theorem. The exact node-count result is already strong and clean; a weaker, implementation-sensitive recursion bound would distract from it.

## Final theorem wording recommended for TC

**Exact DBB Closure Theorem.** Let F,G in Z[x], with F != 0 and G nonconstant. Suppose there exists s >= 1 such that deg(F) < 2^s and G(x)=H(x^(2^s)) for some nonconstant H in Z[x]. Consider canonical normalized weighted *BMDs under exact arithmetic, and perform multiplication in an operand-only manager containing exactly the nodes reachable from the two operands. Then the multiplication creates exactly |V(F)| new nonterminal canonical nodes and the resulting product has exactly |V(F)|+|V(G)| reachable nonterminal nodes:

N_new = |V(F)|,

|V(FG)| = |V(F)| + |V(G)|.

No machine-independent wall-clock complexity claim is part of this theorem.

## Publication consequence

The proof is independent of the empirical 14-case suite and of the CIC application cases. Those experiments serve as exact validation of the theorem predictions, not as evidence needed for the theorem itself.
