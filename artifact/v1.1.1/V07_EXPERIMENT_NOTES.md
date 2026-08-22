# v0.7 Experiment Notes — Operation Closure Map

## Motivation

v0.5 demonstrated a sustained cold-multiplication crossover against sparse for a family whose product remained tiny as a BMD. v0.6 then showed that input compression alone does not predict performance: highly compressed inputs could multiply into a much larger DAG and lose badly to sparse.

The working hypothesis is therefore:

> Cold BMD multiplication is competitive when the **DAG work induced by the operation** is small relative to the explicit work/output, not merely when the operands themselves are compressed.

## Experimental control

Keep both operand cardinalities and BMD operand sizes constant:

- `P = sum_{i=0}^{255} x^i`
- `Q_s = sum_{j=0}^{255} x^(j*2^s)`
- 256 terms in P
- 256 terms in Q
- 65,536 term-pairs in the explicit multiplication
- 8 BMD nodes in P
- 8 BMD nodes in Q

Change only `s`, which moves Q's 8-bit exponent band relative to P's 8-bit band.

## Exact structural reference

The independent reference implementation produced:

| s | overlap bits | result BMD nodes | new BMD nodes | sparse result terms | max coefficient |
|---:|---:|---:|---:|---:|---:|
| 0 | 8 | 510 | 540 | 511 | 256 |
| 1 | 7 | 383 | 410 | 766 | 128 |
| 2 | 6 | 256 | 281 | 1,276 | 64 |
| 3 | 5 | 161 | 184 | 2,296 | 32 |
| 4 | 4 | 98 | 119 | 4,336 | 16 |
| 5 | 3 | 59 | 78 | 8,416 | 8 |
| 6 | 2 | 36 | 53 | 16,576 | 4 |
| 7 | 1 | 23 | 38 | 32,896 | 2 |
| 8+ | 0 | 16 | 8 | 65,536 | 1 |

The endpoint contrast is therefore very strong while the input sizes are fixed.

## Primary question

Does the cold BMD timing fall as `bmd_new_workspace_nodes` / `bmd_result_nodes` fall, and is there a crossover against sparse when the operation becomes sufficiently closed under the BMD representation?

## Candidate empirical selector

If the results are stable, v0.7 will report a bracket rather than prematurely claim a universal formula. Candidate observable features are:

- newly created BMD nodes
- result BMD nodes
- explicit pair-products per newly created BMD node
- sparse result terms per BMD result node

A later experiment should test any threshold on different families before treating it as a general representation-selection criterion.
