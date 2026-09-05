# Open-source CIC operand-union diagnostic erratum

Date: 2026-09-05
Affected script: `application_validation/open_source_cic/run_open_source_cic_dbb.m`.
Affected recorded source revision: `0cee2405bde36486f1fa26c9d9728418b08473e5`.
Correction commit: `b0dfcbfcba3331ff415836603f6331e840cd9503` on `tc-dbb-proof-audit`.

## Defect

`reachableNodeCount` accepts one 1-by-2 weighted reference, not a matrix of roots. Passing `[p;q]` produced an incorrect `V_before=0`. The reported value must not be interpreted as a true operand-union count.

## Correction

The helper `build_bmd_pair_from_sparse_v10` already compacts to the union of the two operand DAGs. After calling that helper, the valid union count is the compact manager's total internal-node count:

```matlab
st_before = mgr.stats();
V_before = double(st_before.total_internal_nodes);
```

No new manager API or change to the frozen implementation is needed. For the fixed R=128 factor case the expected diagnostic is seven nodes. This expectation is a structural deduction, not a new reported MATLAB measurement.

## Effect on historical results

`V_before` is not used in multiplication, in the `nodes_created` difference that defines `N_new`, in either theorem check, in the exact polynomial checks, or in the final conjunction. The recorded values `V_F=3`, `V_G=4`, `N_new=3`, and `V_product=7` and the successful checks are therefore unaffected by this isolated diagnostic calculation. Original run provenance must remain at the recorded source revision rather than being relabeled with the correction commit.

The correction was checked against the frozen helper and stats API; the revised script has not been executed in MATLAB during this revision. A confirming run and preservation of its native output remain tasks for the final artifact assembly. Preserve old records and this erratum instead of silently editing historical result files.

No files under `artifact/v1.1.1` were modified.
