function run_publication_tests_v11()
%RUN_PUBLICATION_TESTS_V11 Static/new-helper tests before timing.
fprintf('Running v1.1 publication helper tests...\n');
% Frozen case count and deterministic build.
for i=[1 6 7 14 15 22 23 30 31 38 39 44 45 50 51 56 57 60]
    c=build_publication_case_v11(i); assert(c.p_terms>0 && c.q_terms>0 && c.pair_products==c.p_terms*c.q_terms);
end
% Alternate sparse baseline exactness on representative coefficient/support cases.
for i=[1 15 31 40 47 52 60]
    c=build_publication_case_v11(i); a=sparse_terms_multiply(c.p_sparse,c.q_sparse); b=sparse_terms_multiply_sortreduce_v11(c.p_sparse,c.q_sparse); assert(sparse_terms_same(a,b));
end
% Packed BMD is concrete and non-empty for nonconstant operands.
c=build_publication_case_v11(15); [m,r]=build_bmd_pair_from_sparse_v10(c.p_sparse,c.q_sparse); p=pack_bmd_roots_v11(m,r); assert(workspace_bytes_v11(p)>0); st=m.stats(); assert(numel(p.levels)==st.total_internal_nodes);
fprintf('v1.1 helper tests PASS.\n');
end
