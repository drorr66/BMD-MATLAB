function ok=run_validation_tests_v10()
%RUN_VALIDATION_TESTS_V10 Contract tests for the frozen v1.0 validation set.
ok=false; fprintf('\nRunning v1.0 validation contract tests...\n');
for i=1:15
    c=build_validation_case_v10(i);
    assert(c.p_terms>0 && c.q_terms>0);
    assert(c.pair_products==double(c.p_terms)*double(c.q_terms));
    assert(all(c.p_sparse.coefficients==round(c.p_sparse.coefficients)));
    assert(all(c.q_sparse.coefficients==round(c.q_sparse.coefficients)));
end
% Exercise one case from each broad construction path.
for i=[1 3 7 9 11 13 14 15]
    c=build_validation_case_v10(i);
    [m,r]=build_bmd_pair_from_sparse_v10(c.p_sparse,c.q_sparse);
    pred=predict_multiply_route_v09(m,r(1,:),r(2,:),c.p_terms,c.q_terms);
    assert(isfield(pred,'recommended_route'));
    assert(any(strcmp(pred.recommended_route,{'BMD','SPARSE','UNCERTAIN'})));
end
fprintf('PASS: 15 frozen cases; v0.9 router reused without threshold changes.\n');
ok=true;
end
