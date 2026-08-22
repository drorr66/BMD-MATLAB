function run_publication_tests_v111()
%RUN_PUBLICATION_TESTS_V111 Validation-only patch tests.
fprintf('Running v1.1.1 validation-only tests...\n');
run_publication_tests_v11();

% Regression for the sole v1.1 validation failure: binomial_diff_04.
% This is (1-x)^18 * (1-x)^14 = (1-x)^32.  v1.1 evaluated near x=1,
% where direct double coefficient summation is catastrophically ill-conditioned.
c=build_publication_case_v11(48);
assert(strcmp(c.name,'binomial_diff_04'));
[m,r]=build_bmd_pair_from_sparse_v10(c.p_sparse,c.q_sparse);
p=m.multiply(r(1,:),r(2,:));
su=sparse_terms_multiply(c.p_sparse,c.q_sparse);
ss=sparse_terms_multiply_sortreduce_v11(c.p_sparse,c.q_sparse);
dp=sparse_terms_to_dense(c.p_sparse,250000); dq=sparse_terms_to_dense(c.q_sparse,250000);
v=validate_publication_result_v111(m,p,su,ss,dp,dq,true,250000);
assert(v.sparse_crosscheck==1);
assert(v.bmd_coefficient_exact==1);
assert(v.dense_crosscheck==1);
assert(v.numeric_check==1);
assert(strcmp(v.validation_mode,'DENSE_COEFFICIENT_EXACT'));

% Non-dense modular-fingerprint path on a separated periodic case.
c=build_publication_case_v11(27); % periodic_mask_05, result degree > dense cap
assert(strcmp(c.name,'periodic_mask_05'));
[m,r]=build_bmd_pair_from_sparse_v10(c.p_sparse,c.q_sparse);
p=m.multiply(r(1,:),r(2,:));
su=sparse_terms_multiply(c.p_sparse,c.q_sparse);
ss=sparse_terms_multiply_sortreduce_v11(c.p_sparse,c.q_sparse);
v=validate_publication_result_v111(m,p,su,ss,[],[],false,250000);
assert(v.sparse_crosscheck==1);
assert(v.bmd_modular_fingerprint==1);
assert(v.numeric_check==1);
assert(strcmp(v.validation_mode,'MODULAR_FINGERPRINT_3X'));

fprintf('v1.1.1 validation-only tests PASS.\n');
end
