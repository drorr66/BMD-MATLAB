function ok=run_predictor_tests_v09()
%RUN_PREDICTOR_TESTS_V09 Static predictor unit/contract tests.
ok=false;
fprintf('\nRunning v0.9 predictor contract tests...\n');

% High-confidence bit-cube, fully disjoint ordered bands.
m=BMDManager();
p=m.indicatorExponents(build_bitcube_support_v08(0:7));
q=m.indicatorExponents(build_bitcube_support_v08([8 10 11 12 13 14 15 16]));
[m,r]=m.compact([p;q]);
x=predict_multiply_route_v09(m,r(1,:),r(2,:),256,256);
assert(strcmp(x.prediction_regime,'ORDERED_DISJOINT_BANDS'));
assert(strcmp(x.recommended_route,'BMD'));
assert(abs(x.predicted_new_nodes-8)<1e-12);

% One collision with a gap: held-out from the v0.8 calibration signatures.
m=BMDManager();
p=m.indicatorExponents(build_bitcube_support_v08(0:7));
q=m.indicatorExponents(build_bitcube_support_v08([7 10 11 12 13 14 15 16]));
[m,r]=m.compact([p;q]);
x=predict_multiply_route_v09(m,r(1,:),r(2,:),256,256);
assert(strcmp(x.prediction_regime,'BITCUBE_RIDGE'));
assert(strcmp(x.recommended_route,'BMD'));
assert(x.predicted_new_nodes>10 && x.predicted_new_nodes<15);

% Near-boundary grid: predictor must abstain rather than force BMD.
m=BMDManager(); p=m.geometricSum(255); q=m.geometricSumShifted(143,11);
[m,r]=m.compact([p;q]);
x=predict_multiply_route_v09(m,r(1,:),r(2,:),256,144);
assert(strcmp(x.prediction_regime,'ORDERED_DISJOINT_BANDS'));
assert(strcmp(x.recommended_route,'UNCERTAIN'));

% General heavy-overlap case: conservative low-confidence sparse route.
bank=make_sharing_template_bank_v06(64,256,128);
[exps,~]=build_sharing_case_v06(bank,2,64,9);
m=BMDManager(); p=m.indicatorExponents(exps); q=m.geometricSum(31);
[m,r]=m.compact([p;q]);
x=predict_multiply_route_v09(m,r(1,:),r(2,:),8192,32);
assert(strcmp(x.prediction_confidence,'LOW'));
assert(strcmp(x.recommended_route,'SPARSE'));
assert(x.predicted_new_nodes>64);

fprintf('PASS: predictor performs no multiply and honors conservative routing zones.\n');
ok=true;
end
