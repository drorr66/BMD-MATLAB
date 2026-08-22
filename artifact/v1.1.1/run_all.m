% BMD-MATLAB v1.1.1 publication validation-only patch
fprintf('\n=== BMD-MATLAB v1.1.1: validation-only patch ===\n');
run_tests();
run_predictor_tests_v09();
run_publication_tests_v111();
[publication_results_v111,publication_trials_v111]=run_publication_validation_v111('Trials',5,'BuildTrials',3,'PredictTrials',5,'DenseMaxCoefficients',250000,'SaveResults',true); %#ok<NASGU,ASGLU>
[publication_summary_v111,publication_family_summary_v111]=run_publication_summary_v111(publication_results_v111,'SaveResults',true); %#ok<NASGU,ASGLU>
write_run_metadata_v111('Trials',5,'BuildTrials',3,'DenseMaxCoefficients',250000);
fprintf('\n=== v1.1.1 complete ===\nResults are in the results folder.\n');
disp(publication_summary_v111);
