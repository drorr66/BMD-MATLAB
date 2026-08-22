function [summary,familySummary]=run_publication_summary_v111(results,varargin)
%RUN_PUBLICATION_SUMMARY_V111 Aggregate v1.1.1 validation-only rerun.
p=inputParser; addParameter(p,'SaveResults',true,@(x)islogical(x)||ismember(x,[0 1])); parse(p,varargin{:});
root=fileparts(mfilename('fullpath'));
ok=strcmp(results.status,'OK'); R=results(ok,:);
if height(R)==0, error('BMD:V11Summary','No successful validation rows.'); end
routed=~strcmp(R.predicted_route,'UNCERTAIN');
app=strcmp(R.source_basis,'APPLICATION_FORMULA');
finiteBE=isfinite(R.break_even_reuses_route_once);
packedSmaller=R.bmd_packed_result_bytes<R.sparse_result_bytes;
sortWins=strcmp(R.sparse_best_algorithm,'SORT_REDUCE');
denseEligible=R.dense_eligible==1;
denseWins=strcmp(R.operation_oracle_winner,'DENSE');
bmdWins=strcmp(R.actual_bs_winner,'BMD');
sparseWins=strcmp(R.actual_bs_winner,'SPARSE');
S=struct();
S.validation_version='v1.1.1-validation-only-20260821'; S.predictor_frozen='v0.9-static-20260820';
S.cases_total=height(results); S.cases_ok=height(R); S.application_formula_cases=sum(app); S.control_cases=sum(~app);
S.families=numel(unique(R.family)); S.routed=sum(routed); S.uncertain=sum(~routed); S.coverage=mean(routed);
if any(routed), S.routed_accuracy=mean(R.routed_correct(routed)); else, S.routed_accuracy=NaN; end
S.false_bmd=sum(R.false_bmd); S.false_sparse=sum(R.false_sparse);
S.bmd_vs_best_sparse_operation_wins=sum(bmdWins); S.best_sparse_operation_wins=sum(sparseWins); S.bmd_robust_wins=sum(strcmp(R.actual_bs_robust,'BMD'));
S.sortreduce_selected_cases=sum(sortWins); S.unique_accumarray_selected_cases=sum(~sortWins);
S.predictor_median_s=median(R.predictor_median_s); S.bmd_conversion_median_s=median(R.bmd_build_median_s);
S.new_node_median_ape=median(R.new_node_abs_pct_error); S.new_node_p90_ape=pct(R.new_node_abs_pct_error,.90);
S.bmd_single_shot_wins_from_sparse=sum(R.bmd_single_shot_speedup_vs_sparse_best>1);
S.router_from_sparse_wins=sum(R.router_from_sparse_speedup_vs_sparse_best>1);
S.median_routed_operation_regret=median(R.routed_operation_regret(isfinite(R.routed_operation_regret)));
S.cases_with_finite_bmd_break_even=sum(finiteBE); if any(finiteBE), S.median_break_even_reuses=median(R.break_even_reuses_route_once(finiteBE)); else, S.median_break_even_reuses=NaN; end
S.dense_eligible_cases=sum(denseEligible); S.dense_operation_oracle_wins=sum(denseWins); S.dense_single_shot_oracle_wins=sum(strcmp(R.single_shot_oracle_winner,'DENSE'));
S.bmd_operation_oracle_wins=sum(strcmp(R.operation_oracle_winner,'BMD')); S.sparse_operation_oracle_wins=sum(strcmp(R.operation_oracle_winner,'SPARSE'));
S.packed_bmd_result_smaller_than_sparse_cases=sum(packedSmaller); S.packed_bmd_result_smaller_fraction=mean(packedSmaller);
S.median_sparse_result_over_bmd_packed_bytes=median(R.sparse_result_over_bmd_packed_bytes);
S.max_sparse_result_over_bmd_packed_bytes=max(R.sparse_result_over_bmd_packed_bytes);
S.exact_coefficient_validation_cases=sum(strcmp(R.validation_mode,'DENSE_COEFFICIENT_EXACT'));
S.modular_fingerprint_validation_cases=sum(strcmp(R.validation_mode,'MODULAR_FINGERPRINT_3X'));
S.sparse_crosscheck_passes=sum(R.sparse_crosscheck==1);
S.bmd_coefficient_exact_passes=sum(R.bmd_coefficient_exact==1);
S.bmd_modular_fingerprint_passes=sum(R.bmd_modular_fingerprint==1);
S.dense_crosscheck_passes=sum(R.dense_crosscheck==1);
S.stable_eval_diagnostic_failures=sum(R.stable_eval_diagnostic~=1);
S.validation_patch_note='v1.1.1 changes validation only; predictor/workloads/timing are frozen from v1.1';
S.memory_note='packed BMD is a concrete MATLAB numeric DAG representation; live manager map/cache overhead is excluded';
summary=struct2table(S);

fams=unique(R.family,'stable'); F=repmat(emptyFamily(),numel(fams),1);
for i=1:numel(fams)
    f=fams{i}; X=R(strcmp(R.family,f),:); q=emptyFamily(); q.family=f; q.cases=height(X);
    rr=~strcmp(X.predicted_route,'UNCERTAIN'); q.routed=sum(rr); q.coverage=mean(rr); if any(rr), q.routed_accuracy=mean(X.routed_correct(rr)); end
    q.false_bmd=sum(X.false_bmd); q.false_sparse=sum(X.false_sparse); q.bmd_operation_wins=sum(strcmp(X.actual_bs_winner,'BMD'));
    q.dense_eligible=sum(X.dense_eligible==1); q.dense_oracle_wins=sum(strcmp(X.operation_oracle_winner,'DENSE'));
    q.sortreduce_selected=sum(strcmp(X.sparse_best_algorithm,'SORT_REDUCE')); q.median_bmd_conversion_s=median(X.bmd_build_median_s);
    q.median_new_node_ape=median(X.new_node_abs_pct_error); q.median_sparse_result_over_bmd_packed_bytes=median(X.sparse_result_over_bmd_packed_bytes);
    be=X.break_even_reuses_route_once(isfinite(X.break_even_reuses_route_once)); if ~isempty(be), q.median_break_even_reuses=median(be); end
    F(i)=q;
end
familySummary=struct2table(F);
if p.Results.SaveResults
    outDir=fullfile(root,'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
    writetable(summary,fullfile(outDir,'publication_validation_summary_v111.csv'));
    writetable(familySummary,fullfile(outDir,'publication_family_summary_v111.csv'));
end
end
function q=pct(x,p), x=sort(x(:)); pos=1+(numel(x)-1)*p; lo=floor(pos); hi=ceil(pos); if lo==hi, q=x(lo); else, q=x(lo)+(pos-lo)*(x(hi)-x(lo)); end, end
function r=emptyFamily()
r=struct('family','','cases',0,'routed',0,'coverage',NaN,'routed_accuracy',NaN,'false_bmd',0,'false_sparse',0,'bmd_operation_wins',0, ...
'dense_eligible',0,'dense_oracle_wins',0,'sortreduce_selected',0,'median_bmd_conversion_s',NaN,'median_new_node_ape',NaN, ...
'median_sparse_result_over_bmd_packed_bytes',NaN,'median_break_even_reuses',NaN);
end
