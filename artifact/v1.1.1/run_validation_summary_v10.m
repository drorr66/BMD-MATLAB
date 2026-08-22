function [summary,familySummary]=run_validation_summary_v10(results,varargin)
%RUN_VALIDATION_SUMMARY_V10 Summarize v1.0 cross-family validation.
p=inputParser; addParameter(p,'SaveResults',true,@(x)islogical(x)||ismember(x,[0 1])); parse(p,varargin{:});
root=fileparts(mfilename('fullpath'));
T=results(strcmp(results.status,'OK'),:);
routes=T.predicted_route; winners=T.actual_bs_winner;
unc=strcmp(routes,'UNCERTAIN'); routed=~unc; correct=false(height(T),1);
for i=1:height(T), if routed(i), correct(i)=strcmp(routes{i},winners{i}); end, end

s=struct();
s.validation_version='v1.0-cross-family-20260820';
s.predictor_frozen='v0.9-static-20260820';
s.cases_total=height(results); s.cases_ok=height(T); s.families=numel(unique(T.family));
s.routed=sum(routed); s.uncertain=sum(unc); s.coverage=sum(routed)/max(1,height(T));
if any(routed), s.routed_accuracy=sum(correct)/sum(routed); else, s.routed_accuracy=NaN; end
s.false_bmd=sum(strcmp(routes,'BMD') & strcmp(winners,'SPARSE'));
s.false_sparse=sum(strcmp(routes,'SPARSE') & strcmp(winners,'BMD'));
s.bmd_operation_wins=sum(strcmp(T.actual_bs_winner,'BMD'));
s.sparse_operation_wins=sum(strcmp(T.actual_bs_winner,'SPARSE'));
s.bmd_robust_wins=sum(strcmp(T.actual_bs_robust,'BMD'));
s.predictor_median_s=medianFinite(T.predictor_median_s);
s.bmd_conversion_median_s=medianFinite(T.bmd_build_median_s);
s.new_node_median_ape=medianFinite(T.new_node_abs_pct_error);
s.new_node_p90_ape=percentileFinite(T.new_node_abs_pct_error,.90);
s.bmd_single_shot_wins_from_sparse=sum(T.bmd_single_shot_speedup_vs_sparse>1);
s.current_router_from_sparse_wins=sum(T.current_router_from_sparse_speedup_vs_sparse>1);
s.median_routed_operation_regret=medianFinite(T.routed_operation_regret);
be=T.break_even_reuses_route_once(isfinite(T.break_even_reuses_route_once));
s.cases_with_finite_bmd_break_even=numel(be); s.median_break_even_reuses=medianFinite(be);
s.dense_eligible_cases=sum(T.dense_eligible~=0);
s.dense_operation_wins=sum(strcmp(T.operation_oracle_winner,'DENSE'));
s.dense_single_shot_wins=sum(strcmp(T.single_shot_oracle_winner,'DENSE'));
s.sparse_single_shot_wins=sum(strcmp(T.single_shot_oracle_winner,'SPARSE'));
s.bmd_single_shot_oracle_wins=sum(strcmp(T.single_shot_oracle_winner,'BMD'));
s.validation_note='v0.9 thresholds frozen; BMD conversion measured from existing sparse term lists; dense is diagnostic only in v1.0';
summary=struct2table(s);

families=unique(T.family,'stable'); fr=repmat(emptyFamily(),numel(families),1);
for j=1:numel(families)
    F=T(strcmp(T.family,families{j}),:); rr=~strcmp(F.predicted_route,'UNCERTAIN'); cc=false(height(F),1);
    for k=1:height(F), if rr(k), cc(k)=strcmp(F.predicted_route{k},F.actual_bs_winner{k}); end, end
    x=emptyFamily(); x.family=families{j}; x.cases=height(F); x.routed=sum(rr); x.coverage=sum(rr)/height(F);
    if any(rr), x.routed_accuracy=sum(cc)/sum(rr); end
    x.false_bmd=sum(strcmp(F.predicted_route,'BMD') & strcmp(F.actual_bs_winner,'SPARSE'));
    x.false_sparse=sum(strcmp(F.predicted_route,'SPARSE') & strcmp(F.actual_bs_winner,'BMD'));
    x.bmd_operation_wins=sum(strcmp(F.actual_bs_winner,'BMD')); x.dense_operation_wins=sum(strcmp(F.operation_oracle_winner,'DENSE'));
    x.median_bmd_conversion_s=medianFinite(F.bmd_build_median_s); x.median_new_node_ape=medianFinite(F.new_node_abs_pct_error);
    b=F.break_even_reuses_route_once(isfinite(F.break_even_reuses_route_once)); x.median_break_even_reuses=medianFinite(b);
    fr(j)=x;
end
familySummary=struct2table(fr);

fprintf('\nv1.0 validation summary\n'); fprintf('=======================\n');
fprintf('Routed %d/%d cases (%.1f%%); routed accuracy %.1f%%; false-BMD=%d; false-sparse=%d\n',s.routed,s.cases_ok,100*s.coverage,100*s.routed_accuracy,s.false_bmd,s.false_sparse);
fprintf('BMD vs sparse operation wins: BMD=%d, sparse=%d; robust BMD wins=%d\n',s.bmd_operation_wins,s.sparse_operation_wins,s.bmd_robust_wins);
fprintf('Median BMD conversion %.3f ms; median predictor %.3f ms\n',1000*s.bmd_conversion_median_s,1000*s.predictor_median_s);
fprintf('BMD single-shot wins from sparse input: %d; current pre-route-from-sparse wins: %d; finite break-even cases: %d\n',s.bmd_single_shot_wins_from_sparse,s.current_router_from_sparse_wins,s.cases_with_finite_bmd_break_even);
fprintf('Dense diagnostic: eligible=%d, operation wins=%d, single-shot wins=%d\n',s.dense_eligible_cases,s.dense_operation_wins,s.dense_single_shot_wins);

if p.Results.SaveResults
    outDir=fullfile(root,'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
    writetable(summary,fullfile(outDir,'validation_summary_v10.csv'));
    writetable(familySummary,fullfile(outDir,'validation_family_summary_v10.csv'));
    writeMetadata(fullfile(outDir,'run_metadata_v10.txt'),s);
end
end

function r=emptyFamily()
r=struct('family','','cases',0,'routed',0,'coverage',NaN,'routed_accuracy',NaN,'false_bmd',0,'false_sparse',0, ...
'bmd_operation_wins',0,'dense_operation_wins',0,'median_bmd_conversion_s',NaN,'median_new_node_ape',NaN,'median_break_even_reuses',NaN);
end
function m=medianFinite(x)
x=x(isfinite(x)); if isempty(x), m=NaN; else, m=median(x); end
end
function q=percentileFinite(x,p)
x=sort(x(isfinite(x))); if isempty(x), q=NaN; return; end
if numel(x)==1, q=x(1); return; end
pos=1+(numel(x)-1)*p; lo=floor(pos); hi=ceil(pos); if lo==hi, q=x(lo); else, q=x(lo)+(pos-lo)*(x(hi)-x(lo)); end
end
function writeMetadata(path,s)
fid=fopen(path,'w'); if fid<0, return; end
c=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'BMD-MATLAB v1.0 cross-family validation\n');
fprintf(fid,'Generated: %s\n',datestr(now,31)); fprintf(fid,'MATLAB version: %s\n',version); fprintf(fid,'Computer: %s\n',computer);
fprintf(fid,'Frozen predictor: v0.9-static-20260820 (no threshold retuning in v1.0).\n');
fprintf(fid,'Validation cases: %d across %d families.\n',s.cases_total,s.families);
fprintf(fid,'BMD conversion timing starts from existing sparse term lists and includes compact().\n');
fprintf(fid,'Break-even assumes conversion + predictor once, followed by repeated BMD multiply versus repeated sparse multiply.\n');
fprintf(fid,'Dense baseline is measured only when result coefficient count <= 250000 and is diagnostic, not routed by v0.9.\n');
fprintf(fid,'UNCERTAIN remains an intentional predictor outcome.\n');
end
