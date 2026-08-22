function summary=run_predictor_summary_v09(replay,holdout,varargin)
%RUN_PREDICTOR_SUMMARY_V09 Summarize conservative routing quality.
p=inputParser; addParameter(p,'SaveResults',true,@(x)islogical(x)||ismember(x,[0 1])); parse(p,varargin{:});
root=fileparts(mfilename('fullpath'));
rs=scoreTable(replay); hs=scoreTable(holdout);

s=struct();
s.predictor='v0.9-static-20260820';
s.replay_cases=height(replay); s.replay_routed=rs.routed; s.replay_coverage=rs.coverage; s.replay_routed_accuracy=rs.accuracy;
s.replay_false_bmd=rs.false_bmd; s.replay_false_sparse=rs.false_sparse; s.replay_uncertain=rs.uncertain;
s.replay_new_node_median_ape=medianFinite(replay.new_node_abs_pct_error);
s.holdout_cases=height(holdout); s.holdout_routed=hs.routed; s.holdout_coverage=hs.coverage; s.holdout_routed_accuracy=hs.accuracy;
s.holdout_false_bmd=hs.false_bmd; s.holdout_false_sparse=hs.false_sparse; s.holdout_uncertain=hs.uncertain;
s.holdout_bmd_recommendations=hs.bmd_recs; s.holdout_bmd_precision=hs.bmd_precision;
s.holdout_sparse_recommendations=hs.sparse_recs; s.holdout_sparse_precision=hs.sparse_precision;
s.holdout_new_node_median_ape=medianFinite(holdout.new_node_abs_pct_error);
s.holdout_new_node_p90_ape=percentileFinite(holdout.new_node_abs_pct_error,.90);
s.holdout_predictor_median_s=medianFinite(holdout.predictor_median_s);
ops=min([holdout.bmd_median_s holdout.sparse_median_s],[],2);
rat=holdout.predictor_median_s./max(realmin,ops);
s.holdout_predictor_over_faster_operation_median=medianFinite(rat);
s.synthetic_exact_training_cases=800; s.synthetic_exact_test_cases=200;
s.synthetic_exact_test_median_ape=0.0535381137144821;
s.synthetic_exact_test_p90_ape=0.143389516806096;
s.synthetic_exact_test_max_factor=1.47826496628548;
s.fresh_exact_reference_cases=200; s.fresh_exact_reference_median_ape=0.0508; s.fresh_exact_reference_p90_ape=0.1364; s.fresh_exact_reference_max_ape=0.2604;
s.routing_policy_note='conservative; UNCERTAIN is an allowed outcome; thresholds are MATLAB-v0.5-v0.8 calibration guides';
summary=struct2table(s);

fprintf('\nv0.9 predictor summary\n'); fprintf('======================\n');
fprintf('Calibration replay: routed %d/%d (%.1f%%), accuracy %.1f%%, false-BMD=%d\n',rs.routed,height(replay),100*rs.coverage,100*rs.accuracy,rs.false_bmd);
fprintf('Fresh holdout:      routed %d/%d (%.1f%%), accuracy %.1f%%, false-BMD=%d, false-sparse=%d\n',hs.routed,height(holdout),100*hs.coverage,100*hs.accuracy,hs.false_bmd,hs.false_sparse);
fprintf('Holdout new-node median APE %.1f%%; p90 %.1f%%\n',100*s.holdout_new_node_median_ape,100*s.holdout_new_node_p90_ape);
fprintf('Median predictor overhead %.3f ms\n',1000*s.holdout_predictor_median_s);

if p.Results.SaveResults
    outDir=fullfile(root,'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
    writetable(summary,fullfile(outDir,'predictor_summary_v09.csv'));
    writeMetadata(fullfile(outDir,'run_metadata_v09.txt'),summary);
end
end

function z=scoreTable(T)
routes=T.predicted_route; winners=T.actual_median_winner;
unc=strcmp(routes,'UNCERTAIN'); routed=~unc; correct=false(height(T),1);
for i=1:height(T), if routed(i), correct(i)=strcmp(routes{i},winners{i}); end, end
z=struct(); z.routed=sum(routed); z.uncertain=sum(unc); z.coverage=z.routed/max(1,height(T));
if z.routed>0, z.accuracy=sum(correct)/z.routed; else, z.accuracy=NaN; end
z.false_bmd=sum(strcmp(routes,'BMD') & strcmp(winners,'SPARSE'));
z.false_sparse=sum(strcmp(routes,'SPARSE') & strcmp(winners,'BMD'));
z.bmd_recs=sum(strcmp(routes,'BMD')); z.sparse_recs=sum(strcmp(routes,'SPARSE'));
if z.bmd_recs>0, z.bmd_precision=sum(strcmp(routes,'BMD') & strcmp(winners,'BMD'))/z.bmd_recs; else, z.bmd_precision=NaN; end
if z.sparse_recs>0, z.sparse_precision=sum(strcmp(routes,'SPARSE') & strcmp(winners,'SPARSE'))/z.sparse_recs; else, z.sparse_precision=NaN; end
end
function m=medianFinite(x)
x=x(isfinite(x)); if isempty(x), m=NaN; else, m=median(x); end
end
function q=percentileFinite(x,p)
x=sort(x(isfinite(x))); if isempty(x), q=NaN; return; end
if numel(x)==1, q=x(1); return; end
pos=1+(numel(x)-1)*p; lo=floor(pos); hi=ceil(pos); if lo==hi, q=x(lo); else, q=x(lo)+(pos-lo)*(x(hi)-x(lo)); end
end
function writeMetadata(path,summary)
fid=fopen(path,'w'); if fid<0, return; end
c=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'BMD-MATLAB v0.9 static pre-multiply predictor\n');
fprintf(fid,'Generated: %s\n',datestr(now,31));
fprintf(fid,'MATLAB version: %s\n',version); fprintf(fid,'Computer: %s\n',computer);
fprintf(fid,'Predictor never calls BMDManager.multiply.\n');
fprintf(fid,'Term counts are caller-supplied metadata in v0.9.\n');
fprintf(fid,'Bitcube ridge training: 800 exact synthetic cases; held-out exact test: 200 cases.\n');
fprintf(fid,'Exact synthetic test median APE: %.6f\n',summary.synthetic_exact_test_median_ape(1));
fprintf(fid,'Exact synthetic test p90 APE: %.6f\n',summary.synthetic_exact_test_p90_ape(1));
fprintf(fid,'Routing calibration basis: measured MATLAB v0.5-v0.8 results.\n');
fprintf(fid,'UNCERTAIN is intentional; v0.9 prioritizes avoiding false BMD routes.\n');
end
