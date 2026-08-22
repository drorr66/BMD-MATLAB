function replay=run_predictor_replay_v09(varargin)
%RUN_PREDICTOR_REPLAY_V09 Replay static predictions against v0.5-v0.8 data.
%
% This is intentionally labeled IN-SAMPLE / calibration replay. It performs
% no multiplication. The independent holdout benchmark is the scientific
% test of routing generalization.
p=inputParser; addParameter(p,'SaveResults',true,@(x)islogical(x)||ismember(x,[0 1])); parse(p,varargin{:});
root=fileparts(mfilename('fullpath'));
rows=repmat(emptyRow(),38,1); pos=0;

% v0.5
T=readtable(fullfile(root,'baseline_v05_matlab','cold_crossover_results_v05.csv'));
for i=1:height(T)
    blocks=T.blocks(i); m=BMDManager(); p0=m.geometricSum(255); q0=m.geometricSumShifted(blocks-1,11); [m,r]=m.compact([p0;q0]);
    pr=predict_multiply_route_v09(m,r(1,:),r(2,:),256,blocks);
    pos=pos+1; rows(pos)=fillRow('v0.5',sprintf('blocks_%d',blocks),pr,T.bmd_new_workspace_nodes(i),T.ratio_median(i),cellText(T.robust_winner,i));
end

% v0.6
T=readtable(fullfile(root,'baseline_v06_matlab','sharing_map_results_v06.csv'));
bank=make_sharing_template_bank_v06(64,256,128);
for i=1:height(T)
    k=T.unique_templates(i); [exps,~]=build_sharing_case_v06(bank,k,64,9);
    m=BMDManager(); p0=m.indicatorExponents(exps); q0=m.geometricSum(31); [m,r]=m.compact([p0;q0]);
    pr=predict_multiply_route_v09(m,r(1,:),r(2,:),8192,32);
    pos=pos+1; rows(pos)=fillRow('v0.6',sprintf('templates_%d',k),pr,T.bmd_new_workspace_nodes(i),T.ratio_median(i),cellText(T.robust_winner,i));
end

% v0.7
T=readtable(fullfile(root,'baseline_v07_matlab','operation_closure_results_v07.csv'));
pExps=build_bitcube_support_v08(0:7);
for i=1:height(T)
    s=T.shift_power(i); qExps=build_bitcube_support_v08(s:(s+7));
    m=BMDManager(); p0=m.indicatorExponents(pExps); q0=m.indicatorExponents(qExps); [m,r]=m.compact([p0;q0]);
    pr=predict_multiply_route_v09(m,r(1,:),r(2,:),256,256);
    pos=pos+1; rows(pos)=fillRow('v0.7',sprintf('shift_%d',s),pr,T.bmd_new_workspace_nodes(i),T.ratio_median(i),cellText(T.robust_winner,i));
end

% v0.8
T=readtable(fullfile(root,'baseline_v08_matlab','threshold_calibration_results_v08.csv'));
for i=1:height(T)
    ov=T.overlap_bit(i);
    if ov<0, qBits=9:16; else, qBits=[ov 9:15]; end
    qExps=build_bitcube_support_v08(qBits);
    m=BMDManager(); p0=m.indicatorExponents(pExps); q0=m.indicatorExponents(qExps); [m,r]=m.compact([p0;q0]);
    pr=predict_multiply_route_v09(m,r(1,:),r(2,:),256,256);
    pos=pos+1; rows(pos)=fillRow('v0.8',sprintf('new_%d',T.target_new_nodes(i)),pr,T.bmd_new_workspace_nodes(i),T.ratio_median(i),cellText(T.robust_winner,i));
end

replay=struct2table(rows(1:pos));
if p.Results.SaveResults
    outDir=fullfile(root,'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
    writetable(replay,fullfile(outDir,'predictor_replay_v09.csv'));
end
end

function r=fillRow(family,caseName,pred,actualNew,ratio,robustWinner)
r=emptyRow(); r.family=family; r.case_name=caseName; r.calibration_replay=1;
r.predicted_new_nodes=pred.predicted_new_nodes; r.actual_new_nodes=actualNew;
r.new_node_abs_pct_error=abs(pred.predicted_new_nodes-actualNew)/max(1,actualNew);
r.pair_products=pred.pair_products; r.pairs_per_predicted_new_node=pred.pairs_per_predicted_new_node;
r.prediction_regime=pred.prediction_regime; r.prediction_confidence=pred.prediction_confidence; r.predicted_route=pred.recommended_route;
r.actual_ratio_median=ratio; if ratio>1, r.actual_median_winner='BMD'; else, r.actual_median_winner='SPARSE'; end
r.actual_robust_winner=robustWinner;
if strcmp(r.predicted_route,'UNCERTAIN'), r.routed_correct=NaN; else, r.routed_correct=double(strcmp(r.predicted_route,r.actual_median_winner)); end
end

function r=emptyRow()
r=struct('family','','case_name','','calibration_replay',1,'predicted_new_nodes',NaN,'actual_new_nodes',NaN, ...
    'new_node_abs_pct_error',NaN,'pair_products',NaN,'pairs_per_predicted_new_node',NaN, ...
    'prediction_regime','','prediction_confidence','','predicted_route','', ...
    'actual_ratio_median',NaN,'actual_median_winner','','actual_robust_winner','','routed_correct',NaN);
end
function s=cellText(v,i)
if iscell(v), s=v{i}; elseif isstring(v), s=char(v(i)); elseif iscategorical(v), s=char(v(i)); else, s=char(string(v(i))); end
end
