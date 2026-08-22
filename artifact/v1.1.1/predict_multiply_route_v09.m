function out = predict_multiply_route_v09(mgr,pRef,qRef,pTerms,qTerms)
%PREDICT_MULTIPLY_ROUTE_V09 Conservative pre-multiply BMD/sparse predictor.
%
% IMPORTANT: this function never calls BMDManager.multiply. It uses only
% existing DAG topology plus caller-supplied sparse term counts.
%
% v0.9 is a research predictor, not a universal routing theorem. It may
% return UNCERTAIN. The bit-cube new-node model was trained offline on 800
% deterministic exact-reference cases and evaluated on 200 held-out cases.

f=bmd_predictor_features_v09(mgr,pRef,qRef,pTerms,qTerms);
out=f;
out.predictor_version='v0.9-static-20260820';
out.predicted_new_nodes=NaN;
out.prediction_regime='';
out.prediction_confidence='LOW';

if f.ordered_disjoint_bands
    % Ordered disjoint variable bands cannot create a same-level carry. In
    % the measured v0.5/v0.7/v0.8 families, the generic multiply creates one
    % new layer per node in the earlier (lower-level) band.
    if maxReachableLevel(mgr,pRef) < minReachableLevel(mgr,qRef)
        out.predicted_new_nodes=max(1,f.p_nodes);
    else
        out.predicted_new_nodes=max(1,f.q_nodes);
    end
    out.prediction_regime='ORDERED_DISJOINT_BANDS';
    out.prediction_confidence='HIGH';
elseif f.both_bitcube_chain
    % Ridge model in log1p space trained only on exact structural counts, not
    % MATLAB timings. Features are cheap topology summaries.
    beta=[ ...
        1.8716285233860908, ...
       -0.49888659949909009, -0.49550136457640398, ...
        0.025865401111516365, 0.30484129754149636, ...
        0.28101502349890534, -0.3836297310014834, ...
        0.046664149634361947, -0.089931009118816949, ...
        0.95963681138947232, -0.058268032408828399, ...
        0.050629935773811183, 0.1887658610473012];
    x=[f.p_nodes,f.q_nodes,f.common_levels,f.carry_runsum,f.carry_run2, ...
       f.carry_max,f.carry_massrun,f.adjacent_mass,f.invdist, ...
       f.suffix_product_sum,f.suffix_carry_sum,f.level_span];
    z=beta(1)+sum(beta(2:end).*log1p(double(x)));
    out.predicted_new_nodes=max(1,expm1(z));
    out.prediction_regime='BITCUBE_RIDGE';
    out.prediction_confidence='HIGH';
elseif f.union_nodes > 64 && f.common_levels > 0
    % Conservative general-overlap estimate. The coefficient 0.82 is the
    % central structural ratio observed in the v0.6 family. We intentionally
    % keep confidence LOW and never route such a case to BMD automatically.
    out.predicted_new_nodes=max([1,0.82*f.suffix_product_sum,f.union_nodes/4]);
    out.prediction_regime='GENERAL_HEAVY_OVERLAP';
    out.prediction_confidence='LOW';
else
    out.predicted_new_nodes=max([1,f.suffix_product_sum,f.union_nodes/2]);
    out.prediction_regime='GENERAL_UNCALIBRATED';
    out.prediction_confidence='LOW';
end

out.pairs_per_predicted_new_node=f.pair_products/max(1,out.predicted_new_nodes);

% Conservative routing policy distilled from v0.5-v0.8 MATLAB measurements.
% The gray zone is intentional. These constants are implementation/session
% calibration guides, NOT universal mathematical thresholds.
if strcmp(out.prediction_confidence,'LOW')
    if out.predicted_new_nodes >= 64 || out.pairs_per_predicted_new_node <= 1800
        out.recommended_route='SPARSE';
        out.route_reason='low-confidence overlap with high predicted DAG work';
    else
        out.recommended_route='UNCERTAIN';
        out.route_reason='structure outside calibrated high-confidence families';
    end
elseif out.predicted_new_nodes >= 64 || out.pairs_per_predicted_new_node <= 1800
    out.recommended_route='SPARSE';
    out.route_reason='predicted BMD closure cost too high';
elseif f.pair_products >= 40000 && out.pairs_per_predicted_new_node >= 3000
    out.recommended_route='BMD';
    out.route_reason='high explicit work per predicted BMD node';
elseif f.pair_products < 30000
    out.recommended_route='SPARSE';
    out.route_reason='explicit sparse workload too small to amortize BMD overhead';
else
    out.recommended_route='UNCERTAIN';
    out.route_reason='near calibrated crossover / insufficient margin';
end
end

function lev=maxReachableLevel(mgr,ref)
ids=reachable(mgr,ref);
if isempty(ids), lev=-Inf; else, lev=max(double(mgr.levels(ids))); end
end
function lev=minReachableLevel(mgr,ref)
ids=reachable(mgr,ref);
if isempty(ids), lev=Inf; else, lev=min(double(mgr.levels(ids))); end
end
function ids=reachable(mgr,ref)
if ref(1)==0 || ref(2)==1, ids=zeros(1,0); return; end
seen=false(1,numel(mgr.levels)); stack=double(ref(2)); ids=zeros(1,0);
while ~isempty(stack)
    n=stack(end); stack(end)=[];
    if n==1 || seen(n), continue; end
    seen(n)=true; ids(end+1)=n; %#ok<AGROW>
    lo=double(mgr.lowNode(n)); hi=double(mgr.highNode(n));
    if lo~=1, stack(end+1)=lo; end %#ok<AGROW>
    if hi~=1, stack(end+1)=hi; end %#ok<AGROW>
end
end
