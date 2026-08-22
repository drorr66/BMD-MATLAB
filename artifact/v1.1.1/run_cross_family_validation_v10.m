function [results,trials]=run_cross_family_validation_v10(varargin)
%RUN_CROSS_FAMILY_VALIDATION_V10 Cross-family + conversion-aware validation.
%
% The v0.9 predictor is frozen.  Prediction is made before any multiplication
% for each case.  v1.0 adds sparse->BMD conversion cost, dense coefficient
% vectors when degree is bounded, and break-even reuse accounting.
p=inputParser;
addParameter(p,'Trials',5,@(x)isnumeric(x)&&isscalar(x)&&x>=5);
addParameter(p,'BuildTrials',3,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'DenseMaxCoefficients',250000,@(x)isnumeric(x)&&isscalar(x)&&x>=1000);
addParameter(p,'SaveResults',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});
nTrials=round(p.Results.Trials); nBuild=round(p.Results.BuildTrials); denseMax=round(p.Results.DenseMaxCoefficients);
root=fileparts(mfilename('fullpath'));
nCases=15; rows=repmat(emptyRow(),nCases,1); raw=repmat(emptyTrialRow(),nCases*nTrials,1); rawPos=0;

fprintf('\nBMD-MATLAB v1.0 cross-family validation\n');
fprintf('=========================================\n');
fprintf('%d frozen cases; %d operation trials/case; %d conversion trials/case.\n\n',nCases,nTrials,nBuild);

for i=1:nCases
    c=build_validation_case_v10(i); row=emptyRow();
    row.case_name=c.name; row.family=c.family; row.application=c.application;
    row.p_terms=c.p_terms; row.q_terms=c.q_terms; row.pair_products=c.pair_products;
    row.p_degree=c.p_degree; row.q_degree=c.q_degree; row.result_degree_bound=c.result_degree_bound;
    fprintf('  %-27s ... ',c.name);
    try
        % Conversion is measured independently and occurs before prediction.
        buildT=zeros(1,nBuild); baseMgr=[]; refs=[];
        for k=1:nBuild
            tic; [mm,rr]=build_bmd_pair_from_sparse_v10(c.p_sparse,c.q_sparse); buildT(k)=toc;
            if k==1, baseMgr=mm; refs=rr; end
        end
        row.bmd_build_median_s=median(buildT);
        bs=baseMgr.stats(); row.bmd_operand_union_nodes=bs.total_internal_nodes;

        % Freeze v0.9 route before any multiplication or warm-up.
        pt=zeros(1,5); pred=[];
        for k=1:5
            tic; pp=predict_multiply_route_v09(baseMgr,refs(1,:),refs(2,:),c.p_terms,c.q_terms); pt(k)=toc;
            if k==1, pred=pp; end
        end
        row.predictor_median_s=median(pt); row.predicted_new_nodes=pred.predicted_new_nodes;
        row.prediction_regime=pred.prediction_regime; row.prediction_confidence=pred.prediction_confidence;
        row.predicted_route=pred.recommended_route; row.pairs_per_predicted_new_node=pred.pairs_per_predicted_new_node;

        % Dense baseline is optional and bounded by result coefficient count.
        resultCoeffs=c.result_degree_bound+1;
        row.dense_eligible=double(resultCoeffs<=denseMax);
        dp=[]; dq=[];
        if row.dense_eligible
            db=zeros(1,nBuild);
            for k=1:nBuild
                tic; aa=sparse_terms_to_dense(c.p_sparse,denseMax); bb=sparse_terms_to_dense(c.q_sparse,denseMax); db(k)=toc;
                if k==1, dp=aa; dq=bb; end
            end
            row.dense_build_median_s=median(db);
        end

        % Warm-up only after prediction is frozen.
        [wm,wr]=baseMgr.compact(refs); wm.multiply(wr(1,:),wr(2,:));
        sparse_terms_multiply(c.p_sparse,c.q_sparse);
        if row.dense_eligible, tmp=conv(dp,dq); tmp(1)=tmp(1); end %#ok<NASGU>

        bt=zeros(1,nTrials); st=zeros(1,nTrials); dt=NaN(1,nTrials);
        firstMgr=[]; firstProd=[]; firstSparse=[];
        for k=1:nTrials
            ord=1+mod(k-1,3);
            if ~row.dense_eligible, ord=1+mod(k-1,2); end
            if row.dense_eligible
                if ord==1
                    [bt(k),mt,prod]=measureBmd(baseMgr,refs); [st(k),spr]=measureSparse(c.p_sparse,c.q_sparse); dt(k)=measureDense(dp,dq);
                    order='BMD_SPARSE_DENSE';
                elseif ord==2
                    [st(k),spr]=measureSparse(c.p_sparse,c.q_sparse); dt(k)=measureDense(dp,dq); [bt(k),mt,prod]=measureBmd(baseMgr,refs);
                    order='SPARSE_DENSE_BMD';
                else
                    dt(k)=measureDense(dp,dq); [bt(k),mt,prod]=measureBmd(baseMgr,refs); [st(k),spr]=measureSparse(c.p_sparse,c.q_sparse);
                    order='DENSE_BMD_SPARSE';
                end
            else
                if ord==1
                    [bt(k),mt,prod]=measureBmd(baseMgr,refs); [st(k),spr]=measureSparse(c.p_sparse,c.q_sparse); order='BMD_SPARSE';
                else
                    [st(k),spr]=measureSparse(c.p_sparse,c.q_sparse); [bt(k),mt,prod]=measureBmd(baseMgr,refs); order='SPARSE_BMD';
                end
            end
            if k==1, firstMgr=mt; firstProd=prod; firstSparse=spr; end
            rawPos=rawPos+1; rr=emptyTrialRow(); rr.case_name=c.name; rr.family=c.family; rr.trial=k; rr.order=order;
            rr.bmd_cold_s=bt(k); rr.sparse_s=st(k); rr.dense_s=dt(k); rr.sparse_over_bmd=st(k)/max(realmin,bt(k)); raw(rawPos)=rr;
        end

        after=firstMgr.stats(firstProd); before=baseMgr.stats();
        row.actual_new_nodes=after.total_internal_nodes-before.total_internal_nodes; row.actual_result_nodes=after.reachable_nodes;
        row.new_node_abs_pct_error=abs(row.predicted_new_nodes-row.actual_new_nodes)/max(1,row.actual_new_nodes);
        row.bmd_median_s=median(bt); row.bmd_q25_s=percentileLinear(bt,.25); row.bmd_q75_s=percentileLinear(bt,.75);
        row.sparse_median_s=median(st); row.sparse_q25_s=percentileLinear(st,.25); row.sparse_q75_s=percentileLinear(st,.75);
        if row.dense_eligible, row.dense_median_s=median(dt); else, row.dense_median_s=NaN; end
        row.sparse_over_bmd=row.sparse_median_s/max(realmin,row.bmd_median_s);
        if row.sparse_over_bmd>1, row.actual_bs_winner='BMD'; else, row.actual_bs_winner='SPARSE'; end
        robustLow=row.sparse_q25_s/max(realmin,row.bmd_q75_s); robustHigh=row.sparse_q75_s/max(realmin,row.bmd_q25_s);
        if robustLow>1, row.actual_bs_robust='BMD'; elseif robustHigh<1, row.actual_bs_robust='SPARSE'; else, row.actual_bs_robust='OVERLAP'; end
        if strcmp(row.predicted_route,'UNCERTAIN'), row.routed_correct=NaN; else, row.routed_correct=double(strcmp(row.predicted_route,row.actual_bs_winner)); end
        row.false_bmd=double(strcmp(row.predicted_route,'BMD') && strcmp(row.actual_bs_winner,'SPARSE'));
        row.false_sparse=double(strcmp(row.predicted_route,'SPARSE') && strcmp(row.actual_bs_winner,'BMD'));

        % End-to-end economics from an existing sparse term list.
        row.bmd_single_shot_from_sparse_s=row.bmd_build_median_s+row.predictor_median_s+row.bmd_median_s;
        row.bmd_single_shot_speedup_vs_sparse=row.sparse_median_s/max(realmin,row.bmd_single_shot_from_sparse_s);
        if strcmp(row.predicted_route,'BMD')
            routeOp=row.bmd_median_s;
        elseif strcmp(row.predicted_route,'SPARSE')
            routeOp=row.sparse_median_s;
        else
            routeOp=NaN;
        end
        if isfinite(routeOp)
            row.bmd_resident_router_total_s=row.predictor_median_s+routeOp;
            row.current_router_from_sparse_s=row.bmd_build_median_s+row.predictor_median_s+routeOp;
            row.current_router_from_sparse_speedup_vs_sparse=row.sparse_median_s/max(realmin,row.current_router_from_sparse_s);
            row.routed_operation_regret=routeOp/max(realmin,min(row.bmd_median_s,row.sparse_median_s));
        end
        denom=row.sparse_median_s-row.bmd_median_s;
        if denom>0
            row.break_even_reuses_route_once=ceil((row.bmd_build_median_s+row.predictor_median_s)/denom);
        else
            row.break_even_reuses_route_once=Inf;
        end
        if row.dense_eligible
            row.dense_single_shot_from_sparse_s=row.dense_build_median_s+row.dense_median_s;
        end
        [row.operation_oracle_winner,row.operation_oracle_s]=winner3(row.bmd_median_s,row.sparse_median_s,row.dense_median_s,row.dense_eligible);
        [row.single_shot_oracle_winner,row.single_shot_oracle_s]=winner3(row.bmd_single_shot_from_sparse_s,row.sparse_median_s,row.dense_single_shot_from_sparse_s,row.dense_eligible);

        row.numeric_check=checkNumeric(firstMgr,firstProd,firstSparse,dp,dq,row.dense_eligible);
        if ~row.numeric_check, error('BMD:V10Numeric','BMD/sparse/dense validation mismatch.'); end
        row.status='OK';
    catch ME
        row.status=['ERROR:' ME.identifier];
    end
    rows(i)=row;
    fprintf('%s',row.status);
    if strcmp(row.status,'OK')
        fprintf(' route=%s actual=%s build=%.2fms op=%.2f/%.2fms',row.predicted_route,row.actual_bs_winner,1000*row.bmd_build_median_s,1000*row.bmd_median_s,1000*row.sparse_median_s);
        if isfinite(row.break_even_reuses_route_once), fprintf(' BE=%g',row.break_even_reuses_route_once); end
    end
    fprintf('\n');
end
results=struct2table(rows); trials=struct2table(raw(1:rawPos));
if p.Results.SaveResults
    outDir=fullfile(root,'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
    writetable(results,fullfile(outDir,'validation_results_v10.csv'));
    writetable(trials,fullfile(outDir,'validation_trials_v10.csv'));
end
end

function [t,m,out]=measureBmd(base,refs)
[m,r]=base.compact(refs); tic; out=m.multiply(r(1,:),r(2,:)); t=toc;
end
function [t,out]=measureSparse(a,b)
tic; out=sparse_terms_multiply(a,b); t=toc;
end
function t=measureDense(a,b)
tic; z=conv(a,b); t=toc; if isempty(z), error('BMD:V10Dense','Unexpected empty dense product.'); end
end
function tf=checkNumeric(m,p,s,dp,dq,denseEligible)
tf=true;
for x=[0.5 0.999]
    a=m.evaluate(p,x); b=sparse_terms_evaluate(s,x); tf=tf && near(a,b);
end
if denseEligible
    z=conv(dp,dq); ss=sparse_terms_to_dense(s,1e7);
    z=trimDense(z); ss=trimDense(ss); tf=tf && numel(z)==numel(ss) && all(abs(z-ss)<=1e-9*max(1,max(abs(ss))));
end
end
function p=trimDense(p)
idx=find(p~=0,1,'first'); if isempty(idx), p=0; else, p=p(idx:end); end
end
function tf=near(a,b)
scale=max([1 abs(a) abs(b)]); tf=abs(a-b)<=1e-8*scale;
end
function q=percentileLinear(x,p)
x=sort(x(:)); pos=1+(numel(x)-1)*p; lo=floor(pos); hi=ceil(pos); if lo==hi, q=x(lo); else, q=x(lo)+(pos-lo)*(x(hi)-x(lo)); end
end
function [name,t]=winner3(b,s,d,denseEligible)
vals=[b s]; names={'BMD','SPARSE'}; if denseEligible && isfinite(d), vals(end+1)=d; names{end+1}='DENSE'; end
[t,idx]=min(vals); name=names{idx};
end
function r=emptyRow()
r=struct('case_name','','family','','application','','status','', ...
'p_terms',NaN,'q_terms',NaN,'pair_products',NaN,'p_degree',NaN,'q_degree',NaN,'result_degree_bound',NaN, ...
'bmd_build_median_s',NaN,'bmd_operand_union_nodes',NaN,'predictor_median_s',NaN,'predicted_new_nodes',NaN, ...
'prediction_regime','','prediction_confidence','','predicted_route','','pairs_per_predicted_new_node',NaN, ...
'actual_new_nodes',NaN,'actual_result_nodes',NaN,'new_node_abs_pct_error',NaN, ...
'bmd_median_s',NaN,'bmd_q25_s',NaN,'bmd_q75_s',NaN,'sparse_median_s',NaN,'sparse_q25_s',NaN,'sparse_q75_s',NaN, ...
'dense_eligible',0,'dense_build_median_s',NaN,'dense_median_s',NaN,'sparse_over_bmd',NaN,'actual_bs_winner','','actual_bs_robust','', ...
'routed_correct',NaN,'false_bmd',0,'false_sparse',0,'bmd_single_shot_from_sparse_s',NaN,'bmd_single_shot_speedup_vs_sparse',NaN, ...
'bmd_resident_router_total_s',NaN,'current_router_from_sparse_s',NaN,'current_router_from_sparse_speedup_vs_sparse',NaN,'routed_operation_regret',NaN, ...
'break_even_reuses_route_once',NaN,'dense_single_shot_from_sparse_s',NaN,'operation_oracle_winner','','operation_oracle_s',NaN, ...
'single_shot_oracle_winner','','single_shot_oracle_s',NaN,'numeric_check',0);
end
function r=emptyTrialRow()
r=struct('case_name','','family','','trial',NaN,'order','','bmd_cold_s',NaN,'sparse_s',NaN,'dense_s',NaN,'sparse_over_bmd',NaN);
end
