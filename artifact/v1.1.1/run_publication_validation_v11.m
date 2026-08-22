function [results,trials]=run_publication_validation_v11(varargin)
%RUN_PUBLICATION_VALIDATION_V11 Frozen publication-readiness validation.
%
% v1.1 does not retune v0.9. It expands to 60 cases, adds an independent
% sort-reduce sparse baseline, conversion-aware economics, and concrete
% MATLAB representation-byte measurements for sparse/dense/packed-BMD data.
p=inputParser;
addParameter(p,'Trials',5,@(x)isnumeric(x)&&isscalar(x)&&x>=3);
addParameter(p,'BuildTrials',3,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'PredictTrials',5,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'DenseMaxCoefficients',250000,@(x)isnumeric(x)&&isscalar(x)&&x>=1000);
addParameter(p,'SaveResults',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});
nTrials=round(p.Results.Trials); nBuild=round(p.Results.BuildTrials); nPred=round(p.Results.PredictTrials); denseMax=round(p.Results.DenseMaxCoefficients);
root=fileparts(mfilename('fullpath')); nCases=60;
rows=repmat(emptyRow(),nCases,1); raw=repmat(emptyTrialRow(),nCases*nTrials,1); rawPos=0;

fprintf('\nBMD-MATLAB v1.1 publication validation\n');
fprintf('=======================================\n');
fprintf('%d frozen cases; %d operation trials/case; v0.9 predictor unchanged.\n\n',nCases,nTrials);

for i=1:nCases
    c=build_publication_case_v11(i); row=emptyRow();
    row.case_name=c.name; row.family=c.family; row.application=c.application; row.source_basis=c.source_basis;
    row.p_terms=c.p_terms; row.q_terms=c.q_terms; row.pair_products=c.pair_products;
    row.p_degree=c.p_degree; row.q_degree=c.q_degree; row.result_degree_bound=c.result_degree_bound;
    row.sparse_operand_bytes=workspace_bytes_v11(c.p_sparse)+workspace_bytes_v11(c.q_sparse);
    fprintf('  %02d/%02d %-29s ... ',i,nCases,c.name);
    try
        % Sparse -> compact BMD conversion (outside operation timing).
        buildT=zeros(1,nBuild); baseMgr=[]; refs=[];
        for k=1:nBuild
            tic; [mm,rr]=build_bmd_pair_from_sparse_v10(c.p_sparse,c.q_sparse); buildT(k)=toc;
            if k==1, baseMgr=mm; refs=rr; end
        end
        row.bmd_build_median_s=median(buildT); bs=baseMgr.stats(); row.bmd_operand_union_nodes=bs.total_internal_nodes;
        row.bmd_packed_operand_bytes=workspace_bytes_v11(pack_bmd_roots_v11(baseMgr,refs));

        % Prediction is frozen before any multiply/warm-up.
        pt=zeros(1,nPred); pred=[];
        for k=1:nPred
            tic; pp=predict_multiply_route_v09(baseMgr,refs(1,:),refs(2,:),c.p_terms,c.q_terms); pt(k)=toc;
            if k==1, pred=pp; end
        end
        row.predictor_median_s=median(pt); row.predicted_new_nodes=pred.predicted_new_nodes;
        row.prediction_regime=pred.prediction_regime; row.prediction_confidence=pred.prediction_confidence;
        row.predicted_route=pred.recommended_route; row.pairs_per_predicted_new_node=pred.pairs_per_predicted_new_node;

        % Optional dense baseline. Eligibility is frozen by output coefficient count only.
        resultCoeffs=c.result_degree_bound+1; row.dense_eligible=double(resultCoeffs<=denseMax); dp=[]; dq=[];
        if row.dense_eligible
            db=zeros(1,nBuild);
            for k=1:nBuild
                tic; aa=sparse_terms_to_dense(c.p_sparse,denseMax); bb=sparse_terms_to_dense(c.q_sparse,denseMax); db(k)=toc;
                if k==1, dp=aa; dq=bb; end
            end
            row.dense_build_median_s=median(db); row.dense_operand_bytes=workspace_bytes_v11(dp)+workspace_bytes_v11(dq);
        end

        % Warm-up after prediction freeze only.
        [wm,wr]=baseMgr.compact(refs); wm.multiply(wr(1,:),wr(2,:));
        sparse_terms_multiply(c.p_sparse,c.q_sparse);
        sparse_terms_multiply_sortreduce_v11(c.p_sparse,c.q_sparse);
        if row.dense_eligible, z=conv(dp,dq); z(1)=z(1); end %#ok<NASGU>

        bt=zeros(1,nTrials); su=zeros(1,nTrials); ss=zeros(1,nTrials); dt=NaN(1,nTrials);
        firstMgr=[]; firstProd=[]; firstSU=[]; firstSS=[]; firstDense=[];
        for k=1:nTrials
            if row.dense_eligible, order=order4(k); else, order=order3(k); end
            bmdOut=[]; suOut=[]; ssOut=[]; denseOut=[]; mt=[];
            for op=order
                switch op
                    case 1, [bt(k),mt,bmdOut]=measureBmd(baseMgr,refs);
                    case 2, [su(k),suOut]=measureSparseUnique(c.p_sparse,c.q_sparse);
                    case 3, [ss(k),ssOut]=measureSparseSort(c.p_sparse,c.q_sparse);
                    case 4, [dt(k),denseOut]=measureDense(dp,dq);
                end
            end
            if k==1, firstMgr=mt; firstProd=bmdOut; firstSU=suOut; firstSS=ssOut; firstDense=denseOut; end
            rawPos=rawPos+1; rr=emptyTrialRow(); rr.case_name=c.name; rr.family=c.family; rr.trial=k; rr.order=orderText(order);
            rr.bmd_cold_s=bt(k); rr.sparse_unique_s=su(k); rr.sparse_sortreduce_s=ss(k); rr.dense_s=dt(k); raw(rawPos)=rr;
        end

        % Exact structural counts and prediction error.
        after=firstMgr.stats(firstProd); before=baseMgr.stats();
        row.actual_new_nodes=after.total_internal_nodes-before.total_internal_nodes; row.actual_result_nodes=after.reachable_nodes;
        row.new_node_abs_pct_error=abs(row.predicted_new_nodes-row.actual_new_nodes)/max(1,row.actual_new_nodes);
        row.bmd_median_s=median(bt); row.bmd_q25_s=pct(bt,.25); row.bmd_q75_s=pct(bt,.75);
        row.sparse_unique_median_s=median(su); row.sparse_sortreduce_median_s=median(ss);
        if row.sparse_unique_median_s<=row.sparse_sortreduce_median_s
            row.sparse_best_algorithm='UNIQUE_ACCUMARRAY'; sb=su; firstSparse=firstSU;
        else
            row.sparse_best_algorithm='SORT_REDUCE'; sb=ss; firstSparse=firstSS;
        end
        row.sparse_best_median_s=median(sb); row.sparse_best_q25_s=pct(sb,.25); row.sparse_best_q75_s=pct(sb,.75);
        row.sortreduce_speedup_vs_unique=row.sparse_unique_median_s/max(realmin,row.sparse_sortreduce_median_s);
        if row.dense_eligible, row.dense_median_s=median(dt); else, row.dense_median_s=NaN; end

        % BMD vs strongest explicit sparse baseline.
        row.sparse_best_over_bmd=row.sparse_best_median_s/max(realmin,row.bmd_median_s);
        if row.sparse_best_over_bmd>1, row.actual_bs_winner='BMD'; else, row.actual_bs_winner='SPARSE'; end
        robustLow=row.sparse_best_q25_s/max(realmin,row.bmd_q75_s); robustHigh=row.sparse_best_q75_s/max(realmin,row.bmd_q25_s);
        if robustLow>1, row.actual_bs_robust='BMD'; elseif robustHigh<1, row.actual_bs_robust='SPARSE'; else, row.actual_bs_robust='OVERLAP'; end
        if strcmp(row.predicted_route,'UNCERTAIN'), row.routed_correct=NaN; else, row.routed_correct=double(strcmp(row.predicted_route,row.actual_bs_winner)); end
        row.false_bmd=double(strcmp(row.predicted_route,'BMD') && strcmp(row.actual_bs_winner,'SPARSE'));
        row.false_sparse=double(strcmp(row.predicted_route,'SPARSE') && strcmp(row.actual_bs_winner,'BMD'));

        % Concrete representation bytes. Packed BMD excludes live map/cache overhead by definition.
        row.bmd_packed_result_bytes=workspace_bytes_v11(pack_bmd_roots_v11(firstMgr,firstProd));
        row.sparse_result_bytes=workspace_bytes_v11(firstSparse);
        row.sparse_result_terms=numel(firstSparse.exponents);
        row.sparse_result_over_bmd_packed_bytes=row.sparse_result_bytes/max(1,row.bmd_packed_result_bytes);
        if row.dense_eligible
            if isempty(firstDense), firstDense=conv(dp,dq); end
            row.dense_result_bytes=workspace_bytes_v11(firstDense);
        end

        % End-to-end economics from existing sparse term lists.
        row.bmd_single_shot_from_sparse_s=row.bmd_build_median_s+row.predictor_median_s+row.bmd_median_s;
        row.bmd_single_shot_speedup_vs_sparse_best=row.sparse_best_median_s/max(realmin,row.bmd_single_shot_from_sparse_s);
        if strcmp(row.predicted_route,'BMD'), routeOp=row.bmd_median_s; elseif strcmp(row.predicted_route,'SPARSE'), routeOp=row.sparse_best_median_s; else, routeOp=NaN; end
        if isfinite(routeOp)
            row.bmd_resident_router_total_s=row.predictor_median_s+routeOp;
            row.router_from_sparse_total_s=row.bmd_build_median_s+row.predictor_median_s+routeOp;
            row.router_from_sparse_speedup_vs_sparse_best=row.sparse_best_median_s/max(realmin,row.router_from_sparse_total_s);
            row.routed_operation_regret=routeOp/max(realmin,min(row.bmd_median_s,row.sparse_best_median_s));
        end
        denom=row.sparse_best_median_s-row.bmd_median_s;
        if denom>0, row.break_even_reuses_route_once=ceil((row.bmd_build_median_s+row.predictor_median_s)/denom); else, row.break_even_reuses_route_once=Inf; end
        if row.dense_eligible, row.dense_single_shot_from_sparse_s=row.dense_build_median_s+row.dense_median_s; end
        [row.operation_oracle_winner,row.operation_oracle_s]=winner3(row.bmd_median_s,row.sparse_best_median_s,row.dense_median_s,row.dense_eligible);
        [row.single_shot_oracle_winner,row.single_shot_oracle_s]=winner3(row.bmd_single_shot_from_sparse_s,row.sparse_best_median_s,row.dense_single_shot_from_sparse_s,row.dense_eligible);

        row.numeric_check=checkNumeric(firstMgr,firstProd,firstSU,firstSS,dp,dq,row.dense_eligible);
        if ~row.numeric_check, error('BMD:V11Numeric','BMD/sparse/dense validation mismatch.'); end
        row.status='OK';
    catch ME
        row.status=['ERROR:' ME.identifier];
    end
    rows(i)=row;
    fprintf('%s',row.status);
    if strcmp(row.status,'OK')
        fprintf(' route=%s actual=%s BMD=%.3fms sparse*=%.3fms',row.predicted_route,row.actual_bs_winner,1000*row.bmd_median_s,1000*row.sparse_best_median_s);
        if row.dense_eligible, fprintf(' dense=%.3fms',1000*row.dense_median_s); end
    end
    fprintf('\n');
end
results=struct2table(rows); trials=struct2table(raw(1:rawPos));
if p.Results.SaveResults
    outDir=fullfile(root,'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
    writetable(results,fullfile(outDir,'publication_validation_results_v11.csv'));
    writetable(trials,fullfile(outDir,'publication_validation_trials_v11.csv'));
end
end

function [t,m,out]=measureBmd(base,refs), [m,r]=base.compact(refs); tic; out=m.multiply(r(1,:),r(2,:)); t=toc; end
function [t,out]=measureSparseUnique(a,b), tic; out=sparse_terms_multiply(a,b); t=toc; end
function [t,out]=measureSparseSort(a,b), tic; out=sparse_terms_multiply_sortreduce_v11(a,b); t=toc; end
function [t,out]=measureDense(a,b), tic; out=conv(a,b); t=toc; end
function ord=order3(k), P={[1 2 3],[2 3 1],[3 1 2],[1 3 2],[2 1 3]}; ord=P{1+mod(k-1,numel(P))}; end
function ord=order4(k), P={[1 2 3 4],[2 3 4 1],[3 4 1 2],[4 1 2 3],[1 3 2 4]}; ord=P{1+mod(k-1,numel(P))}; end
function s=orderText(ord), names={'BMD','SPARSE_UNIQUE','SPARSE_SORT','DENSE'}; s=strjoin(names(ord),'_'); end
function tf=checkNumeric(m,p,su,ss,dp,dq,denseEligible)
tf=sparse_terms_same(su,ss);
for x=[0.5 0.999], a=m.evaluate(p,x); b=sparse_terms_evaluate(su,x); tf=tf && near(a,b); end
if denseEligible
    z=conv(dp,dq); sd=sparse_terms_to_dense(su,1e7); z=trimDense(z); sd=trimDense(sd);
    tf=tf && numel(z)==numel(sd) && all(abs(z-sd)<=1e-9*max(1,max(abs(sd))));
end
end
function p=trimDense(p), idx=find(p~=0,1,'first'); if isempty(idx), p=0; else, p=p(idx:end); end, end
function tf=near(a,b), scale=max([1 abs(a) abs(b)]); tf=abs(a-b)<=1e-8*scale; end
function q=pct(x,p), x=sort(x(:)); pos=1+(numel(x)-1)*p; lo=floor(pos); hi=ceil(pos); if lo==hi, q=x(lo); else, q=x(lo)+(pos-lo)*(x(hi)-x(lo)); end, end
function [name,t]=winner3(b,s,d,denseEligible), vals=[b s]; names={'BMD','SPARSE'}; if denseEligible && isfinite(d), vals(end+1)=d; names{end+1}='DENSE'; end; [t,idx]=min(vals); name=names{idx}; end
function r=emptyRow()
r=struct('case_name','','family','','application','','source_basis','','status','', ...
'p_terms',NaN,'q_terms',NaN,'pair_products',NaN,'p_degree',NaN,'q_degree',NaN,'result_degree_bound',NaN, ...
'sparse_operand_bytes',NaN,'bmd_build_median_s',NaN,'bmd_operand_union_nodes',NaN,'bmd_packed_operand_bytes',NaN, ...
'predictor_median_s',NaN,'predicted_new_nodes',NaN,'prediction_regime','','prediction_confidence','','predicted_route','','pairs_per_predicted_new_node',NaN, ...
'actual_new_nodes',NaN,'actual_result_nodes',NaN,'new_node_abs_pct_error',NaN, ...
'bmd_median_s',NaN,'bmd_q25_s',NaN,'bmd_q75_s',NaN, ...
'sparse_unique_median_s',NaN,'sparse_sortreduce_median_s',NaN,'sparse_best_algorithm','','sparse_best_median_s',NaN,'sparse_best_q25_s',NaN,'sparse_best_q75_s',NaN,'sortreduce_speedup_vs_unique',NaN, ...
'dense_eligible',0,'dense_build_median_s',NaN,'dense_median_s',NaN,'dense_operand_bytes',NaN, ...
'sparse_best_over_bmd',NaN,'actual_bs_winner','','actual_bs_robust','','routed_correct',NaN,'false_bmd',0,'false_sparse',0, ...
'bmd_packed_result_bytes',NaN,'sparse_result_bytes',NaN,'dense_result_bytes',NaN,'sparse_result_terms',NaN,'sparse_result_over_bmd_packed_bytes',NaN, ...
'bmd_single_shot_from_sparse_s',NaN,'bmd_single_shot_speedup_vs_sparse_best',NaN,'bmd_resident_router_total_s',NaN,'router_from_sparse_total_s',NaN,'router_from_sparse_speedup_vs_sparse_best',NaN,'routed_operation_regret',NaN, ...
'break_even_reuses_route_once',NaN,'dense_single_shot_from_sparse_s',NaN,'operation_oracle_winner','','operation_oracle_s',NaN,'single_shot_oracle_winner','','single_shot_oracle_s',NaN,'numeric_check',0);
end
function r=emptyTrialRow()
r=struct('case_name','','family','','trial',NaN,'order','','bmd_cold_s',NaN,'sparse_unique_s',NaN,'sparse_sortreduce_s',NaN,'dense_s',NaN);
end
