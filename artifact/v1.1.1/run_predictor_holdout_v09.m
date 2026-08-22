function [results,trials]=run_predictor_holdout_v09(varargin)
%RUN_PREDICTOR_HOLDOUT_V09 Fresh timing holdout for the v0.9 static router.
%
% Prediction is computed and frozen BEFORE any multiplication/warm-up for
% each case. The subsequent BMD/sparse timings exist only to score it.
p=inputParser;
addParameter(p,'Trials',7,@(x)isnumeric(x)&&isscalar(x)&&x>=5);
addParameter(p,'SaveResults',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:}); nTrials=round(p.Results.Trials);
root=fileparts(mfilename('fullpath'));
cases=holdoutCases();
rows=repmat(emptyRow(),numel(cases),1);
raw=repmat(emptyTrialRow(),numel(cases)*nTrials,1); rawPos=0;
bank=make_sharing_template_bank_v06(64,256,128);

fprintf('\nBMD-MATLAB v0.9 predictor holdout\n');
fprintf('=================================\n');
fprintf('%d fresh cases, %d timing trials/case. Prediction is frozen before warm-up.\n\n',numel(cases),nTrials);

for i=1:numel(cases)
    c=cases(i); row=emptyRow(); row.case_name=c.name; row.family=c.family;
    fprintf('  %-25s ... ',c.name);
    try
        [baseMgr,refs,sp,sq,pTerms,qTerms]=buildCase(c,bank);
        % Predictor first. Repeat only the predictor to get a stable overhead
        % number; no multiplication occurs in this block.
        pt=zeros(1,5); pred=[];
        for k=1:5
            tic; pp=predict_multiply_route_v09(baseMgr,refs(1,:),refs(2,:),pTerms,qTerms); pt(k)=toc;
            if k==1, pred=pp; end
        end
        row.predictor_median_s=median(pt);
        row.predicted_new_nodes=pred.predicted_new_nodes;
        row.prediction_regime=pred.prediction_regime;
        row.prediction_confidence=pred.prediction_confidence;
        row.predicted_route=pred.recommended_route;
        row.pair_products=pred.pair_products;
        row.pairs_per_predicted_new_node=pred.pairs_per_predicted_new_node;

        % Warm-up only AFTER the route has been frozen.
        [wm,wr]=baseMgr.compact(refs); wm.multiply(wr(1,:),wr(2,:));
        sparse_terms_multiply(sp,sq);

        bt=zeros(1,nTrials); st=zeros(1,nTrials); firstMgr=[]; firstProd=[]; firstSparse=[];
        for k=1:nTrials
            if mod(k,2)==1
                [bt(k),mt,prod]=measureBmd(baseMgr,refs); [st(k),spr]=measureSparse(sp,sq);
            else
                [st(k),spr]=measureSparse(sp,sq); [bt(k),mt,prod]=measureBmd(baseMgr,refs);
            end
            if k==1, firstMgr=mt; firstProd=prod; firstSparse=spr; end
            rawPos=rawPos+1; rr=emptyTrialRow(); rr.case_name=c.name; rr.trial=k;
            if mod(k,2)==1, rr.order='BMD_FIRST'; else, rr.order='SPARSE_FIRST'; end
            rr.bmd_cold_s=bt(k); rr.sparse_s=st(k); rr.sparse_over_bmd=st(k)/max(realmin,bt(k)); raw(rawPos)=rr;
        end

        after=firstMgr.stats(firstProd); before=baseMgr.stats();
        row.actual_new_nodes=after.total_internal_nodes-before.total_internal_nodes;
        row.actual_result_nodes=after.reachable_nodes;
        row.new_node_abs_pct_error=abs(row.predicted_new_nodes-row.actual_new_nodes)/max(1,row.actual_new_nodes);
        row.bmd_median_s=median(bt); row.bmd_q25_s=percentileLinear(bt,.25); row.bmd_q75_s=percentileLinear(bt,.75);
        row.sparse_median_s=median(st); row.sparse_q25_s=percentileLinear(st,.25); row.sparse_q75_s=percentileLinear(st,.75);
        row.actual_ratio_median=row.sparse_median_s/max(realmin,row.bmd_median_s);
        if row.actual_ratio_median>1, row.actual_median_winner='BMD'; else, row.actual_median_winner='SPARSE'; end
        robustLow=row.sparse_q25_s/max(realmin,row.bmd_q75_s); robustHigh=row.sparse_q75_s/max(realmin,row.bmd_q25_s);
        if robustLow>1, row.actual_robust_winner='BMD'; elseif robustHigh<1, row.actual_robust_winner='SPARSE'; else, row.actual_robust_winner='OVERLAP'; end
        if strcmp(row.predicted_route,'UNCERTAIN'), row.routed_correct=NaN; else, row.routed_correct=double(strcmp(row.predicted_route,row.actual_median_winner)); end
        row.false_bmd=double(strcmp(row.predicted_route,'BMD') && strcmp(row.actual_median_winner,'SPARSE'));
        row.false_sparse=double(strcmp(row.predicted_route,'SPARSE') && strcmp(row.actual_median_winner,'BMD'));
        row.numeric_check=checkNumeric(firstMgr,firstProd,firstSparse);
        if ~row.numeric_check, error('BMD:V09Numeric','BMD/sparse holdout evaluation mismatch.'); end
        row.status='OK';
    catch ME
        row.status=['ERROR:' ME.identifier];
    end
    rows(i)=row;
    fprintf('%s',row.status);
    if strcmp(row.status,'OK')
        fprintf(' pred=%s new=%.1f/%d ratio=%.3fx actual=%s',row.predicted_route,row.predicted_new_nodes,row.actual_new_nodes,row.actual_ratio_median,row.actual_median_winner);
    end
    fprintf('\n');
end
results=struct2table(rows); trials=struct2table(raw(1:rawPos));
if p.Results.SaveResults
    outDir=fullfile(root,'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
    writetable(results,fullfile(outDir,'predictor_holdout_results_v09.csv'));
    writetable(trials,fullfile(outDir,'predictor_holdout_trials_v09.csv'));
end
end

function C=holdoutCases()
C=[ ...
 struct('name','bitcube_disjoint_adj','family','BITCUBE','param',{[8 10 11 12 13 14 15 16]}), ...
 struct('name','bitcube_single_gap','family','BITCUBE','param',{[7 10 11 12 13 14 15 16]}), ...
 struct('name','bitcube_carry_chain','family','BITCUBE','param',{[5 8 10 12 14 16 18 20]}), ...
 struct('name','bitcube_two_overlap','family','BITCUBE','param',{[6 7 9 11 13 15 17 19]}), ...
 struct('name','bitcube_two_separate','family','BITCUBE','param',{[1 4 9 11 13 15 17 19]}), ...
 struct('name','grid_blocks_144','family','GRID','param',144), ...
 struct('name','grid_blocks_224','family','GRID','param',224), ...
 struct('name','grid_blocks_320','family','GRID','param',320), ...
 struct('name','sharing_templates_3','family','SHARING','param',3), ...
 struct('name','sharing_templates_6','family','SHARING','param',6), ...
 struct('name','sharing_templates_12','family','SHARING','param',12)];
end

function [m,refs,sp,sq,np,nq]=buildCase(c,bank)
switch c.family
    case 'BITCUBE'
        pe=build_bitcube_support_v08(0:7); qe=build_bitcube_support_v08(c.param);
        m=BMDManager(); p=m.indicatorExponents(pe); q=m.indicatorExponents(qe); [m,refs]=m.compact([p;q]);
        sp=sparse_terms(pe,ones(1,numel(pe))); sq=sparse_terms(qe,ones(1,numel(qe))); np=numel(pe); nq=numel(qe);
    case 'GRID'
        blocks=c.param; m=BMDManager(); p=m.geometricSum(255); q=m.geometricSumShifted(blocks-1,11); [m,refs]=m.compact([p;q]);
        sp=build_sparse_family('geometric_sum',255,0); sq=build_sparse_shifted_geometric(blocks,10); np=256; nq=blocks;
    case 'SHARING'
        T=c.param; [exps,~]=build_sharing_case_v06(bank,T,64,9);
        m=BMDManager(); p=m.indicatorExponents(exps); q=m.geometricSum(31); [m,refs]=m.compact([p;q]);
        sp=sparse_terms(exps,ones(1,numel(exps))); sq=build_sparse_family('geometric_sum',31,0); np=numel(exps); nq=32;
    otherwise
        error('BMD:V09Family','Unknown holdout family.');
end
end

function [t,m,out]=measureBmd(base,refs)
[m,r]=base.compact(refs); tic; out=m.multiply(r(1,:),r(2,:)); t=toc;
end
function [t,out]=measureSparse(a,b)
tic; out=sparse_terms_multiply(a,b); t=toc;
end
function tf=checkNumeric(m,p,s)
tf=true; for x=[1.0 0.999999], tf=tf && near(m.evaluate(p,x),sparse_terms_evaluate(s,x)); end
end
function tf=near(a,b)
scale=max([1 abs(a) abs(b)]); tf=abs(a-b)<=1e-8*scale;
end
function q=percentileLinear(x,p)
x=sort(x(:)); pos=1+(numel(x)-1)*p; lo=floor(pos); hi=ceil(pos); if lo==hi, q=x(lo); else, q=x(lo)+(pos-lo)*(x(hi)-x(lo)); end
end
function r=emptyRow()
r=struct('case_name','','family','','status','','predictor_median_s',NaN,'predicted_new_nodes',NaN, ...
    'actual_new_nodes',NaN,'actual_result_nodes',NaN,'new_node_abs_pct_error',NaN,'pair_products',NaN, ...
    'pairs_per_predicted_new_node',NaN,'prediction_regime','','prediction_confidence','','predicted_route','', ...
    'bmd_median_s',NaN,'bmd_q25_s',NaN,'bmd_q75_s',NaN,'sparse_median_s',NaN,'sparse_q25_s',NaN,'sparse_q75_s',NaN, ...
    'actual_ratio_median',NaN,'actual_median_winner','','actual_robust_winner','','routed_correct',NaN, ...
    'false_bmd',0,'false_sparse',0,'numeric_check',0);
end
function r=emptyTrialRow()
r=struct('case_name','','trial',NaN,'order','','bmd_cold_s',NaN,'sparse_s',NaN,'sparse_over_bmd',NaN);
end
