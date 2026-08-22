function results = run_structural_operations_v04(varargin)
%RUN_STRUCTURAL_OPERATIONS_V04 First/repeated operations with shared blocks.
%
% Measures generic BMD multiplication of two compact operands with disjoint
% bit ranges, exact-repeat cache benefit, related-repeat cache reuse, and
% canonical equality. Dense and sparse baselines use already-built operands.
p=inputParser;
addParameter(p,'Quick',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'ColdRepeats',3,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'SaveResults',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:}); quick=logical(p.Results.Quick); reps=round(p.Results.ColdRepeats);
if quick
    cases=[31 7 32; 63 8 64; 255 10 64; 255 10 128];
else
    cases=[15 6 8; 31 7 16; 31 7 32; 31 7 64; 63 8 32; 63 8 64; 63 8 128; 255 10 32; 255 10 64; 255 10 128];
end
rows=repmat(emptyRow(),size(cases,1),1);
fprintf('\nBMD-MATLAB v0.4 structural operations\n');
fprintf('=====================================\n');
for ii=1:size(cases,1)
    innerN=cases(ii,1); stridePower=cases(ii,2); blocks=cases(ii,3);
    row=emptyRow(); row.inner_n=innerN; row.stride_power=stridePower; row.blocks=blocks;
    fprintf('  inner=%-4g blocks=%-4g stride=2^%-2g ... ',innerN,blocks,stridePower);
    try
        [baseMgr,baseRefs,denseA,denseB,sparseA,sparseB]=fixture(innerN,blocks,stridePower);
        [row.bmd_cold_multiply_s,detail]=coldMultiply(baseMgr,baseRefs,reps);
        row.bmd_operand_nodes=detail.operand_nodes; row.bmd_result_nodes=detail.result_nodes;
        row.bmd_new_workspace_nodes=detail.new_workspace_nodes; row.bmd_mul_cache_entries=detail.mul_cache_entries;

        [wm,wr]=baseMgr.compact(baseRefs); wm.multiply(wr(1,:),wr(2,:));
        row.bmd_warm_exact_multiply_s=timeit(@() bmdMulSink(wm,wr));
        row.dense_multiply_s=timeit(@() denseMulSink(denseA,denseB));
        row.sparse_multiply_s=timeit(@() sparseMulSink(sparseA,sparseB));
        row.dense_over_bmd_cold=row.dense_multiply_s/max(realmin,row.bmd_cold_multiply_s);
        row.sparse_over_bmd_cold=row.sparse_multiply_s/max(realmin,row.bmd_cold_multiply_s);
        row.dense_over_bmd_warm=row.dense_multiply_s/max(realmin,row.bmd_warm_exact_multiply_s);
        row.sparse_over_bmd_warm=row.sparse_multiply_s/max(realmin,row.bmd_warm_exact_multiply_s);

        row.bmd_related_cold_s=relatedMultiply(innerN,blocks,stridePower,reps,false);
        row.bmd_related_reuse_s=relatedMultiply(innerN,blocks,stridePower,reps,true);
        row.related_cache_speedup=row.bmd_related_cold_s/max(realmin,row.bmd_related_reuse_s);

        % Correctness of first multiplication against both baselines.
        [cm,cr]=baseMgr.compact(baseRefs); prod=cm.multiply(cr(1,:),cr(2,:));
        direct=cm.geometricGrid(innerN,blocks,stridePower);
        if ~cm.same(prod,direct), error('BMD:GridCanonical','Generic product != direct canonical grid.'); end
        dprod=conv(trimLeading(denseA),trimLeading(denseB));
        sprod=sparse_terms_multiply(sparseA,sparseB);
        xv=0.999; yb=cm.evaluate(prod,xv); yd=polyval(dprod,xv); ys=sparse_terms_evaluate(sprod,xv);
        row.numeric_check=near(yb,yd)&&near(yb,ys);

        % Equality: two independent baseline materializations vs canonical refs.
        eqm=BMDManager(); ea=eqm.geometricSum(innerN); eb=eqm.geometricSumShifted(blocks-1,stridePower+1);
        ep=eqm.multiply(ea,eb); ed=eqm.geometricGrid(innerN,blocks,stridePower);
        assert(eqm.same(ep,ed));
        denseEqA=build_dense_grid(innerN,blocks,stridePower); denseEqB=build_dense_grid(innerN,blocks,stridePower);
        sparseEqA=build_sparse_grid(innerN,blocks,stridePower); sparseEqB=build_sparse_grid(innerN,blocks,stridePower);
        row.bmd_equality_s=timeit(@() eqm.same(ep,ed));
        row.dense_equality_s=timeit(@() isequal(denseEqA,denseEqB));
        row.sparse_equality_s=timeit(@() sparse_terms_same(sparseEqA,sparseEqB));
        row.dense_over_bmd_equality=row.dense_equality_s/max(realmin,row.bmd_equality_s);
        row.sparse_over_bmd_equality=row.sparse_equality_s/max(realmin,row.bmd_equality_s);
        if ~row.numeric_check, error('BMD:GridNumeric','Structural multiply validation failed.'); end
        row.status='OK';
    catch ME
        row.status=['ERROR:' ME.identifier];
    end
    rows(ii)=row; fprintf('%s\n',row.status);
end
results=struct2table(rows); disp(results);
if p.Results.SaveResults
    outDir=fullfile(fileparts(mfilename('fullpath')),'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
    writetable(results,fullfile(outDir,'structural_operations_v04.csv'));
end
end

function [m,refs,dA,dB,sA,sB]=fixture(innerN,blocks,stridePower)
m=BMDManager(); a=m.geometricSum(innerN); b=m.geometricSumShifted(blocks-1,stridePower+1); refs=[a;b];
dA=ones(1,innerN+1); dB=build_dense_shifted_geometric(blocks,stridePower);
sA=build_sparse_family('geometric_sum',innerN,0); sB=build_sparse_shifted_geometric(blocks,stridePower);
end

function [t,detail]=coldMultiply(baseMgr,baseRefs,reps)
times=zeros(1,reps); detail=struct('operand_nodes',NaN,'result_nodes',NaN,'new_workspace_nodes',NaN,'mul_cache_entries',NaN);
for kk=1:reps
    [m,r]=baseMgr.compact(baseRefs); before=m.stats();
    tic; out=m.multiply(r(1,:),r(2,:)); times(kk)=toc;
    if kk==1
        after=m.stats(out); detail.operand_nodes=before.total_internal_nodes;
        detail.result_nodes=after.reachable_nodes; detail.new_workspace_nodes=after.total_internal_nodes-before.total_internal_nodes;
        detail.mul_cache_entries=after.mul_cache_entries;
    end
end
t=median(times);
end

function t=relatedMultiply(innerN,blocks,stridePower,reps,usePriorCache)
times=zeros(1,reps);
for kk=1:reps
    m=BMDManager(); a=m.geometricSum(innerN); b=m.geometricSumShifted(blocks-1,stridePower+1);
    if usePriorCache, m.multiply(a,b); end
    bplus=m.geometricSumShifted(blocks,stridePower+1);
    tic; m.multiply(a,bplus); times(kk)=toc;
end
t=median(times);
end

function y=bmdMulSink(m,r), out=m.multiply(r(1,:),r(2,:)); y=double(out(1))+double(out(2)); end
function y=denseMulSink(a,b)
persistent out
out=conv(trimLeading(a),trimLeading(b)); y=double(numel(out))+out(1)+out(end);
end
function y=sparseMulSink(a,b)
persistent out
out=sparse_terms_multiply(a,b); y=double(numel(out.exponents))+out.coefficients(1)+out.coefficients(end);
end
function p=trimLeading(p)
idx=find(p~=0,1,'first');
if isempty(idx)
    p=0;
else
    p=p(idx:end);
end
end
function tf=near(a,b), tf=abs(a-b)<=1e-8*max(1,max(abs(a),abs(b))); end
function r=emptyRow()
r=struct('inner_n',NaN,'stride_power',NaN,'blocks',NaN,'status','', ...
    'bmd_cold_multiply_s',NaN,'bmd_warm_exact_multiply_s',NaN,'dense_multiply_s',NaN,'sparse_multiply_s',NaN, ...
    'dense_over_bmd_cold',NaN,'sparse_over_bmd_cold',NaN,'dense_over_bmd_warm',NaN,'sparse_over_bmd_warm',NaN, ...
    'bmd_related_cold_s',NaN,'bmd_related_reuse_s',NaN,'related_cache_speedup',NaN, ...
    'bmd_operand_nodes',NaN,'bmd_result_nodes',NaN,'bmd_new_workspace_nodes',NaN,'bmd_mul_cache_entries',NaN, ...
    'bmd_equality_s',NaN,'dense_equality_s',NaN,'sparse_equality_s',NaN, ...
    'dense_over_bmd_equality',NaN,'sparse_over_bmd_equality',NaN,'numeric_check',false);
end
