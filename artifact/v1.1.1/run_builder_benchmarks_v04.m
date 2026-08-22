function results = run_builder_benchmarks_v04(varargin)
%RUN_BUILDER_BENCHMARKS_V04 Direct bit-recursive vs naive geometric builder.
%
% The represented polynomial is identical in every row: 1+x+...+x^n.
% The experiment isolates construction strategy.  v0.3 built by adding
% n+1 monomials sequentially; v0.4 recursively partitions exponents by
% binary digits and should create only reachable nodes.
p=inputParser;
addParameter(p,'Quick',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'SaveResults',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:}); quick=logical(p.Results.Quick);
if quick, ns=[100 1000 5000 10000]; else, ns=[100 500 1000 2500 5000 10000 20000]; end
rows=repmat(emptyRow(),numel(ns),1);
fprintf('\nBMD-MATLAB v0.4 builder benchmarks\n');
fprintf('===================================\n');
for ii=1:numel(ns)
    n=ns(ii); row=emptyRow(); row.n=n;
    fprintf('  geometric n=%-8g ... ',n);
    try
        row.direct_time_s=timeit(@() builderSink('geometric_sum',n));
        [md,rd]=build_bmd_family('geometric_sum',n,0); sd=md.stats(rd);
        row.direct_reachable_nodes=sd.reachable_nodes;
        row.direct_workspace_nodes=sd.total_internal_nodes;
        row.direct_garbage_nodes=sd.garbage_internal_nodes;

        row.naive_time_s=timeit(@() builderSink('geometric_sum_naive',n));
        [mn,rn]=build_bmd_family('geometric_sum_naive',n,0); sn=mn.stats(rn);
        row.naive_reachable_nodes=sn.reachable_nodes;
        row.naive_workspace_nodes=sn.total_internal_nodes;
        row.naive_garbage_nodes=sn.garbage_internal_nodes;

        row.naive_over_direct_time=row.naive_time_s/max(realmin,row.direct_time_s);
        row.workspace_reduction=sn.total_internal_nodes/max(1,sd.total_internal_nodes);
        row.garbage_eliminated=sn.garbage_internal_nodes-sd.garbage_internal_nodes;
        xv=0.9375;
        row.numeric_check=near(md.evaluate(rd,xv),mn.evaluate(rn,xv));
        row.node_check=(sd.reachable_nodes==sn.reachable_nodes);
        if ~(row.numeric_check && row.node_check && sd.garbage_internal_nodes==0)
            error('BMD:BuilderCheck','Direct/naive builder validation failed.');
        end
        row.status='OK';
    catch ME
        row.status=['ERROR:' ME.identifier];
    end
    rows(ii)=row; fprintf('%s\n',row.status);
end
results=struct2table(rows); disp(results);
if p.Results.SaveResults
    outDir=fullfile(fileparts(mfilename('fullpath')),'results');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    writetable(results,fullfile(outDir,'builder_results_v04.csv'));
end
end

function y=builderSink(family,n)
persistent heldMgr heldRoot
[heldMgr,heldRoot]=build_bmd_family(family,n,0);
st=heldMgr.stats(heldRoot);
y=double(st.total_internal_nodes)+double(heldRoot(1))+double(heldRoot(2));
end
function tf=near(a,b), tf=abs(a-b)<=1e-12*max(1,max(abs(a),abs(b))); end
function r=emptyRow()
r=struct('n',NaN,'status','', ...
    'direct_time_s',NaN,'direct_reachable_nodes',NaN,'direct_workspace_nodes',NaN,'direct_garbage_nodes',NaN, ...
    'naive_time_s',NaN,'naive_reachable_nodes',NaN,'naive_workspace_nodes',NaN,'naive_garbage_nodes',NaN, ...
    'naive_over_direct_time',NaN,'workspace_reduction',NaN,'garbage_eliminated',NaN, ...
    'numeric_check',false,'node_check',false);
end
