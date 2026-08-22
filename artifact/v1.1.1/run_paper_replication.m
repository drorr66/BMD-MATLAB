function results = run_paper_replication(varargin)
%RUN_PAPER_REPLICATION Reproduce the structural benchmark behind paper Table 2.
%
% IMPORTANT: Table 2 uses x_k (subscript; distinct variables), not x^k.
% It reports prod_{k=1}^n (x_k+1)^r for r=1,4,8, with final *BMD node
% counts r*n+1. This script checks that identity and measures construction
% time in the MATLAB prototype.
%
% The v0.4 exact-weight guard (+/-flintmax) limits large n for r=4 and r=8.
% r=1 has unit coefficients and can scale much further.

p = inputParser;
addParameter(p,'Quick',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'Repeats',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'SaveResults',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});
cfg = p.Results;

if cfg.Quick
    plan = {
        1, 10; 1, 25; 1, 50; 1, 100;
        4, 2;  4, 4;  4, 6;  4, 8;
        8, 1;  8, 2;  8, 3;  8, 5;
        };
else
    plan = {
        1, 50; 1,100; 1,200; 1,300; 1,500;
        4, 2; 4,4; 4,6; 4,8; 4,10; 4,12;
        8, 1; 8,2; 8,3; 8,4; 8,5; 8,6; 8,7;
        };
end

rows = repmat(emptyRow(),0,1);
fprintf('Paper Table 2 structural replication\n');
fprintf('MATLAB/engine: %s\n\n',version);

for i=1:size(plan,1)
    rpow = plan{i,1};
    n = plan{i,2};
    row = emptyRow(); row.r = rpow; row.n = n;
    row.expected_final_nodes = rpow*n + 1;
    fprintf('[%2d/%2d] r=%d n=%d ... ',i,size(plan,1),rpow,n);
    try
        [mgr,root] = buildCase(rpow,n);
        st = mgr.stats(root);
        row.status = 'OK';
        row.final_nodes = st.reachable_nodes + 1;
        row.total_internal_nodes = st.total_internal_nodes;
        row.node_count_matches_paper = (row.final_nodes == row.expected_final_nodes);
        row.mul_cache_entries = st.mul_cache_entries;
        row.time_s = medianTime(@() sink(rpow,n),cfg.Repeats);
        fprintf('%.4gs, final nodes=%d (expected %d)\n',row.time_s,row.final_nodes,row.expected_final_nodes);
    catch ME
        row.status = [ME.identifier ': ' ME.message];
        fprintf('STOP: %s\n',row.status);
    end
    rows(end+1,1)=row; %#ok<AGROW>
end

results=struct2table(rows);
disp(results);

if cfg.SaveResults
    outDir=fullfile(fileparts(mfilename('fullpath')),'results');
    if ~exist(outDir,'dir'),mkdir(outDir);end
    writetable(results,fullfile(outDir,'paper_table2_replication.csv'));
    save(fullfile(outDir,'paper_table2_replication.mat'),'results','cfg');
end
end

function [mgr,root]=buildCase(rpow,n)
mgr=BMDManager();
root=mgr.one();
for k=1:n
    base=mgr.add(mgr.variable(k),mgr.one());
    factor=mgr.power(base,rpow);
    root=mgr.multiply(root,factor);
end
end

function y=sink(rpow,n)
[m,r]=buildCase(rpow,n);
st=m.stats(r);
y=double(r(1))+double(r(2))+st.reachable_nodes;
end

function t=medianTime(f,repeats)
times=zeros(1,repeats);
for k=1:repeats
    tic; z=f(); %#ok<NASGU>
    times(k)=toc;
end
t=median(times);
end

function r=emptyRow()
r=struct('r',NaN,'n',NaN,'status','', ...
    'expected_final_nodes',NaN,'final_nodes',NaN, ...
    'node_count_matches_paper',false,'total_internal_nodes',NaN, ...
    'mul_cache_entries',NaN,'time_s',NaN);
end
