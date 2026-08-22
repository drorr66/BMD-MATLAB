function [results, summary, trialsTable] = run_cold_crossover_v05(varargin)
%RUN_COLD_CROSSOVER_V05 Focused cold-multiply crossover: BMD vs sparse.
%
% Fixes the structural family used in v0.4 at:
%   inner_n = 255, stride = 2^10
% and varies only the number of repeated blocks:
%   96, 128, 160, 192, 256, 384, 512, 1024.
%
% The timed BMD operation is always genuinely cold: each trial starts from
% a fresh compact manager containing only the two operands, so the multiply
% cache is empty. Sparse operands are pre-built and only multiplication is
% timed.  Multiple independent trials are retained to distinguish a real
% crossover from timer noise.

p = inputParser;
addParameter(p,'Trials',9,@(x)isnumeric(x)&&isscalar(x)&&x>=5);
addParameter(p,'SaveResults',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});
numTrials = round(p.Results.Trials);

innerN = 255;
stridePower = 10;
blocksList = [96 128 160 192 256 384 512 1024];

rows = repmat(emptyRow(),numel(blocksList),1);
raw = repmat(emptyTrialRow(),numel(blocksList)*numTrials,1);
rawPos = 0;

fprintf('\nBMD-MATLAB v0.5 focused cold-multiply crossover\n');
fprintf('=================================================\n');
fprintf('Fixed family: inner_n=%d, stride=2^%d, trials=%d\n\n',innerN,stridePower,numTrials);

for ii = 1:numel(blocksList)
    blocks = blocksList(ii);
    row = emptyRow();
    row.inner_n = innerN;
    row.stride_power = stridePower;
    row.blocks = blocks;
    row.logical_terms = (innerN+1)*blocks;
    fprintf('  blocks=%-4d logical_terms=%-7d ... ',blocks,row.logical_terms);

    try
        % Build operands once, then compact them.  Compaction is deliberately
        % outside the multiplication timing.
        baseMgr = BMDManager();
        a = baseMgr.geometricSum(innerN);
        b = baseMgr.geometricSumShifted(blocks-1,stridePower+1);
        [baseMgr,refs] = baseMgr.compact([a;b]);
        baseStats = baseMgr.stats();
        row.bmd_operand_nodes = baseStats.total_internal_nodes;

        sparseA = build_sparse_family('geometric_sum',innerN,0);
        sparseB = build_sparse_shifted_geometric(blocks,stridePower);
        row.sparse_terms_a = numel(sparseA.exponents);
        row.sparse_terms_b = numel(sparseB.exponents);

        % Untimed JIT/allocation warm-up on throwaway state.
        [warmMgr,warmRefs] = baseMgr.compact(refs);
        warmMgr.multiply(warmRefs(1,:),warmRefs(2,:));
        sparse_terms_multiply(sparseA,sparseB);

        bmdTimes = zeros(1,numTrials);
        sparseTimes = zeros(1,numTrials);
        firstProd = [];
        firstMgr = [];

        % Alternate measurement order to reduce systematic order bias.
        for kk = 1:numTrials
            if mod(kk,2)==1
                [bmdTimes(kk),mTmp,pTmp] = measureBmdCold(baseMgr,refs);
                [sparseTimes(kk),sTmp] = measureSparse(sparseA,sparseB); %#ok<NASGU>
            else
                [sparseTimes(kk),sTmp] = measureSparse(sparseA,sparseB); %#ok<NASGU>
                [bmdTimes(kk),mTmp,pTmp] = measureBmdCold(baseMgr,refs);
            end
            if kk==1
                firstMgr = mTmp;
                firstProd = pTmp;
            end

            rawPos = rawPos + 1;
            rr = emptyTrialRow();
            rr.inner_n = innerN;
            rr.stride_power = stridePower;
            rr.blocks = blocks;
            rr.trial = kk;
            if mod(kk,2)==1, rr.order = 'BMD_FIRST'; else, rr.order = 'SPARSE_FIRST'; end
            rr.bmd_cold_s = bmdTimes(kk);
            rr.sparse_s = sparseTimes(kk);
            rr.sparse_over_bmd = sparseTimes(kk)/max(realmin,bmdTimes(kk));
            raw(rawPos) = rr;
        end

        % Correctness and structural diagnostics from the first cold trial.
        direct = firstMgr.geometricGrid(innerN,blocks,stridePower);
        if ~firstMgr.same(firstProd,direct)
            error('BMD:V05Canonical','Cold generic multiply != direct canonical grid.');
        end
        after = firstMgr.stats(firstProd);
        row.bmd_result_nodes = after.reachable_nodes;
        row.bmd_new_workspace_nodes = after.total_internal_nodes - row.bmd_operand_nodes;
        row.bmd_mul_cache_entries = after.mul_cache_entries;

        sparseProd = sparse_terms_multiply(sparseA,sparseB);
        expectedTerms = (innerN+1)*blocks;
        if numel(sparseProd.exponents) ~= expectedTerms || any(sparseProd.coefficients ~= 1)
            error('BMD:V05SparseStructure','Sparse product structure differs from expected disjoint grid.');
        end
        xv = 0.999;
        yb = firstMgr.evaluate(firstProd,xv);
        ys = sparse_terms_evaluate(sparseProd,xv);
        row.numeric_check = near(yb,ys);
        if ~row.numeric_check
            error('BMD:V05Numeric','BMD and sparse evaluations differ.');
        end

        % Robust timing statistics.
        row.bmd_median_s = median(bmdTimes);
        row.bmd_q25_s = percentileLinear(bmdTimes,0.25);
        row.bmd_q75_s = percentileLinear(bmdTimes,0.75);
        row.bmd_min_s = min(bmdTimes);
        row.bmd_max_s = max(bmdTimes);
        row.bmd_cv = std(bmdTimes)/max(realmin,mean(bmdTimes));

        row.sparse_median_s = median(sparseTimes);
        row.sparse_q25_s = percentileLinear(sparseTimes,0.25);
        row.sparse_q75_s = percentileLinear(sparseTimes,0.75);
        row.sparse_min_s = min(sparseTimes);
        row.sparse_max_s = max(sparseTimes);
        row.sparse_cv = std(sparseTimes)/max(realmin,mean(sparseTimes));

        ratios = sparseTimes./max(realmin,bmdTimes);
        row.ratio_median = row.sparse_median_s/max(realmin,row.bmd_median_s);
        row.paired_ratio_median = median(ratios);
        % Conservative interval: compare the faster-looking quartile of the
        % denominator against the slower-looking quartile of the numerator.
        row.robust_ratio_low = row.sparse_q25_s/max(realmin,row.bmd_q75_s);
        row.robust_ratio_high = row.sparse_q75_s/max(realmin,row.bmd_q25_s);

        if row.robust_ratio_low > 1
            row.robust_winner = 'BMD';
        elseif row.robust_ratio_high < 1
            row.robust_winner = 'SPARSE';
        else
            row.robust_winner = 'OVERLAP';
        end
        row.status = 'OK';
    catch ME
        row.status = ['ERROR:' ME.identifier];
    end
    rows(ii) = row;
    fprintf('%s',row.status);
    if strcmp(row.status,'OK')
        fprintf(' ratio=%.3fx robust=[%.3f, %.3f] %s', ...
            row.ratio_median,row.robust_ratio_low,row.robust_ratio_high,row.robust_winner);
    end
    fprintf('\n');
end

results = struct2table(rows);
trialsTable = struct2table(raw(1:rawPos));
summary = makeSummary(results);

disp(results(:,{'blocks','bmd_median_s','sparse_median_s','ratio_median','robust_ratio_low','robust_ratio_high','robust_winner','status'}));
fprintf('\nCrossover summary:\n');
disp(summary);

if p.Results.SaveResults
    outDir = fullfile(fileparts(mfilename('fullpath')),'results');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    writetable(results,fullfile(outDir,'cold_crossover_results_v05.csv'));
    writetable(trialsTable,fullfile(outDir,'cold_crossover_trials_v05.csv'));
    writetable(summary,fullfile(outDir,'cold_crossover_summary_v05.csv'));
    writeMetadata(fullfile(outDir,'run_metadata_v05.txt'),numTrials,innerN,stridePower,blocksList);
end
end

function [t,m,out] = measureBmdCold(baseMgr,refs)
[m,r] = baseMgr.compact(refs);
tic;
out = m.multiply(r(1,:),r(2,:));
t = toc;
% Force observable access after timing without contaminating the timed body.
if out(2) < 1, error('BMD:V05Sink','Invalid BMD result.'); end
end

function [t,out] = measureSparse(a,b)
tic;
out = sparse_terms_multiply(a,b);
t = toc;
if isempty(out.exponents), error('BMD:V05Sink','Unexpected empty sparse result.'); end
end

function q = percentileLinear(x,p)
x = sort(x(:));
if isempty(x), q = NaN; return; end
if numel(x)==1, q = x(1); return; end
pos = 1 + (numel(x)-1)*p;
lo = floor(pos); hi = ceil(pos);
if lo==hi
    q = x(lo);
else
    q = x(lo) + (pos-lo)*(x(hi)-x(lo));
end
end

function s = makeSummary(results)
s = struct('comparison','sparse_over_bmd_cold','median_status','', ...
    'median_n_low',NaN,'median_n_high',NaN,'median_ratio_low',NaN,'median_ratio_high',NaN,'median_crossover_estimate',NaN, ...
    'robust_status','', 'robust_n_low',NaN,'robust_n_high',NaN,'robust_low_at_low',NaN,'robust_low_at_high',NaN,'robust_crossover_estimate',NaN, ...
    'sustained_median_bmd_from',NaN,'sustained_robust_bmd_from',NaN);

ok = strcmp(results.status,'OK');
n = results.blocks(ok); r = results.ratio_median(ok); rl = results.robust_ratio_low(ok);
[s.median_status,s.median_n_low,s.median_n_high,s.median_ratio_low,s.median_ratio_high,s.median_crossover_estimate] = bracket(n,r);
[s.robust_status,s.robust_n_low,s.robust_n_high,s.robust_low_at_low,s.robust_low_at_high,s.robust_crossover_estimate] = bracket(n,rl);
s.sustained_median_bmd_from = sustainedFrom(n,r>1);
s.sustained_robust_bmd_from = sustainedFrom(n,rl>1);
s = struct2table(s);
end

function [status,nlo,nhi,rlo,rhi,est] = bracket(n,r)
status='NO_CROSSING'; nlo=NaN; nhi=NaN; rlo=NaN; rhi=NaN; est=NaN;
for k=2:numel(n)
    if r(k-1)<1 && r(k)>=1
        status='BRACKETED'; nlo=n(k-1); nhi=n(k); rlo=r(k-1); rhi=r(k);
        if rlo>0 && rhi>0 && rhi~=rlo
            % Log-linear interpolation in block count and ratio.
            a = (0-log(rlo))/(log(rhi)-log(rlo));
            est = exp(log(nlo)+a*(log(nhi)-log(nlo)));
        end
        return;
    end
end
if all(r>1), status='BMD_FASTER_ALL'; elseif all(r<1), status='SPARSE_FASTER_ALL'; end
end

function n0 = sustainedFrom(n,win)
n0=NaN;
for k=1:numel(n)
    if all(win(k:end))
        n0=n(k); return;
    end
end
end

function writeMetadata(path,numTrials,innerN,stridePower,blocksList)
fid=fopen(path,'w');
if fid<0, return; end
c=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'BMD-MATLAB v0.5 focused cold crossover\n');
fprintf(fid,'Generated: %s\n',datestr(now,31));
fprintf(fid,'MATLAB version: %s\n',version);
fprintf(fid,'Computer: %s\n',computer);
fprintf(fid,'Trials per case: %d\n',numTrials);
fprintf(fid,'inner_n: %d\n',innerN);
fprintf(fid,'stride_power: %d\n',stridePower);
fprintf(fid,'blocks: %s\n',mat2str(blocksList));
fprintf(fid,'Winner rule: ratio=sparse_time/BMD_time; >1 means BMD faster.\n');
fprintf(fid,'Robust BMD win: sparse Q25 / BMD Q75 > 1.\n');
end

function tf = near(a,b)
tf = abs(a-b) <= 1e-8*max(1,max(abs(a),abs(b)));
end

function r = emptyRow()
r = struct('inner_n',NaN,'stride_power',NaN,'blocks',NaN,'logical_terms',NaN,'status','', ...
    'bmd_median_s',NaN,'bmd_q25_s',NaN,'bmd_q75_s',NaN,'bmd_min_s',NaN,'bmd_max_s',NaN,'bmd_cv',NaN, ...
    'sparse_median_s',NaN,'sparse_q25_s',NaN,'sparse_q75_s',NaN,'sparse_min_s',NaN,'sparse_max_s',NaN,'sparse_cv',NaN, ...
    'ratio_median',NaN,'paired_ratio_median',NaN,'robust_ratio_low',NaN,'robust_ratio_high',NaN,'robust_winner','', ...
    'bmd_operand_nodes',NaN,'bmd_result_nodes',NaN,'bmd_new_workspace_nodes',NaN,'bmd_mul_cache_entries',NaN, ...
    'sparse_terms_a',NaN,'sparse_terms_b',NaN,'numeric_check',false);
end

function r = emptyTrialRow()
r = struct('inner_n',NaN,'stride_power',NaN,'blocks',NaN,'trial',NaN,'order','', ...
    'bmd_cold_s',NaN,'sparse_s',NaN,'sparse_over_bmd',NaN);
end
