function [results, summary, trialsTable] = run_threshold_calibration_v08(varargin)
%RUN_THRESHOLD_CALIBRATION_V08 Calibrate cold BMD-vs-sparse closure threshold.
%
% Fixed inputs in every case:
%   P support uses exponent bits 0..7             -> 256 terms, 8 BMD nodes
%   Q support uses eight independent bit choices  -> 256 terms, 8 BMD nodes
%   sparse pair-products                          -> 65,536
%
% The Q bit pattern is chosen so the exact cold BMD multiplication creates
% 8,12,16,...,40 new workspace nodes.  This densely samples the v0.7
% crossover bracket (8 new nodes = BMD win, 38 = sparse win) without
% changing input term counts or individual BMD operand sizes.

p = inputParser;
addParameter(p,'Trials',17,@(x)isnumeric(x)&&isscalar(x)&&x>=9);
addParameter(p,'SaveResults',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});
numTrials = round(p.Results.Trials);

pBits = 0:7;
highSeven = 9:15;
% NaN means fully disjoint eight-bit Q support (bits 9..16).
overlapBits = [NaN 7 6 5 4 3 2 1 0];
expectedNew = [8 12 16 20 24 28 32 36 40];
expectedResult = [16 16 17 18 19 20 21 22 23];
nTerms = 256;
pairProducts = nTerms*nTerms;

pExponents = build_bitcube_support_v08(pBits);
assert(numel(pExponents)==nTerms);

rows = repmat(emptyRow(),numel(expectedNew),1);
raw = repmat(emptyTrialRow(),numel(expectedNew)*numTrials,1);
rawPos = 0;

fprintf('\nBMD-MATLAB v0.8 threshold calibration\n');
fprintf('======================================\n');
fprintf('P terms=%d, Q terms=%d, pair-products=%d, trials=%d\n', ...
    nTerms,nTerms,pairProducts,numTrials);
fprintf('Target new-node grid: [%s]\n\n',strtrim(sprintf('%d ',expectedNew)));

for ii = 1:numel(expectedNew)
    row = emptyRow();
    row.case_index = ii;
    row.target_new_nodes = expectedNew(ii);
    row.expected_result_nodes = expectedResult(ii);
    row.p_terms = nTerms;
    row.q_terms = nTerms;
    row.sparse_pair_products = pairProducts;

    if isnan(overlapBits(ii))
        qBits = 9:16;
        row.overlap_bit = -1;
        row.overlap_count = 0;
    else
        qBits = [overlapBits(ii) highSeven];
        row.overlap_bit = overlapBits(ii);
        row.overlap_count = 1;
    end
    fprintf('  target_new=%-2d overlap_bit=%2d ... ',row.target_new_nodes,row.overlap_bit);

    try
        qExponents = build_bitcube_support_v08(qBits);
        if numel(qExponents) ~= nTerms || numel(unique(qExponents)) ~= nTerms
            error('BMD:V08QCardinality','Q support must contain exactly 256 unique exponents.');
        end

        baseMgr = BMDManager();
        pRef = baseMgr.indicatorExponents(pExponents);
        qRef = baseMgr.indicatorExponents(qExponents);
        pStats = baseMgr.stats(pRef);
        qStats = baseMgr.stats(qRef);
        [baseMgr,refs] = baseMgr.compact([pRef;qRef]);
        unionStats = baseMgr.stats();
        row.bmd_p_nodes = pStats.reachable_nodes;
        row.bmd_q_nodes = qStats.reachable_nodes;
        row.bmd_operand_union_nodes = unionStats.total_internal_nodes;
        if row.bmd_p_nodes ~= 8 || row.bmd_q_nodes ~= 8 || row.bmd_operand_union_nodes ~= 16
            error('BMD:V08OperandSize','Expected P=8 nodes, Q=8 nodes, union=16 nodes.');
        end

        sparseP = sparse_terms(pExponents,ones(1,nTerms));
        sparseQ = sparse_terms(qExponents,ones(1,nTerms));
        if numel(sparseP.exponents) ~= nTerms || numel(sparseQ.exponents) ~= nTerms
            error('BMD:V08SparseCardinality','Sparse input cardinality changed unexpectedly.');
        end

        % Untimed warm-up on throwaway state to avoid first-call/JIT effects.
        [wm,wr] = baseMgr.compact(refs);
        wm.multiply(wr(1,:),wr(2,:));
        sparse_terms_multiply(sparseP,sparseQ);

        bmdTimes = zeros(1,numTrials);
        sparseTimes = zeros(1,numTrials);
        firstMgr = [];
        firstProd = [];
        firstSparseProd = [];

        for kk = 1:numTrials
            if mod(kk,2)==1
                [bmdTimes(kk),mTmp,pTmp] = measureBmdCold(baseMgr,refs);
                [sparseTimes(kk),sTmp] = measureSparse(sparseP,sparseQ);
            else
                [sparseTimes(kk),sTmp] = measureSparse(sparseP,sparseQ);
                [bmdTimes(kk),mTmp,pTmp] = measureBmdCold(baseMgr,refs);
            end
            if kk==1
                firstMgr = mTmp;
                firstProd = pTmp;
                firstSparseProd = sTmp;
            end

            rawPos = rawPos + 1;
            rr = emptyTrialRow();
            rr.case_index = ii;
            rr.target_new_nodes = row.target_new_nodes;
            rr.overlap_bit = row.overlap_bit;
            rr.trial = kk;
            if mod(kk,2)==1, rr.order='BMD_FIRST'; else, rr.order='SPARSE_FIRST'; end
            rr.bmd_cold_s = bmdTimes(kk);
            rr.sparse_s = sparseTimes(kk);
            rr.sparse_over_bmd = sparseTimes(kk)/max(realmin,bmdTimes(kk));
            raw(rawPos) = rr;
        end

        after = firstMgr.stats(firstProd);
        row.bmd_result_nodes = after.reachable_nodes;
        row.bmd_new_workspace_nodes = after.total_internal_nodes - row.bmd_operand_union_nodes;
        row.bmd_mul_cache_entries = after.mul_cache_entries;
        row.sparse_pairs_per_new_bmd_node = pairProducts/max(1,row.bmd_new_workspace_nodes);
        row.sparse_result_terms = numel(firstSparseProd.exponents);
        row.sparse_result_max_coeff = max(abs(firstSparseProd.coefficients));
        row.sparse_result_terms_per_bmd_result_node = row.sparse_result_terms/max(1,row.bmd_result_nodes);

        % Hard structural calibration guard. If these fail, this is no longer
        % the v0.8 experiment that was independently validated.
        if row.bmd_new_workspace_nodes ~= row.target_new_nodes
            error('BMD:V08NewNodes','Expected %d new nodes, observed %d.', ...
                row.target_new_nodes,row.bmd_new_workspace_nodes);
        end
        if row.bmd_result_nodes ~= row.expected_result_nodes
            error('BMD:V08ResultNodes','Expected %d result nodes, observed %d.', ...
                row.expected_result_nodes,row.bmd_result_nodes);
        end

        % Exact coefficient-level check, untimed. Degrees are below 2^16.
        denseBmd = firstMgr.toDense(firstProd,1e6);
        denseSparse = sparse_terms_to_dense(firstSparseProd,1e6);
        row.numeric_check = isequal(denseBmd,denseSparse);
        if ~row.numeric_check
            error('BMD:V08Numeric','Exact BMD and sparse product coefficients differ.');
        end

        row.bmd_median_s = median(bmdTimes);
        row.bmd_q25_s = percentileLinear(bmdTimes,0.25);
        row.bmd_q75_s = percentileLinear(bmdTimes,0.75);
        row.bmd_cv = std(bmdTimes)/max(realmin,mean(bmdTimes));
        row.sparse_median_s = median(sparseTimes);
        row.sparse_q25_s = percentileLinear(sparseTimes,0.25);
        row.sparse_q75_s = percentileLinear(sparseTimes,0.75);
        row.sparse_cv = std(sparseTimes)/max(realmin,mean(sparseTimes));
        row.ratio_median = row.sparse_median_s/max(realmin,row.bmd_median_s);
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

    rows(ii)=row;
    fprintf('%s',row.status);
    if strcmp(row.status,'OK')
        fprintf(' result=%-2d ratio=%6.3fx robust=%s', ...
            row.bmd_result_nodes,row.ratio_median,row.robust_winner);
    end
    fprintf('\n');
end

results = struct2table(rows);
trialsTable = struct2table(raw(1:rawPos));
summary = makeSummary(results);

disp(results(:,{'target_new_nodes','overlap_bit','bmd_result_nodes','sparse_result_terms', ...
    'bmd_median_s','sparse_median_s','ratio_median','robust_winner','status'}));
fprintf('\nThreshold calibration summary:\n');
disp(summary);

if p.Results.SaveResults
    outDir = fullfile(fileparts(mfilename('fullpath')),'results');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    writetable(results,fullfile(outDir,'threshold_calibration_results_v08.csv'));
    writetable(trialsTable,fullfile(outDir,'threshold_calibration_trials_v08.csv'));
    writetable(summary,fullfile(outDir,'threshold_calibration_summary_v08.csv'));
    writeMetadata(fullfile(outDir,'run_metadata_v08.txt'),numTrials,pBits,overlapBits,expectedNew,expectedResult,pairProducts);
end
end

function [t,m,out] = measureBmdCold(baseMgr,refs)
[m,r] = baseMgr.compact(refs);
tic;
out = m.multiply(r(1,:),r(2,:));
t = toc;
if out(2) < 1, error('BMD:V08Sink','Invalid BMD result.'); end
end

function [t,out] = measureSparse(a,b)
tic;
out = sparse_terms_multiply(a,b);
t = toc;
if isempty(out.exponents), error('BMD:V08Sink','Unexpected empty sparse result.'); end
end

function q = percentileLinear(x,p)
x=sort(x(:));
if isempty(x), q=NaN; return; end
if numel(x)==1, q=x(1); return; end
pos=1+(numel(x)-1)*p; lo=floor(pos); hi=ceil(pos);
if lo==hi, q=x(lo); else, q=x(lo)+(pos-lo)*(x(hi)-x(lo)); end
end

function s = makeSummary(results)
s = struct('comparison','sparse_over_bmd_cold', ...
    'median_transition_status','', ...
    'median_last_bmd_new_nodes',NaN,'median_first_sparse_new_nodes',NaN, ...
    'median_ratio_bmd_side',NaN,'median_ratio_sparse_side',NaN,'median_new_nodes_estimate',NaN, ...
    'robust_transition_status','', ...
    'robust_last_bmd_new_nodes',NaN,'robust_first_sparse_new_nodes',NaN, ...
    'robust_ratio_bmd_side',NaN,'robust_ratio_sparse_side',NaN,'robust_new_nodes_estimate',NaN, ...
    'median_last_bmd_pairs_per_new_node',NaN,'median_first_sparse_pairs_per_new_node',NaN, ...
    'robust_last_bmd_pairs_per_new_node',NaN,'robust_first_sparse_pairs_per_new_node',NaN, ...
    'corr_log_bmd_time_vs_log_new_nodes',NaN, ...
    'corr_log_ratio_vs_log_pairs_per_new_node',NaN, ...
    'bmd_time_exponent_vs_new_nodes',NaN, ...
    'bmd_linear_intercept_s',NaN,'bmd_linear_s_per_new_node',NaN,'bmd_linear_r2',NaN, ...
    'median_sparse_time_s',NaN,'linear_predicted_crossover_new_nodes',NaN, ...
    'max_bmd_speedup',NaN,'max_bmd_speedup_new_nodes',NaN);

ok=strcmp(results.status,'OK');
n=double(results.bmd_new_workspace_nodes(ok));
ratio=double(results.ratio_median(ok));
robust=double(results.robust_ratio_low(ok));
pairratio=double(results.sparse_pairs_per_new_bmd_node(ok));
bmdt=double(results.bmd_median_s(ok));
spt=double(results.sparse_median_s(ok));
[n,ord]=sort(n); ratio=ratio(ord); robust=robust(ord); pairratio=pairratio(ord); bmdt=bmdt(ord); spt=spt(ord);

[s.median_transition_status,s.median_last_bmd_new_nodes,s.median_first_sparse_new_nodes, ...
    s.median_ratio_bmd_side,s.median_ratio_sparse_side,s.median_new_nodes_estimate,idxMed] = fallingBracket(n,ratio);
[s.robust_transition_status,s.robust_last_bmd_new_nodes,s.robust_first_sparse_new_nodes, ...
    s.robust_ratio_bmd_side,s.robust_ratio_sparse_side,s.robust_new_nodes_estimate,idxRob] = fallingBracket(n,robust);
if ~isnan(idxMed)
    s.median_last_bmd_pairs_per_new_node=pairratio(idxMed-1);
    s.median_first_sparse_pairs_per_new_node=pairratio(idxMed);
end
if ~isnan(idxRob)
    s.robust_last_bmd_pairs_per_new_node=pairratio(idxRob-1);
    s.robust_first_sparse_pairs_per_new_node=pairratio(idxRob);
end

s.corr_log_bmd_time_vs_log_new_nodes=safeLogCorr(n,bmdt);
s.corr_log_ratio_vs_log_pairs_per_new_node=safeLogCorr(pairratio,ratio);
s.bmd_time_exponent_vs_new_nodes=logSlope(n,bmdt);
if numel(n)>=3
    pp=polyfit(n,bmdt,1);
    s.bmd_linear_s_per_new_node=pp(1); s.bmd_linear_intercept_s=pp(2);
    pred=polyval(pp,n);
    ssres=sum((bmdt-pred).^2); sstot=sum((bmdt-mean(bmdt)).^2);
    if sstot>0, s.bmd_linear_r2=1-ssres/sstot; end
end
s.median_sparse_time_s=median(spt);
if isfinite(s.bmd_linear_s_per_new_node) && s.bmd_linear_s_per_new_node>0
    s.linear_predicted_crossover_new_nodes=(s.median_sparse_time_s-s.bmd_linear_intercept_s)/s.bmd_linear_s_per_new_node;
end
if ~isempty(ratio)
    [s.max_bmd_speedup,ix]=max(ratio); s.max_bmd_speedup_new_nodes=n(ix);
end
s=struct2table(s);
end

function [status,nlo,nhi,rlo,rhi,est,idx] = fallingBracket(n,r)
status='NO_CROSSING'; nlo=NaN; nhi=NaN; rlo=NaN; rhi=NaN; est=NaN; idx=NaN;
for k=2:numel(n)
    if r(k-1)>=1 && r(k)<1
        status='BRACKETED'; nlo=n(k-1); nhi=n(k); rlo=r(k-1); rhi=r(k); idx=k;
        if rlo>0 && rhi>0 && rhi~=rlo
            a=(0-log(rlo))/(log(rhi)-log(rlo));
            est=nlo+a*(nhi-nlo);
        end
        return;
    end
end
if all(r>=1), status='BMD_FASTER_ALL'; elseif all(r<1), status='SPARSE_FASTER_ALL'; end
end

function c=safeLogCorr(x,y)
mask=isfinite(x)&isfinite(y)&x>0&y>0; x=x(mask); y=y(mask);
if numel(x)<3 || numel(unique(x))<2 || numel(unique(y))<2, c=NaN; return; end
C=corrcoef(log(x),log(y)); c=C(1,2);
end

function a=logSlope(x,y)
mask=isfinite(x)&isfinite(y)&x>0&y>0; x=x(mask); y=y(mask);
if numel(x)<3 || numel(unique(x))<2, a=NaN; return; end
pp=polyfit(log(x),log(y),1); a=pp(1);
end

function writeMetadata(path,numTrials,pBits,overlapBits,expectedNew,expectedResult,pairProducts)
fid=fopen(path,'w'); if fid<0, return; end
c=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'BMD-MATLAB v0.8 cold threshold calibration\n');
fprintf(fid,'Generated: %s\n',datestr(now,31));
fprintf(fid,'MATLAB version: %s\n',version);
fprintf(fid,'Computer: %s\n',computer);
fprintf(fid,'Trials per case: %d\n',numTrials);
fprintf(fid,'P bits: [%s]\n',strtrim(sprintf('%d ',pBits)));
fprintf(fid,'P terms: 256\nQ terms: 256\n');
fprintf(fid,'pair products fixed: %d\n',pairProducts);
fprintf(fid,'Q high bits for one-overlap cases: [9 10 11 12 13 14 15]\n');
fprintf(fid,'Overlap bits (NaN=fully disjoint Q bits 9..16): [%s]\n',strtrim(sprintf('%g ',overlapBits)));
fprintf(fid,'Expected new nodes: [%s]\n',strtrim(sprintf('%d ',expectedNew)));
fprintf(fid,'Expected result nodes: [%s]\n',strtrim(sprintf('%d ',expectedResult)));
fprintf(fid,'Winner ratio = sparse_time/BMD_time; >1 means BMD faster.\n');
fprintf(fid,'Robust BMD win: sparse Q25 / BMD Q75 > 1.\n');
fprintf(fid,'Robust sparse win: sparse Q75 / BMD Q25 < 1.\n');
end

function r=emptyRow()
r=struct('case_index',NaN,'target_new_nodes',NaN,'expected_result_nodes',NaN, ...
    'overlap_bit',NaN,'overlap_count',NaN, ...
    'p_terms',NaN,'q_terms',NaN,'sparse_pair_products',NaN,'status','', ...
    'bmd_p_nodes',NaN,'bmd_q_nodes',NaN,'bmd_operand_union_nodes',NaN, ...
    'bmd_result_nodes',NaN,'bmd_new_workspace_nodes',NaN,'bmd_mul_cache_entries',NaN, ...
    'sparse_pairs_per_new_bmd_node',NaN,'sparse_result_terms',NaN,'sparse_result_max_coeff',NaN, ...
    'sparse_result_terms_per_bmd_result_node',NaN, ...
    'bmd_median_s',NaN,'bmd_q25_s',NaN,'bmd_q75_s',NaN,'bmd_cv',NaN, ...
    'sparse_median_s',NaN,'sparse_q25_s',NaN,'sparse_q75_s',NaN,'sparse_cv',NaN, ...
    'ratio_median',NaN,'robust_ratio_low',NaN,'robust_ratio_high',NaN,'robust_winner','', ...
    'numeric_check',false);
end

function r=emptyTrialRow()
r=struct('case_index',NaN,'target_new_nodes',NaN,'overlap_bit',NaN,'trial',NaN,'order','', ...
    'bmd_cold_s',NaN,'sparse_s',NaN,'sparse_over_bmd',NaN);
end
