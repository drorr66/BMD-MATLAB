function [results, summary, trialsTable] = run_operation_closure_v07(varargin)
%RUN_OPERATION_CLOSURE_V07 Controlled operation-closure map: BMD vs sparse.
%
% Both operands have exactly 256 sparse terms and exactly 8 BMD nodes.
%
%   P   = sum_{i=0}^{255} x^i
%   Q_s = sum_{j=0}^{255} x^(j * 2^s)
%
% Only s changes. P occupies exponent bits 0..7 and Q_s occupies bits
% s..s+7.  Thus s controls how much the operand bit bands overlap, while
% explicit input cardinalities and individual BMD operand sizes stay fixed.
%
% When s>=8 the bit bands are disjoint and the product is a clean shared
% grid.  For smaller s, carries/collisions make the canonical product and
% the intermediate BMD workspace progressively larger.  The experiment is
% designed to test whether cold-multiply time tracks operation closure
% (result/new DAG nodes) more strongly than input representation size.

p = inputParser;
addParameter(p,'Trials',9,@(x)isnumeric(x)&&isscalar(x)&&x>=5);
addParameter(p,'SaveResults',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});
numTrials = round(p.Results.Trials);

nTerms = 256;
maxExponent = nTerms-1;
shiftList = [0 1 2 3 4 5 6 7 8 9 10 12];
bitWidth = 8; % 256 terms == 2^8
pairProducts = nTerms*nTerms;

rows = repmat(emptyRow(),numel(shiftList),1);
raw = repmat(emptyTrialRow(),numel(shiftList)*numTrials,1);
rawPos = 0;

fprintf('\nBMD-MATLAB v0.7 operation-closure map\n');
fprintf('======================================\n');
fprintf('P terms=%d, Q terms=%d, pair-products=%d, trials=%d\n\n', ...
    nTerms,nTerms,pairProducts,numTrials);

for ii = 1:numel(shiftList)
    shiftPower = shiftList(ii);
    row = emptyRow();
    row.shift_power = shiftPower;
    row.q_step = 2^shiftPower;
    row.bit_band_overlap = max(0,bitWidth-shiftPower);
    row.normalized_band_overlap = row.bit_band_overlap/bitWidth;
    row.normalized_band_separation = 1-row.normalized_band_overlap;
    row.p_terms = nTerms;
    row.q_terms = nTerms;
    row.sparse_pair_products = pairProducts;

    fprintf('  shift=%-2d overlap_bits=%d ... ',shiftPower,row.bit_band_overlap);

    try
        % Build operands outside timed region, then compact the union so each
        % cold BMD trial begins from only the operand DAG and empty caches.
        baseMgr = BMDManager();
        pRef = baseMgr.geometricSum(maxExponent);
        qRef = baseMgr.geometricSumShifted(maxExponent,shiftPower+1);
        pStats = baseMgr.stats(pRef);
        qStats = baseMgr.stats(qRef);
        [baseMgr,refs] = baseMgr.compact([pRef;qRef]);
        unionStats = baseMgr.stats();
        row.bmd_p_nodes = pStats.reachable_nodes;
        row.bmd_q_nodes = qStats.reachable_nodes;
        row.bmd_operand_union_nodes = unionStats.total_internal_nodes;
        if row.bmd_p_nodes ~= bitWidth || row.bmd_q_nodes ~= bitWidth
            error('BMD:V07OperandSize','Expected both BMD operands to have %d nodes.',bitWidth);
        end

        sparseP = build_sparse_family('geometric_sum',maxExponent,0);
        sparseQ = build_sparse_shifted_geometric(nTerms,shiftPower);
        if numel(sparseP.exponents) ~= nTerms || numel(sparseQ.exponents) ~= nTerms
            error('BMD:V07SparseCardinality','Sparse input cardinality changed unexpectedly.');
        end

        % Untimed warm-up on throwaway state.
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
            rr.shift_power = shiftPower;
            rr.bit_band_overlap = row.bit_band_overlap;
            rr.trial = kk;
            if mod(kk,2)==1, rr.order='BMD_FIRST'; else, rr.order='SPARSE_FIRST'; end
            rr.bmd_cold_s = bmdTimes(kk);
            rr.sparse_s = sparseTimes(kk);
            rr.sparse_over_bmd = sparseTimes(kk)/max(realmin,bmdTimes(kk));
            raw(rawPos) = rr;
        end

        % Structural closure diagnostics from the first genuinely cold run.
        after = firstMgr.stats(firstProd);
        row.bmd_result_nodes = after.reachable_nodes;
        row.bmd_new_workspace_nodes = after.total_internal_nodes - row.bmd_operand_union_nodes;
        row.bmd_mul_cache_entries = after.mul_cache_entries;
        row.bmd_result_over_operand_union = row.bmd_result_nodes/max(1,row.bmd_operand_union_nodes);
        row.bmd_new_over_operand_union = row.bmd_new_workspace_nodes/max(1,row.bmd_operand_union_nodes);
        row.sparse_pairs_per_new_bmd_node = pairProducts/max(1,row.bmd_new_workspace_nodes);

        row.sparse_result_terms = numel(firstSparseProd.exponents);
        row.sparse_result_max_coeff = max(abs(firstSparseProd.coefficients));
        row.sparse_result_terms_per_bmd_result_node = row.sparse_result_terms/max(1,row.bmd_result_nodes);

        checks = true;
        for xv = [0.9995 0.9999]
            yb = firstMgr.evaluate(firstProd,xv);
            ys = sparse_terms_evaluate(firstSparseProd,xv);
            checks = checks && near(yb,ys);
        end
        row.numeric_check = checks;
        if ~checks
            error('BMD:V07Numeric','BMD and sparse product evaluations differ.');
        end

        % For fully disjoint bands the product is the direct geometric grid.
        if shiftPower >= bitWidth
            direct = firstMgr.geometricGrid(maxExponent,nTerms,shiftPower);
            if ~firstMgr.same(firstProd,direct)
                error('BMD:V07Closure','Disjoint-band product != direct canonical grid.');
            end
            if row.sparse_result_terms ~= pairProducts || any(firstSparseProd.coefficients ~= 1)
                error('BMD:V07SparseStructure','Disjoint-band sparse product is not collision-free.');
            end
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
        fprintf(' new=%-4d result=%-4d out=%-6d ratio=%6.3fx robust=%s', ...
            row.bmd_new_workspace_nodes,row.bmd_result_nodes,row.sparse_result_terms, ...
            row.ratio_median,row.robust_winner);
    end
    fprintf('\n');
end

results = struct2table(rows);
trialsTable = struct2table(raw(1:rawPos));
summary = makeSummary(results);

disp(results(:,{'shift_power','bit_band_overlap','bmd_result_nodes','bmd_new_workspace_nodes', ...
    'sparse_result_terms','bmd_median_s','sparse_median_s','ratio_median','robust_winner','status'}));
fprintf('\nOperation-closure summary:\n');
disp(summary);

if p.Results.SaveResults
    outDir = fullfile(fileparts(mfilename('fullpath')),'results');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    writetable(results,fullfile(outDir,'operation_closure_results_v07.csv'));
    writetable(trialsTable,fullfile(outDir,'operation_closure_trials_v07.csv'));
    writetable(summary,fullfile(outDir,'operation_closure_summary_v07.csv'));
    writeMetadata(fullfile(outDir,'run_metadata_v07.txt'),numTrials,nTerms,shiftList,bitWidth,pairProducts);
end
end

function [t,m,out] = measureBmdCold(baseMgr,refs)
[m,r] = baseMgr.compact(refs);
tic;
out = m.multiply(r(1,:),r(2,:));
t = toc;
if out(2) < 1, error('BMD:V07Sink','Invalid BMD result.'); end
end

function [t,out] = measureSparse(a,b)
tic;
out = sparse_terms_multiply(a,b);
t = toc;
if isempty(out.exponents), error('BMD:V07Sink','Unexpected empty sparse result.'); end
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
    'median_transition_status','', 'median_shift_low',NaN,'median_shift_high',NaN, ...
    'median_ratio_low',NaN,'median_ratio_high',NaN,'median_shift_estimate',NaN, ...
    'robust_transition_status','', 'robust_shift_low',NaN,'robust_shift_high',NaN, ...
    'robust_ratio_low_side',NaN,'robust_ratio_high_side',NaN,'robust_shift_estimate',NaN, ...
    'last_sparse_new_nodes',NaN,'first_bmd_new_nodes',NaN, ...
    'last_sparse_result_nodes',NaN,'first_bmd_result_nodes',NaN, ...
    'last_sparse_pairs_per_new_node',NaN,'first_bmd_pairs_per_new_node',NaN, ...
    'sustained_median_bmd_from_shift',NaN,'sustained_robust_bmd_from_shift',NaN, ...
    'corr_log_bmd_time_vs_log_new_nodes',NaN,'corr_log_bmd_time_vs_log_result_nodes',NaN, ...
    'corr_log_sparse_time_vs_log_result_terms',NaN,'corr_log_ratio_vs_log_pairs_per_new_node',NaN, ...
    'bmd_time_exponent_vs_new_nodes',NaN,'sparse_time_exponent_vs_result_terms',NaN, ...
    'max_bmd_speedup',NaN,'max_bmd_speedup_shift',NaN);

ok=strcmp(results.status,'OK');
sh=results.shift_power(ok); ratio=results.ratio_median(ok); robust=results.robust_ratio_low(ok);
newn=results.bmd_new_workspace_nodes(ok); resn=results.bmd_result_nodes(ok);
pairratio=results.sparse_pairs_per_new_bmd_node(ok);
bmdt=results.bmd_median_s(ok); spt=results.sparse_median_s(ok); outterms=results.sparse_result_terms(ok);

[s.median_transition_status,s.median_shift_low,s.median_shift_high,s.median_ratio_low,s.median_ratio_high,s.median_shift_estimate] = risingBracket(sh,ratio);
[s.robust_transition_status,s.robust_shift_low,s.robust_shift_high,s.robust_ratio_low_side,s.robust_ratio_high_side,s.robust_shift_estimate] = risingBracket(sh,robust);
s.sustained_median_bmd_from_shift=sustainedFrom(sh,ratio>1);
s.sustained_robust_bmd_from_shift=sustainedFrom(sh,robust>1);

idx=find(robust>1,1,'first');
if ~isempty(idx)
    s.first_bmd_new_nodes=newn(idx); s.first_bmd_result_nodes=resn(idx); s.first_bmd_pairs_per_new_node=pairratio(idx);
    if idx>1
        s.last_sparse_new_nodes=newn(idx-1); s.last_sparse_result_nodes=resn(idx-1); s.last_sparse_pairs_per_new_node=pairratio(idx-1);
    end
end

s.corr_log_bmd_time_vs_log_new_nodes=safeLogCorr(newn,bmdt);
s.corr_log_bmd_time_vs_log_result_nodes=safeLogCorr(resn,bmdt);
s.corr_log_sparse_time_vs_log_result_terms=safeLogCorr(outterms,spt);
s.corr_log_ratio_vs_log_pairs_per_new_node=safeLogCorr(pairratio,ratio);
s.bmd_time_exponent_vs_new_nodes=logSlope(newn,bmdt);
s.sparse_time_exponent_vs_result_terms=logSlope(outterms,spt);
if ~isempty(ratio)
    [s.max_bmd_speedup,ix]=max(ratio); s.max_bmd_speedup_shift=sh(ix);
end
s=struct2table(s);
end

function [status,xlo,xhi,rlo,rhi,est] = risingBracket(x,r)
status='NO_CROSSING'; xlo=NaN; xhi=NaN; rlo=NaN; rhi=NaN; est=NaN;
for k=2:numel(x)
    if r(k-1)<1 && r(k)>=1
        status='BRACKETED'; xlo=x(k-1); xhi=x(k); rlo=r(k-1); rhi=r(k);
        if rlo>0 && rhi>0 && rhi~=rlo
            a=(0-log(rlo))/(log(rhi)-log(rlo));
            est=xlo+a*(xhi-xlo);
        end
        return;
    end
end
if all(r>=1), status='BMD_FASTER_ALL'; elseif all(r<1), status='SPARSE_FASTER_ALL'; end
end

function x0=sustainedFrom(x,win)
x0=NaN;
for k=1:numel(x)
    if all(win(k:end)), x0=x(k); return; end
end
end

function c=safeLogCorr(x,y)
mask=isfinite(x)&isfinite(y)&x>0&y>0;
x=x(mask); y=y(mask);
if numel(x)<3 || numel(unique(x))<2 || numel(unique(y))<2, c=NaN; return; end
C=corrcoef(log(double(x)),log(double(y))); c=C(1,2);
end

function a=logSlope(x,y)
mask=isfinite(x)&isfinite(y)&x>0&y>0;
x=x(mask); y=y(mask);
if numel(x)<3 || numel(unique(x))<2, a=NaN; return; end
p=polyfit(log(double(x)),log(double(y)),1); a=p(1);
end

function writeMetadata(path,numTrials,nTerms,shiftList,bitWidth,pairProducts)
fid=fopen(path,'w'); if fid<0, return; end
c=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'BMD-MATLAB v0.7 operation-closure map\n');
fprintf(fid,'Generated: %s\n',datestr(now,31));
fprintf(fid,'MATLAB version: %s\n',version);
fprintf(fid,'Computer: %s\n',computer);
fprintf(fid,'Trials per case: %d\n',numTrials);
fprintf(fid,'P terms: %d\n',nTerms);
fprintf(fid,'Q terms: %d\n',nTerms);
fprintf(fid,'P BMD nodes expected: %d\n',bitWidth);
fprintf(fid,'Q BMD nodes expected: %d\n',bitWidth);
fprintf(fid,'pair products fixed: %d\n',pairProducts);
fprintf(fid,'shift_power: [%s]\n',strtrim(sprintf('%d ',shiftList)));
fprintf(fid,'P bit band: [0,%d]\n',bitWidth-1);
fprintf(fid,'Q_s bit band: [s,s+%d]\n',bitWidth-1);
fprintf(fid,'Winner ratio = sparse_time/BMD_time; >1 means BMD faster.\n');
fprintf(fid,'Robust BMD win: sparse Q25 / BMD Q75 > 1.\n');
fprintf(fid,'Robust sparse win: sparse Q75 / BMD Q25 < 1.\n');
end

function tf=near(a,b)
tol=1e-8*max([1 abs(a) abs(b)]); tf=abs(a-b)<=tol;
end

function r=emptyRow()
r=struct('shift_power',NaN,'q_step',NaN,'bit_band_overlap',NaN,'normalized_band_overlap',NaN,'normalized_band_separation',NaN, ...
    'p_terms',NaN,'q_terms',NaN,'sparse_pair_products',NaN,'status','', ...
    'bmd_p_nodes',NaN,'bmd_q_nodes',NaN,'bmd_operand_union_nodes',NaN, ...
    'bmd_result_nodes',NaN,'bmd_new_workspace_nodes',NaN,'bmd_mul_cache_entries',NaN, ...
    'bmd_result_over_operand_union',NaN,'bmd_new_over_operand_union',NaN,'sparse_pairs_per_new_bmd_node',NaN, ...
    'sparse_result_terms',NaN,'sparse_result_max_coeff',NaN,'sparse_result_terms_per_bmd_result_node',NaN, ...
    'bmd_median_s',NaN,'bmd_q25_s',NaN,'bmd_q75_s',NaN,'bmd_cv',NaN, ...
    'sparse_median_s',NaN,'sparse_q25_s',NaN,'sparse_q75_s',NaN,'sparse_cv',NaN, ...
    'ratio_median',NaN,'robust_ratio_low',NaN,'robust_ratio_high',NaN,'robust_winner','', ...
    'numeric_check',false);
end

function r=emptyTrialRow()
r=struct('shift_power',NaN,'bit_band_overlap',NaN,'trial',NaN,'order','', ...
    'bmd_cold_s',NaN,'sparse_s',NaN,'sparse_over_bmd',NaN);
end
