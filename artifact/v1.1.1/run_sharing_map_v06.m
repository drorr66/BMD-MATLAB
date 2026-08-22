function [results,summary,trialsTable] = run_sharing_map_v06(varargin)
%RUN_SHARING_MAP_V06 Map BMD usefulness versus controlled structural sharing.
%
% Fixed problem size:
%   blocks            = 64
%   inner universe    = 256 low exponents per block
%   terms/template    = 128
%   stride            = 2^9 = 512
%   total P terms     = 8192 for EVERY sharing level
%   multiplier Q      = 1+x+...+x^31 (32 terms)
%
% Only the number of distinct block templates changes:
%   1,2,4,8,16,24,32,48,64.
% With one template, all 64 blocks repeat the same low support. With 64
% templates, every block has its own deterministic support pattern. Sparse
% operand cardinality therefore stays fixed while BMD structural sharing is
% progressively removed.

p = inputParser;
addParameter(p,'Trials',7,@(x)isnumeric(x)&&isscalar(x)&&x>=5);
addParameter(p,'SaveResults',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});
numTrials = round(p.Results.Trials);

blocks = 64;
innerUniverse = 256;
termsPerTemplate = 128;
stridePower = 9;
multiplierN = 31;
templateCounts = [1 2 4 8 16 24 32 48 64];

bank = make_sharing_template_bank_v06(blocks,innerUniverse,termsPerTemplate);
rows = repmat(emptyRow(),numel(templateCounts),1);
raw = repmat(emptyTrialRow(),numel(templateCounts)*numTrials,1);
rawPos = 0;

fprintf('\nBMD-MATLAB v0.6 structural-sharing map\n');
fprintf('========================================\n');
fprintf('Fixed P terms=%d, blocks=%d, K=%d of U=%d, stride=2^%d\n', ...
    blocks*termsPerTemplate,blocks,termsPerTemplate,innerUniverse,stridePower);
fprintf('Fixed Q = 1+x+...+x^%d (%d terms), trials=%d\n\n',multiplierN,multiplierN+1,numTrials);

for ii = 1:numel(templateCounts)
    T = templateCounts(ii);
    row = emptyRow();
    row.unique_templates = T;
    row.blocks = blocks;
    row.template_diversity = T/blocks;
    row.normalized_reuse = (blocks-T)/(blocks-1);
    row.mean_blocks_per_template = blocks/T;
    row.inner_universe = innerUniverse;
    row.terms_per_template = termsPerTemplate;
    row.stride_power = stridePower;
    row.multiplier_terms = multiplierN+1;
    fprintf('  templates=%-3d reuse=%6.2fx ... ',T,row.mean_blocks_per_template);

    try
        [exponents,templateIds] = build_sharing_case_v06(bank,T,blocks,stridePower);
        row.logical_terms = numel(exponents);
        row.unique_templates_used = numel(unique(templateIds));
        if row.logical_terms ~= blocks*termsPerTemplate || row.unique_templates_used ~= T
            error('BMD:V06Invariant','Fixed-size/template-count invariant failed.');
        end

        % Build BMD operand P with the direct support-set builder. Build time
        % is reported but deliberately excluded from multiplication timing.
        mgr = BMDManager();
        tic;
        pRef = mgr.indicatorExponents(exponents);
        row.bmd_build_time_s = toc;
        qRef = mgr.geometricSum(multiplierN);
        [baseMgr,refs] = mgr.compact([pRef;qRef]);
        ps = baseMgr.stats(refs(1,:));
        qs = baseMgr.stats(refs(2,:));
        unionStats = baseMgr.stats();
        row.bmd_p_nodes = ps.reachable_nodes;
        row.bmd_q_nodes = qs.reachable_nodes;
        row.bmd_operand_union_nodes = unionStats.total_internal_nodes;
        row.bmd_p_payload_bytes_lb = ps.reachable_payload_bytes_lower_bound;
        row.terms_per_bmd_p_node = row.logical_terms/max(1,row.bmd_p_nodes);

        tic;
        sparseP = sparse_terms(exponents,ones(1,numel(exponents)));
        row.sparse_build_time_s = toc;
        sparseQ = build_sparse_family('geometric_sum',multiplierN,0);
        row.sparse_p_terms = numel(sparseP.exponents);
        row.sparse_q_terms = numel(sparseQ.exponents);
        row.sparse_p_bytes = sparse_terms_bytes(sparseP);
        row.sparse_bytes_over_bmd_payload = row.sparse_p_bytes/max(1,row.bmd_p_payload_bytes_lb);

        % JIT/allocation warm-up on throwaway state. The measured BMD trials
        % still start from a fresh compact manager with empty computed caches.
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
            rr.unique_templates = T;
            rr.trial = kk;
            if mod(kk,2)==1, rr.order='BMD_FIRST'; else, rr.order='SPARSE_FIRST'; end
            rr.bmd_cold_s = bmdTimes(kk);
            rr.sparse_s = sparseTimes(kk);
            rr.sparse_over_bmd = sparseTimes(kk)/max(realmin,bmdTimes(kk));
            raw(rawPos) = rr;
        end

        % Exact structural diagnostics and floating evaluation check.
        after = firstMgr.stats(firstProd);
        row.bmd_result_nodes = after.reachable_nodes;
        row.bmd_new_workspace_nodes = after.total_internal_nodes - row.bmd_operand_union_nodes;
        row.bmd_mul_cache_entries = after.mul_cache_entries;
        row.sparse_result_terms = numel(firstSparseProd.exponents);
        checks = true;
        for xv = [0.997 0.9993]
            yb = firstMgr.evaluate(firstProd,xv);
            ys = sparse_terms_evaluate(firstSparseProd,xv);
            checks = checks && near(yb,ys);
        end
        row.numeric_check = checks;
        if ~checks
            error('BMD:V06Numeric','BMD and sparse product evaluations differ.');
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
        fprintf(' Pnodes=%-5d compression=%7.2fx ratio=%6.3fx robust=[%.3f, %.3f] %s', ...
            row.bmd_p_nodes,row.sparse_bytes_over_bmd_payload,row.ratio_median, ...
            row.robust_ratio_low,row.robust_ratio_high,row.robust_winner);
    end
    fprintf('\n');
end

results = struct2table(rows);
trialsTable = struct2table(raw(1:rawPos));
summary = makeSummary(results);

disp(results(:,{'unique_templates','mean_blocks_per_template','bmd_p_nodes','terms_per_bmd_p_node', ...
    'sparse_bytes_over_bmd_payload','bmd_median_s','sparse_median_s','ratio_median','robust_winner','status'}));
fprintf('\nSharing-map summary:\n');
disp(summary);

if p.Results.SaveResults
    outDir = fullfile(fileparts(mfilename('fullpath')),'results');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    writetable(results,fullfile(outDir,'sharing_map_results_v06.csv'));
    writetable(trialsTable,fullfile(outDir,'sharing_map_trials_v06.csv'));
    writetable(summary,fullfile(outDir,'sharing_map_summary_v06.csv'));
    writeMetadata(fullfile(outDir,'run_metadata_v06.txt'),numTrials,blocks,innerUniverse,termsPerTemplate,stridePower,multiplierN,templateCounts);
end
end

function [t,m,out] = measureBmdCold(baseMgr,refs)
[m,r] = baseMgr.compact(refs);
tic;
out = m.multiply(r(1,:),r(2,:));
t = toc;
if out(2) < 1, error('BMD:V06Sink','Invalid BMD result.'); end
end

function [t,out] = measureSparse(a,b)
tic;
out = sparse_terms_multiply(a,b);
t = toc;
if isempty(out.exponents), error('BMD:V06Sink','Unexpected empty sparse result.'); end
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
    'median_transition_status','', 'median_templates_low',NaN,'median_templates_high',NaN, ...
    'median_ratio_low',NaN,'median_ratio_high',NaN,'median_templates_estimate',NaN,'median_bmd_nodes_estimate',NaN, ...
    'last_robust_bmd_templates',NaN,'first_robust_sparse_templates',NaN, ...
    'sustained_sparse_median_from',NaN,'sustained_sparse_robust_from',NaN, ...
    'max_compression_sparse_over_bmd',NaN,'min_compression_sparse_over_bmd',NaN, ...
    'max_bmd_speedup',NaN,'max_bmd_speedup_templates',NaN);
ok=strcmp(results.status,'OK');
T=results.unique_templates(ok); r=results.ratio_median(ok); nodes=results.bmd_p_nodes(ok);
w=results.robust_winner(ok); comp=results.sparse_bytes_over_bmd_payload(ok);
[s.median_transition_status,s.median_templates_low,s.median_templates_high,s.median_ratio_low,s.median_ratio_high,s.median_templates_estimate,s.median_bmd_nodes_estimate] = fallingBracket(T,r,nodes);
bmdIdx=find(strcmp(w,'BMD')); sparseIdx=find(strcmp(w,'SPARSE'));
if ~isempty(bmdIdx), s.last_robust_bmd_templates=T(bmdIdx(end)); end
if ~isempty(sparseIdx), s.first_robust_sparse_templates=T(sparseIdx(1)); end
s.sustained_sparse_median_from=sustainedFrom(T,r<1);
s.sustained_sparse_robust_from=sustainedFrom(T,strcmp(w,'SPARSE'));
if ~isempty(comp), s.max_compression_sparse_over_bmd=max(comp); s.min_compression_sparse_over_bmd=min(comp); end
if ~isempty(r)
    [s.max_bmd_speedup,ix]=max(r); s.max_bmd_speedup_templates=T(ix);
end
s=struct2table(s);
end

function [status,tlo,thi,rlo,rhi,estT,estNodes] = fallingBracket(T,r,nodes)
status='NO_TRANSITION'; tlo=NaN; thi=NaN; rlo=NaN; rhi=NaN; estT=NaN; estNodes=NaN;
for k=2:numel(T)
    if r(k-1)>=1 && r(k)<1
        status='BRACKETED'; tlo=T(k-1); thi=T(k); rlo=r(k-1); rhi=r(k);
        if rlo>0 && rhi>0 && rhi~=rlo
            a=(0-log(rlo))/(log(rhi)-log(rlo));
            estT=exp(log(tlo)+a*(log(thi)-log(tlo)));
            estNodes=exp(log(max(1,nodes(k-1)))+a*(log(max(1,nodes(k)))-log(max(1,nodes(k-1)))));
        end
        return;
    end
end
if all(r>=1), status='BMD_FASTER_ALL'; elseif all(r<1), status='SPARSE_FASTER_ALL'; end
end

function t0 = sustainedFrom(T,win)
t0=NaN;
for k=1:numel(T)
    if all(win(k:end)), t0=T(k); return; end
end
end

function writeMetadata(path,numTrials,blocks,U,K,stridePower,multiplierN,templateCounts)
fid=fopen(path,'w'); if fid<0, return; end
c=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'BMD-MATLAB v0.6 structural-sharing map\n');
fprintf(fid,'Generated: %s\n',datestr(now,31));
fprintf(fid,'MATLAB version: %s\n',version);
fprintf(fid,'Computer: %s\n',computer);
fprintf(fid,'Trials per case: %d\n',numTrials);
fprintf(fid,'blocks: %d\n',blocks);
fprintf(fid,'inner_universe: %d\n',U);
fprintf(fid,'terms_per_template: %d\n',K);
fprintf(fid,'stride_power: %d\n',stridePower);
fprintf(fid,'multiplier_n: %d\n',multiplierN);
fprintf(fid,'unique_templates: [%s]\n',strtrim(sprintf('%d ',templateCounts)));
fprintf(fid,'All P cases have exactly blocks*terms_per_template terms.\n');
fprintf(fid,'Winner ratio = sparse_time/BMD_time; >1 means BMD faster.\n');
fprintf(fid,'Robust BMD win: sparse Q25 / BMD Q75 > 1.\n');
fprintf(fid,'Robust sparse win: sparse Q75 / BMD Q25 < 1.\n');
end

function tf = near(a,b)
tol=1e-9*max([1 abs(a) abs(b)]);
tf=abs(a-b)<=tol;
end

function r=emptyRow()
r=struct('unique_templates',NaN,'blocks',NaN,'template_diversity',NaN,'normalized_reuse',NaN,'mean_blocks_per_template',NaN, ...
    'inner_universe',NaN,'terms_per_template',NaN,'stride_power',NaN,'multiplier_terms',NaN,'logical_terms',NaN,'unique_templates_used',NaN, ...
    'status','','bmd_build_time_s',NaN,'bmd_p_nodes',NaN,'bmd_q_nodes',NaN,'bmd_operand_union_nodes',NaN,'bmd_p_payload_bytes_lb',NaN, ...
    'terms_per_bmd_p_node',NaN,'sparse_build_time_s',NaN,'sparse_p_terms',NaN,'sparse_q_terms',NaN,'sparse_p_bytes',NaN,'sparse_bytes_over_bmd_payload',NaN, ...
    'bmd_median_s',NaN,'bmd_q25_s',NaN,'bmd_q75_s',NaN,'bmd_cv',NaN,'sparse_median_s',NaN,'sparse_q25_s',NaN,'sparse_q75_s',NaN,'sparse_cv',NaN, ...
    'ratio_median',NaN,'robust_ratio_low',NaN,'robust_ratio_high',NaN,'robust_winner','', ...
    'bmd_result_nodes',NaN,'bmd_new_workspace_nodes',NaN,'bmd_mul_cache_entries',NaN,'sparse_result_terms',NaN,'numeric_check',false);
end

function r=emptyTrialRow()
r=struct('unique_templates',NaN,'trial',NaN,'order','','bmd_cold_s',NaN,'sparse_s',NaN,'sparse_over_bmd',NaN);
end
