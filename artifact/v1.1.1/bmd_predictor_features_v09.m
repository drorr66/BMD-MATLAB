function f = bmd_predictor_features_v09(mgr,pRef,qRef,pTerms,qTerms)
%BMD_PREDICTOR_FEATURES_V09 Pre-multiply structural features; never multiplies.
%
% The feature extractor traverses only the two existing canonical DAGs. It
% does NOT call BMDManager.multiply and does not create predicted result
% nodes. pTerms/qTerms are explicit term-count metadata supplied by the
% caller; v0.9 does not yet attempt to infer them from an arbitrary BMD.

if nargin < 5
    error('BMD:V09Terms','pTerms and qTerms are required in v0.9.');
end
if pTerms < 1 || qTerms < 1 || pTerms ~= floor(pTerms) || qTerms ~= floor(qTerms)
    error('BMD:V09Terms','pTerms and qTerms must be positive integer counts.');
end

[pIds,pHist,pLevels,pChain] = inspectRoot(mgr,pRef);
[qIds,qHist,qLevels,qChain] = inspectRoot(mgr,qRef);
allIds = unique([pIds qIds]);
commonLevels = intersect(pLevels,qLevels);
unionLevels = union(pLevels,qLevels);

f = struct();
f.p_nodes = numel(pIds);
f.q_nodes = numel(qIds);
f.union_nodes = numel(allIds);
f.p_levels = numel(pLevels);
f.q_levels = numel(qLevels);
f.union_levels = numel(unionLevels);
f.common_levels = numel(commonLevels);
f.p_is_bitcube_chain = pChain;
f.q_is_bitcube_chain = qChain;
f.both_bitcube_chain = pChain && qChain;
f.ordered_disjoint_bands = false;
if ~isempty(pLevels) && ~isempty(qLevels) && isempty(commonLevels)
    f.ordered_disjoint_bands = max(pLevels) < min(qLevels) || max(qLevels) < min(pLevels);
end

f.carry_runsum = 0;
f.carry_run2 = 0;
f.carry_max = 0;
f.carry_massrun = 0;
f.adjacent_mass = 0;
f.invdist = 0;
f.suffix_product_sum = 0;
f.suffix_carry_sum = 0;

for ii = 1:numel(commonLevels)
    lev = commonLevels(ii);
    run = consecutiveRunAbove(lev,unionLevels);
    f.carry_runsum = f.carry_runsum + run;
    f.carry_run2 = f.carry_run2 + (run+1)^2;
    f.carry_max = max(f.carry_max,run);
    ca = histCount(pHist,lev); cb = histCount(qHist,lev);
    f.carry_massrun = f.carry_massrun + ca*cb*(run+1);
    sa = suffixMass(pHist,lev); sb = suffixMass(qHist,lev);
    f.suffix_product_sum = f.suffix_product_sum + sa*sb;
    f.suffix_carry_sum = f.suffix_carry_sum + sa*sb*(run+1);
end

for ia = 1:size(pHist,1)
    la = pHist(ia,1); ca = pHist(ia,2);
    for ib = 1:size(qHist,1)
        lb = qHist(ib,1); cb = qHist(ib,2);
        d = abs(la-lb); w = ca*cb;
        f.invdist = f.invdist + w/(1+d);
        if d == 1
            f.adjacent_mass = f.adjacent_mass + w;
        end
    end
end

if isempty(unionLevels)
    f.level_span = 0;
else
    f.level_span = max(unionLevels)-min(unionLevels)+1;
end
f.p_terms = double(pTerms);
f.q_terms = double(qTerms);
f.pair_products = double(pTerms)*double(qTerms);
end

function [ids,hist,levels,isChain] = inspectRoot(mgr,ref)
ids = reachableIds(mgr,ref);
levels = zeros(1,numel(ids));
for k=1:numel(ids)
    levels(k)=double(mgr.levels(ids(k)));
end
if isempty(levels)
    hist=zeros(0,2);
else
    u=unique(levels);
    hist=zeros(numel(u),2);
    for k=1:numel(u)
        hist(k,:)=[u(k),sum(levels==u(k))];
    end
end
levels=unique(levels);
isChain=detectBitcubeChain(mgr,ref,ids);
end

function ids = reachableIds(mgr,ref)
if ref(1)==0 || ref(2)==1
    ids=zeros(1,0);
    return;
end
seen=false(1,numel(mgr.levels));
stack=double(ref(2));
ids=zeros(1,0);
while ~isempty(stack)
    n=stack(end); stack(end)=[];
    if n==1 || seen(n), continue; end
    seen(n)=true; ids(end+1)=n; %#ok<AGROW>
    lo=double(mgr.lowNode(n)); hi=double(mgr.highNode(n));
    if lo~=1, stack(end+1)=lo; end %#ok<AGROW>
    if hi~=1, stack(end+1)=hi; end %#ok<AGROW>
end
ids=sort(ids);
end

function tf = detectBitcubeChain(mgr,ref,ids)
% A bit-cube indicator is a single canonical chain whose low/high branches
% point to the same child with the same nonzero normalized edge weight.
if ref(1)==0
    tf=false; return;
end
if ref(2)==1
    tf=true; return;
end
if isempty(ids)
    tf=false; return;
end
seen=false(1,numel(mgr.levels));
n=double(ref(2)); count=0; tf=true;
while n~=1
    if seen(n), tf=false; return; end
    seen(n)=true; count=count+1;
    lo=double(mgr.lowNode(n)); hi=double(mgr.highNode(n));
    lw=mgr.lowWeight(n); hw=mgr.highWeight(n);
    if lo~=hi || lw~=hw || lw==0
        tf=false; return;
    end
    n=lo;
end
tf = tf && count==numel(ids);
end

function c = histCount(hist,lev)
idx=find(hist(:,1)==lev,1);
if isempty(idx), c=0; else, c=hist(idx,2); end
end

function s = suffixMass(hist,lev)
if isempty(hist), s=0; else, s=sum(hist(hist(:,1)>=lev,2)); end
end

function r = consecutiveRunAbove(lev,levels)
r=0; k=lev+1;
while any(levels==k)
    r=r+1; k=k+1;
end
end
