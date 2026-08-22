function packed = pack_bmd_roots_v11(mgr,roots)
%PACK_BMD_ROOTS_V11 Portable packed representation for memory accounting.
%
% This is NOT the live BMDManager footprint. It stores exactly the numeric
% DAG payload required by the selected roots, with child IDs remapped to a
% compact local numbering. workspace_bytes_v11(packed) therefore measures
% an actual MATLAB packed representation, while deliberately excluding
% containers.Map unique/computed-table overhead.
roots=int64(roots);
if isvector(roots) && numel(roots)==2, roots=reshape(roots,1,2); end
if size(roots,2)~=2, error('BMD:V11Roots','roots must be Nx2 refs.'); end
seen=false(1,numel(mgr.levels)); stack=double(roots(roots(:,2)~=1,2)).';
while ~isempty(stack)
    n=stack(end); stack(end)=[];
    if n==1 || seen(n), continue; end
    seen(n)=true;
    lo=double(mgr.lowNode(n)); hi=double(mgr.highNode(n));
    if lo~=1, stack(end+1)=lo; end %#ok<AGROW>
    if hi~=1, stack(end+1)=hi; end %#ok<AGROW>
end
ids=find(seen);
map=zeros(1,numel(mgr.levels),'uint32');
map(1)=uint32(0); % packed terminal id
for k=1:numel(ids), map(ids(k))=uint32(k); end
n=numel(ids);
levels=zeros(1,n,'uint32'); lowNode=zeros(1,n,'uint32'); highNode=zeros(1,n,'uint32');
lowWeight=zeros(1,n,'int64'); highWeight=zeros(1,n,'int64');
for k=1:n
    old=ids(k); levels(k)=mgr.levels(old);
    lo=double(mgr.lowNode(old)); hi=double(mgr.highNode(old));
    lowNode(k)=map(lo); highNode(k)=map(hi);
    lowWeight(k)=mgr.lowWeight(old); highWeight(k)=mgr.highWeight(old);
end
proots=zeros(size(roots),'int64');
proots(:,1)=roots(:,1);
for k=1:size(roots,1)
    if roots(k,2)==1, proots(k,2)=0; else, proots(k,2)=int64(map(double(roots(k,2)))); end
end
packed=struct('levels',levels,'lowNode',lowNode,'highNode',highNode, ...
    'lowWeight',lowWeight,'highWeight',highWeight,'roots',proots);
end
