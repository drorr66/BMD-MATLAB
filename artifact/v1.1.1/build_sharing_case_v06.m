function [exponents,templateIds] = build_sharing_case_v06(bank,uniqueTemplates,blocks,stridePower)
%BUILD_SHARING_CASE_V06 Fixed-size polynomial with controlled template reuse.
%
% The polynomial is a union of `blocks` disjoint low-exponent blocks.  Each
% block contains the same number of terms; only the number of distinct block
% templates changes.  uniqueTemplates=1 maximizes block-level reuse, while
% uniqueTemplates=blocks removes deliberate block-template reuse.
if uniqueTemplates < 1 || uniqueTemplates > blocks || uniqueTemplates ~= floor(uniqueTemplates)
    error('BMD:V06Case','uniqueTemplates must be an integer in 1..blocks.');
end
if size(bank,1) < blocks
    error('BMD:V06Case','Template bank has fewer rows than blocks.');
end
stride = uint64(2^double(stridePower));
innerUniverse = size(bank,2);
if innerUniverse > double(stride)
    error('BMD:V06Case','innerUniverse must be <= 2^stridePower.');
end
termsPerTemplate = sum(bank(1,:));
if any(sum(bank(1:blocks,:),2) ~= termsPerTemplate)
    error('BMD:V06Case','All templates must have equal cardinality.');
end

templateIds = 1 + mod(0:blocks-1,uniqueTemplates);
exponents = zeros(1,blocks*termsPerTemplate,'uint64');
pos = 0;
for jj = 1:blocks
    low = uint64(find(bank(templateIds(jj),:))-1);
    idx = pos + (1:numel(low));
    exponents(idx) = uint64(jj-1)*stride + low;
    pos = pos + numel(low);
end
exponents = exponents(1:pos);
if numel(unique(exponents)) ~= numel(exponents)
    error('BMD:V06Case','Disjoint-block support invariant failed.');
end
end
