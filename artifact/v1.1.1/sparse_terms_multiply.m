function r = sparse_terms_multiply(a,b)
%SPARSE_TERMS_MULTIPLY Multiply sparse term lists and combine exponents.
ea=a.exponents; eb=b.exponents; ca=a.coefficients; cb=b.coefficients;
na=numel(ea); nb=numel(eb);
if na==0 || nb==0
    r=sparse_terms(uint64([]),[]); return;
end
% Important fast path: the sparse baseline should not pay generic hashing/
% sorting overhead for the monomial case.
if na==1 && nb==1
    r=struct('exponents',ea+eb,'coefficients',ca.*cb); return;
end
% Vectorized outer product, then canonicalize equal exponents.
E = ea(:) + eb(:).';
C = ca(:) .* cb(:).';
e = E(:); c = C(:);
[ue,~,idx] = unique(e,'sorted');
cc = accumarray(double(idx),c,[],@sum);
nz = (cc~=0);
ue=ue(nz); cc=cc(nz);
[ue,ord]=sort(ue,'descend'); cc=cc(ord);
r=struct('exponents',ue(:).','coefficients',cc(:).');
end
