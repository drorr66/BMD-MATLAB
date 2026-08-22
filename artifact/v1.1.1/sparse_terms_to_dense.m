function p = sparse_terms_to_dense(s,maxDegree)
%SPARSE_TERMS_TO_DENSE Validation helper; descending dense coefficients.
if nargin<2, maxDegree=1e6; end
if isempty(s.exponents), p=0; return; end
deg=double(max(s.exponents));
if deg>maxDegree, error('BMD:SparseDenseTooLarge','Degree exceeds maxDegree.'); end
p=zeros(1,deg+1);
for k=1:numel(s.exponents)
    idx=deg-double(s.exponents(k))+1;
    p(idx)=p(idx)+s.coefficients(k);
end
idx=find(p~=0,1,'first'); if isempty(idx), p=0; else, p=p(idx:end); end
end
