function s = sparse_terms(exponents,coefficients)
%SPARSE_TERMS Construct normalized sparse polynomial term list.
e = uint64(exponents(:).');
c = double(coefficients(:).');
if numel(e) ~= numel(c)
    error('BMD:SparseShape','Exponent and coefficient arrays must have equal length.');
end
if isempty(e)
    s = struct('exponents',uint64([]),'coefficients',double([]));
    return;
end
nz = (c~=0);
e=e(nz); c=c(nz);
if isempty(e)
    s = struct('exponents',uint64([]),'coefficients',double([]));
    return;
end
[eSorted,ord] = sort(e,'descend');
cSorted = c(ord);
% Combine accidental duplicate exponents while preserving descending order.
outE = zeros(1,numel(eSorted),'uint64');
outC = zeros(1,numel(eSorted));
w=0; i=1;
while i<=numel(eSorted)
    ee=eSorted(i); cc=0;
    while i<=numel(eSorted) && eSorted(i)==ee
        cc=cc+cSorted(i); i=i+1;
    end
    if cc~=0
        w=w+1; outE(w)=ee; outC(w)=cc;
    end
end
s=struct('exponents',outE(1:w),'coefficients',outC(1:w));
end
