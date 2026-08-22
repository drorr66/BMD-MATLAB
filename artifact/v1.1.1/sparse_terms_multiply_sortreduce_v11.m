function r = sparse_terms_multiply_sortreduce_v11(a,b)
%SPARSE_TERMS_MULTIPLY_SORTREDUCE_V11 Alternate vectorized sparse baseline.
%
% Forms all explicit term pairs, sorts pair exponents once, and performs a
% linear run reduction. This is intentionally independent of the v0.x
% unique+accumarray sparse baseline so publication validation can compare
% BMD against the faster of two explicit sparse implementations.
ea=a.exponents; eb=b.exponents; ca=a.coefficients; cb=b.coefficients;
na=numel(ea); nb=numel(eb);
if na==0 || nb==0
    r=sparse_terms(uint64([]),[]); return;
end
if na==1 && nb==1
    r=struct('exponents',ea+eb,'coefficients',ca.*cb); return;
end
E=ea(:)+eb(:).';
C=ca(:).*cb(:).';
e=E(:); c=C(:);
[e,ord]=sort(e,'ascend'); c=c(ord);
newRun=[true; e(2:end)~=e(1:end-1)];
starts=find(newRun);
stops=[starts(2:end)-1; numel(e)];
outE=zeros(numel(starts),1,'uint64'); outC=zeros(numel(starts),1);
for k=1:numel(starts)
    outE(k)=e(starts(k));
    outC(k)=sum(c(starts(k):stops(k)));
end
nz=(outC~=0); outE=outE(nz); outC=outC(nz);
[outE,ord]=sort(outE,'descend'); outC=outC(ord);
r=struct('exponents',outE(:).','coefficients',outC(:).');
end
