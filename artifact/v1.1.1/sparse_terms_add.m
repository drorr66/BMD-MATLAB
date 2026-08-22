function r = sparse_terms_add(a,b)
%SPARSE_TERMS_ADD Add two normalized descending sparse term lists.
ea=a.exponents; eb=b.exponents; ca=a.coefficients; cb=b.coefficients;
na=numel(ea); nb=numel(eb);
if na==0, r=b; return; end
if nb==0, r=a; return; end
outE=zeros(1,na+nb,'uint64'); outC=zeros(1,na+nb);
i=1; j=1; w=0;
while i<=na && j<=nb
    if ea(i)>eb(j)
        w=w+1; outE(w)=ea(i); outC(w)=ca(i); i=i+1;
    elseif eb(j)>ea(i)
        w=w+1; outE(w)=eb(j); outC(w)=cb(j); j=j+1;
    else
        cc=ca(i)+cb(j);
        if cc~=0
            w=w+1; outE(w)=ea(i); outC(w)=cc;
        end
        i=i+1; j=j+1;
    end
end
while i<=na, w=w+1; outE(w)=ea(i); outC(w)=ca(i); i=i+1; end
while j<=nb, w=w+1; outE(w)=eb(j); outC(w)=cb(j); j=j+1; end
r=struct('exponents',outE(1:w),'coefficients',outC(1:w));
end
