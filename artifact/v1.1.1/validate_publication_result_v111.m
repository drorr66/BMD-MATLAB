function v = validate_publication_result_v111(m,p,su,ss,dp,dq,denseEligible,denseMax)
%VALIDATE_PUBLICATION_RESULT_V111 Validation-only replacement for v1.1 checkNumeric.
%
% This function is deliberately outside every timed region.  It removes the
% ill-conditioned x=0.999 gate from v1.1 and exposes each validation subcheck.
%
% Dense-eligible results:
%   - exact sparse-vs-sparse term-list equality
%   - coefficient-for-coefficient equality of BMD.toDense vs sparse result
%   - existing tolerance-based dense conv vs sparse check
%
% Dense-ineligible results:
%   - exact sparse-vs-sparse term-list equality
%   - three deterministic modular fingerprints of BMD vs sparse result
%
% A stable x=0.5 evaluation is retained as a diagnostic only.  It is not a
% gate because the stronger coefficient/modular checks above are used.
if nargin < 8 || isempty(denseMax), denseMax=250000; end
v=struct('validation_mode','','sparse_crosscheck',0,'bmd_coefficient_exact',NaN, ...
    'bmd_modular_fingerprint',NaN,'dense_crosscheck',NaN,'stable_eval_diagnostic',0, ...
    'numeric_check',0,'failure_component','');

v.sparse_crosscheck=double(sparse_terms_same(su,ss));
try
    a=m.evaluate(p,0.5); b=sparse_terms_evaluate(su,0.5);
    v.stable_eval_diagnostic=double(nearScalar(a,b));
catch
    v.stable_eval_diagnostic=0;
end

if denseEligible
    v.validation_mode='DENSE_COEFFICIENT_EXACT';
    try
        bd=trimDense(m.toDense(p,denseMax-1));
        sd=trimDense(sparse_terms_to_dense(su,denseMax));
        v.bmd_coefficient_exact=double(isequal(bd,sd));
    catch
        v.bmd_coefficient_exact=0;
    end
    try
        z=trimDense(conv(dp,dq));
        sd2=trimDense(sparse_terms_to_dense(su,denseMax));
        v.dense_crosscheck=double(nearVector(z,sd2));
    catch
        v.dense_crosscheck=0;
    end
    v.numeric_check=double(v.sparse_crosscheck==1 && v.bmd_coefficient_exact==1 && v.dense_crosscheck==1);
else
    v.validation_mode='MODULAR_FINGERPRINT_3X';
    try
        primes=[1000003 1000033 1000037];
        xs=[2 123457 765431];
        ok=true;
        for k=1:numel(primes)
            bv=bmdModEval(m,p,xs(k),primes(k));
            sv=sparseModEval(su,xs(k),primes(k));
            ok=ok && (bv==sv);
        end
        v.bmd_modular_fingerprint=double(ok);
    catch
        v.bmd_modular_fingerprint=0;
    end
    v.numeric_check=double(v.sparse_crosscheck==1 && v.bmd_modular_fingerprint==1);
end

if v.numeric_check~=1
    parts={};
    if v.sparse_crosscheck~=1, parts{end+1}='SPARSE_CROSSCHECK'; end %#ok<AGROW>
    if denseEligible
        if v.bmd_coefficient_exact~=1, parts{end+1}='BMD_COEFFICIENT_EXACT'; end %#ok<AGROW>
        if v.dense_crosscheck~=1, parts{end+1}='DENSE_CROSSCHECK'; end %#ok<AGROW>
    else
        if v.bmd_modular_fingerprint~=1, parts{end+1}='BMD_MODULAR_FINGERPRINT'; end %#ok<AGROW>
    end
    v.failure_component=strjoin(parts,'+');
end
end

function tf=nearScalar(a,b)
scale=max([1 abs(a) abs(b)]); tf=abs(a-b)<=1e-8*scale;
end

function tf=nearVector(a,b)
a=trimDense(a); b=trimDense(b);
if numel(a)~=numel(b), tf=false; return; end
scale=max(1,max(abs(b))); tf=all(abs(a-b)<=1e-9*scale);
end

function p=trimDense(p)
idx=find(p~=0,1,'first'); if isempty(idx), p=0; else, p=p(idx:end); end
end

function y=bmdModEval(m,r,x,prime)
% Integer-exact modular arithmetic in double: prime~1e6 keeps products <2^53.
n=numel(m.levels); memo=NaN(1,n); memo(1)=1; basisMemo=NaN(1,128);
y=mod(mod(double(r(1)),prime)*evalNode(r(2)),prime);
    function val=evalNode(nodeId)
        idx=double(nodeId);
        if ~isnan(memo(idx)), val=memo(idx); return; end
        lev=double(m.levels(idx));
        if lev<1 || lev>127, error('BMD:V111ModLevel','Unexpected non-univariate level.'); end
        if isnan(basisMemo(lev))
            base=mod(double(x),prime);
            for jj=2:lev, base=mod(base*base,prime); end
            basisMemo(lev)=base;
        end
        lw=mod(double(m.lowWeight(idx)),prime); hw=mod(double(m.highWeight(idx)),prime);
        lo=mod(lw*evalNode(m.lowNode(idx)),prime);
        hi=mod(hw*evalNode(m.highNode(idx)),prime);
        val=mod(lo + mod(basisMemo(lev)*hi,prime),prime);
        memo(idx)=val;
    end
end

function y=sparseModEval(s,x,prime)
y=0;
for k=1:numel(s.exponents)
    c=mod(s.coefficients(k),prime);
    pw=powMod(x,s.exponents(k),prime);
    y=mod(y + mod(c*pw,prime),prime);
end
end

function y=powMod(x,e,prime)
y=1; b=mod(double(x),prime); e=uint64(e);
while e>0
    if bitand(e,uint64(1))~=0, y=mod(y*b,prime); end
    e=bitshift(e,-1);
    if e>0, b=mod(b*b,prime); end
end
end
