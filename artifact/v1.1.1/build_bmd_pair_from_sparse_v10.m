function [mgr,refs] = build_bmd_pair_from_sparse_v10(pSparse,qSparse)
%BUILD_BMD_PAIR_FROM_SPARSE_V10 Convert two sparse term lists to one BMD manager.
%
% Unit-coefficient supports use the direct indicator builder.  General exact
% integer coefficients use the existing canonical add/monomial operations.
% The returned manager is compacted so conversion cost includes removal of
% construction garbage and operation timings start from reachable operands.

mgr=BMDManager();
p=onePolynomial(mgr,pSparse);
q=onePolynomial(mgr,qSparse);
[mgr,refs]=mgr.compact([p;q]);
end

function r=onePolynomial(mgr,s)
if isempty(s.exponents), r=mgr.zero(); return; end
c=s.coefficients;
if any(~isfinite(c)) || any(c~=round(c))
    error('BMD:V10Coefficient','v1.0 BMD conversion requires exact integer coefficients.');
end
if all(c==1)
    r=mgr.indicatorExponents(s.exponents);
    return;
end
if all(c==c(1))
    r=mgr.scale(mgr.indicatorExponents(s.exponents),c(1));
    return;
end
r=mgr.zero();
% Build low exponents first to keep the conversion deterministic.
[e,ord]=sort(s.exponents,'ascend'); cc=c(ord);
for k=1:numel(e)
    r=mgr.add(r,mgr.monomial(double(e(k)),cc(k)));
end
end
