function y = sparse_terms_evaluate(s,x)
%SPARSE_TERMS_EVALUATE Evaluate sparse term list at scalar x.
if isempty(s.exponents), y=0; return; end
y = sum(s.coefficients .* (double(x) .^ double(s.exponents)));
end
