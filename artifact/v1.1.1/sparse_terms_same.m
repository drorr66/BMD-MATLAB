function tf = sparse_terms_same(a,b)
%SPARSE_TERMS_SAME Exact structural equality of normalized term lists.
tf = isequal(a.exponents,b.exponents) && isequal(a.coefficients,b.coefficients);
end
