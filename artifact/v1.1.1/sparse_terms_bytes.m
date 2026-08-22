function bytes = sparse_terms_bytes(s)
%SPARSE_TERMS_BYTES Numeric payload bytes of exponent+coefficient arrays.
e=s.exponents; c=s.coefficients; %#ok<NASGU>
i1=whos('e'); i2=whos('c');
bytes=i1.bytes+i2.bytes;
end
