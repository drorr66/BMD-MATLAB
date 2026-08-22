function s = build_sparse_shifted_geometric(blocks,stridePower)
%BUILD_SPARSE_SHIFTED_GEOMETRIC Sparse shifted geometric factor.
if blocks<1, error('BMD:BadGrid','blocks must be positive.'); end
stride=uint64(2^double(stridePower));
exps=stride*uint64(0:blocks-1);
s=sparse_terms(exps,ones(1,blocks));
end
