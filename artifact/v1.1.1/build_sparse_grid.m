function s = build_sparse_grid(innerN,blocks,stridePower)
%BUILD_SPARSE_GRID Sparse term-list for v0.4 structural-sharing grid.
stride=uint64(2^double(stridePower));
if innerN<0 || blocks<1 || double(innerN)>=double(stride)
    error('BMD:BadGrid','Require innerN>=0, blocks>=1, innerN<2^stridePower.');
end
inner=uint64((0:innerN).');
blockOffsets=stride*uint64(0:blocks-1);
E=inner+blockOffsets;
exps=reshape(E,1,[]);
s=sparse_terms(exps,ones(1,numel(exps)));
end
