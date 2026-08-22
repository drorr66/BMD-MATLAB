function p = build_dense_grid(innerN,blocks,stridePower)
%BUILD_DENSE_GRID Dense coefficient vector for v0.4 structural-sharing grid.
stride=2^double(stridePower);
if innerN<0 || blocks<1 || innerN>=stride
    error('BMD:BadGrid','Require innerN>=0, blocks>=1, innerN<2^stridePower.');
end
degree=(blocks-1)*stride+innerN;
p=zeros(1,degree+1);
for jj=0:blocks-1
    exps=jj*stride+(0:innerN);
    p(degree-exps+1)=1;
end
end
