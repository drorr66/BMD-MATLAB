function p = build_dense_shifted_geometric(blocks,stridePower)
%BUILD_DENSE_SHIFTED_GEOMETRIC Sum_{j=0}^{blocks-1} x^(j*2^stridePower).
if blocks<1, error('BMD:BadGrid','blocks must be positive.'); end
stride=2^double(stridePower);
degree=(blocks-1)*stride;
p=zeros(1,degree+1);
exps=(0:blocks-1)*stride;
p(degree-exps+1)=1;
end
