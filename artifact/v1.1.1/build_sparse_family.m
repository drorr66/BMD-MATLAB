function s = build_sparse_family(family,n,param)
%BUILD_SPARSE_FAMILY Build a sparse term-list baseline.
%
% Representation:
%   s.exponents    uint64 row vector, strictly descending
%   s.coefficients double row vector, same length
%
% This is intentionally a simple MATLAB-native term-list baseline rather
% than Symbolic Math Toolbox. It answers the key methodological question:
% does BMD beat dense MATLAB only because dense vectors store zeros?
if nargin < 3, param = 0; end
switch family
    case 'monomial'
        s = sparse_terms(uint64(n),1);
    case 'geometric_sum'
        s = sparse_terms(uint64(n:-1:0),ones(1,n+1));
    case 'weighted_sum'
        if n==0
            s = sparse_terms(uint64([]),[]);
        else
            s = sparse_terms(uint64(n:-1:1),double(n:-1:1));
        end
    case 'binomial'
        s = sparse_terms(uint64(0),1);
        factor = sparse_terms(uint64([1 0]),[1 1]);
        for k=1:n
            s = sparse_terms_multiply(s,factor);
        end
    case 'product_powers_plus_1'
        r = param;
        s = sparse_terms(uint64(0),1);
        for k=1:n
            base = sparse_terms(uint64([k 0]),[1 1]);
            factor = sparse_terms(uint64(0),1);
            for j=1:r
                factor = sparse_terms_multiply(factor,base);
            end
            s = sparse_terms_multiply(s,factor);
        end
    otherwise
        error('BMD:UnknownBenchmark','Unknown family %s.',family);
end
end
