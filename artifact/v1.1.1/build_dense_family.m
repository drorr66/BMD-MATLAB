function dense = build_dense_family(family,n,param)
%BUILD_DENSE_FAMILY Construct the same univariate polynomial as a dense vector.
if nargin < 3, param = 0; end
switch family
    case 'monomial'
        dense = [1 zeros(1,n)];
    case 'geometric_sum'
        dense = ones(1,n+1);
    case 'weighted_sum'
        dense = [n:-1:1 0];
    case 'binomial'
        dense = 1;
        for k = 1:n
            dense = conv(dense,[1 1]);
        end
    case 'product_powers_plus_1'
        r = param;
        dense = 1;
        for k = 1:n
            base = [1 zeros(1,k-1) 1];
            factor = 1;
            for j = 1:r
                factor = conv(factor,base);
            end
            dense = conv(dense,factor);
        end
    otherwise
        error('BMD:UnknownBenchmark','Unknown family %s.',family);
end
end
