function [mgr,root] = build_bmd_family(family,n,param)
%BUILD_BMD_FAMILY Construct benchmark polynomials in a fresh BMD manager.
% v0.4 uses the direct bit-recursive builder for geometric_sum.  The old
% sequential construction remains available as geometric_sum_naive so the
% builder experiment can quantify eliminated transient/garbage nodes.
if nargin < 3, param = 0; end
mgr = BMDManager();
switch family
    case 'monomial'
        root = mgr.monomial(n);
    case 'geometric_sum'
        root = mgr.geometricSum(n);
    case 'geometric_sum_naive'
        root = mgr.zero();
        for k = 0:n
            root = mgr.add(root,mgr.monomial(k));
        end
    case 'weighted_sum'
        root = mgr.zero();
        for k = 1:n
            root = mgr.add(root,mgr.monomial(k,k));
        end
    case 'binomial'
        factor = mgr.add(mgr.monomial(1),mgr.one());
        root = mgr.one();
        for k = 1:n
            root = mgr.multiply(root,factor);
        end
    case 'product_powers_plus_1'
        r = param;
        root = mgr.one();
        for k = 1:n
            base = mgr.add(mgr.monomial(k),mgr.one());
            factor = mgr.power(base,r);
            root = mgr.multiply(root,factor);
        end
    otherwise
        error('BMD:UnknownBenchmark','Unknown family %s.',family);
end
end
