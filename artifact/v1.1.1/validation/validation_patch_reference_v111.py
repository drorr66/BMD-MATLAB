#!/usr/bin/env python3
"""Independent rationale check for the v1.1.1 validation-only patch."""
import math

x=0.999
coeff=[math.comb(32,k)*((-1)**k) for k in range(33)]
naive=sum(c*(x**k) for k,c in enumerate(coeff))
closed=(1-x)**32
err=abs(naive-closed)
print('v1.1.1 validation patch rationale')
print(f'binomial_diff_04 direct coefficient evaluation at x=0.999: {naive:.17g}')
print(f'closed-form (1-x)^32: {closed:.17g}')
print(f'absolute error: {err:.17g}')
print(f'old 1e-8 gate multiples: {err/1e-8:.9g}')
assert err > 1e-8
print('PASS: x=0.999 is unsuitable as a correctness gate for this alternating polynomial.')
