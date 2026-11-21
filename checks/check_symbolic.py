#!/usr/bin/env python3
"""Compare every x-polynomial coefficient of the two formal t-series.

Compares [t^(k+r) y^k] for all k, r >= 0 with k + r <= --order, y = (x-1)/2,
as exact polynomials over QQ in the remaining parameters (half-parameter
factors cancelled by exact division). The default order is sized so that
checks/run_all.py finishes quickly; --order 40 reproduces the saved log
checks/symbolic.txt.
"""

import argparse
from functools import cache
from math import comb, perm

import sympy
from sympy.polys.rings import ring


polynomials, s, z = ring('S, Z', sympy.QQ)


@cache
def rising(base, shift, length, step=1):
    value = polynomials.one
    for index in range(length):
        value *= base + shift + step * index
    return value


@cache
def half_parameter_ratio(k, n):
    # Before division the factors are 2^k(S/2)_k and 2^k((S+1)/2)_k.
    numerator = rising(s, 0, k, 2) * rising(s, 1, k, 2)
    return numerator.exquo(rising(s, 0, k + n))


def left_polynomial(k, r):
    # Z stands for A=alpha+1 here, with B=S-A+1. The independent
    # gamma factor (gamma+k)_r is multiplied in after dropping A.
    a, b = z, s - z + 1
    coefficient = polynomials.zero
    for i in range(k + 1):
        coefficient += comb(k, i) * rising(a, i, k - i) * rising(b, k + r - i, i)
    if coefficient.degree(z) > 0:
        raise ArithmeticError(f'unexpected dependence on A at k={k}, r={r}')
    return coefficient * rising(z, k, r)


def right_polynomial(k, r):
    # Z stands for gamma here. The t^r factor left after removing t^k
    # comes from (-t)^n and the binomial series of (1-t)^(-gamma-2k).
    coefficient = polynomials.zero
    for n in range(min(k, r) + 1):
        multiplier = (-1)**n * comb(k, n) * perm(r, n)
        coefficient += (multiplier * half_parameter_ratio(k, n)
                        * rising(s - z, 0, n) * rising(z, 2 * k, r - n))
    return coefficient


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--order', type=int, default=30)
    args = parser.parse_args()
    if args.order < 0:
        parser.error('--order must be nonnegative')

    print(f'SymPy {sympy.__version__}, exact polynomial arithmetic over QQ', flush=True)
    print('Both series expanded in t and y=(x-1)/2. No parameter values assigned.', flush=True)
    print('For each [t^(k+r) y^k] multiply both sides by (A)_k k! r! / (gamma)_k.', flush=True)
    count = 0
    for degree in range(args.order + 1):
        for k in range(degree + 1):
            r = degree - k
            residual = left_polynomial(k, r) - right_polynomial(k, r)
            if residual:
                raise ArithmeticError(f'[t^{degree} y^{k}] mismatch: {residual}')
            count += 1
        print(f't^{degree}: {degree + 1} zero residuals', flush=True)
    print(f'PASS through t^{args.order} inclusive, {count} exact coefficient comparisons.')


if __name__ == '__main__':
    main()
