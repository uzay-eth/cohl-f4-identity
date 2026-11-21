#!/usr/bin/env python3
"""Sum the two double series as printed, with independent degree cutoffs.

Each case is evaluated at 85 and 115 decimal working digits with two
different total-degree cutoffs; a case fails unless all three values agree
to better than 1e-65 relative to their own magnitudes. Random cases are
seeded, so every run tests the same points.

By default a reduced number of seeded cases per family runs, so that
checks/run_all.py finishes quickly; the nine deliberate special cases
always run. --full runs all 36 seeded cases per family and reproduces the
saved log checks/numeric.txt.
"""

import argparse
from random import Random

import mpmath as mp


def products(start, length, step=1):
    values = [mp.mpf(1)]
    for index in range(length):
        values.append(values[-1] * (start + step * index))
    return values


def powers(base, degree):
    values = [mp.mpf(1)]
    for _ in range(degree):
        values.append(values[-1] * base)
    return values


def evaluate(alpha, beta, gamma, x, t, degree):
    a, b, s = alpha + 1, beta + 1, alpha + beta + 1
    u, v = t * (x - 1) / 2, t * (x + 1) / 2
    first = 2 * t * (x - 1) / (t - 1) ** 2
    second = -2 * t**2 * (x - 1) / (t - 1) ** 2
    factorial = products(mp.mpf(1), degree)
    ag = products(gamma, degree)
    aa = products(a, degree)
    bb = products(b, degree)
    sg = products(s - gamma, degree)
    half_s = products(s / 2, degree)
    half_s_plus = products((s + 1) / 2, degree)
    ss = products(s, 2 * degree)

    left_u = [power / (aa[i] * factorial[i]) for i, power in enumerate(powers(u, degree))]
    left_v = [power / (bb[j] * factorial[j]) for j, power in enumerate(powers(v, degree))]
    right_u = [power / factorial[i] for i, power in enumerate(powers(first, degree))]
    right_v = [power * sg[j] / factorial[j] for j, power in enumerate(powers(second, degree))]

    left, right = [], []
    for k in range(degree + 1):
        if ag[k] == 0:
            break
        left.append(ag[k] * bb[k] * mp.fsum(left_u[k - j] * left_v[j] for j in range(k + 1)))
        diagonal = mp.fsum(right_u[k - j] * right_v[j] / ss[k + j] for j in range(k + 1))
        right.append(ag[k] * half_s[k] * half_s_plus[k] * diagonal / aa[k])
    return mp.fsum(left), mp.fsum(right) * (1 - t) ** (-gamma)


def growth_rates(x, t):
    u, v = t * (x - 1) / 2, t * (x + 1) / 2
    return (mp.sqrt(abs(u)) + mp.sqrt(abs(v))) ** 2, 4 * abs(u / (1 - t) ** 2)


def cases(seed, per_family):
    rng = Random(seed)

    def grid(real_min, real_max, imag_limit=0):
        real = mp.mpf(rng.randint(real_min, real_max)) / 1000
        if imag_limit:
            imag = mp.mpf(rng.randint(-imag_limit, imag_limit)) / 1000
            return mp.mpc(real, imag)
        return real

    for imag_limit in (0, 1050):
        for index in range(per_family):
            limit = mp.mpf('.83' if index >= 30 else '.65')
            lower = mp.mpf('.68' if index >= 30 else '.05')
            while True:
                x = grid(-2350, 2950, imag_limit)
                t = grid(-650, 650, min(imag_limit, 650))
                if lower < max(growth_rates(x, t)) < limit:
                    break
            while True:
                alpha = grid(-2800, 3600, imag_limit)
                beta = grid(-2200, 4250, imag_limit)
                gamma = grid(-3600, 4800, imag_limit)
                # Keep randomized cases away from accidental singularities.
                if all(abs(z - mp.nint(z.real)) > mp.mpf('.03') or z.real > 0
                       for z in (alpha + 1, beta + 1, alpha + beta + 1)):
                    break
            yield ('complex' if imag_limit else 'real'), (alpha, beta, gamma, x, t)

    specials = (
        ('t=0', ('0.3', '0.5', '1.7', '2.2', '0')),
        ('x=1', ('-1.8', '0.6', '1.4', '1', '0.8')),
        ('x=-1', ('0.2', '-1.1', '-0.7', '-1', '-0.55')),
        ('gamma=-5', ('0.8', '1.2', '-5', '-0.3', '-0.21')),
        ('s-gamma=0', ('1.25', '0.5', '2.75', '0.4', '0.14')),
        ('s-gamma=-3', ('-0.625', '1.5', '4.875', '-0.3', '-0.2')),
        ('near A pole', ('-0.999999999999', '0.2', '1.5', '0.4', '0.08')),
        ('near B pole', ('0.3', '-0.999999999999', '-0.4', '0.4', '0.08')),
        ('near S pole', ('0.4', '-1.400000000001', '-0.7', '0.4', '-0.1')),
    )
    for label, values in specials:
        yield label, tuple(map(mp.mpf, values))


def distance(first, second):
    return abs(first - second) / max(mp.mpf(1), abs(first), abs(second))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--per-family', type=int, default=None,
                        help='seeded cases per family (default 2, or 36 with --full)')
    parser.add_argument('--full', action='store_true',
                        help='run all 81 cases of the saved log checks/numeric.txt')
    args = parser.parse_args()
    per_family = args.per_family if args.per_family is not None \
        else (36 if args.full else 2)
    seed = 1700261948
    with mp.workdps(115):
        print(f'mpmath {mp.__version__}, seed={seed}, 85 and 115 decimal working digits')
        print('Cutoffs include all terms of total degree <= N.')
        print('Errors are scaled by max(1, abs(first value), abs(second value)).')
        maxima = [mp.mpf(0)] * 3
        for index, (label, parameters) in enumerate(cases(seed, per_family), 1):
            left_rate, right_rate = growth_rates(*parameters[-2:])
            rate = max(left_rate, right_rate)
            assert rate < 1
            cutoff = 24 if rate == 0 else max(24, int(mp.ceil(-80 / mp.log10(rate))) + 12)
            last_cutoff = cutoff + max(35, cutoff // 5)
            with mp.workdps(85):
                early = evaluate(*parameters, cutoff)
            left, right = evaluate(*parameters, last_cutoff)
            errors = (distance(left, right), distance(left, early[0]), distance(right, early[1]))
            maxima = list(map(max, zip(maxima, errors)))
            if max(errors) >= mp.mpf('1e-65'):
                raise ArithmeticError(f'case {index} failed: {parameters!r}, errors={errors!r}')
            print(f'{index:02d} {label}  (alpha,beta,gamma,x,t) = '
                  + ', '.join(mp.nstr(parameter, 16) for parameter in parameters))
            print(f'   rates {mp.nstr(left_rate, 6)}, {mp.nstr(right_rate, 6)}'
                  f' | N={cutoff},{last_cutoff} | errors '
                  + ', '.join(mp.nstr(error, 6) for error in errors))
            print('   LHS', mp.nstr(left, 72))
            print('   RHS', mp.nstr(right, 72), flush=True)
        print(f'PASS {index} cases; maximum errors (LHS vs RHS, LHS stability, RHS stability) = '
              + ', '.join(mp.nstr(error, 8) for error in maxima))


if __name__ == '__main__':
    main()
