"""Independent check of Cohl's Appell F4 identity.

Written from the posted image and DLMF 16.13.4 only; it does not reuse the
proof's own checks/check_numeric.py or checks/check_symbolic.py.

  F4(g, b+1; a+1, b+1; t(x-1)/2, t(x+1)/2)
    = (1-t)^(-g) sum_{m,n} (a+b-g+1)_n (g)_{m+n} ((a+b+1)/2)_{m+n} ((a+b+2)/2)_{m+n}
                         / ((a+1)_{m+n} (a+b+1)_{m+2n} m! n!)
                 * (2t(x-1)/(t-1)^2)^m * (-2t^2(x-1)/(t-1)^2)^n

with F4(p,q;r,s;X,Y) = sum_{i,j} (p)_{i+j}(q)_{i+j}/((r)_i (s)_j i! j!) X^i Y^j.

Part 1: exact rational comparison of [t^N] of both sides, N = 0..K.
Part 2: fully symbolic comparison (a, b, g, x indeterminates) for small N.
Part 3: mpmath comparison at points inside the convergence domain, with the
        left side also evaluated as sum_n (g)_n/(a+1)_n P_n^{(a,b)}(x) t^n.
Part 4: negative controls: perturbed right sides must fail.

Default arguments are sized so checks/run_all.py finishes quickly; the file
checks/opus-output.txt is the complete log of the full-size run
`python3 checks/check_opus.py --full` (same output as `check_opus.py 16 6`).
"""
import argparse
import sys
from fractions import Fraction as Fr

import mpmath as mp
import sympy as sp


def poch(a, n):
    r = a * 0 + 1  # exact 1 of the same type as a (Fraction or SymPy)
    for i in range(n):
        r = r * (a + i)
    return r


def fact(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return Fr(r)  # exact; int / int would give a float


def binom_series(c, r):
    """[t^r] (1-t)^(-c) = (c)_r / r!"""
    return poch(c, r) / fact(r)


def lhs_coeffs(a, b, g, x, K):
    A, B = a + 1, b + 1
    u0, v0 = (x - 1) / 2, (x + 1) / 2
    out = []
    for N in range(K + 1):
        s = 0
        for i in range(N + 1):
            j = N - i
            s += (poch(g, N) * poch(B, N) / (poch(A, i) * poch(B, j) * fact(i) * fact(j))
                  * u0 ** i * v0 ** j)
        out.append(s)
    return out


def rhs_coeffs(a, b, g, x, K, sign=-1, shift=1, den=0):
    """Expand the printed right side literally: each (m,n) term is
    c_{mn} (2(x-1))^m (-2(x-1))^n t^{m+2n} (1-t)^{-g-2(m+n)}."""
    out = [0] * (K + 1)
    for m in range(K + 1):
        for n in range(K + 1):
            deg = m + 2 * n
            if deg > K:
                continue
            k = m + n
            c = (poch(a + b - g + shift, n) * poch(g, k) * poch((a + b + 1) / 2, k)
                 * poch((a + b + 2) / 2, k)
                 / (poch(a + 1, k) * poch(a + b + 1, m + 2 * n + den) * fact(m) * fact(n)))
            c = c * (2 * (x - 1)) ** m * (sign * 2 * (x - 1)) ** n
            for r in range(K - deg + 1):
                out[deg + r] += c * binom_series(g + 2 * k, r)
    return out


def part1(K):
    F = Fr
    cases = [
        # (alpha, beta, gamma, x)
        (F(1, 3), F(2, 7), F(5, 11), F(3, 5)),
        (F(-1, 2), F(1, 3), F(-7, 4), F(-2, 3)),
        (F(-5, 3), F(-2, 5), F(9, 7), F(7, 2)),
        (F(2), F(-1, 7), F(-13, 6), F(0)),
        (F(-7, 2), F(8, 3), F(1, 9), F(-5, 4)),
        (F(3, 4), F(-9, 4), F(-1, 2), F(1)),      # x = 1: right side is (1-t)^-g
        (F(1, 5), F(4, 5), F(-3), F(11, 3)),      # terminating gamma
        (F(-1, 3), F(-1, 3), F(2, 3), F(-1)),     # x = -1
        (F(5, 2), F(7, 3), F(-11, 5), F(13, 7)),
        (F(-4, 3), F(5, 6), F(17, 3), F(-9, 8)),
    ]
    ok = True
    for (a, b, g, x) in cases:
        A, B, S = a + 1, b + 1, a + b + 1
        for p in (A, B, S):  # literal denominators must not vanish
            assert not (p.denominator == 1 and p <= 0), (a, b, p)
        L = lhs_coeffs(a, b, g, x, K)
        R = rhs_coeffs(a, b, g, x, K)
        diff = [L[i] - R[i] for i in range(K + 1)]
        good = all(d == 0 for d in diff)
        ok &= good
        print(f"alpha={a}, beta={b}, gamma={g}, x={x}: [t^0..t^{K}] equal: {good}")
        print(f"    [t^{K}] left  = {L[K]}")
        print(f"    [t^{K}] right = {R[K]}")
    return ok


def part4(K):
    """Negative controls: perturbed right sides must fail."""
    F = Fr
    a, b, g, x = F(1, 3), F(2, 7), F(5, 11), F(3, 5)
    L = lhs_coeffs(a, b, g, x, K)
    ok = True
    for name, kw in [("second base with + sign", dict(sign=1)),
                     ("(alpha+beta-gamma)_n instead of (alpha+beta-gamma+1)_n", dict(shift=0)),
                     ("(alpha+beta+1)_{m+2n+1} denominator", dict(den=1))]:
        R = rhs_coeffs(a, b, g, x, K, **kw)
        bad = [N for N in range(K + 1) if L[N] != R[N]]
        print(f"perturbation '{name}': first differing [t^N] at N = {bad[0] if bad else None}")
        ok &= bool(bad)
    return ok


def part2(K):
    a, b, g, x = sp.symbols("alpha beta gamma x")
    L = lhs_coeffs(a, b, g, x, K)
    R = rhs_coeffs(a, b, g, x, K)
    ok = True
    for N in range(K + 1):
        d = sp.cancel(sp.together(L[N] - R[N]))
        print(f"symbolic [t^{N}] left - right = {d}")
        ok &= (d == 0)
    return ok


def mp_lhs(a, b, g, x, t, K):
    A, B = a + 1, b + 1
    u, v = t * (x - 1) / 2, t * (x + 1) / 2
    s = mp.mpf(0)
    for i in range(K + 1):
        ci = mp.rf(g, i) * mp.rf(B, i) / (mp.rf(A, i) * mp.factorial(i)) * u ** i
        # inner over j, with (g)_{i+j}(B)_{i+j}/(B)_j/j! v^j built by recurrence
        term = ci
        for j in range(K + 1 - i):
            s += term
            term = term * (g + i + j) * (B + i + j) / ((B + j) * (j + 1)) * v
    return s


def mp_lhs_jacobi(a, b, g, x, t, K):
    s = mp.mpf(0)
    for n in range(K + 1):
        s += mp.rf(g, n) / mp.rf(a + 1, n) * mp.jacobi(n, a, b, x) * t ** n
    return s


def rf_table(a, n):
    """[(a)_0, ..., (a)_n]"""
    out = [mp.mpf(1)]
    for i in range(n):
        out.append(out[-1] * (a + i))
    return out


def mp_rhs(a, b, g, x, t, K):
    S = a + b + 1
    X = 2 * t * (x - 1) / (t - 1) ** 2
    Y = -2 * t ** 2 * (x - 1) / (t - 1) ** 2
    Sg, G, S2, S3 = (rf_table(c, K) for c in (S - g, g, S / 2, (S + 1) / 2))
    Ap, Sp, fa = rf_table(a + 1, K), rf_table(S, 2 * K), rf_table(1, K)
    s = mp.mpf(0)
    for k in range(K + 1):
        for n in range(k + 1):
            m = k - n
            c = (Sg[n] * G[k] * S2[k] * S3[k]
                 / (Ap[k] * Sp[m + 2 * n] * fa[m] * fa[n]))
            s += c * X ** m * Y ** n
    return (1 - t) ** (-g) * s


def mp_h3(a, b, g, x, t, K):
    """(1-t)^-g H3(g, S-g; A; w, (1-t) w), w = t(x-1)/(2(1-t)^2),
    H3(p, q; r; X, Y) = sum (p)_{2m+n} (q)_n / ((r)_{m+n} m! n!) X^m Y^n."""
    S, A = a + b + 1, a + 1
    w = t * (x - 1) / (2 * (1 - t) ** 2)
    y = (1 - t) * w
    G, Sg, Ap, fa = rf_table(g, 2 * K), rf_table(S - g, K), rf_table(A, K), rf_table(1, K)
    s = mp.mpf(0)
    for k in range(K + 1):
        for n in range(k + 1):
            m = k - n
            s += G[2 * m + n] * Sg[n] / (Ap[k] * fa[m] * fa[n]) * w ** m * y ** n
    return (1 - t) ** (-g) * s


KS = (250, 350)


def part3():
    mp.mp.dps = 60
    pts = [
        (mp.mpf("0.3"), mp.mpf("-0.45"), mp.mpf("1.7"), mp.mpf("0.5"), mp.mpf("0.3")),
        (mp.mpf("-1.25"), mp.mpf("0.6"), mp.mpf("-2.35"), mp.mpf("-0.8"), mp.mpf("-0.25")),
        (mp.mpf("2.5"), mp.mpf("-0.7"), mp.mpf("0.42"), mp.mpf("2.2"), mp.mpf("0.12")),
        (mp.mpc("0.3", "0.2"), mp.mpc("-0.4", "0.1"), mp.mpc("1.1", "-0.6"),
         mp.mpc("0.2", "0.7"), mp.mpc("0.15", "0.2")),
        (mp.mpf("-0.6"), mp.mpf("1.3"), mp.mpf("-0.75"), mp.mpf("4"), mp.mpf("-0.08")),
    ]
    worst = mp.mpf(0)
    for (a, b, g, x, t) in pts:
        u, v = t * (x - 1) / 2, t * (x + 1) / 2
        Lq = mp.sqrt(abs(u)) + mp.sqrt(abs(v))
        Rq = 4 * abs(u) / abs(1 - t) ** 2
        res = {}
        for K in KS:
            res[K] = (mp_lhs(a, b, g, x, t, K), mp_lhs_jacobi(a, b, g, x, t, K),
                      mp_rhs(a, b, g, x, t, K), mp_h3(a, b, g, x, t, K))
        l, lj, r, h = res[KS[1]]
        scale = max(1, abs(l))
        e1, e2, e3 = abs(l - r) / scale, abs(lj - r) / scale, abs(h - r) / scale
        tail = max(abs(res[KS[1]][i] - res[KS[0]][i]) for i in range(4)) / scale
        worst = max(worst, e1, e2, e3)
        print(f"alpha={mp.nstr(a, 6)}, beta={mp.nstr(b, 6)}, gamma={mp.nstr(g, 6)}, "
              f"x={mp.nstr(x, 6)}, t={mp.nstr(t, 6)}; L={mp.nstr(Lq, 5)}, R={mp.nstr(Rq, 5)}")
        print(f"    F4 double sum   = {mp.nstr(l, 40)}")
        print(f"    Jacobi gen. fn. = {mp.nstr(lj, 40)}")
        print(f"    right side (1)  = {mp.nstr(r, 40)}")
        print(f"    H3 form         = {mp.nstr(h, 40)}")
        print(f"    |F4-RHS|={mp.nstr(e1, 3)}, |Jac-RHS|={mp.nstr(e2, 3)}, "
              f"|H3-RHS|={mp.nstr(e3, 3)}, change K={KS[0]}->{KS[1]}: {mp.nstr(tail, 3)}")
    print(f"worst relative discrepancy: {mp.nstr(worst, 3)}")
    return worst < mp.mpf(10) ** -40


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('exact_order', nargs='?', type=int, default=None,
                        help='Part 1/Part 4 expansion order K (default 12)')
    parser.add_argument('symbolic_order', nargs='?', type=int, default=None,
                        help='Part 2 symbolic order (default 4)')
    parser.add_argument('--full', action='store_true',
                        help='full-size run of the saved log: K=16, symbolic order 6')
    args = parser.parse_args()
    K = args.exact_order if args.exact_order is not None else (16 if args.full else 12)
    Ks = args.symbolic_order if args.symbolic_order is not None else (6 if args.full else 4)
    print(f"mpmath {mp.__version__}, sympy {sp.__version__}\n")
    print(f"== Part 1: exact rational coefficients, N = 0..{K} ==")
    ok1 = part1(K)
    print(f"\n== Part 2: symbolic in alpha, beta, gamma, x, N = 0..{Ks} ==")
    ok2 = part2(Ks)
    print("\n== Part 4: negative controls ==")
    ok4 = part4(K)
    print("\n== Part 3: numeric, 60 digits ==")
    ok3 = part3()
    print(f"\nPART1 {'PASS' if ok1 else 'FAIL'}, PART2 {'PASS' if ok2 else 'FAIL'}, "
          f"PART3 {'PASS' if ok3 else 'FAIL'}, PART4 {'PASS' if ok4 else 'FAIL'}")
    return 0 if ok1 and ok2 and ok3 and ok4 else 1


if __name__ == "__main__":
    raise SystemExit(main())