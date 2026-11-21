#!/usr/bin/env python3
"""Independent re-verification of Cohl's Appell F4 identity.

Written from equation (1) of proof.md and the F4 convention of DLMF 16.13.4
only, NOT from the proof's own check programs.

  S1  Equation (2) of the proof as a polynomial identity in S, k,n <= 10.
  S2  The Vandermonde step of the proof, symbolically in A,B,r for k <= 6.
  S3  Full statement, coefficientwise in TWO formal variables t and
      y = (x-1)/2: [t^r y^s] for all 0 <= s <= r <= 12, exact Fractions,
      for several rational (alpha,beta,gamma) including non-integer and
      negative non-integer values.
  S4  Full statement at fixed rational x: coefficients [t^r], r <= 12, as
      exact Fractions, for a table of rational (alpha, beta, gamma, x).
  S5  High-precision numeric summation of both literal series at points
      inside the common convergence region (mpmath dps=60, two cutoffs).
  S6  Numeric check of the Horn H3 identity used in the proof's "known
      reductions" section:
        F4(g,B;A,B;u,v) = (1-t)^(-g) H3(g,S-g;A; w,(1-t)w),
      with H3(a,b;c;x,y) = sum_{m,n}(a)_{2m+n}(b)_n/((c)_{m+n} m! n!) x^m y^n.
  S7  Proof steps (5) and (6) as identities of formal t-series (exact
      Fractions through t^12, several k and rational parameters).

By default S5 and S6 run a reduced set of cases so that checks/run_all.py
finishes quickly; --full runs the complete set whose log is saved in
checks/glm-output.txt.
"""

import argparse
import math
import sys
from fractions import Fraction

import mpmath
import sympy

K = 12  # all exact comparisons go through degree 12 in t


def poch(a, n):
    r = Fraction(1) if isinstance(a, Fraction) else 1
    for i in range(n):
        r *= a + i
    return r


def fact(n):
    return math.factorial(n)


# ---------------------------------------------------------------- S1: eq (2)
def check_eq2():
    S = sympy.symbols("S")
    for k in range(0, 11):
        for n in range(0, k + 1):
            # 4^k (S/2)_k ((S+1)/2)_k = (S)_{2k}  and
            # (S)_{2k} = (S)_{k+n} (S+k+n)_{k-n}  (so the ratio in (2) equals
            # the polynomial (S+k+n)_{k-n}; poles cancel because n <= k).
            d1 = sympy.expand(
                4 ** k * poch(S / 2, k) * poch((S + 1) / 2, k) - poch(S, 2 * k)
            )
            d2 = sympy.expand(poch(S, k + n) * poch(S + k + n, k - n) - poch(S, 2 * k))
            if d1 != 0 or d2 != 0:
                print(f"S1 FAIL: equation (2) wrong at (k,n)=({k},{n})")
                return False
    print("S1  equation (2) holds as a polynomial identity in S for all k<=10, 0<=n<=k ... OK")
    return True


# ------------------------------------------------------- S2: Vandermonde step
def check_vandermonde():
    a, b, r = sympy.symbols("a b r")
    for k in range(0, 7):
        def binom(x, j):  # generalized C(x,j) = (x-j+1)_j / j!
            return poch(x - j + 1, j) / sympy.factorial(j)

        tot = 0
        for i in range(0, k + 1):
            tot += binom(a + k - 1, k - i) * binom(b + r + k - 1, i)
        c = a + b + r + 2 * k - 2  # = S+k+r+k-1 with S=a+b-1
        target = sympy.expand(poch(c - k + 1, k) / sympy.factorial(k))  # = C(c,k)
        if sympy.expand(tot - target) != 0:
            print(f"S2 FAIL: Vandermonde step wrong at k={k}")
            return False
    print("S2  sum_i C(A+k-1,k-i) C(B+r+k-1,i) = C(A+B+r+2k-2,k) = (S+k+r)_k/k!, k<=6 ... OK")
    return True


# ------------------------------- statement series, exact, both formal in t,y
def lhs_coeffs(alpha, beta, gamma, K):
    """[t^r y^s] of F4(gamma,B;A,B; t*y, t*(y+1)), DLMF 16.13.4 convention,
    with A=alpha+1, B=beta+1 and y=(x-1)/2.  Returns {(r,s): Fraction}."""
    A, B = alpha + 1, beta + 1
    out = {}
    for i in range(0, K + 1):
        for j in range(0, K + 1 - i):
            coef = (
                Fraction(1)
                * poch(gamma, i + j)
                * poch(B, i + j)
                / (poch(A, i) * poch(B, j) * fact(i) * fact(j))
            )
            # u^i v^j = t^{i+j} y^i (1+y)^j
            for h in range(0, j + 1):  # (1+y)^j = sum_h C(j,h) y^h
                s = i + h
                key = (i + j, s)
                add = coef * fact(j) / (fact(h) * fact(j - h))
                out[key] = out.get(key, Fraction(0)) + add
    return out


def rhs_coeffs(alpha, beta, gamma, K):
    """[t^r y^s] of the right side of (1) exactly as printed (x-1=2y).

    Term = c_{m,n} 4^m (-4)^n y^{m+n} t^{m+2n} (1-t)^{-gamma-2m-2n}, with
    (1-t)^{-p} = sum_r (p)_r/r! t^r,  and c_{m,n} the literal coefficient:
      (S-gamma)_n (gamma)_{m+n} (S/2)_{m+n} ((S+1)/2)_{m+n}
        / ((A)_{m+n} (S)_{m+2n} m! n!).
    """
    A, S = alpha + 1, alpha + beta + 1
    out = {}
    for m in range(0, K + 1):
        for n in range(0, K + 1 - m):
            if m + 2 * n > K:
                continue
            c = (
                Fraction(1)
                * poch(S - gamma, n)
                * poch(gamma, m + n)
                * poch(Fraction(S, 2), m + n)
                * poch(Fraction(S + 1, 2), m + n)
                / (poch(A, m + n) * poch(S, m + 2 * n) * fact(m) * fact(n))
            )
            base = c * Fraction(4) ** m * Fraction(-4) ** n
            p = gamma + 2 * (m + n)
            for extra in range(0, K + 1 - (m + 2 * n)):
                key = (m + 2 * n + extra, m + n)
                out[key] = out.get(key, Fraction(0)) + base * poch(p, extra) / fact(extra)
    return out


PARAMS = [
    (Fraction(1, 2), Fraction(3, 2), Fraction(1, 3)),
    (Fraction(-1, 2), Fraction(5, 3), Fraction(-7, 5)),
    (Fraction(2), Fraction(-3, 2), Fraction(5, 2)),
    (Fraction(-5, 3), Fraction(-2, 7), Fraction(-11, 6)),
    (Fraction(0), Fraction(3, 5), Fraction(-4, 3)),
    (Fraction(3), Fraction(1), Fraction(-2)),  # gamma=-2: both sides terminate
    (Fraction(7, 4), Fraction(-9, 4), Fraction(3, 7)),
]

X_VALUES = [
    Fraction(2, 5),
    Fraction(-3, 4),
    Fraction(1, 3),
    Fraction(7, 2),
    Fraction(-7, 5),
    Fraction(-2),
    Fraction(-5, 9),
    Fraction(1, 2),
    Fraction(1),   # x=1: u=0
    Fraction(-1),  # x=-1: v=0
    Fraction(0),
]


def check_statement():
    """S3 and S4 together: two-variable coefficient comparison plus fixed x."""
    ok = True
    for (al, be, ga) in PARAMS:
        A, B, S = al + 1, be + 1, al + be + 1
        if any(q <= 0 and q.denominator == 1 for q in (A, B, S)):
            print(f"S3 SKIP (alpha,beta,gamma)=({al},{be},{ga}): A/B/S nonpositive integer")
            continue
        lhs = lhs_coeffs(al, be, ga, K)
        rhs = rhs_coeffs(al, be, ga, K)
        keys = sorted(set(lhs) | set(rhs))
        bad = [ky for ky in keys if lhs.get(ky, 0) != rhs.get(ky, 0)]
        if bad:
            print(
                f"S3 FAIL (alpha,beta,gamma)=({al},{be},{ga}): {len(bad)} differing [t^r y^s],"
                f" first {bad[:5]}"
            )
            ok = False
            continue
        nk = sum(1 for (r, s) in keys if s <= r)
        print(
            f"S3  (alpha,beta,gamma)=({al},{be},{ga}): all {nk} coefficients "
            f"[t^r y^s], 0<=s<=r<={K}, agree exactly ... OK"
        )
        # S4: evaluate y=(x-1)/2 at the rational x-values, compare [t^r]
        for x in X_VALUES:
            y = Fraction(x - 1, 2)
            lpoly = [Fraction(0)] * (K + 1)
            rpoly = [Fraction(0)] * (K + 1)
            for (r, s), cf in lhs.items():
                lpoly[r] += cf * y ** s
            for (r, s), cf in rhs.items():
                rpoly[r] += cf * y ** s
            if lpoly != rpoly:
                diff = [r for r in range(K + 1) if lpoly[r] != rpoly[r]]
                print(
                    f"S4 FAIL (alpha,beta,gamma,x)=({al},{be},{ga},{x}): "
                    f"differing t-degrees {diff[:5]}"
                )
                ok = False
            else:
                print(
                    f"S4  (alpha,beta,gamma,x)=({al},{be},{ga},{x}): [t^r], r=0..{K}, agree ... OK"
                )
    return ok


# ------------------------------------------------------ S5/S6: numeric checks
def f4_sum(gamma, B, A, u, v, cutoff):
    """Literal F4 double series (DLMF 16.13.4), total degree <= cutoff."""
    s = mpmath.mpf(0)
    for i in range(cutoff + 1):
        for j in range(cutoff + 1 - i):
            k = i + j
            s += (
                mpmath.rf(gamma, k)
                * mpmath.rf(B, k)
                / (mpmath.rf(A, i) * mpmath.rf(B, j))
                / mpmath.factorial(i)
                / mpmath.factorial(j)
                * u ** i
                * v ** j
            )
    return s


def rhs_sum(alpha, beta, gamma, x, t, cutoff):
    """Right side of (1) exactly as printed, cutoff on m+2n."""
    A, S = alpha + 1, alpha + beta + 1
    first = 2 * t * (x - 1) / (t - 1) ** 2
    second = -2 * t * t * (x - 1) / (t - 1) ** 2
    pf = [mpmath.mpf(1)]
    ps = [mpmath.mpf(1)]
    sg = [mpmath.mpf(1)]  # (S-gamma)_n
    for _ in range(cutoff + 1):
        pf.append(pf[-1] * first)
        ps.append(ps[-1] * second)
    for n in range(1, cutoff + 1):
        sg.append(sg[-1] * (S - gamma + n - 1))
    s = mpmath.mpf(0)
    for n in range(0, cutoff // 2 + 1):
        for m in range(0, cutoff - 2 * n + 1):
            k = m + n
            if m + 2 * n > cutoff:
                continue
            coef = (
                sg[n]
                * mpmath.rf(gamma, k)
                * mpmath.rf(mpmath.mpf(S) / 2, k)
                * mpmath.rf(mpmath.mpf(S + 1) / 2, k)
                / (mpmath.rf(A, k) * mpmath.rf(S, m + 2 * n))
                / mpmath.factorial(m)
                / mpmath.factorial(n)
            )
            s += coef * pf[m] * ps[n]
    return (1 - t) ** (-gamma) * s


def h3_sum(a, b, c, x, y, cutoff):
    """Horn H3(a,b;c;x,y) = sum_{m,n} (a)_{2m+n}(b)_n/((c)_{m+n} m! n!) x^m y^n."""
    s = mpmath.mpf(0)
    for n in range(0, cutoff + 1):
        for m in range(0, cutoff + 1 - n):
            s += (
                mpmath.rf(a, 2 * m + n)
                * mpmath.rf(b, n)
                / mpmath.rf(c, m + n)
                / mpmath.factorial(m)
                / mpmath.factorial(n)
                * x ** m
                * y ** n
            )
    return s


# (alpha, beta, gamma, x, t, cutoff1, cutoff2): all inside the common region
# sqrt|u|+sqrt|v|<1 and 4|u|<|1-t|^2, gamma nonterminating.
# Cutoffs are sized so the truncation tail at the larger cutoff is below
# 1e-40 relative: the left series decays like L^{2k}, the right like R^k.
NUM_CASES = [
    (mpmath.mpf(1) / 2, mpmath.mpf(3) / 2, mpmath.mpf(1) / 3, mpmath.mpf("0.6"), mpmath.mpf("0.2"), 110, 180),
    (mpmath.mpf(-1) / 2, mpmath.mpf(5) / 3, mpmath.mpf(-7) / 5, mpmath.mpf("-0.75"), mpmath.mpf("0.12"), 220, 300),  # R=0.54: slow
    (mpmath.mpf(2), mpmath.mpf(-3) / 2, mpmath.mpf(5) / 2, mpmath.mpf(3), mpmath.mpf("-0.08"), 160, 240),  # L=0.68: slow
    (mpmath.mpf(-5) / 3, mpmath.mpf(-2) / 7, mpmath.mpf(-11) / 6, mpmath.mpf(1) / 3, mpmath.mpf("0.12"), 110, 180),
    (mpmath.mpf(1) / 2, mpmath.mpf(3) / 2, mpmath.mpf(1) / 3, mpmath.mpf(2), mpmath.mpf("0.15"), 200, 280),  # L=0.75: slow
    (mpmath.mpf(9) / 2, mpmath.mpf(-7) / 3, mpmath.mpf(8) / 5, mpmath.mpf(0), mpmath.mpf("-0.3"), 210, 280),  # L=0.77: slow
]


# --full restores the six cases and cutoff 300 of the saved full run.
ACTIVE_CASES = NUM_CASES
S6_CUTOFF = 300


def check_numeric():
    mpmath.mp.dps = 60
    worst = mpmath.mpf(0)
    for (al, be, ga, x, t, n1, n2) in ACTIVE_CASES:
        A, B = al + 1, be + 1
        u, v = t * (x - 1) / 2, t * (x + 1) / 2
        L = mpmath.sqrt(abs(u)) + mpmath.sqrt(abs(v))
        R = 4 * abs(u) / abs(1 - t) ** 2
        assert L < 1 and R < 1, "case outside the convergence region"
        left1, left2 = f4_sum(ga, B, A, u, v, n1), f4_sum(ga, B, A, u, v, n2)
        right1, right2 = rhs_sum(al, be, ga, x, t, n1), rhs_sum(al, be, ga, x, t, n2)
        d2 = abs(left2 - right2) / max(1, abs(left2), abs(right2))
        drift = max(abs(left1 - left2), abs(right1 - right2)) / max(1, abs(left2), abs(right2))
        worst = max(worst, d2)
        print(
            f"S5  alpha={al} beta={be} gamma={ga} x={x} t={t}  (L={mpmath.nstr(L,4)}, R={mpmath.nstr(R,4)})\n"
            f"     left  = {mpmath.nstr(left2, 40)}\n"
            f"     right = {mpmath.nstr(right2, 40)}\n"
            f"     scaled diff at cutoff {n2} = {mpmath.nstr(d2, 5)}, "
            f"cutoff-doubling drift = {mpmath.nstr(drift, 5)}"
        )
        if d2 > mpmath.mpf(10) ** (-40) or drift > mpmath.mpf(10) ** (-40):
            print("S5 FAIL")
            return False
    print(f"S5  maximum scaled difference over all cases: {mpmath.nstr(worst, 5)} ... OK")
    return True


def check_h3():
    mpmath.mp.dps = 60
    for (al, be, ga, x, t, _n1, n2) in ACTIVE_CASES[:4]:
        A, B, S = al + 1, be + 1, al + be + 1
        u, v = t * (x - 1) / 2, t * (x + 1) / 2
        w = u / (1 - t) ** 2
        left = f4_sum(ga, B, A, u, v, S6_CUTOFF)
        right = (1 - t) ** (-ga) * h3_sum(ga, S - ga, A, w, (1 - t) * w, S6_CUTOFF)
        d = abs(left - right) / max(1, abs(left), abs(right))
        print(
            f"S6  alpha={al} beta={be} gamma={ga} x={x} t={t}: "
            f"F4 vs (1-t)^-g H3(g,S-g;A;w,(1-t)w) scaled diff = {mpmath.nstr(d, 5)}"
        )
        if d > mpmath.mpf(10) ** (-40):
            print("S6 FAIL")
            return False
    print("S6  Horn H3 reduction of the report's novelty section verified numerically ... OK")
    return True


# ------------------------------------------------- S7: proof steps (5) and (6)
def check_steps_5_6():
    cases = [
        (Fraction(1, 4), Fraction(7, 3), 3),
        (Fraction(-7, 5), Fraction(-2, 9), 5),
        (Fraction(13, 6), Fraction(-11, 7), 0),
    ]
    T = K
    for (g, S, k) in cases:
        # (5) left:  sum_r (g+k)_r/r! t^r (S+k+r)_k/k!
        a = [Fraction(0)] * (T + 1)
        for r in range(0, T + 1):
            a[r] = poch(g + k, r) / fact(r) * poch(S + k + r, k) / fact(k)
        # (5) right: [z^k] (1-z)^{g-S} (1-t-z)^{-g-k} expanded in t:
        # (1-t-z)^{-g-k} = sum_n (g+k)_n/n! t^n (1-z)^{-g-k-n},
        # [z^k](1-z)^{-S-k-n} = (S+k+n)_k/k!
        b = [Fraction(0)] * (T + 1)
        for n in range(0, T + 1):
            b[n] = poch(g + k, n) / fact(n) * poch(S + k + n, k) / fact(k)
        # (6): (1-t)^{-g-2k} sum_{n<=k} (S-g)_n (S+k+n)_{k-n}/(n!(k-n)!) (-t)^n
        c = [Fraction(0)] * (T + 1)
        pp = g + 2 * k
        for r in range(0, T + 1):
            acc = Fraction(0)
            for n in range(0, min(k, r) + 1):
                poly = poch(S - g, n) * poch(S + k + n, k - n) / (fact(n) * fact(k - n))
                acc += ((-1) ** n) * poly * poch(pp, r - n) / fact(r - n)
            c[r] = acc
        bad = [r for r in range(T + 1) if not (a[r] == b[r] == c[r])]
        if bad:
            print(f"S7 FAIL gamma={g} S={S} k={k}: differing t-degrees {bad[:5]}")
            return False
        print(f"S7  (5) and (6) hold as exact t-series through t^{T} for gamma={g}, S={S}, k={k} ... OK")
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true',
                        help='run the complete case lists of checks/glm-output.txt')
    args = parser.parse_args()
    global ACTIVE_CASES, S6_CUTOFF
    ACTIVE_CASES = NUM_CASES if args.full else NUM_CASES[:3]
    S6_CUTOFF = 300 if args.full else 150
    print(f"python {sys.version.split()[0]}, sympy {sympy.__version__}, mpmath {mpmath.__version__}\n")
    res = {
        "S1": check_eq2(),
        "S2": check_vandermonde(),
        "S3/S4": check_statement(),
        "S5": check_numeric(),
        "S6": check_h3(),
        "S7": check_steps_5_6(),
    }
    print()
    for name in sorted(res):
        print(f"  {name}: {'PASS' if res[name] else 'FAIL'}")
    allok = all(res.values())
    print("ALL CHECKS PASSED" if allok else "SOME CHECKS FAILED")
    return 0 if allok else 1


if __name__ == "__main__":
    raise SystemExit(main())
