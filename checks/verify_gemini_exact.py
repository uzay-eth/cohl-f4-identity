#!/usr/bin/env python3
"""
Independent exact rational power-series verification of Howard Cohl's Appell F4 identity.

Written directly from Cohl's posted image (cohl-f4.jpg / DLMF 16.13.4):
  LHS = F4(gamma, beta+1; alpha+1, beta+1; t(x-1)/2, t(x+1)/2)
  RHS = (1-t)^(-gamma) * sum_{m,n>=0}
          [(alpha+beta-gamma+1)_n (gamma)_{m+n} ((alpha+beta+1)/2)_{m+n} ((alpha+beta+2)/2)_{m+n}]
          / [(alpha+1)_{m+n} (alpha+beta+1)_{m+2n} m! n!]
          * (2 t (x-1) / (t-1)^2)^m * (-2 t^2 (x-1) / (t-1)^2)^n

Pure Python 3 standard library (fractions.Fraction); no external dependencies.
"""

from fractions import Fraction
import math


def poch(a: Fraction, k: int) -> Fraction:
    """Exact Pochhammer symbol (a)_k = a (a+1) ... (a+k-1) for integer k >= 0."""
    res = Fraction(1)
    for j in range(k):
        res *= (a + j)
    return res


def lhs_t_coeffs(alpha: Fraction, beta: Fraction, gamma: Fraction, x: Fraction, order: int):
    """
    Exact rational coefficients [t^r] (r = 0..order) of the literal F4 series (DLMF 16.13.4):
      sum_{i,j>=0} [(gamma)_{i+j} (beta+1)_{i+j} / ((alpha+1)_i (beta+1)_j i! j!)] * u^i v^j
    where u = t*(x-1)/2, v = t*(x+1)/2.
    """
    A = alpha + 1
    B = beta + 1
    u_coeff = (x - 1) / 2
    v_coeff = (x + 1) / 2

    # Precompute factorials and Pochhammer tables
    fact = [Fraction(math.factorial(k)) for k in range(order + 1)]
    p_gamma = [poch(gamma, k) for k in range(order + 1)]
    p_B = [poch(B, k) for k in range(order + 1)]
    p_A = [poch(A, k) for k in range(order + 1)]
    u_pow = [u_coeff ** k for k in range(order + 1)]
    v_pow = [v_coeff ** k for k in range(order + 1)]

    coeffs = []
    for r in range(order + 1):
        total = Fraction(0)
        for i in range(r + 1):
            j = r - i
            term = (p_gamma[r] * p_B[r]) / (p_A[i] * p_B[j] * fact[i] * fact[j])
            total += term * u_pow[i] * v_pow[j]
        coeffs.append(total)
    return coeffs


def rhs_t_coeffs(
    alpha: Fraction,
    beta: Fraction,
    gamma: Fraction,
    x: Fraction,
    order: int,
    *,
    error_mode: str = "none",
):
    """
    Exact rational coefficients [t^r] (r = 0..order) of the literal RHS in Cohl's image:
      (1-t)^(-gamma) * sum_{m,n>=0} C(m,n) * (2 t (x-1) / (t-1)^2)^m * (-2 t^2 (x-1) / (t-1)^2)^n
    Note (t-1)^2 = (1-t)^2, so each (m,n) summand contributes
      C(m,n) * (2(x-1))^(m+n) * (-1)^n * t^(m+2n) * (1-t)^(-gamma - 2(m+n)).
    Using (1-t)^(-p) = sum_{s>=0} (p)_s / s! * t^s, only m + 2n <= order contribute.

    Optional error_mode plants deliberate bugs to test fault detection:
      - "sign_flip": drops the minus sign in the second base (-2t^2(x-1) -> +2t^2(x-1))
      - "denom_m_plus_n": replaces (alpha+beta+1)_{m+2n} by (alpha+beta+1)_{m+n}
      - "shift_param": replaces (alpha+beta-gamma+1)_n by (alpha+beta-gamma)_n
      - "half_param_shift": replaces ((alpha+beta+2)/2)_{m+n} by ((alpha+beta+3)/2)_{m+n}
    """
    A = alpha + 1
    S = alpha + beta + 1
    half1 = S / 2
    half2 = (S + 1) / 2 if error_mode != "half_param_shift" else (S + 2) / 2
    c_n_param = (S - gamma) if error_mode != "shift_param" else (S - gamma - 1)
    base_x = 2 * (x - 1)

    fact = [Fraction(math.factorial(k)) for k in range(order + 1)]
    coeffs = [Fraction(0) for _ in range(order + 1)]

    for n in range(order // 2 + 1):
        for m in range(order - 2 * n + 1):
            k = m + n
            deg_t = m + 2 * n
            denom_idx = (m + 2 * n) if error_mode != "denom_m_plus_n" else (m + n)

            num = poch(c_n_param, n) * poch(gamma, k) * poch(half1, k) * poch(half2, k)
            den = poch(A, k) * poch(S, denom_idx) * fact[m] * fact[n]
            sign = 1 if (error_mode == "sign_flip" or n % 2 == 0) else -1
            pref = (num / den) * (base_x ** k) * sign

            # Multiply by (1-t)^(-gamma - 2k) = sum_{s>=0} (gamma + 2k)_s / s! * t^s
            p = gamma + 2 * k
            binom_term = Fraction(1)
            for s in range(order - deg_t + 1):
                if s > 0:
                    binom_term *= (p + s - 1) / s
                coeffs[deg_t + s] += pref * binom_term

    return coeffs


# --- Polynomial arithmetic in x over Fraction for symbolic-x verification ---

def poly_add(p, q):
    n = max(len(p), len(q))
    out = [Fraction(0)] * n
    for i, c in enumerate(p):
        out[i] += c
    for i, c in enumerate(q):
        out[i] += c
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_mul(p, q):
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0:
            continue
        for j, b in enumerate(q):
            out[i + j] += a * b
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_scale(p, c: Fraction):
    if c == 0:
        return [Fraction(0)]
    return [a * c for a in p]


def poly_pow(p, e: int):
    res = [Fraction(1)]
    base = list(p)
    while e > 0:
        if e & 1:
            res = poly_mul(res, base)
        base = poly_mul(base, base)
        e >>= 1
    return res


def check_polynomial_in_x(alpha: Fraction, beta: Fraction, gamma: Fraction, order: int):
    """Compare [t^r] LHS and [t^r] RHS as exact polynomials in x over Q."""
    A = alpha + 1
    B = beta + 1
    S = alpha + beta + 1
    fact = [Fraction(math.factorial(k)) for k in range(order + 1)]

    # u_poly = (x-1)/2, v_poly = (x+1)/2
    u_poly = [Fraction(-1, 2), Fraction(1, 2)]
    v_poly = [Fraction(1, 2), Fraction(1, 2)]
    u_pows = [poly_pow(u_poly, k) for k in range(order + 1)]
    v_pows = [poly_pow(v_poly, k) for k in range(order + 1)]

    lhs_polys = []
    for r in range(order + 1):
        acc = [Fraction(0)]
        for i in range(r + 1):
            j = r - i
            c = (poch(gamma, r) * poch(B, r)) / (poch(A, i) * poch(B, j) * fact[i] * fact[j])
            acc = poly_add(acc, poly_scale(poly_mul(u_pows[i], v_pows[j]), c))
        lhs_polys.append(acc)

    # RHS polynomial expansion in x: 2(x-1) = [-2, 2]
    bx_poly = [Fraction(-2), Fraction(2)]
    bx_pows = [poly_pow(bx_poly, k) for k in range(order + 1)]
    rhs_polys = [[Fraction(0)] for _ in range(order + 1)]

    for n in range(order // 2 + 1):
        for m in range(order - 2 * n + 1):
            k = m + n
            deg_t = m + 2 * n
            num = poch(S - gamma, n) * poch(gamma, k) * poch(S / 2, k) * poch((S + 1) / 2, k)
            den = poch(A, k) * poch(S, m + 2 * n) * fact[m] * fact[n]
            sign = 1 if n % 2 == 0 else -1
            term_poly = poly_scale(bx_pows[k], (num / den) * sign)

            p = gamma + 2 * k
            binom_term = Fraction(1)
            for s in range(order - deg_t + 1):
                if s > 0:
                    binom_term *= (p + s - 1) / s
                rhs_polys[deg_t + s] = poly_add(rhs_polys[deg_t + s], poly_scale(term_poly, binom_term))

    mismatches = []
    for r in range(order + 1):
        if lhs_polys[r] != rhs_polys[r]:
            mismatches.append(r)
    return mismatches


def h3_t_coeffs(alpha: Fraction, beta: Fraction, gamma: Fraction, x: Fraction, order: int):
    """
    Exact rational coefficients [t^r] (r = 0..order) of the Horn H3 form in proof.md lines 360-363:
      (1-t)^(-gamma) * H3(gamma, S-gamma; A; w, (1-t)w)
    where w = t(x-1)/(2(1-t)^2) and H3(a,b;c;z1,z2) = sum_{p,q>=0} (a)_{2p+q} (b)_q / ((c)_{p+q} p! q!) z1^p z2^q.
    Each (p,q) term is:
      [(gamma)_{2p+q} (S-gamma)_q / ((A)_{p+q} p! q!)] * ((x-1)/2)^(p+q) * t^(p+q) * (1-t)^(-gamma - 2p - q).
    """
    A = alpha + 1
    S = alpha + beta + 1
    u_coeff = (x - 1) / 2
    fact = [Fraction(math.factorial(k)) for k in range(order + 1)]
    coeffs = [Fraction(0) for _ in range(order + 1)]

    for p in range(order + 1):
        for q in range(order - p + 1):
            k = p + q
            c = (poch(gamma, 2 * p + q) * poch(S - gamma, q)) / (poch(A, k) * fact[p] * fact[q])
            pref = c * (u_coeff ** k)
            exp_pow = gamma + 2 * p + q
            binom_term = Fraction(1)
            for s in range(order - k + 1):
                if s > 0:
                    binom_term *= (exp_pow + s - 1) / s
                coeffs[k + s] += pref * binom_term
    return coeffs


def main():
    print("=== EXACT RATIONAL POWER-SERIES VERIFICATION (COHL F4 IDENTITY) ===")
    cases = [
        ("Generic positive rationals", Fraction(1, 3), Fraction(2, 5), Fraction(7, 4), Fraction(3, 7), 20),
        ("Awkward negative non-integers", Fraction(-7, 3), Fraction(-5, 2), Fraction(-11, 5), Fraction(-2, 3), 20),
        ("Negative non-integers, x = +1", Fraction(-3, 2), Fraction(-7, 4), Fraction(-1, 3), Fraction(1), 20),
        ("Negative non-integers, x = -1", Fraction(-3, 2), Fraction(-7, 4), Fraction(-1, 3), Fraction(-1), 20),
        ("Negative non-integers, x = 0", Fraction(-3, 2), Fraction(-7, 4), Fraction(-1, 3), Fraction(0), 20),
        ("Large negative non-integers, x = 5", Fraction(-17, 5), Fraction(-13, 3), Fraction(-9, 2), Fraction(5), 20),
        ("Terminating gamma = -5, x = 0", Fraction(-5, 3), Fraction(-11, 4), Fraction(-5), Fraction(0), 20),
        ("Terminating gamma = -5, x = -1", Fraction(-5, 3), Fraction(-11, 4), Fraction(-5), Fraction(-1), 20),
        ("Terminating gamma = -5, x = 13/4", Fraction(-5, 3), Fraction(-11, 4), Fraction(-5), Fraction(13, 4), 20),
        ("Specialization gamma = alpha+beta+1", Fraction(-4, 3), Fraction(-7, 5), Fraction(-26, 15), Fraction(2, 5), 20),
        ("Near denominator poles", Fraction(-1999, 1000), Fraction(-2990, 997), Fraction(-1, 7), Fraction(1, 2), 18),
    ]

    all_ok = True
    for name, a, b, g, x, order in cases:
        lc = lhs_t_coeffs(a, b, g, x, order)
        rc = rhs_t_coeffs(a, b, g, x, order)
        diffs = [(r, lc[r] - rc[r]) for r in range(order + 1) if lc[r] != rc[r]]
        status = "PASS" if not diffs else f"FAIL ({len(diffs)} mismatches)"
        if diffs:
            all_ok = False
        extra = ""
        if g == Fraction(-5):
            tail_nonzero = [r for r in range(6, order + 1) if lc[r] != 0 or rc[r] != 0]
            extra = f" | tail(r>5) zero: {len(tail_nonzero) == 0}"
            if tail_nonzero:
                all_ok = False
        print(f"[{status}] {name:38s} (a={a}, b={b}, g={g}, x={x}, ord={order}){extra}")

    print("\n=== EXACT POLYNOMIAL-IN-x CHECK (ALL x SIMULTANEOUSLY) ===")
    mismatches = check_polynomial_in_x(Fraction(-7, 3), Fraction(-5, 2), Fraction(-11, 5), 16)
    poly_status = "PASS" if not mismatches else f"FAIL at orders {mismatches}"
    if mismatches:
        all_ok = False
    print(f"[{poly_status}] Symbolic polynomial in x for a=-7/3, b=-5/2, g=-11/5 through t^16 (153 (r,d) coeffs)")

    print("\n=== EXACT HORN H3 REDUCTION CHECK (proof.md eq. lines 360-363) ===")
    a, b, g, x = Fraction(-7, 3), Fraction(-5, 2), Fraction(-11, 5), Fraction(-2, 3)
    lc = lhs_t_coeffs(a, b, g, x, 16)
    hc = h3_t_coeffs(a, b, g, x, 16)
    h3_diffs = [r for r in range(17) if lc[r] != hc[r]]
    h3_status = "PASS" if not h3_diffs else f"FAIL at orders {h3_diffs}"
    if h3_diffs:
        all_ok = False
    print(f"[{h3_status}] F4 == (1-t)^(-gamma) H3(gamma, S-gamma; A; w, (1-t)w) through t^16")

    print("\n=== PLANTED DELIBERATE ERROR DETECTION TESTS ===")
    pa, pb, pg, px = Fraction(-7, 3), Fraction(-5, 2), Fraction(-11, 5), Fraction(-2, 3)
    ref_lhs = lhs_t_coeffs(pa, pb, pg, px, 10)
    planted = [
        ("sign_flip", "Drop minus sign in second base: -2t^2(x-1) -> +2t^2(x-1)"),
        ("denom_m_plus_n", "Denominator subscript typo: (alpha+beta+1)_{m+2n} -> (alpha+beta+1)_{m+n}"),
        ("shift_param", "Parameter shift typo: (alpha+beta-gamma+1)_n -> (alpha+beta-gamma)_n"),
        ("half_param_shift", "Half-parameter typo: ((alpha+beta+2)/2)_{m+n} -> ((alpha+beta+3)/2)_{m+n}"),
    ]
    all_caught = True
    for mode, desc in planted:
        bad_rhs = rhs_t_coeffs(pa, pb, pg, px, 10, error_mode=mode)
        first_bad = next(((r, ref_lhs[r] - bad_rhs[r]) for r in range(11) if ref_lhs[r] != bad_rhs[r]), None)
        if first_bad is None:
            all_caught = False
            print(f"[MISSED] {mode}: {desc}")
        else:
            r, diff = first_bad
            print(f"[CAUGHT at t^{r}] {mode:16s} | {desc}")
            print(f"                 exact residual LHS - RHS = {diff}")

    print("\n=== SUMMARY ===")
    print(f"All exact identity checks passed: {all_ok}")
    print(f"All {len(planted)} planted errors caught:     {all_caught}")
    if not (all_ok and all_caught):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
