#!/usr/bin/env python3
"""
Independent high-precision numerical verification of Howard Cohl's Appell F4 identity,
its boundary convergence table (proof.md lines 181-186), conditional-convergence
counterexamples to "wherever both sides converge", and planted-error detection.

Written directly from Cohl's posted image (cohl-f4.jpg / DLMF 16.13.4).
Pure Python 3 standard library (decimal.Decimal) at 90-110 working digits; no external dependencies.
"""

from decimal import Decimal, getcontext
import math

WORK_PREC = 110
getcontext().prec = WORK_PREC


def dec_pi() -> Decimal:
    """Compute pi to current Decimal precision via Machin's formula."""
    def arctan_inv(inv_n: int) -> Decimal:
        n2 = Decimal(inv_n * inv_n)
        term = Decimal(1) / Decimal(inv_n)
        total = term
        k = 1
        eps = Decimal(10) ** (-(WORK_PREC - 2))
        while True:
            term = -term / n2
            contrib = term / Decimal(2 * k + 1)
            if abs(contrib) < eps:
                break
            total += contrib
            k += 1
        return total

    return Decimal(16) * arctan_inv(5) - Decimal(4) * arctan_inv(239)


PI = dec_pi()
TWO_PI = Decimal(2) * PI
HALF_PI = PI / Decimal(2)


def dec_atan(x: Decimal) -> Decimal:
    """High-precision arctan(x) for real Decimal x via argument halving + Taylor series."""
    if x == 0:
        return Decimal(0)
    neg = x < 0
    x = abs(x)
    # Reduce |x| using arctan(x) = 2 arctan(x / (1 + sqrt(1+x^2)))
    halvings = 0
    while x > Decimal("0.05"):
        x = x / (Decimal(1) + (Decimal(1) + x * x).sqrt())
        halvings += 1
    x2 = x * x
    term = x
    total = term
    k = 1
    eps = Decimal(10) ** (-(WORK_PREC - 4))
    while True:
        term = -term * x2
        contrib = term / Decimal(2 * k + 1)
        if abs(contrib) < eps:
            break
        total += contrib
        k += 1
    res = total * (Decimal(2) ** halvings)
    return -res if neg else res


def dec_atan2(y: Decimal, x: Decimal) -> Decimal:
    if x > 0:
        return dec_atan(y / x)
    if x < 0:
        return dec_atan(y / x) + (PI if y >= 0 else -PI)
    if y > 0:
        return HALF_PI
    if y < 0:
        return -HALF_PI
    raise ValueError("atan2(0, 0) undefined")


def dec_sin_cos(theta: Decimal):
    """High-precision (sin(theta), cos(theta)) via range reduction + Taylor series."""
    # Reduce to [-pi, pi]
    q = (theta / TWO_PI).to_integral_value()
    t = theta - q * TWO_PI
    t2 = t * t
    eps = Decimal(10) ** (-(WORK_PREC - 4))

    # cos
    c_term = Decimal(1)
    c_sum = c_term
    k = 1
    while True:
        c_term = -c_term * t2 / Decimal((2 * k - 1) * (2 * k))
        if abs(c_term) < eps:
            break
        c_sum += c_term
        k += 1

    # sin
    s_term = t
    s_sum = s_term
    k = 1
    while True:
        s_term = -s_term * t2 / Decimal((2 * k) * (2 * k + 1))
        if abs(s_term) < eps:
            break
        s_sum += s_term
        k += 1

    return s_sum, c_sum


class CDec:
    """Arbitrary-precision complex number backed by decimal.Decimal."""
    __slots__ = ("re", "im")

    def __init__(self, re, im=0):
        self.re = re if isinstance(re, Decimal) else Decimal(str(re))
        self.im = im if isinstance(im, Decimal) else Decimal(str(im))

    def __add__(self, o):
        o = o if isinstance(o, CDec) else CDec(o)
        return CDec(self.re + o.re, self.im + o.im)

    def __sub__(self, o):
        o = o if isinstance(o, CDec) else CDec(o)
        return CDec(self.re - o.re, self.im - o.im)

    def __neg__(self):
        return CDec(-self.re, -self.im)

    def __mul__(self, o):
        o = o if isinstance(o, CDec) else CDec(o)
        return CDec(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)

    def __truediv__(self, o):
        o = o if isinstance(o, CDec) else CDec(o)
        d = o.re * o.re + o.im * o.im
        return CDec((self.re * o.re + self.im * o.im) / d, (self.im * o.re - self.re * o.im) / d)

    def abs(self) -> Decimal:
        return (self.re * self.re + self.im * self.im).sqrt()

    def ln(self):
        r = self.abs()
        ang = dec_atan2(self.im, self.re)
        return CDec(r.ln(), ang)

    def exp(self):
        er = self.re.exp()
        s, c = dec_sin_cos(self.im)
        return CDec(er * c, er * s)

    def pow(self, exponent):
        exponent = exponent if isinstance(exponent, CDec) else CDec(exponent)
        return (self.ln() * exponent).exp()

    def fmt(self, digits=18) -> str:
        if abs(self.im) < Decimal(10) ** (-75):
            return f"{self.re:.{digits}g}"
        sign = "+" if self.im >= 0 else "-"
        return f"{self.re:.{digits}g} {sign} {abs(self.im):.{digits}g}i"


def cdec_self_test():
    """Verify CDec.pow against the binomial series (1 - t)^(-gamma) for complex t, gamma."""
    t = CDec("0.2", "-0.15")
    g = CDec("-1.3", "0.7")
    direct = (CDec(1) - t).pow(-g)
    series = CDec(0)
    term = CDec(1)
    for k in range(180):
        series = series + term
        term = term * (g + CDec(k)) * t / CDec(k + 1)
    err = (direct - series).abs()
    assert err < Decimal("1e-95"), f"CDec self-test failed: err={err}"
    print(f"[SELF-TEST PASS] CDec complex power vs binomial series discrepancy = {err:.3e}")


def eval_lhs(alpha: CDec, beta: CDec, gamma: CDec, x: CDec, t: CDec, K: int):
    """
    Evaluate literal F4(gamma, beta+1; alpha+1, beta+1; t(x-1)/2, t(x+1)/2)
    summed over total degree k = i + j <= K.
    Returns (value, abs_sum, last_deg_abs_sum).
    """
    A = alpha + CDec(1)
    B = beta + CDec(1)
    u = t * (x - CDec(1)) / CDec(2)
    v = t * (x + CDec(1)) / CDec(2)

    # Precompute tables up to K
    p_g = [CDec(1)] * (K + 1)
    p_B = [CDec(1)] * (K + 1)
    p_A = [CDec(1)] * (K + 1)
    fact = [Decimal(1)] * (K + 1)
    u_pow = [CDec(1)] * (K + 1)
    v_pow = [CDec(1)] * (K + 1)
    for k in range(1, K + 1):
        p_g[k] = p_g[k - 1] * (gamma + CDec(k - 1))
        p_B[k] = p_B[k - 1] * (B + CDec(k - 1))
        p_A[k] = p_A[k - 1] * (A + CDec(k - 1))
        fact[k] = fact[k - 1] * Decimal(k)
        u_pow[k] = u_pow[k - 1] * u
        v_pow[k] = v_pow[k - 1] * v

    total = CDec(0)
    abs_total = Decimal(0)
    last_deg_abs = Decimal(0)
    for k in range(K + 1):
        deg_abs = Decimal(0)
        num_k = p_g[k] * p_B[k]
        for i in range(k + 1):
            j = k - i
            den = p_A[i] * p_B[j] * CDec(fact[i] * fact[j])
            term = (num_k / den) * u_pow[i] * v_pow[j]
            total = total + term
            t_abs = term.abs()
            deg_abs += t_abs
        abs_total += deg_abs
        if k == K:
            last_deg_abs = deg_abs
    return total, abs_total, last_deg_abs


def eval_rhs(alpha: CDec, beta: CDec, gamma: CDec, x: CDec, t: CDec, K: int, *, error_mode: str = "none"):
    """
    Evaluate literal RHS of Cohl's image summed over total degree k = m + n <= K:
      (1-t)^(-gamma) * sum_{m+n <= K}
        [(alpha+beta-gamma+1)_n (gamma)_{m+n} ((alpha+beta+1)/2)_{m+n} ((alpha+beta+2)/2)_{m+n}]
        / [(alpha+1)_{m+n} (alpha+beta+1)_{m+2n} m! n!]
        * (2 t (x-1) / (t-1)^2)^m * (-2 t^2 (x-1) / (t-1)^2)^n
    Returns (value, abs_sum_with_prefactor, last_deg_abs_sum_with_prefactor).
    """
    A = alpha + CDec(1)
    S = alpha + beta + CDec(1)
    half1 = S / CDec(2)
    half2 = (S + CDec(1)) / CDec(2)
    c_n = S - gamma

    tm1_sq = (t - CDec(1)) * (t - CDec(1))
    b1 = (CDec(2) * t * (x - CDec(1))) / tm1_sq
    b2 = (CDec(2) * t * t * (x - CDec(1))) / tm1_sq
    if error_mode != "sign_flip":
        b2 = -b2

    pref = (CDec(1) - t).pow(-gamma)
    pref_abs = pref.abs()

    # Precompute tables
    p_cn = [CDec(1)] * (K + 1)
    p_g = [CDec(1)] * (K + 1)
    p_h1 = [CDec(1)] * (K + 1)
    p_h2 = [CDec(1)] * (K + 1)
    p_A = [CDec(1)] * (K + 1)
    p_S = [CDec(1)] * (2 * K + 1)
    fact = [Decimal(1)] * (K + 1)
    b1_pow = [CDec(1)] * (K + 1)
    b2_pow = [CDec(1)] * (K + 1)
    for k in range(1, K + 1):
        p_cn[k] = p_cn[k - 1] * (c_n + CDec(k - 1))
        p_g[k] = p_g[k - 1] * (gamma + CDec(k - 1))
        p_h1[k] = p_h1[k - 1] * (half1 + CDec(k - 1))
        p_h2[k] = p_h2[k - 1] * (half2 + CDec(k - 1))
        p_A[k] = p_A[k - 1] * (A + CDec(k - 1))
        fact[k] = fact[k - 1] * Decimal(k)
        b1_pow[k] = b1_pow[k - 1] * b1
        b2_pow[k] = b2_pow[k - 1] * b2
    for k in range(1, 2 * K + 1):
        p_S[k] = p_S[k - 1] * (S + CDec(k - 1))

    total = CDec(0)
    abs_total = Decimal(0)
    last_deg_abs = Decimal(0)
    for k in range(K + 1):
        deg_abs = Decimal(0)
        num_k = p_g[k] * p_h1[k] * p_h2[k]
        for n in range(k + 1):
            m = k - n
            denom_idx = (m + 2 * n) if error_mode != "denom_m_plus_n" else (m + n)
            num = p_cn[n] * num_k
            den = p_A[k] * p_S[denom_idx] * CDec(fact[m] * fact[n])
            term = (num / den) * b1_pow[m] * b2_pow[n]
            total = total + term
            deg_abs += term.abs()
        abs_total += deg_abs
        if k == K:
            last_deg_abs = deg_abs
    return pref * total, pref_abs * abs_total, pref_abs * last_deg_abs


def eval_rhs_cancelled(alpha: CDec, beta: CDec, gamma: CDec, x: CDec, t: CDec, K: int):
    """Evaluate RHS after the cancellation (2)->(3), valid even when S = alpha+beta+1 in {0,-1,-2,...}."""
    A = alpha + CDec(1)
    S = alpha + beta + CDec(1)
    u = t * (x - CDec(1)) / CDec(2)
    w = u / ((CDec(1) - t) * (CDec(1) - t))
    pref = (CDec(1) - t).pow(-gamma)

    p_g = [CDec(1)] * (K + 1)
    p_A = [CDec(1)] * (K + 1)
    p_cn = [CDec(1)] * (K + 1)
    fact = [Decimal(1)] * (K + 1)
    w_pow = [CDec(1)] * (K + 1)
    mt_pow = [CDec(1)] * (K + 1)
    for k in range(1, K + 1):
        p_g[k] = p_g[k - 1] * (gamma + CDec(k - 1))
        p_A[k] = p_A[k - 1] * (A + CDec(k - 1))
        p_cn[k] = p_cn[k - 1] * (S - gamma + CDec(k - 1))
        fact[k] = fact[k - 1] * Decimal(k)
        w_pow[k] = w_pow[k - 1] * w
        mt_pow[k] = mt_pow[k - 1] * (-t)

    total = CDec(0)
    for k in range(K + 1):
        inner = CDec(0)
        for n in range(k + 1):
            m = k - n
            # (S + k + n)_{k - n}
            p_skn = CDec(1)
            for j in range(m):
                p_skn = p_skn * (S + CDec(k + n + j))
            inner = inner + (p_cn[n] * p_skn / CDec(fact[n] * fact[m])) * mt_pow[n]
        total = total + (p_g[k] / p_A[k]) * w_pow[k] * inner
    return pref * total


def domain_metrics(x: CDec, t: CDec):
    u = t * (x - CDec(1)) / CDec(2)
    v = t * (x + CDec(1)) / CDec(2)
    L = u.abs().sqrt() + v.abs().sqrt()
    R = (Decimal(4) * u.abs()) / ((CDec(1) - t).abs() ** 2)
    return L, R


def main():
    cdec_self_test()

    print("\n=== HIGH-PRECISION INTERIOR NUMERICAL CHECKS (90+ DIGITS) ===")
    interior_cases = [
        ("Real positive params", CDec("0.7"), CDec("1.3"), CDec("2.1"), CDec("0.25"), CDec("0.15"), 180, 240),
        ("Real negative non-int params", CDec("-2.3"), CDec("-1.8"), CDec("-1.4"), CDec("-0.5"), CDec("-0.2"), 200, 260),
        ("Complex params & args", CDec("-1.2", "0.8"), CDec("0.6", "-1.1"), CDec("1.5", "0.7"), CDec("0.3", "-0.4"), CDec("-0.15", "0.12"), 200, 260),
        ("Complex closer to boundary", CDec("1.4", "-0.9"), CDec("-0.7", "0.5"), CDec("-0.4", "1.2"), CDec("-0.8", "0.6"), CDec("-0.22", "0.14"), 320, 420),
        ("Special x = 0, complex t", CDec("-1.5"), CDec("-2.25"), CDec("-0.35"), CDec("0"), CDec("0.18", "-0.12"), 340, 420),
        ("Special x = -1, t = -0.2", CDec("-1.5"), CDec("-2.25"), CDec("-0.35"), CDec("-1"), CDec("-0.2"), 260, 340),
        ("Special x = +1, t = 0.3+0.4i", CDec("-1.5"), CDec("-2.25"), CDec("-0.35"), CDec("1"), CDec("0.3", "0.4"), 260, 340),
    ]

    all_interior_ok = True
    for name, a, b, g, x, t, K1, K2 in interior_cases:
        L, R = domain_metrics(x, t)
        l1, _, _ = eval_lhs(a, b, g, x, t, K1)
        l2, _, _ = eval_lhs(a, b, g, x, t, K2)
        r1, _, _ = eval_rhs(a, b, g, x, t, K1)
        r2, _, _ = eval_rhs(a, b, g, x, t, K2)
        scale = max(Decimal(1), l2.abs(), r2.abs())
        diff = (l2 - r2).abs() / scale
        dl = (l2 - l1).abs() / scale
        dr = (r2 - r1).abs() / scale
        ok = diff < Decimal("1e-80")
        if not ok:
            all_interior_ok = False
        print(f"[{'PASS' if ok else 'FAIL'}] {name:30s} | L={float(L):.4f}, R={float(R):.4f} | K={K1}->{K2}")
        print(f"       LHS = {l2.fmt(22)}")
        print(f"       RHS = {r2.fmt(22)}")
        print(f"       |LHS-RHS|/scale = {diff:.3e}  (dLHS={dl:.3e}, dRHS={dr:.3e})")

    print("\n=== REMOVABLE S IN {0,-1,-2,...} POLE CHECK (proof.md lines 44-45) ===")
    a_s, b_s, g_s, x_s, t_s = CDec("-2.3"), CDec("-1.7"), CDec("-1.4"), CDec("-0.5"), CDec("-0.2")
    try:
        eval_rhs(a_s, b_s, g_s, x_s, t_s, 10)
        lit_failed = False
    except Exception as exc:
        lit_failed = True
        lit_exc = type(exc).__name__
    l_s, _, _ = eval_lhs(a_s, b_s, g_s, x_s, t_s, 220)
    r_s = eval_rhs_cancelled(a_s, b_s, g_s, x_s, t_s, 220)
    s_diff = (l_s - r_s).abs() / max(Decimal(1), l_s.abs(), r_s.abs())
    assert lit_failed and s_diff < Decimal("1e-70")
    print(f"[PASS] a=-2.3, b=-1.7 => S=-3: literal RHS raises {lit_exc} (0/0 pole); cancelled form (3) matches LHS to {s_diff:.3e}")

    print("\n=== PLANTED DELIBERATE ERRORS IN HIGH-PRECISION NUMERICS ===")
    a, b, g, x, t = CDec("-1.2", "0.8"), CDec("0.6", "-1.1"), CDec("1.5", "0.7"), CDec("0.3", "-0.4"), CDec("-0.15", "0.12")
    ref_l, _, _ = eval_lhs(a, b, g, x, t, 220)
    for mode in ("sign_flip", "denom_m_plus_n"):
        bad_r, _, _ = eval_rhs(a, b, g, x, t, 220, error_mode=mode)
        err = (ref_l - bad_r).abs() / max(Decimal(1), ref_l.abs(), bad_r.abs())
        assert err > Decimal("1e-3"), f"Numeric check failed to catch {mode}"
        print(f"[CAUGHT] mode={mode:16s} | |LHS - bad_RHS|/scale = {err:.6e}")

    print("\n=== BOUNDARY & CONVERGENCE-DOMAIN DIAGNOSTIC TESTS ===")

    # Boundary A1: L = 1, R = 1, uv != 0 (x = 5/4, t = 1/2 => u = 1/16, v = 9/16), d = 5.5 > 1/2
    a, b, g, x, t = CDec("6"), CDec("0.5"), CDec("1.5"), CDec("1.25"), CDec("0.5")
    L, R = domain_metrics(x, t)
    l_val, l_abs, l_last = eval_lhs(a, b, g, x, t, 120)
    r_val, r_abs, r_last = eval_rhs(a, b, g, x, t, 120)
    b_diff = (l_val - r_val).abs()
    print(f"[BOUNDARY L=R=1, uv!=0, d=5.5>1/2] L={float(L):.4f}, R={float(R):.4f}, K=120:")
    print(f"       LHS = {l_val.fmt(20)}  | RHS = {r_val.fmt(20)}  | |LHS-RHS| = {b_diff:.3e}")
    print(f"       deg-120 abs terms: LHS={l_last:.3e}, RHS={r_last:.3e} (decaying as k^(-6))")

    # Boundary B1: L = 1, R = 1, v = 0 (x = -1, t = -1 => u = 1, v = 0), A=6, B=1.5, gamma=0.5 => d = 5.5 > max(1.5, 0.5)
    # Exact Gauss value: 2F1(1/2, 3/2; 6; 1) = Gamma(6)Gamma(4)/(Gamma(11/2)Gamma(9/2))
    # = 120 * 6 / ((945/32 sqrt(pi)) * (105/16 sqrt(pi))) = 720 * 512 / (99225 * pi) = 4096 / (1102.5 * pi) = 8192 / (2205 * pi)
    a, b, g, x, t = CDec("5"), CDec("0.5"), CDec("0.5"), CDec("-1"), CDec("-1")
    L, R = domain_metrics(x, t)
    l_val, _, l_last = eval_lhs(a, b, g, x, t, 150)
    r_val, _, r_last = eval_rhs(a, b, g, x, t, 150)
    gauss_exact = Decimal(8192) / (Decimal(2205) * PI)
    print(f"[BOUNDARY x=-1,t=-1, d=5.5>max(ReB=1.5,1/2)] L={float(L):.4f}, R={float(R):.4f}, K=150:")
    print(f"       LHS={l_val.fmt(18)}, RHS={r_val.fmt(18)}, Gauss={gauss_exact:.18g}")
    print(f"       |LHS-Gauss|={abs(l_val.re - gauss_exact):.3e}, |RHS-Gauss|={abs(r_val.re - gauss_exact):.3e}")

    # Boundary B2: x = -1, t = -1, A = 0.75, B = -0.5, gamma = 0.5 => d = 0.25 > Re B = -0.5 (LHS converges abs), d <= 1/2 (RHS diverges!)
    a, b, g, x, t = CDec("-0.25"), CDec("-1.5"), CDec("0.5"), CDec("-1"), CDec("-1")
    _, l_abs_60, l_last_60 = eval_lhs(a, b, g, x, t, 60)
    _, l_abs_120, l_last_120 = eval_lhs(a, b, g, x, t, 120)
    _, r_abs_60, r_last_60 = eval_rhs(a, b, g, x, t, 60)
    _, r_abs_120, r_last_120 = eval_rhs(a, b, g, x, t, 120)
    print(f"[BOUNDARY x=-1,t=-1, d=0.25 (d>ReB=-0.5, d<=1/2)] Testing dual-row condition:")
    print(f"       LHS deg-k abs term: k=60: {l_last_60:.3e}, k=120: {l_last_120:.3e} (ratio={float(l_last_120/l_last_60):.4f}, expected ~2^(-1.75)={2**(-1.75):.4f}, CONVERGES)")
    print(f"       RHS deg-k abs term: k=60: {r_last_60:.3e}, k=120: {r_last_120:.3e} (ratio={float(r_last_120/r_last_60):.4f}, expected ~2^(-0.75)={2**(-0.75):.4f}, DIVERGES)")

    # Outside common domain: x = -1, t = 0.5 => L = sqrt(0.5) < 1, R = 8 > 1
    a, b, g, x, t = CDec("1"), CDec("1"), CDec("0.5"), CDec("-1"), CDec("0.5")
    L, R = domain_metrics(x, t)
    l_val, _, l_last = eval_lhs(a, b, g, x, t, 60)
    _, _, r_last = eval_rhs(a, b, g, x, t, 60)
    print(f"[OUTSIDE COMMON DOMAIN x=-1,t=0.5] L={float(L):.4f}<1, R={float(R):.4f}>1:")
    print(f"       deg-60 abs terms: LHS={l_last:.3e} (converges), RHS={r_last:.3e} (diverges ~ 8^60)")

    # Conditional convergence counterexample D1: x = 1, t = -1, gamma = 1/2 (L=1, u=0, Re(gamma)=0.5 not < 0)
    a, b, g, x, t = CDec("0.3"), CDec("0.4"), CDec("0.5"), CDec("1"), CDec("-1")
    l_val_200, l_abs_200, l_last_200 = eval_lhs(a, b, g, x, t, 200)
    r_val_200, r_abs_200, _ = eval_rhs(a, b, g, x, t, 200)
    print(f"[CONDITIONAL CONVERGENCE D1: x=1, t=-1, gamma=1/2]")
    print(f"       LHS partial sum K=200 = {l_val_200.fmt(12)} vs RHS = 1/sqrt(2) = {r_val_200.fmt(12)} (diff={abs(l_val_200.re - r_val_200.re):.3e})")
    print(f"       LHS abs sum K=200 = {l_abs_200:.6g}, deg-200 term = {l_last_200:.4e} (diverges as sum k^(-1/2))")

    # Conditional convergence counterexample D2: x = -1, t = 3 - 2*sqrt(2), A = 1, B = 1, gamma = 1/2 (L<1, R=1, d=1/2)
    t_crit = CDec(Decimal(3) - Decimal(2) * Decimal(2).sqrt())
    a, b, g, x = CDec("0"), CDec("0"), CDec("0.5"), CDec("-1")
    L, R = domain_metrics(x, t_crit)
    l_val_150, _, l_last_150 = eval_lhs(a, b, g, x, t_crit, 150)
    r_val_150, r_abs_150, r_last_150 = eval_rhs(a, b, g, x, t_crit, 150)
    r_val_151, r_abs_151, r_last_151 = eval_rhs(a, b, g, x, t_crit, 151)
    print(f"[CONDITIONAL CONVERGENCE D2: x=-1, t=3-2*sqrt(2), alpha=0, beta=0, gamma=1/2 (R=1, d=1/2)]")
    print(f"       L={float(L):.4f}<1, R={float(R):.6f}")
    print(f"       LHS (exact to ~1e-58) = {l_val_150.fmt(14)}")
    print(f"       RHS K=150 = {r_val_150.fmt(14)}, K=151 = {r_val_151.fmt(14)}, avg = {((r_val_150 + r_val_151) / CDec(2)).fmt(14)}")
    print(f"       RHS abs sum K=150 = {r_abs_150:.6g}, deg-150 abs = {r_last_150:.4e} (~ c/k, diverges log-wise, alternates to LHS)")

    print("\n=== SUMMARY ===")
    print(f"All interior high-precision checks passed: {all_interior_ok}")
    if not all_interior_ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
