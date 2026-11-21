*This review is republished as written, lightly redacted for the
repository-internal file paths of the working tree it reviewed.*

VERDICT: CORRECT
NOVELTY: FOLLOWS FROM DLMF 16.16.5_5 + 16.16.8 combined with the standard quadratic Horn F1→H3 reduction (Ernst–Karlsson eq. (42); print source not openable by me)

Independent review of `proof.md` (commit 6f70adb), by a reviewer from a
different model family. I did not read `check_numeric.py` or `check_symbolic.py` before
writing my own program; my script is `checks/verify_glm.py` (output in
`checks/glm-output.txt`), written from equation (1) and the DLMF 16.13.4 convention as
transcribed in the report only.

- Statement section: correct, with one unavoidable caveat about the image (below).
- Proof (eqs. (2)–(6)), "Passage to convergent series", "Absolute convergence, including
  the boundary": every step checked by hand and by exact computation; **no error found**.
- Recomputation: full exact agreement, coefficientwise in two variables through degree 12
  and numerically to ~1e-60 at points inside the convergence region.
- Novelty: the identity is the published repeated-parameter F4→F1 reduction plus a
  classical quadratic F1→H3 reduction, re-collected by the duplication formula and one
  Vandermonde convolution. The specific double series (1) is not displayed in any source
  I could find, but it is a known result in disguise.

## 1. Statement

**Image.** I cannot view images (no vision capability in this environment, no OCR tool, no
`pip`/`apt` access to install one). So per the task instructions I checked the
transcription against the report's text only, plus everything checkable around the image:

- Downloaded <https://pbs.twimg.com/media/F5iLOskbAAA0o8A.jpg?name=orig> (opened): a
  52,711-byte progressive JPEG, 2048×362 — same dimensions the report states.
- The post metadata via <https://api.fxtwitter.com/hcohl001/status/1700261948808700126>
  (opened): post by @hcohl001 dated Fri Sep 08 21:36:52 UTC 2023, text
  "Here's a conjecture. Can anyone prove it? The F4 function is the Appell F4 double
  hypergeometric series https://dlmf.nist.gov/16.13.E4", media = the above image,
  0 replies, 0 quotes. This matches the report's description of the post, its date, and
  the DLMF link in the post (which also corroborates that the F4 convention was meant to
  be DLMF 16.13.4). No platform provides alt text of the formula, so the pixels themselves
  remain unchecked by me. What I *can* say: the transcribed equation is a true identity
  (Section 3 below proves the transcribed statement exactly in characteristic 0), and my
  mutation test shows the verification is sensitive to single-sign, single-shift changes
  in every part of the formula, so a "repairing" transcription error is essentially ruled
  out; only an error that maps the posted formula to an equally true one could survive.

**F4 convention.** DLMF 16.13.4, opened, TeX encoding fetched verbatim from
<https://dlmf.nist.gov/16.13.E4.tex>:

```
F4(α,β;γ,γ′;x,y) = Σ_{m,n≥0} (α)_{m+n}(β)_{m+n}/((γ)_m(γ′)_n m!n!) x^m y^n
```

This is exactly the convention used in the report's equation (1) and in my script. ✓

**Internal consistency of the transcription.** Verified by hand: with A=α+1, B=β+1,
S=α+β+1, u=t(x−1)/2, v=t(x+1)/2, the report's claims u+t=v, w=u/(1−t)², bases 4w and
−4tw, and the open-domain reformulation √|u|+√|v|<1, 4|u|<|1−t|² ⟺ R<1 are all mutually
consistent; the implied |t|<1 on the open domain is correct since
|t| ≤ |u|+|v| ≤ (√|u|+√|v|)² < 1. The parameter exclusions (A,B,S ∉ {0,−1,…} for the
literal series) are the right ones: A,B appear alone in F4 denominators; S only through
the collection (2).

## 2. Proof, step by step (line numbers from the report)

**(2), lines 62–64** — Verified symbolically (my check S1: `4^k(S/2)_k((S+1)/2)_k =
(S)_{2k} = (S)_{k+n}(S+k+n)_{k−n}` as polynomial identities for all 0≤n≤k≤10) and by hand.
The report's remark that the S-exclusion is removable is correct and the reason is
elementary: in the literal coefficient the denominator is (S)_{m+2n} = (S)_{k+n} with
n ≤ k, so its zeros all lie among the zeros of the numerator 4^k(S/2)_k((S+1)/2)_k. This
one-line justification is not spelled out in the report; it is a fixable presentational
omission, not a mathematical gap (the report says "The last expression is polynomial in
S", which is true).

**(3), lines 67–77** — Regrouping by k=m+n exactly as claimed; correct by inspection.

**(4), lines 84–88** — [u^k] extraction on F4(γ,B;A,B;u,t+u): re-derived by hand; the
binomial coefficient claim "(t+u)^j, divided by j!, gives 1/(r!(k−i)!)" and
"(B)_{k+r}/(B)_{r+k−i} = (B+r+k−i)_i" are both right; each [u^k t^r] receives finitely
many contributions so the formal manipulation is legitimate. ✓

**Vandermonde step, lines 90–97** — Verified symbolically in A,B,r for all k≤6 (my S2),
and by hand: C(A+k−1,k−i) = (A)_k/((A)_i(k−i)!), C(B+r+k−1,i) = (B+r+k−i)_i/i!, the middle
equality is the generalized binomial theorem, and C(A+B+r+2k−2,k) = (S+k+r)_k/k! since
A+B+r+2k−2 = (S+k+r)+k−1. ✓

**(5), lines 100–112** — Verified exactly as a t-series through t^12 for three
parameter choices and k∈{0,3,5} (my S7), matching both the series form and the
coefficient-extraction form. ✓

**(6), lines 115–134** — The report flags the rescaling z=(1−t)q as its weakest algebraic
step; it is in fact trivially safe in the formal t-ring, because (1−t) is a unit and
[q^k] f((1−t)q) = (1−t)^k [z^k] f needs only the first k Taylor coefficients of f. The
exponent bookkeeping (−γ−2k = −k from the rescaling and −γ−k from 1−t−z=(1−t)(1−q)) is
right, as is the final expansion of (1+tq/(1−q))^{γ−S} (my S7 confirms the last line
through t^12). The residue "(1−z)^{−S−k−n} has [q^{k−n}]-coefficient (S+k+n)_{k−n}/(k−n)!"
is right. ✓

**End of Proof, lines 136–146** — The claim that specializing u=t(x−1)/2 stays formal is
correct (each [t^r] receives contributions with u-degree ≤ r only). My two-variable
coefficient check S3 is the direct verification of this. My exact checks include γ=−2
(the terminating case asserted at lines 50–54 and lines 169–171) and the axis cases
x=±1 (u=0, v=0), which also agree exactly. ✓

**Passage to convergent series, lines 148–172** — Checked by hand:
- The formal identity gives equality near t=0 at fixed x because both sides are analytic
  there (the two defining series converge locally uniformly by the estimates of the next
  section) and their Taylor coefficients agree. The report compresses this to one
  sentence; standard and acceptable.
- The λ-segment construction: I re-derived d/dλ[λ/|1−λt|²] = (1−λ²|t|²)/|1−λt|⁴ ≠ 0 for
  λ<1, |t|≤1. With the endpoint in the open domain (|t|<1 there), |w_λ| is monotone and
  stays ≤ R<1; √(λ|u|)+√(λ|v|) = √λ·L < 1. Identity theorem along the segment: both sides
  analytic in λ∈(0,1] (1−λt≠0 because |λt|≤|t|<1 open, or |t|≤1, t≠1 at the boundary with
  Re(1−λt)>0 — a branch claim I verified: for |t|≤1, t≠1, Re(1−t)>0). ✓
- The boundary case (lines 160–169): |t| ≤ L² ≤ 1 forces the endpoint bound; dominated
  convergence along the segment is justified because |w_λ|, |λu|, |λv|, |λt| are
  monotone/monotone-dominated by the endpoint values, so the endpoint absolute sums
  dominate; the phase factors (1−λt)^{−γ−2(m+n)} are uniformly bounded on the segment.
  The parenthetical claim at |t|=1, L≤1 (uv=0; then u=0⟹R=0, or v=0⟹t=−1, R=1) is
  correct: 4|u| ≤ |1−t|² ≤ 4 and equality forces t=−1, |u|=1. ✓

**Absolute convergence, including the boundary, lines 173–278** — I checked every
inequality and asymptotic claim by hand:
- The degree-k term decomposition into C(k,i)²|u|^i|v|^{k−j} times the four ratios
  (lines 186–193) is an identity; the ≍ bounds on (a)_j/j! follow from DLMF 5.11.E12
  (opened, fetched its TeX: Γ(z+a)/Γ(z+b) ∼ z^{a−b}), and are indeed uniform in i,j for
  fixed parameters (small indices are finitely many cases).
- L^{2k}/(k+1) ≤ Σ C(k,i)²|u|^i|v|^{k−i} ≤ L^{2k} (line 201): the upper bound is
  (Σx)² ≥ Σx² applied to x_i = C(k,i)(√|u|)^i(√|v|)^{k−i}, the lower is Cauchy–Schwarz.
  ✓
- L=1, uv≠0 (lines 204–223): the binomial-probability setup is right; binomial tails are
  exponentially small (the Chernoff-type bound with derivative ε−p<0 is correct);
  Σ_i P_i² ≍ k^{−1/2} via max P_i = O(k^{−1/2}) for the upper bound and a Chebyshev
  3/4-probability window of O(√k) indices with Cauchy–Schwarz for the lower bound. The
  four-ratio product on central indices is ≍ k^{−d} with d=Re(A−γ) — I re-derived the
  exponent. Test d>1/2: correct.
- Axis cases (lines 224–226): orders k^{Re γ−1} (u=0) and k^{Re(γ+B−A)−1} (v=0),
  d>Re B ⟺ Re(γ+B−A)<0: correct.
- Right side: |(S+k+n)_{k−n}|/(k−n)! ≍ C(2k,k−n) (lines 229–236): the quotient is a
  product of (1+(S−1)/(j+1)) over k+n ≤ j < 2k; Σ 1/(j+1) ≤ 1 there, so two-sided
  uniform bounds exist — plus the nonpositive-integer caveat the report handles.
  (7), (k+1)^{−d} prefactor, and the n=0 lower bound R^k k^{−d−1/2}: correct.
- |t|=1 analysis (lines 257–278): C(2k,k−n)/C(2k,k) = Π_{j=1}^n (k−j+1)/(k+j) =
  Π(1−(2j−1)/(k+j)) ≤ exp(−Σ(2j−1)/(2k)) = exp(−n²/2k) — the report's bound is right;
  the band-splitting bound at 2^j√k with decay e^{−4^j/2} is right; δ = d+Re B−1 is
  right; summability exponents: −d−1/2+δ/2 < −1 ⟺ d > Re B (δ>0), and O(log k)·k^{−d−1/2}
  summable iff d>1/2 (δ≤0) — both check out, and the combined display (t=−1, x=−1)
  requires both tests as the table says.
- The four-row boundary table (lines 181–186) is consistent with all of the above. I
  could not falsify it with any parameter/point probe; the derivation supports it.

No error in any of these sections. The only fixable omissions are presentational (the
one-line pole-cancellation argument behind "polynomial in S"; the branch/analyticity
facts near t=0 asserted rather than restated), and the report itself flags and checks
its weakest step.

## 3. My recomputation (from the statement, not from their scripts)

Script: `checks/verify_glm.py`; full log: `checks/glm-output.txt`. Exact parts
use only `fractions.Fraction`; numeric parts use mpmath 1.3.0 (dps=60) and sympy 1.14.0
for the symbolic step checks (packages installed from PyPI).

- **S3 — exact, two variables.** Both sides expanded in t and y=(x−1)/2 and compared at
  every degree pair [t^r y^s], 0≤s≤r≤12 (91 coefficients per parameter set), for 7 rational
  (α,β,γ) sets: (1/2,3/2,1/3), (−1/2,5/3,−7/5), (2,−3/2,5/2), (−5/3,−2/7,−11/6),
  (0,3/5,−4/3), (7/4,−9/4,3/7) and the terminating γ=−2 case. **All coefficients agree
  exactly in every set.** This is a stronger check than t^r alone: it fixes the polynomial
  dependence on x too.
- **S4 — exact, fixed rational x.** For each of the 7 parameter sets and 11 rational x
  values (2/5, −3/4, 1/3, 7/2, −7/5, −2, −5/9, 1/2, 0, 1, −1), including non-integer and
  negative non-integer α, β, γ, and the degenerate axes x=1 (u=0) and x=−1 (v=0): the
  coefficients of t^0..t^12 agree exactly (77 comparisons). **All pass.**
- **S1/S2/S7 — exact proof-step checks.** Equation (2) as polynomial identity in S
  (k≤10); the Vandermonde convolution symbolically in (A,B,r) (k≤6); steps (5)–(6) as
  exact t-series through t^12. **All pass.**
- **S5 — numeric, inside the convergence region.** Six (α,β,γ,x,t) points with
  γ non-terminating, all satisfying √|u|+√|v|<1 and 4|u|<|1−t|², evaluated at 60 decimal
  digits with two cutoffs each (cutoffs sized so the truncation tail is below the
  comparison threshold). Maximum scaled |left−right| = **3.93e−60**; cutoff-doubling drift
  at every point ≤ 4.3e−55 and consistent with the estimated truncation tails. **Pass.**
- **S6 — numeric, the novelty chain.** F4(γ,B;A,B;u,v) vs (1−t)^{−γ}H_{a,b}^{c}(x,y) with
  H3(a,b;c;x,y)=Σ (a)_{2m+n}(b)_n/((c)_{m+n} m!n!) x^m y^n at w=u/(1−t)², (1−t)w, four
  points: scaled differences ≤ **8.7e−62**. Pass. This independently confirms the key
  identity of the report's final section, which matters for the novelty verdict.
- **Mutation test.** Five separate single-edit perturbations of the statement (flip the sign of
  the (−4)^n base; (S+1)/2→(S+2)/2; (S)_{m+2n}→(S)_{m+n}; (S−γ)_n→(S−γ)_{m+n}; swap the
  F4 denominator parameters (A)_i(B)_j) each break 66–78 of the 91 compared coefficients:
  my verification is genuinely sensitive to the exact form of (1).
- **Extra cross-check.** The left side equals Σ_j (γ)_j/(α+1)_j P_j^{(α,β)}(x)t^j, the
  Jacobi generating function the report attributes to Chan–Chen–Chyan–Srivastava (71):
  confirmed exactly through t^8 in two parameter sets.

## 4. Novelty: "How the known reductions imply (1)" (lines 343–end)

The section's chain is:

1. ζ(1−ζ)=w, X=(1−t)ζ, Y=X+t ⟹ X(1−Y)=u, Y(1−X)=v, (1−Y)(1+b)=1−t for b=ζ/(1−ζ) —
   I verified all four by direct algebra. ✓
2. "Combining DLMF 16.16.5_5 and 16.16.8 yields F4(γ,B;A,B;u,v) = (1−Y)^{−γ}
   F1(γ;S−γ,γ+1−A;A;X,X/(1−Y))". I fetched both DLMF equations (TeX above and at
   <https://dlmf.nist.gov/16.16.E5_5.tex>, <https://dlmf.nist.gov/16.16.E8.tex>) and
   re-derived the combination by hand: applying the second F1 transformation of 16.16.8
   to the output of 16.16.5_5 gives exactly the report's reduced form, with prefactor
   (1−X)^{−γ}(1−Y)^{−γ}·(1−X)^{γ} = (1−Y)^{−γ}, new upper parameters (S−γ, γ+1−A), and
   arguments (X, X/(1−Y)). **Correct, and both ingredients are published: 16.16.5_5 was
   added to the DLMF effective release 1.1.12 (2023-12-15, i.e., after Cohl's post), but
   the same repeated-parameter F4→F1 reduction is classical** (the report cites Schlosser
   eq. (43), Chan et al. eq. (75), Ernst–Karlsson eq. (22)).
3. The F1→H3 quadratic step (attributed to Ernst–Karlsson eq. (42), b/(1+b)²=w etc.):
   the algebra of the arguments checks by hand (b/(1+b)²=ζ(1−ζ)=w,
   X/(1+b)=(1−t)w), and my S6 confirms the resulting identity numerically to 1e−60 from
   the standard H3 series. I could not open the Ernst–Karlsson preprint itself: the
   cloudfront PDF link and the Authorea copy of "Multiple (q-)hypergeometric quadratic
   transformations..." (Ernst & Karlsson) are both behind Cloudflare bot checks as of
   2026-11-21; zbMATH's record for the 2025 paper also did not resolve for me. Search
   snippets of that paper (opened) confirm it contains numbered theorems on F1/F4
   quadratic transformations and cross-references to its eq. (42), consistent with the
   report's citation, but I did not read eq. (42) itself. The Srivastava–Karlsson
   treatise (cited via Shpot for the F4–H3 relation, p. 57) exists on archive.org but the
   download timed out for me.
4. The re-collection of the H3 series into (1): the coefficient formula
   (γ)_{2k−j}(S−γ)_j… → (3) via one Vandermonde/binomial convolution, then the
   duplication identity (2) to reach the literal printed form — verified by hand and by
   my S1/S3 (the collected form (3) and the literal (1) agree exactly, so the whole chain
   lands on the posted formula).

**Judgment.** (1) does follow from published transformations in a few lines: the
repeated-parameter F4→F1 reduction (DLMF 16.16.5_5 + 16.16.8, or the older references),
one quadratic change of variables to Horn's H3 function (a classical Bailey-type
transformation; published in the sources the report cites, which I could only verify
numerically, not read), the duplication formula, and one Vandermonde convolution. So the
identity is a known result in disguise: nothing about the F4→H3 side is essentially new.
What does appear to be new — and I found nothing contradicting the report's search here —
is the explicit collected double series with the (α+β+1)_{m+2n} denominator as displayed;
nobody seems to have written it out. But "writing out the re-collected H3 series" is a
mechanical step, so I put this as FOLLOWS FROM rather than NEW. The γ=S specialization
reducing the n-sum to n=0 and matching Cohl–MacKenzie eq. (7) is trivially consistent
((S−γ)_n = 0 for n≥1).

## 5. Sources

| Source | Opened? | Use |
| --- | --- | --- |
| <https://pbs.twimg.com/media/F5iLOskbAAA0o8A.jpg?name=orig> | downloaded, not viewable (no vision/OCR in this environment) | dimensions/bytes only |
| <https://api.fxtwitter.com/hcohl001/status/1700261948808700126> | yes | post text, date, media link, reply/quote counts |
| <https://dlmf.nist.gov/16.13.E4.tex> | yes | F4 definition matches the report's convention |
| <https://dlmf.nist.gov/16.16.E5_5.tex>, <https://dlmf.nist.gov/16.16.E8.tex> | yes | the two transformations behind the F4→F1 reduction; 5_5 annotated "Addition (effective with 1.1.12)" |
| <https://dlmf.nist.gov/15.8.E1.tex> | yes | Euler transformation, as cited |
| <https://dlmf.nist.gov/5.11.E12.tex> | yes | Γ(z+a)/Γ(z+b) ∼ z^{a−b}, supports the ≍ bounds |
| Ernst–Karlsson preprint (cloudfront link in the report) and its Authorea copy | no — Cloudflare "Just a moment" for both | F1→H3 step verified numerically instead (S6) |
| <https://archive.org/details/in.ernet.dli.2015.134480> (Srivastava–Karlsson book) | metadata yes, download timed out | F4–H3 relation per Shpot p. 57 not read |
| Search snippets of "Multiple (q-)hypergeometric quadratic transformations …" (Ernst & Karlsson) via web search | yes (snippets only) | consistent with the report's citation of eqs. (22), (42)–(43) |

## 6. Left unchecked

- The pixels of the posted image: no image-viewing or OCR capability in this sandbox.
- The exact statement of Ernst–Karlsson eqs. (42)–(43) and the Srivastava–Karlsson
  book's H3–F4 list: downloads blocked/timed out; the step was verified numerically
  instead (S6), so correctness of the report is unaffected; only the *attribution* rests
  on the report's word (and on corroborating snippets).
- The report's own numeric/symbolic logs (`numeric.txt`, `symbolic.txt`): not re-run
  (they are superseded by my independent S3–S5).
- Several rows of the literature table (Srivastava 1973, Niukkanen 1999, Vyas–Fatawat,
  Qureshi et al., Shpot, Cohl–MacKenzie, Schlosser, Chan et al.): not opened; the two
  claims of theirs I could test without opening them (γ=S reduction; Jacobi generating
  function (71)) both check out.

## 7. Reproducing my checks

```sh
python3 checks/verify_glm.py   # prints S1..S7 and PASS/FAIL
```

`output.txt` is the saved output of that run. All checks report PASS: equation (1) as
transcribed is an identity of formal power series in t (with polynomial dependence on x,
rational in the parameters) inside the stated parameter region, extends to the
terminating case γ=−N, and both defining series converge absolutely and agree on the
common open domain; the report's own boundary analysis also survived scrutiny.
