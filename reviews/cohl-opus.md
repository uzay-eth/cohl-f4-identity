*This review is republished as written, lightly redacted for the
repository-internal file paths of the working tree it reviewed.*

VERDICT: CORRECT
NOVELTY: FOLLOWS FROM DLMF 18.5.7–18.5.8 with Euler's transformation DLMF 15.8.1 and the duplication DLMF 5.2.8 (a five-line textbook derivation; the case gamma = alpha+beta+1 is DLMF 18.12.3)

# Independent review: Cohl's Appell $F_4$ identity (`proof.md`, commit 6f70adb)

Reviewer: Claude Opus 5.5, 2026-11-21. I did not edit any of the reviewed
files. My script and its output are in
[`checks/`](../checks/). I wrote the script from the posted image and
DLMF 16.13.4. I did not look at the proof's own check programs.

I found no mathematical error in the proof. The minor points below concern
exposition and citations; none changes the verdict.

## 1. Statement

I downloaded the original image
(<https://pbs.twimg.com/media/F5iLOskbAAA0o8A.jpg?name=orig>, opened: 2048×362
JPEG, sha256 `24d9544d…3dcd6`) and viewed it. I compared it symbol by symbol with
equation (1) of the report (REPORT.md lines 16–30):

| Item | Image | Report (1) |
| --- | --- | --- |
| $F_4$ parameters | $\gamma,\beta+1;\alpha+1,\beta+1$ | same |
| $F_4$ arguments | $t(x-1)/2,\ t(x+1)/2$ | same |
| prefactor | $1/(1-t)^\gamma$ | same |
| summation | $\sum_{m,n=0}^\infty$ | same |
| numerator | $(\alpha+\beta-\gamma+1)_n(\gamma)_{m+n}(\tfrac12(\alpha+\beta+1))_{m+n}(\tfrac12(\alpha+\beta+2))_{m+n}$ | same |
| denominator | $(\alpha+1)_{m+n}(\alpha+\beta+1)_{m+2n}\,m!\,n!$ | same |
| powers | $\big(\frac{2t(x-1)}{(t-1)^2}\big)^m\big(-\frac{2t^2(x-1)}{(t-1)^2}\big)^n$ | same |

The transcription is exact. I did not open the X post itself
(<https://x.com/hcohl001/status/1700261948808700126>), only the image.

$F_4$ convention: DLMF 16.13.4 (<https://dlmf.nist.gov/16.13.E4>; I fetched
its TeX source at <https://dlmf.nist.gov/16.13.E4.tex>) reads
$F_4(\alpha,\beta;\gamma,\gamma';x,y)=\sum_{m,n}\frac{(\alpha)_{m+n}(\beta)_{m+n}}{(\gamma)_m(\gamma')_n m!n!}x^my^n$.
This is the convention in REPORT.md lines 34–37. I did not check the
convergence condition printed on the DLMF page because the report derives its
own.

## 2. Proof, step by step

I used the report's notation: $A=\alpha+1$, $B=\beta+1$, $S=A+B-1$,
$u=t(x-1)/2$, $v=u+t$, $w=u/(1-t)^2$.

### "Proof" (lines 55–146)

- **(2), lines 57–65.** $(S)_{2k}=4^k(S/2)_k((S+1)/2)_k$ is DLMF 5.2.8. With
  $m+2n=k+n$, $(S)_{2k}/(S)_{k+n}=(S+k+n)_{k-n}$. Correct.
- **(3), lines 67–76.** $2t(x-1)/(t-1)^2=4w$ and $-2t^2(x-1)/(t-1)^2=-4tw$. The
  $4^k$ cancels against (2), leaving $w^k(-t)^n$, and $S-\gamma=\alpha+\beta-\gamma+1$.
  Correct.
- **(4), lines 78–88.** I rederived the index bookkeeping: $j=r+k-i$, the
  binomial factor divided by $j!$ is $1/(r!(k-i)!)$, and
  $(B)_{k+r}/(B)_{r+k-i}=(B+r+k-i)_i$. Correct.
- **Vandermonde step, lines 91–101.** $(A)_k/((A)_i(k-i)!)=\binom{A+k-1}{k-i}$ and
  $(B+r+k-i)_i/i!=\binom{B+r+k-1}{i}$. Chu–Vandermonde holds as a polynomial
  identity in $A$ and $B$, and $\binom{S+r+2k-1}{k}=(S+k+r)_k/k!$. Correct.
- **(5), lines 102–115.** $(\gamma)_{k+r}=(\gamma)_k(\gamma+k)_r$ and
  $(S+k+r)_k/k!=[z^k](1-z)^{-S-k-r}$. The sum over $r$ gives
  $(1-z)^{-S-k}(1-t/(1-z))^{-\gamma-k}=(1-z)^{\gamma-S}(1-t-z)^{-\gamma-k}$. This
  is valid in $\mathbb{Q}(\alpha,\beta,\gamma)[[t,z]]$ because every base has
  constant term 1. Correct.
- **(6), lines 116–133.** In $R[[t]][[q]]$, $1-t$ is a unit, so the substitution
  $z=(1-t)q$ is legitimate and rescales $[z^k]$ by $(1-t)^{-k}$. Also
  $1-t-z=(1-t)(1-q)$ and $1-(1-t)q=(1-q)(1+tq/(1-q))$. The binomial expansion
  $(1+y)^{\gamma-S}=\sum_n(S-\gamma)_n(-y)^n/n!$ yields the sum in (6). I checked
  the exponent $-\gamma-2k$ (lines 139–142) and the sign $(-t)^n$. Correct.
- **Formal conclusion, lines 133–137.** Each coefficient of $u^kt^N$ receives
  only finitely many contributions. Substituting $u=t(x-1)/2$, which has no
  constant term, is legitimate. Correct.
- **Euler remark, lines 144–146.** I confirmed it independently. The $r$-sum in
  (5) equals $\frac{(S+k)_k}{k!}\,{}_2F_1(\gamma+k,S+2k;S+k;t)$. The last form of DLMF
  15.8.1 (TeX fetched from <https://dlmf.nist.gov/15.8.E1.tex>) turns this into
  $(1-t)^{-\gamma-2k}{}_2F_1(S-\gamma,-k;S+k;t)$, which is the $n$-sum in (6).

### "Passage to convergent series" (lines 148–171)

- **Lines 150–151, local uniform convergence (minor gap, cosmetic).** The
  report asserts this without proof. It does follow from the estimates in the
  next section. The constants in those bounds depend only on the parameters,
  and the terms are bounded by $\mathrm{poly}(k)\,L^{2k}$ and
  $\mathrm{poly}(k)\,R^k$. The Weierstrass M-test on compact sets with
  $L\le L_0<1$ and $R\le R_0<1$ then applies. Weierstrass's theorem lets the
  Taylor coefficients at $t=0$ be computed term by term. On the right side,
  term $(m,n)$ has $t$-order $\ge m+2n$, so each Taylor coefficient is a finite
  sum and equals the formal coefficient. One sentence would close this.
- **Lines 152–159, star-shapedness.** $L(\lambda)=\sqrt\lambda\,L$. I verified
  the derivative
  $\frac{d}{d\lambda}\frac{\lambda}{|1-\lambda t|^2}=\frac{1-\lambda^2|t|^2}{|1-\lambda t|^4}$
  by expanding $|1-\lambda t|^2=1-2\lambda\operatorname{Re}t+\lambda^2|t|^2$.
  So the $t$-domain at fixed $x$ is star-shaped about 0, and the identity
  theorem applies. Correct.
- **Lines 161–168, boundary.** For $\lambda<1$ we have $L(\lambda)<1$, and
  $R(\lambda)<R(1)\le1$ by strict monotonicity, since $\lambda|t|<1$. Every term
  at $\lambda$ is bounded in modulus by the term at the endpoint, and the
  prefactor is continuous. Dominated convergence therefore applies, and the
  branch is continuous because $\operatorname{Re}(1-t)>0$. Correct.
- **Lines 169–171, terminating $\gamma=-N$.** Here a polynomial equals a
  rational function with its only possible pole at $t=1$. Correct.

### "Absolute convergence, including the boundary" (lines 173–277)

- **Lines 195–202.** I checked the rewriting
  $|(\gamma)_k(B)_k|/(|(A)_i(B)_j|\,i!\,j!)=\binom{k}{i}^2\times$ (four ratios)
  algebraically. Cauchy–Schwarz bounds the weight sum between $L^{2k}/(k+1)$
  and $L^{2k}$. Correct.
- **Lines 204–227, $L=1$.** When $L=1$ the weights are exactly $P_i^2$. The
  Chernoff tail bound (whose derivative at $h=0$ is $\epsilon-p$) is correct.
  The Chebyshev and Cauchy–Schwarz lower bound $\sum_IP_i^2\ge(9/16)/|I|$ is
  correct. On central indices the four ratios give
  $k^{\operatorname{Re}\gamma-1+\operatorname{Re}B-1+1-\operatorname{Re}A+1-\operatorname{Re}B}=k^{-d}$.
  The axis cases give $k^{\operatorname{Re}\gamma-1}$ and
  $k^{\operatorname{Re}(\gamma+B-A)-1}$, hence the tests $\operatorname{Re}\gamma<0$
  and $d>\operatorname{Re}B$. Correct.
- **Lines 229–253, right side.** The quotient of $(S+k+n)_{k-n}/(k-n)!$ by
  $\binom{2k}{k-n}$ is $\prod_{j=k+n}^{2k-1}(1+(S-1)/(j+1))$. The sum of
  $1/(j+1)$ over this range is at most $\ln2+o(1)$. No factor vanishes for
  large $k$, even when $S$ is a nonpositive integer. Correct. The $n=0$ lower
  bound is $R^kk^{-d-1/2}$ and the upper bound uses $\binom{2k}{k-n}\le\binom{2k}{k}$.
  Correct.
- **Lines 255–277, $|t|=1$.** The chain $1\le|u|+|v|\le L^2\le1$ forces
  $uv=0$. With $v=0$ we have $x=-1$, $u=-t$ and $R=4/|1-t|^2$, so only $t=-1$
  qualifies. The bound $\binom{2k}{k-n}/\binom{2k}{k}\le e^{-n^2/(2k)}$ uses
  $k+j\le2k$. I checked the exponent $-d-\tfrac12+\tfrac\delta2=-\tfrac{d-\operatorname{Re}B}{2}-1$.
  It is consistent with the necessity of $d>1/2$ because, for $\delta>0$,
  $d>\max(\operatorname{Re}B,1-\operatorname{Re}B)\ge\tfrac12$. Correct.

No step fails. The only gap is the unproved "local uniform convergence"
sentence (lines 150–151), and it can be fixed in one sentence as described above.

## 3. Independent recomputation

[`checks/check_opus.py`](../checks/check_opus.py) has four parts.
[`checks/opus-output.txt`](../checks/opus-output.txt) is the complete log of
`python3 checks/check_opus.py 16 6` (mpmath 1.3.0, sympy 1.14.0, run time
about 80 s).

- **Part 1, exact rational arithmetic.** Python `Fraction` throughout. For
  each case I expanded the left side from the defining $F_4$ series. I expanded
  the right side literally: every $(m,n)$ term with the unsimplified
  $(\tfrac12(\alpha+\beta+1))_{m+n}(\tfrac12(\alpha+\beta+2))_{m+n}/(\alpha+\beta+1)_{m+2n}$,
  times $t^{m+2n}(1-t)^{-\gamma-2(m+n)}$ expanded binomially. I compared
  $[t^N]$ for $N=0,\dots,16$ at 10 parameter choices. They include negative
  non-integer $\alpha,\beta,\gamma$, the points $x=1$, $x=-1$, $x=0$ and
  $x=7/2$, and the terminating case $\gamma=-3$. **All 170 coefficients agree
  exactly.** For example, at $(\alpha,\beta,\gamma,x)=(1/3,2/7,5/11,3/5)$ both
  sides give $[t^{16}]=-55728719782763788582730910032720712101703/35165635373667978288010892398634356689453125$.
- **Part 2, symbolic.** $\alpha,\beta,\gamma,x$ are SymPy indeterminates and
  $N\le6$. Every left-minus-right difference cancels to 0.
- **Part 4, negative controls.** These show the comparison is sensitive.
  Flipping the sign of the second base, using $(\alpha+\beta-\gamma)_n$, or
  using $(\alpha+\beta+1)_{m+2n+1}$ makes the two sides differ already at
  $t^2$, $t^2$ and $t^0$ respectively.
- **Part 3, numerical, 60 digits.** I computed four independent quantities at
  5 points inside the convergence domain. Four points are real, with $L$ from
  0.63 to 0.79 and $R$ from 0.37 to 0.61. One point has all of
  $\alpha,\beta,\gamma,x,t$ complex.
  - the $F_4$ double series;
  - the Jacobi series $\sum_n(\gamma)_n P_n^{(\alpha,\beta)}(x)t^n/(\alpha+1)_n$, using mpmath's `jacobi`;
  - the printed right side;
  - $(1-t)^{-\gamma}H_3(\gamma,S-\gamma;A;w,(1-t)w)$, the form in the report's novelty section.

  Each was truncated at total degree 250 and 350. The worst discrepancy
  between any of them and the right side was $9.7\times10^{-58}$ in relative
  terms. This came from the complex point, where the change between degree
  250 and degree 350 was $5.8\times10^{-42}$. At every point all four values
  agree to the 40 printed digits.

Final line of the log: `PART1 PASS, PART2 PASS, PART3 PASS, PART4 PASS`.

## 4. Novelty: "How the known reductions imply (1)" (lines 343–381)

**Checking the section.**

- I verified the $F_1$ step algebraically. In DLMF 16.16.5_5 put
  $(\alpha,\beta,\gamma)\to(\gamma,B,A)$ with $X(1-Y)=u$ and $Y(1-X)=v$. Then
  apply the second form of 16.16.8, whose first argument is $X/(X-1)$. This
  gives exactly
  $F_4(\gamma,B;A,B;u,v)=(1-Y)^{-\gamma}F_1(\gamma;S-\gamma,\gamma+1-A;A;X,X/(1-Y))$,
  because $A-(A-B)-(\gamma-A+1)=S-\gamma$ and
  $(y-x)/(1-x)=X/(1-Y)$. Both equations came from the DLMF TeX files.
- The parametrisation $X=(1-t)\zeta$ with $\zeta(1-\zeta)=w$, together with
  $b/(1+b)^2=w$, $X/(1+b)=(1-t)w$ and $(1-Y)(1+b)=1-t$, is all correct.
- The $F_1\to H_3$ step cites Ernst and Karlsson (42). I could not open that
  preprint because its host did not resolve from this machine
  (<https://d197for5662m48.cloudfront.net/documents/publicationstatus/245830/preprint_pdf/793746cb6d0de25bbcca2c3cbce48018.pdf>,
  not opened). The resulting $H_3$ identity is nevertheless confirmed in two
  ways. The re-expansion in lines 368–379 is a correct Chu–Vandermonde
  convolution, which I rederived. Part 3 also checks the $H_3$ form
  numerically to 60 digits.

**My judgement: (1) is a routine consequence of textbook formulas.**
The left side of (1) is $\sum_n(\gamma)_nP_n^{(\alpha,\beta)}(x)t^n/(\alpha+1)_n$,
by the binomial-sum form DLMF 18.5.8. Write
$P_n^{(\alpha,\beta)}(x)=\frac{(A)_n}{n!}{}_2F_1(-n,n+S;A;\tfrac{1-x}2)$ (DLMF
18.5.7) and interchange the sums. The coefficient of
$\bigl(t(x-1)/2\bigr)^k$ becomes
$\frac{(\gamma)_k(S+k)_k}{(A)_kk!}{}_2F_1(\gamma+k,S+2k;S+k;t)$.
Euler's transformation (DLMF 15.8.1) turns this into the terminating series
$(1-t)^{-\gamma-2k}{}_2F_1(S-\gamma,-k;S+k;t)$. Expanding it gives (3), and the
duplication formula DLMF 5.2.8 gives (1). That is five lines. The report's own
proof, equations (4)–(6), is the same computation with Vandermonde in place of
18.5.7, and its line 144 says so.

The special case $\gamma=S$ is a classical formula: DLMF 18.12.3 after
$\alpha\leftrightarrow\beta$, $x\to-x$, $z\to-t$ (TeX fetched from
<https://dlmf.nist.gov/18.12.E3.tex>). In that case $(S-\gamma)_n=0$ for
$n\ge1$, and only the $m$-sum survives.

Shpot (<https://arxiv.org/abs/0711.2742>; opened the arXiv HTML v1 through a
fetch agent) notes in his §IV that an $F_4$ with the same repeated-parameter
pattern $F_4(a,b;c,b;\cdot,\cdot)$ "can be alternatively expressed in terms of
the Horn function $H_3$". He cites Srivastava and Manocha, *A Treatise on
Generating Functions* (1984), p. 57, but displays no formula. I could not
consult that book (not opened), so I cannot confirm that the exact $H_3$ form
is printed there. I found no source that displays (1) itself. Cohl's
conjecture is therefore best described as a disguised form of standard
material, not a new identity. "FOLLOWS FROM" is the accurate label, even
though no source writes (1) in Cohl's double-series form.

**Citation problems in the report's literature table (they do not affect correctness).**

- Line 333 says Cohl and MacKenzie's "equation (7)" is the case $\gamma=S$. In
  the arXiv HTML v3 (<https://arxiv.org/html/1210.0039v3>, opened through a
  fetch agent), (7) is the definition of the Legendre function $P^\mu_\nu$, and
  the $\gamma=S$ formula is (13). The PDF numbering may differ: the HTML skips
  numbers, and I could not render the PDF at
  <https://arxiv.org/pdf/1210.0039> (fetched, not readable here). The report
  should either cite DLMF 18.12.3 directly or recheck the number.
- Line 334 implies more than Shpot gives. Shpot displays only Bailey's
  $F_4\to F_1$ formula (his (21)) and mentions $H_3$ in words.

## Sources

| Source | URL | Opened? |
| --- | --- | --- |
| Posted image | <https://pbs.twimg.com/media/F5iLOskbAAA0o8A.jpg?name=orig> | yes (downloaded and viewed) |
| X post | <https://x.com/hcohl001/status/1700261948808700126> | no |
| DLMF 16.13.4 | <https://dlmf.nist.gov/16.13.E4> (TeX: `16.13.E4.tex`) | yes (TeX) |
| DLMF 16.16.5_5, 16.16.8 | <https://dlmf.nist.gov/16.16> | yes, through a fetch agent (HTML and TeX) |
| DLMF 15.8.1 | <https://dlmf.nist.gov/15.8.E1> | yes (TeX) |
| DLMF 5.2.8 | <https://dlmf.nist.gov/5.2.E8> | yes (TeX of 5.2.E8a) |
| DLMF 18.5.7, 18.5.8 | <https://dlmf.nist.gov/18.5> | yes, through a fetch agent (TeX) |
| DLMF 18.12.3 | <https://dlmf.nist.gov/18.12.E3> | yes (TeX) |
| DLMF 5.11.12 (cited by the report) | <https://dlmf.nist.gov/5.11.E12> | no |
| Shpot, J. Math. Phys. 48 (2007) 123512 | <https://arxiv.org/abs/0711.2742> | yes, HTML v1 through a fetch agent |
| Cohl and MacKenzie | <https://arxiv.org/abs/1210.0039> | HTML v3 through a fetch agent; PDF not readable |
| Ernst and Karlsson preprint | link above | no (host unreachable) |
| Srivastava and Manocha, *A Treatise on Generating Functions*, 1984 | no URL | no |

## Left unchecked (the 60-minute limit was not reached; these are access limits)

- Ernst and Karlsson eqs. (22) and (42), and Srivastava and Manocha p. 57. I
  could not open either, so whether the exact $H_3$ form is printed in the
  literature is unconfirmed. The novelty label does not depend on it.
- The numbering of Cohl and MacKenzie's equation in the PDF version.
- The report's "Computational checks" and literature-search sections, beyond
  the rows discussed above. I did not rerun the proof's own checks,
  by design.
