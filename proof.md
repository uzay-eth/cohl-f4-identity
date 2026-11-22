# Cohl’s Appell $`F_4`$ conjecture

**Proved as written.** The identity holds formally in $`t`$ and as an equality of
the two absolutely convergent defining series. Known transformations also
imply the formula via Horn's $`H_3`$ function. The literature search found no
prior source displaying the full double series as posted. The direct X quote
listing was inaccessible.

## Statement from the image

Howard S. Cohl posted the following formula on 2023-09-08. The supplied
2048-by-362 image and the original media URL served byte-identical files.
[Original post](https://x.com/hcohl001/status/1700261948808700126) ·
[Original image](https://pbs.twimg.com/media/F5iLOskbAAA0o8A.jpg?name=orig)

```math
\begin{aligned}
&F_4\left(\gamma,\beta+1;\alpha+1,\beta+1;
      \frac{t(x-1)}{2},\frac{t(x+1)}{2}\right)\\
&=\frac{1}{(1-t)^\gamma}
  \sum_{m,n=0}^\infty
  \frac{(\alpha+\beta-\gamma+1)_n(\gamma)_{m+n}
        \left(\tfrac12(\alpha+\beta+1)\right)_{m+n}
        \left(\tfrac12(\alpha+\beta+2)\right)_{m+n}}
       {(\alpha+1)_{m+n}(\alpha+\beta+1)_{m+2n}\,m!\,n!}
  \left(\frac{2t(x-1)}{(t-1)^2}\right)^m
  \left(-\frac{2t^2(x-1)}{(t-1)^2}\right)^n . \qquad (1)
\end{aligned}
```

Here $`(a)_j=\prod_{i=0}^{j-1}(a+i)`$, with $`(a)_0=1`$. The $`F_4`$ convention is

```math
F_4(a,b;c,d;u,v)=
\sum_{i,j\ge0}\frac{(a)_{i+j}(b)_{i+j}}{(c)_i(d)_j\,i!j!}u^i v^j .
```

[Definition in DLMF 16.13.4](https://dlmf.nist.gov/16.13.E4)

Write $`A=\alpha+1`$, $`B=\beta+1`$, $`S=\alpha+\beta+1=A+B-1`$,
$`u=t(x-1)/2`$, and $`v=t(x+1)/2=u+t`$.
Take $`\alpha,\beta,\gamma,x,t`$ complex and
$`A,B,S\notin\{0,-1,-2,\ldots\}`$ for the literal coefficients in (1).
The exclusion on $`S`$ is removable after the cancellation in (2).
There is no parameter restriction on $`\gamma`$. For nonterminating $`\gamma`$ the
common open domain of absolute convergence is
$`\sqrt{|u|}+\sqrt{|v|}<1`$ and $`4|u|<|1-t|^2`$. These imply $`|t|<1`$.
Use the branch of $`(1-t)^{-\gamma}`$ equal to $`1`$ at $`t=0`$.
The exact boundary conditions are given below. If $`\gamma=-N`$ for an integer
$`N\ge0`$, both sums terminate at total degree $`N`$. In that case (1) holds
for every $`x`$ and every $`t\ne1`$, and the combined right side has a removable
singularity at $`t=1`$.

## Proof

For $`k=m+n`$, split the product $`(S)_{2k}`$ into its even and odd factors.
It follows that

```math
\frac{4^k(S/2)_k((S+1)/2)_k}{(S)_{m+2n}}
  =\frac{(S)_{2k}}{(S)_{k+n}}
  =(S+k+n)_{k-n}. \qquad (2)
```

The last expression is polynomial in $`S`$. Let $`w=u/(1-t)^2`$. The two
powers in (1) have bases $`4w`$ and $`-4tw`$, so its right side becomes

```math
(1-t)^{-\gamma}
\sum_{k=0}^\infty \frac{(\gamma)_k}{(A)_k}w^k
\sum_{n=0}^k
 \frac{(S-\gamma)_n(S+k+n)_{k-n}}{n!(k-n)!}(-t)^n . \qquad (3)
```

First treat $`u,t`$ as independent formal variables. To extract $`[u^k]`$ on
the left, expand $`(t+u)^j`$ in the defining $`F_4`$ series. The term with first
index $`i`$ contributes to $`[u^kt^r]`$ when $`0\le i\le k`$ and $`j=r+k-i`$.
The binomial coefficient from $`(t+u)^j`$, divided by $`j!`$, is
$`1/(r!(k-i)!)`$. Also $`(B)_{k+r}/(B)_{r+k-i}=(B+r+k-i)_i`$. Hence

```math
[u^k] F_4(\gamma,B;A,B;u,t+u)
 =\sum_{r\ge0}\frac{(\gamma)_{k+r}}{r!} t^r
  \sum_{i=0}^k\frac{(B+r+k-i)_i}{(A)_i i!(k-i)!}. \qquad (4)
```

Multiply the inner sum by $`(A)_k`$. With
$`\binom{a}{j}=a(a-1)\cdots(a-j+1)/j!`$, it becomes

```math
\sum_{i=0}^k \binom{A+k-1}{k-i}\binom{B+r+k-1}{i}
 =\binom{A+B+r+2k-2}{k}
 =\frac{(S+k+r)_k}{k!}.
```

The middle equality takes the coefficient of $`z^k`$ in
$`(1+z)^{A+k-1}(1+z)^{B+r+k-1}`$ by the binomial theorem.
Thus (4) is $`(\gamma)_k/(A)_k`$ times

```math
\begin{aligned}
\sum_{r\ge0}\frac{(\gamma+k)_r}{r!}t^r\frac{(S+k+r)_k}{k!}
 &= [z^k](1-z)^{-S-k}\left(1-\frac{t}{1-z}\right)^{-\gamma-k}\\
 &= [z^k](1-z)^{\gamma-S}(1-t-z)^{-\gamma-k}. \qquad (5)
\end{aligned}
```

Every factor here has a formal binomial expansion with constant term $`1`$.
For the first equality, use
$`(S+k+r)_k/k!=[z^k](1-z)^{-S-k-r}`$ and sum in $`t`$.
Now set $`z=(1-t)q`$ inside the coefficient extraction. This multiplies
the $`k`$th coefficient by $`(1-t)^k`$, so the last line of (5) equals

```math
\begin{aligned}
(1-t)^{-\gamma-2k}[q^k](1-(1-t)q)^{\gamma-S}(1-q)^{-\gamma-k}
 &= (1-t)^{-\gamma-2k}[q^k](1-q)^{-S-k}
    \left(1+\frac{tq}{1-q}\right)^{\gamma-S}\\
 &= (1-t)^{-\gamma-2k}
    \sum_{n=0}^k\frac{(S-\gamma)_n(S+k+n)_{k-n}}{n!(k-n)!}(-t)^n . \qquad (6)
\end{aligned}
```

For the last line, the term with index $`n`$ supplies $`q^n`$ and leaves
$`(1-q)^{-S-k-n}`$, whose coefficient of $`q^{k-n}`$ is
$`(S+k+n)_{k-n}/(k-n)!`$. Multiplying (6) by $`(\gamma)_k/(A)_k`$
gives the $`u^k`$ coefficient of (3), which proves the formal identity.
All operations are coefficientwise, with only finite sums contributing
to each coefficient of the formal variables. Substituting $`u=t(x-1)/2`$
preserves that property and proves (1) as a formal $`t`$-series with polynomial
coefficients in $`x`$ and rational dependence on $`\alpha,\beta,\gamma`$.

The weakest algebraic step to check is the scaling $`z=(1-t)q`$ in (6).
The rescaling of a $`k`$th coefficient contributes $`(1-t)^{-k}`$,
and factoring $`1-t-z=(1-t)(1-q)`$ contributes
$`(1-t)^{-\gamma-k}`$. Together they give the exponent $`-\gamma-2k`$.
The term $`(-t)^n`$ in (6) gives the minus sign in the second base of (1).
Equivalently, (5)–(6) apply Euler's Gauss transformation to a Gauss series
whose transformed upper parameter is $`-k`$.
[DLMF 15.8.1](https://dlmf.nist.gov/15.8.E1)

### Passage to convergent series

First we justify the local uniform convergence on the common open domain,
using the estimates derived in the next section; after that, the formal
identity is turned into an equality of analytic functions.

Fix $`x`$ and $`t\ne1`$, keep the branch convention above, and recall
$`L=\sqrt{|u|}+\sqrt{|v|}`$ and $`R=4|w|`$. Let $`K`$ be a compact set of such
$`t`$ with $`L<1`$ and $`R<1`$ throughout, so that $`L\le L_0<1`$ and $`R\le R_0<1`$
on $`K`$ for some $`L_0,R_0`$, and $`|1-t|`$ has a positive minimum on $`K`$.
The estimates of the next section bound the absolute value of the
degree-$`k`$ part of the left series by a constant multiple of
$`(k+1)^P L^{2k}`$, and the absolute
value of the degree-$`k`$ part of (3) — that is, of the right side apart
from the prefactor $`(1-t)^{-\gamma}`$ — by a constant multiple of
$`(k+1)^{P'} R^{k}`$. The constants depend only on the parameters, on $`x`$,
and on the minimum of $`|1-t|`$ on $`K`$: the four gamma-ratio bounds $`\asymp`$
carry parameter-only constants, and the remaining factors are literal
powers of $`u,v,w`$ and $`t`$ whose moduli are maximized on $`K`$ by continuity.
Both bounds are summable in $`k`$. By the Weierstrass
M-test, both defining series converge absolutely and locally uniformly on
the common open domain, so both sides of (1) are analytic in $`t`$ there
and their Taylor coefficients at $`t=0`$ can be computed term by term.
Those Taylor coefficients agree with the formal ones: on the left, the
term indexed $`(i,j)`$ is a constant multiple of $`t^{i+j}`$, so $`[t^r]`$
receives only the terms with $`i+j\le r`$, a finite sum taken literally
from the $`F_4`$ series; on the right, the term indexed $`(m,n)`$ is a
constant multiple of $`t^{m+2n}(1-t)^{-\gamma-2(m+n)}`$, whose binomial
expansion contributes only powers $`t^r`$ with $`r\ge m+2n`$, so $`[t^r]`$ is
again a finite sum taken literally from the printed series. In both cases
the Taylor coefficient is the corresponding coefficient in the formal
identity, which states that these agree. The formal identity of the
preceding section, with $`u=t(x-1)/2`$ substituted, therefore proves that
the two sides are equal as analytic functions of $`t`$ on a neighborhood
of $`t=0`$.

For a given point of the common open domain, set $`t_\lambda=\lambda t`$ for
$`0\le\lambda\le1`$ and keep $`x`$ fixed. The left arguments become
$`\lambda u,\lambda v`$. The modulus of $`w`$ along this segment is
$`|u|\lambda/|1-\lambda t|^2`$. The derivative of
$`\lambda/|1-\lambda t|^2`$ is
$`(1-\lambda^2|t|^2)/|1-\lambda t|^4`$, which is positive for $`0\le\lambda<1`$.
Thus the segment stays in the common open domain. The identity theorem
extends equality to the point.

For a point on the common absolute convergence boundary with $`t\ne1`$, take
the same segment with $`0\le\lambda<1`$. At the endpoint $`|t|\le1`$, so the
same monotonicity holds and the segment before the endpoint is in the open
domain. Both left arguments and the right arguments $`4w,-4tw`$ have moduli
bounded by their respective endpoint moduli. The absolute sums at the
endpoint therefore dominate the series along the segment. Taking
$`\lambda\uparrow1`$ proves equality there by dominated convergence, with the
same branch since $`\mathrm{Re}(1-t)>0`$ when $`|t|\le1`$ and $`t\ne1`$.
When $`\gamma=-N`$, the left side is polynomial in $`t`$ at fixed $`x`$, and the
finite right side is rational. Equality of their Taylor series gives
equality of the rational functions and the removable value at $`t=1`$.

### Absolute convergence, including the boundary

Suppose $`\gamma\notin\{0,-1,-2,\ldots\}`$ and $`t\ne1`$. Put
$`L=\sqrt{|u|}+\sqrt{|v|}`$, $`R=4|w|`$ and $`d=\mathrm{Re}(A-\gamma)`$.
The common domain of absolute convergence is given by $`L\le1`$ and $`R\le1`$,
with the following additional tests whenever the corresponding equality
holds.

| Boundary | Additional requirement |
| --- | --- |
| $`L=1`$ and $`uv\ne0`$ | $`d>1/2`$ |
| $`L=1`$ and $`v=0`$ | $`d>\mathrm{Re}\,B`$ |
| $`L=1`$ and $`u=0`$ | $`\mathrm{Re}\,\gamma<0`$ |
| $`R=1`$ | $`d>1/2`$ |

Here is a derivation. For fixed $`a\notin\{0,-1,-2,\ldots\}`$,
$`|(a)_j|/j!`$ is bounded above and below by positive constant multiples of
$`(j+1)^{\mathrm{Re}\,a-1}`$ by the gamma ratio asymptotic.
We write $`\asymp`$ for such two-sided bounds. If $`a`$ is a nonpositive
integer the sequence is eventually zero and every bound used below for it
is an upper bound. [DLMF 5.11.12](https://dlmf.nist.gov/5.11.E12)

On the left, the term at $`(i,j)`$ with $`i+j=k`$, in absolute value, is
$`\binom{k}{i}^{2}|u|^i|v|^j`$ times the product of four ratios
$`|(\gamma)_k|/k!`$, $`|(B)_k|/k!`$, $`i!/|(A)_i|`$, and $`j!/|(B)_j|`$.
This product lies between positive constant multiples of fixed powers
of $`k+1`$, uniformly in $`i,j`$. The sum of the binomial weights lies between
$`L^{2k}/(k+1)`$ and $`L^{2k}`$ by Cauchy–Schwarz and the binomial theorem
applied to their square roots. Therefore $`L<1`$ gives convergence and
$`L>1`$ gives divergence.

If $`L=1`$ and $`uv\ne0`$, set $`p=\sqrt{|u|}/L\in(0,1)`$. The binomial
weights are the squares of the probabilities
$`P_i=\binom{k}{i}p^i(1-p)^{k-i}`$. For any
$`0<\epsilon<\min(p,1-p)`$, indices with $`i<\epsilon k`$ or $`k-i<\epsilon k`$
have exponentially small total probability. For example, the first
tail is bounded by
$`[e^{h\epsilon}(1-p+pe^{-h})]^k`$ for $`h>0`$, and the base is less than
$`1`$ for small $`h`$ because its derivative at $`0`$ is $`\epsilon-p<0`$.
Use $`1-p`$ for the other tail. Uniformly on the remaining indices the
four ratios have product $`\asymp k^{-d}`$. There
$`\sum_i P_i^2\asymp k^{-1/2}`$. For the upper bound, $`\sum_iP_i=1`$ and
the maximum of $`P_i`$ is $`O(k^{-1/2})`$ by Stirling's formula at a mode
$`\lfloor(k+1)p\rfloor`$, which follows by comparing consecutive
probabilities. For the lower bound, Chebyshev's inequality places at
least $`3/4`$ of the probability within $`2\sqrt{kp(1-p)}`$ of the mean $`kp`$.
For large $`k`$ this interval lies among the remaining indices and contains
$`O(\sqrt{k})`$ integers, so Cauchy–Schwarz gives the lower bound. The
exponentially small tails stay negligible after multiplication by any
fixed power of $`k+1`$. The absolute degree-$`k`$ sum is consequently
$`\asymp k^{-d-1/2}`$, giving the stated test. If $`u=0`$ or $`v=0`$, only
one axis contributes. The surviving coefficients have orders
$`k^{\mathrm{Re}\,\gamma-1}`$ and
$`k^{\mathrm{Re}(\gamma+B-A)-1}`$, respectively, giving the other
two left boundary tests.

On the right, $`L\le1`$ implies $`|t|=|v-u|\le|v|+|u|\le1`$.
Uniformly for large $`k`$ and $`0\le n\le k`$,

```math
\frac{|(S+k+n)_{k-n}|}{(k-n)!}\asymp\binom{2k}{k-n}.
```

Indeed, the quotient by the binomial coefficient is the product of
$`|1+(S-1)/(j+1)|`$ for $`k+n\le j<2k`$. Each factor is within
$`|S-1|/(j+1)`$ of $`1`$, and the sum of these bounds is bounded independently of
$`k,n`$, giving upper and positive lower bounds on the product for large $`k`$.
Consequently the absolute degree-$`k`$ sum in (3), apart from the prefactor,
is comparable to

```math
(k+1)^{-d}|w|^k
\sum_{n=0}^k \frac{|(S-\gamma)_n|}{n!}\binom{2k}{k-n}|t|^n. \qquad (7)
```

The $`n=0`$ term alone gives a lower bound comparable to
$`R^k k^{-d-1/2}`$ by Stirling's formula. Since every binomial factor is
at most $`\binom{2k}{k}`$, the remaining sum adds at most a polynomial factor
for $`|t|\le1`$ and a bounded factor for $`|t|<1`$. This proves the right
tests for $`R<1`$, $`R>1`$, and $`R=1`$ with $`|t|<1`$.

For $`|t|=1`$ and $`L\le1`$, equality in
$`1\le |u|+|v|\le L^2\le1`$ forces $`uv=0`$. If $`u=0`$, then $`R=0`$.
If $`v=0`$, then $`|u|=1`$, and $`R\le1`$ is possible only at $`t=-1`$ with
$`R=1`$. The left test already requires $`d>\mathrm{Re}\,B`$ there.
For $`k\ge1`$ one has

```math
\frac{\binom{2k}{k-n}}{\binom{2k}{k}}
 =\prod_{j=1}^n\left(1-\frac{2j-1}{k+j}\right)
 \le \exp\left(-\frac{n^2}{2k}\right).
```

Write $`\delta=\mathrm{Re}(S-\gamma)=d+\mathrm{Re}\,B-1`$.
The sum in (7), divided by $`\binom{2k}{k}`$, is at most a constant times
$`\sum_{n=0}^k (n+1)^{\delta-1}\exp(-n^2/(2k))`$. The part
$`n\le\sqrt{k}`$ is $`O(k^{\delta/2})`$ for $`\delta>0`$,
$`O(\log(k+2))`$ for $`\delta=0`$, and $`O(1)`$ for $`\delta<0`$.
Splitting the rest at $`2^j\sqrt{k}`$ gives the same bounds, since the
exponential decreases as $`\exp(-4^j/2)`$ on each such band. When $`\delta>0`$,
the resulting degree-$`k`$ bound is $`O(k^{-d-1/2+\delta/2})`$, summable
because $`d>\mathrm{Re}\,B`$. When $`\delta\le0`$, $`d>1/2`$ gives
summability even with the logarithm. Necessity of $`d>1/2`$ again follows
from $`n=0`$. This completes the common boundary test.

## Computational checks

The numerical check uses the two series exactly as transcribed, evaluated
separately in mpmath 1.3.0. Each of 81 cases was evaluated at 85 and 115
decimal working digits with two different total-degree cutoffs. The cases
include 36 seeded real and 36 seeded complex examples, including points
closer to the boundary, and nine deliberate cases covering $`t=0`$, $`x=\pm1`$,
termination and parameters close to denominator poles. The log lists the
parameters, both approximations, and changes under increasing both precision
and cutoff. Errors are scaled by $`\max(1,|a|,|b|)`$ for the two values being
compared. The maximum scaled discrepancy between the 115-digit evaluations was
$`1.6528182\times10^{-91}`$. The maximum changes between the low- and
high-precision runs were $`3.9028162\times10^{-75}`$ and
$`2.0863434\times10^{-75}`$ on the left and right, respectively.
The comparison is numerical evidence.

The exact check uses SymPy 1.14.0 to expand each literal side through
$`t^{40}`$ inclusive and compare all 861 coefficients in $`t`$ and
$`y=(x-1)/2`$. Since $`x=2y+1`$, this tests full symbolic dependence on $`x`$.
It multiplies $`[t^{k+r}y^k]`$ by $`(A)_k k!r!/(\gamma)_k`$ in the rational
function field, cancels the paired half-parameters by exact polynomial
division, and compares expanded polynomials. It checks cancellation of $`A`$
on the left symbolically before using the remaining variables $`S,\gamma`$.
The coefficient expressions are obtained directly from the defining $`F_4`$
and the binomial expansion of the printed right side. All residuals are
zero. [Numerical log](checks/numeric.txt) · [Symbolic log](checks/symbolic.txt)

Run all four checks of this repository from the repository root with
Python 3 and the pinned packages from [requirements.txt](requirements.txt);
each check prints `PASS` on success:

```sh
python3 -m pip install -r requirements.txt
python3 checks/run_all.py
```

The individual checks are `checks/check_numeric.py` and
`checks/check_symbolic.py` (the two runs described above),
`checks/check_opus.py` and `checks/verify_glm.py` (the independent
reviewers' programs). Each program defaults to a reduced size so the whole
suite finishes in about three minutes. The two logs above are the full-size
runs, with 81 cases and expansion through $`t^{40}`$; `--order 40`
regenerates [symbolic.txt](checks/symbolic.txt) and `--full`
regenerates [numeric.txt](checks/numeric.txt).

## Literature search

Search date for every row is 2026-11-21 UTC. Searches included the status
ID and media ID, Cohl's handle and name, and combinations of
`Appell F4`, `Jacobi generating function`, `transformation`,
`conjecture`, `proof`, and `solved`. I looked for the specific
$`F_4(\gamma,B;A,B;u,u+t)`$ with the two rational arguments and the
$`(S)_{m+2n}`$ denominator on the proposed right side.

| Source | Date read | Finding |
| --- | --- | --- |
| [Original X post](https://x.com/hcohl001/status/1700261948808700126), [quotes](https://x.com/hcohl001/status/1700261948808700126/quotes), [author timeline search](https://x.com/search?q=from%3Ahcohl001%20%28Appell%20OR%20F4%29%20since%3A2023-09-08&src=typed_query&f=live), [public post metadata](https://api.fxtwitter.com/hcohl001/status/1700261948808700126) | 2026-11-21 | The post is dated 2023-09-08. The public metadata response matched the supplied snapshot byte for byte and reports zero replies and zero quotes. X's web pages returned access errors or a generic shell with no quote contents. Searches of publicly indexed posts found no follow-up. The metadata is one snapshot, so the live listing could not be independently verified. |
| [Cohl's arXiv author page](https://arxiv.org/a/cohl_h_1.html), [NIST publications](https://www.nist.gov/publications/search_by_author/1159096?page=0), [Costas-Santos's bibliography](https://rscosan.com/papers) | 2026-11-21 | Screened the relevant titles and abstracts after the post, earlier manuscripts with later versions or publication, and full text for the closest classical Jacobi papers below. No match surfaced in those lists. |
| Cohl and Costas-Santos, [Multi-integral representations for Jacobi functions](https://arxiv.org/pdf/2308.13652), [NIST publication record](https://www.nist.gov/publications/multi-integral-representations-jacobi-functions-first-and-second-kind) | 2026-11-21 | Preprint 2023-08-25, NIST publication date 2023-10-05. The paper derives integral and derivative representations for Jacobi functions. Full-text search found no Appell reference or target identity. |
| Cohl, Costas-Santos, Durand, Montoya and Ólafsson, [Double summation addition theorems](https://arxiv.org/abs/2306.03035), [2025 chapter](https://repository.lsu.edu/mathematics_pubs/1608/), [full text](https://rscosan.com/documents/paper_37.pdf) | 2026-11-21 | June 2023 preprint, 2025 publication. The double sums are Jacobi addition theorems. Full-text search found no Appell reference or target identity. |
| Cohl and Durand, [Jacobi product integral](https://arxiv.org/pdf/2508.08085), [journal article](https://www.sciencedirect.com/science/article/abs/pii/S0022247X26000624), [2026 erratum to their related Nicholson paper](https://epubs.siam.org/doi/full/10.1137/25M1794759) | 2026-11-21 | 2025 product-integral preprint, 2026 article. The preprint's Jacobi products and integrals and the erratum's Nicholson formulas address different identities. |
| Cohl and Schlosser, [Askey–Wilson product formulas](https://arxiv.org/html/2411.03571v4), [bilateral basic series](https://arxiv.org/abs/2504.21782), and Cohl and Costas-Santos, [symmetric q-Askey polynomials](https://arxiv.org/abs/2307.04884) | 2026-11-21 | 2024–2026 papers or revisions on basic hypergeometric series. The Askey–Wilson paper's Appell reference is to a $`q`$-Appell $`\Phi^{(1)}`$. The survey found no statement of (1) in these sources. |
| [DLMF Appell definitions](https://dlmf.nist.gov/16.13), [transformations](https://dlmf.nist.gov/16.16), [Jacobi generating functions](https://dlmf.nist.gov/18.12), [release history](https://dlmf.nist.gov/errata/#V1.1.12) | 2026-11-21 | DLMF added [16.16.5_5](https://dlmf.nist.gov/16.16.E5_5) on 2023-12-15. It gives the older $`F_1`$ reduction in Schlosser's eq. (43) and Chan et al.'s eq. (75). That reduction applies to the left side of (1). |
| Schlosser, [Multiple hypergeometric series, eq. (43)](https://www.mat.univie.ac.at/~schlosse/pdf/appell-survey.pdf), and Chan, Chen, Chyan and Srivastava, [Some multiple hypergeometric transformations, eqs. (71), (75)](https://www.academia.edu/94011355/Some_multiple_hypergeometric_transformations_and_associated_reduction_formulas) | 2026-11-21 | Earlier repeated-parameter $`F_4`$ to $`F_1`$ reductions. In the 2004 paper, (71) also gives the left-side Jacobi generating function $`\sum_{j\ge0}(\gamma)_j P_j^{(\alpha,\beta)}(x)t^j/(A)_j`$ by taking one numerator weight $`\gamma`$ and one denominator weight $`A`$. |
| Cohl and MacKenzie, [Generalizations and specializations of generating functions, eq. (7)](https://arxiv.org/pdf/1210.0039) | 2026-11-21 | 2013 revision. In the arXiv PDF and LaTeX source the specialization $`\gamma=S`$ is equation (7); the arXiv HTML rendering numbers the same formula (13). After identifying the left side with the Jacobi generating function above, it is also DLMF 18.12.3 after $`\alpha\leftrightarrow\beta`$, $`x\to-x`$, $`z\to-t`$. The posted $`n`$-sum then reduces to $`n=0`$. |
| Shpot, [A massive Feynman integral and some reduction relations for Appell functions, §IV](https://arxiv.org/pdf/0711.2742#page=9) | 2026-11-21 | 2007 article. Shpot remarks in words that a repeated-parameter $`F_4`$ of his integral can alternatively be expressed in terms of Horn's $`H_3`$, citing p. 57 of Srivastava and Manocha's 1984 treatise, but displays no $`H_3`$ formula. He works out the $`F_1`$ form using the general repeated-parameter reduction, his eq. (21). |
| Ernst and Karlsson, [On multiple hypergeometric quadratic transformations, §§3–5](https://doi.org/10.22541/au.173980497.70111458/v1) | 2026-11-21 | 2025-02-17 preprint. Its (22) gives the repeated-parameter $`F_4`$ to $`F_1`$ reduction in a form equivalent to the DLMF formula, and (42) reduces that $`F_1`$ to $`H_3`$. The combination yields (1) after the re-expansion below. The preprint displays the component transformations. |
| Srivastava and Panda, [Some hypergeometric transformations involving Horn's function H₃](https://zbmath.org/3421169), [zbMATH API record](https://api.zbmath.org/v1/document/_search?search_string=ti%3A%22Some%20hypergeometric%20transformations%20involving%20Horn%22) | 2026-11-21 | Cited by Ernst and Karlsson for transformations between $`H_3`$ and $`F_4`$. I retrieved the bibliographic record for the 1974 article using the zbMATH API. No full text was available in the sources consulted, so its exact formulas are unverified. |
| Vyas and Fatawat, [Transformations and Reductions of Srivastava–Daoust Type Double Hypergeometric Functions](https://link.springer.com/chapter/10.1007/978-3-031-56307-2_15), and Liu and Wang, [Transformation and summation formulae for Kampé de Fériet series](https://www.sciencedirect.com/science/article/pii/S0022247X13006264) | 2026-11-21 | The publicly accessible abstracts describe related double-series methods, including Euler transformations in Liu–Wang (2014). The full formulas in these two papers were unavailable in the public previews, so an overlap could not be ruled out from these records. Vyas–Fatawat appeared online on 2024-03-29. |
| Qureshi et al., [Three General Double-Series Identities](https://www.mdpi.com/2504-3110/7/10/700) and [Notes on Various Implications of Bailey Transformations](https://link.springer.com/article/10.1007/s40819-023-01576-6) | 2026-11-21 | Published on 2023-09-23 and 2023-09-29, respectively. The first paper's identities use terminating Whipple transformations. The second paper's public abstract describes terminating Bailey transformations and reductions for equal arguments. The accessible material yielded no match. |
| Srivastava, [On the reducibility of Appell's function F₄](https://doi.org/10.4153/CMB-1973-049-9), Niukkanen, [Special transformations, eq. (13)](https://www.mathnet.ru/eng/rm242), and Ananthanarayan et al., [New analytic continuations](https://arxiv.org/abs/2005.07170) | 2026-11-21 | 1973, 1999 and 2020 papers. Srivastava's eq. (4) uses the same Vandermonde convolution with equal arguments. Niukkanen's transformation uses a common product denominator for its two variables. The 2020 formulas impose a different parameter relation and use quadratic Gauss transformations. |
| [MathOverflow API search](https://api.stackexchange.com/2.3/search/advanced?site=mathoverflow&q=Appell%20F4), [Math Stack Exchange API search](https://api.stackexchange.com/2.3/search/advanced?site=math&q=Appell%20F4), [indexed MathOverflow result](https://mathoverflow.net/questions/499734/solutions-of-appell-system-of-type-f-4) | 2026-11-21 | The API search found no MathOverflow hit and one Math Stack Exchange question about Laplace transforms. Supplementary indexed searches found the MathOverflow question on local solutions of the Appell system. None of these yielded the conjecture or a solution. |
| [Google Scholar search](https://scholar.google.com/scholar?q=%22Howard+Cohl%22+%22Appell%22+%22F4%22) and indexed public proof searches on [GitHub](https://github.com/search?q=%221700261948808700126%22&type=code) and [Axiom](https://www.axiommath.ai/research) | 2026-11-21 | Direct Scholar access was unavailable. No match surfaced in public indexing searches by post ID, handle or exact keywords. |

### How the known reductions imply (1)

Near $`u=t=0`$, let $`\zeta`$ be the solution of $`\zeta(1-\zeta)=w`$ with
$`\zeta=0`$ at $`w=0`$, and put $`X=(1-t)\zeta`$, $`Y=X+t`$. Then
$`X(1-Y)=u`$ and $`Y(1-X)=v`$. Combining DLMF 16.16.5_5 and 16.16.8 yields

```math
F_4(\gamma,B;A,B;u,v)
 =(1-Y)^{-\gamma}
 F_1\left(\gamma;S-\gamma,\gamma+1-A;A;X,\frac{X}{1-Y}\right).
```

This is also Ernst and Karlsson's (22). In their (42), the second $`F_1`$
argument $`b=X/(1-Y)=\zeta/(1-\zeta)`$ produces the $`H_3`$ arguments
$`b/(1+b)^2=w`$ and $`X/(1+b)=(1-t)w`$. Since $`(1-Y)(1+b)=1-t`$, the
result is

```math
F_4(\gamma,B;A,B;u,v)
   =(1-t)^{-\gamma}H_3(\gamma,S-\gamma;A;w,(1-t)w).
```

[DLMF 16.16.5_5 and 16.16.8](https://dlmf.nist.gov/16.16) ·
[Ernst and Karlsson, eqs. (22), (42)–(43)](https://doi.org/10.22541/au.173980497.70111458/v1)

The coefficient of $`w^k`$ in this $`H_3`$ is the sum, for $`0\le j\le k`$, of
$`(\gamma)_{2k-j}(S-\gamma)_j(1-t)^j/((A)_k(k-j)!j!)`$.
Expanding $`(1-t)^j`$, the coefficient of $`(-t)^n`$ becomes

```math
\frac{(\gamma)_k(S-\gamma)_n}{(A)_k n!}
\sum_{j=n}^k
  \frac{(\gamma+k)_{k-j}(S-\gamma+n)_{j-n}}{(k-j)!(j-n)!}
 =\frac{(\gamma)_k(S-\gamma)_n(S+k+n)_{k-n}}{(A)_k n!(k-n)!}.
```

The equality is the binomial convolution in degree $`k-n`$.
This recovers (3) and explains the connection between the known reductions
and the exact array of coefficients in the post.
