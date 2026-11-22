*This review is republished as written, lightly redacted for the
repository-internal file paths of the working tree it reviewed; the
planned reply's repository link is shown as a placeholder.*

VERDICT: CORRECT
SENSE: CHANGE TEXT (It holds as a formal power series in t, and for t ≠ 1 wherever both sides converge absolutely: one Vandermonde sum, then Euler's transformation (DLMF 15.8.1); also via DLMF 16.16.5_5 and Horn's H3. AI-drafted proof and checks: <link to this repository>)

# Cold Review of `proof.md` and Planned Reply

Review date: 2026-11-21 UTC. Conducted cold without inspecting any other review or branch, and without reusing any of the proof's own check programs.

---

## 1. Is the Proof Correct?

**Yes (`VERDICT: CORRECT`).** Every step of the proof in `proof.md` — image transcription, formal power-series identity, analytic continuation to the common open domain, dominated-convergence passage to the absolute-convergence boundary, terminating case $`\gamma = -N`$, exact boundary classification, and the alternative Horn $`H_3`$ derivation — holds under adversarial line-by-line verification and independent recomputation.

### 1.1 Image Transcription & DLMF 16.13.4 Convention Check

- **Image integrity:** the repository's `cohl-f4.jpg` (2048×362 JPEG) and the original media URL `https://pbs.twimg.com/media/F5iLOskbAAA0o8A.jpg?name=orig` share the exact SHA-256 digest `24d9544de99e7fb324f0879ef166dc2c6be43a180b92e1bded7751ae39c3dcd6`.
- **Symbol-by-symbol check of eq. (1) (`proof.md` lines 16–30 vs `cohl-f4.jpg`):**
  - Left-hand side: $`F_4\bigl(\gamma, \beta+1; \alpha+1, \beta+1; \frac{t(x-1)}{2}, \frac{t(x+1)}{2}\bigr)`$ — exact match on function name, parameter order, semicolons/commas, and both arguments.
  - Prefactor: $`\frac{1}{(1-t)^\gamma}`$ (note $`(1-t)`$ here versus $`(t-1)^2`$ in the power bases) — exact match.
  - Summation limits: $`\sum_{m,n=0}^\infty`$ — exact match.
  - Numerator factors: $`(\alpha+\beta-\gamma+1)_n`$, $`(\gamma)_{m+n}`$, $`(\frac{1}{2}(\alpha+\beta+1))_{m+n}`$, $`(\frac{1}{2}(\alpha+\beta+2))_{m+n}`$ — exact match on signs, constants, and subscripts.
  - Denominator factors: $`(\alpha+1)_{m+n}`$, $`(\alpha+\beta+1)_{m+2n}`$, $`m!\,n!`$ — exact match (specifically confirming $`m+2n`$ in the second Pochhammer symbol).
  - Power factors: $`\bigl(\frac{2t(x-1)}{(t-1)^2}\bigr)^m`$ and $`\bigl(-\frac{2t^2(x-1)}{(t-1)^2}\bigr)^n`$ — exact match on signs, factors of $`2`$, powers of $`t`$, and $`(x-1)/(t-1)^2`$.
- **$`F_4`$ convention (`proof.md` lines 35–37 vs DLMF 16.13.4):** Checked directly against `https://dlmf.nist.gov/16.13.E4.tex`, which defines
  ```math
  F_4(\alpha,\beta;\gamma,\gamma';x,y)=\sum_{m,n=0}^\infty \frac{(\alpha)_{m+n}(\beta)_{m+n}}{(\gamma)_m(\gamma')_n\,m!\,n!}x^m y^n, \qquad \sqrt{|x|}+\sqrt{|y|}<1,
  ```
  with $`\gamma, \gamma' \notin \{0,-1,-2,\dots\}`$. Exact match.

### 1.2 Step-by-Step Audit of the Proof

No errors or logical gaps were found. Each section was checked line by line:

1. **Domain statement and parameters (`proof.md` lines 41–53):**
   - $`A=\alpha+1`$, $`B=\beta+1`$, $`S=\alpha+\beta+1=A+B-1`$, $`u=t(x-1)/2`$, $`v=t(x+1)/2=u+t`$.
   - Excluding $`A, B, S \notin \{0,-1,-2,\dots\}`$ for the literal coefficients in (1) is necessary and accurate; indeed, our numerical script hit `DivisionUndefined` (`0/0`) on the literal RHS when $`\alpha=-2.3, \beta=-1.7`$ (so $`S=-3`$), while the cancelled form (3) agreed with the LHS to $`2.47 \times 10^{-87}`$, confirming line 45. *(Minor remark, not a gap: the $`B`$ exclusion is also algebraically removable on the LHS via $`(B)_{i+j}/(B)_j=(B+j)_i`$, though excluding $`B \in \mathbb{Z}_{\le 0}`$ is required by DLMF 16.13.4 for the literal $`F_4`$ series.)*
   - For nonterminating $`\gamma`$, $`L = \sqrt{|u|}+\sqrt{|v|} < 1`$ alone already implies $`|t| = |v-u| \le |u|+|v| \le L^2 < 1`$, so the open set $`\{L < 1, 4|u| < |1-t|^2\}`$ stays inside $`|t| < 1`$ and $`\mathrm{Re}(1-t) > 0`$, fixing a unique holomorphic branch of $`(1-t)^{-\gamma}`$ with value $`1`$ at $`t=0`$.

2. **Formal power-series proof (`proof.md` lines 55–147):**
   - **Eq. (2) (lines 57–65):** Splitting $`(S)_{2k} = 4^k (S/2)_k ((S+1)/2)_k`$ and dividing by $`(S)_{m+2n} = (S)_{k+n}`$ (with $`k=m+n \ge n`$) yields $`(S+k+n)_{k-n}`$, a polynomial in $`S`$. Exact.
   - **Eq. (3) (lines 67–76):** With $`w=u/(1-t)^2`$, the bases are $`4w`$ and $`-4tw`$, so $`(4w)^m(-4tw)^n = 4^k w^k (-t)^n`$. Re-indexing by $`k=m+n`$ is valid in $`\mathbb{C}[[u,t]]`$ because $`w`$ has $`u`$-valuation $`1`$, so only finitely many $`(m,n)`$ contribute to $`[u^k t^r]`$. Exact.
   - **Eq. (4) & Chu–Vandermonde (lines 78–102):** Expanding $`(t+u)^j`$, extracting $`[u^k t^r]`$ gives $`j = r+k-i`$ ($`0 \le i \le k`$), $`\binom{j}{r}/j! = 1/(r!(k-i)!)`$, and $`(B)_{k+r}/(B)_{r+k-i} = (B+r+k-i)_i`$. Multiplying the inner sum by $`(A)_k`$ converts it to $`\sum_{i=0}^k \binom{A+k-1}{k-i}\binom{B+r+k-1}{i} = \binom{A+B+r+2k-2}{k} = (S+k+r)_k/k!`$, valid identically as polynomials in $`A, B \in \mathbb{C}`$. Exact.
   - **Eqs. (5)–(6) (lines 104–146):** Using $`(S+k+r)_k/k! = [z^k](1-z)^{-S-k-r}`$ and summing over $`r`$ in $`\mathbb{C}[[z,t]]`$ gives $`[z^k](1-z)^{\gamma-S}(1-t-z)^{-\gamma-k}`$. All factors $`(1-z)`$, $`(1-t/(1-z))`$, $`(1-t-z)`$, $`(1-t)`$, $`(1-q)`$, $`(1-(1-t)q)`$, and $`(1+tq/(1-q))`$ have constant term $`1`$ in the formal power series ring, justifying all power laws. Substituting $`z=(1-t)q`$ multiplies $`[z^k]`$ by $`(1-t)^{-k}`$, which combines with $`(1-t)^{-\gamma-k}`$ to give $`(1-t)^{-\gamma-2k}`$; expanding $`(1+tq/(1-q))^{\gamma-S}`$ with $`\binom{\gamma-S}{n} = (-1)^n (S-\gamma)_n/n!`$ recovers the $`u^k`$ coefficient of (3). Exact.
   - *(Minor remark on lines 144–146: the remark noting equivalence to Euler's transformation DLMF 15.8.1 applied to $`{}_2F_1(\gamma+k, S+2k; S+k; t)`$ is accurate for $`S+k \notin \mathbb{Z}_{\le 0}`$; the report's actual proof via (5)–(6) avoids dividing by $`(S+k)_r`$ and hence has no pole issue there.)*

3. **Passage to convergent series (`proof.md` lines 148–172):**
   - **Open domain (lines 150–159):** At fixed $`x`$, both sides converge uniformly on compact subsets of $`\{L<1, R<1\}`$ (using $`\sup_K |t| \le \sup_K L^2 < 1`$ for uniform boundedness of the $`n`$-sum in (7)), defining holomorphic functions of $`t`$ whose Maclaurin coefficients at $`t=0`$ agree termwise with the formal $`t`$-series (each $`(m,n)`$ summand on the RHS vanishes to order $`m+2n`$ at $`t=0`$, so only finitely many terms contribute to any derivative at $`t=0`$). Along $`t_\lambda = \lambda t`$ ($`\lambda \in [0,1]`$), $`\frac{d}{d\lambda}\bigl(\frac{\lambda}{|1-\lambda t|^2}\bigr) = \frac{1-\lambda^2|t|^2}{|1-\lambda t|^4} > 0`$ for $`\lambda \in [0,1)`$ whenever $`|t| \le 1`$, keeping the segment inside the open domain and connecting $`t`$ to $`0`$. The identity theorem applies rigorously.
   - **Absolute-convergence boundary ($`t \ne 1`$, lines 161–168):** At any boundary point with $`L \le 1, R \le 1, t \ne 1`$, we have $`|t| \le L^2 \le 1`$ and $`\mathrm{Re}(1-t) > 0`$. For $`0 \le \lambda < 1`$, $`L(\lambda)=\sqrt{\lambda}L < 1`$ and $`R(\lambda) < R \le 1`$ (or $`R(\lambda)=0`$ when $`u=0`$), so the open segment lies strictly inside the open domain. Moduli of all arguments ($`u, v, 4w, -4tw`$) are bounded by their endpoint values, so the finite absolute endpoint sums dominate the series along the segment; dominated convergence as $`\lambda \uparrow 1`$ gives equality on the boundary with the same principal branch. Rigorous.
   - **Terminating case $`\gamma = -N`$ (lines 50–53, 169–171):** Both sums truncate at total degree $`N`$, giving a polynomial in $`t`$ on the LHS and a rational function with possible pole only at $`t=1`$ on the RHS; identical Taylor series at $`t=0`$ imply identity as rational functions for $`t \ne 1`$ and a removable singularity at $`t=1`$. Rigorous.

4. **Absolute convergence domain & boundary table (`proof.md` lines 173–278):**
   - **Gamma-ratio bounds (lines 188–193):** Verified via DLMF 5.11.12; when $`S-\gamma \in \mathbb{Z}_{\le 0}`$, the $`n=0`$ factor is still $`1`$ (preserving the lower bound) and every $`n \ge 1`$ bound used in (7) is an upper bound.
   - **Left-hand side (lines 195–227):** The rewrite in terms of $`\binom{k}{i}^2 |u|^i|v|^j`$ and four ratios is exact; Cauchy–Schwarz + binomial theorem settle $`L<1`$ and $`L>1`$. At $`L=1, uv \ne 0`$, the Chernoff tail bound $`[e^{h\epsilon}(1-p+pe^{-h})]^k`$, mode at $`\lfloor (k+1)p \rfloor`$, Stirling upper bound $`O(k^{-1/2})`$, and Chebyshev ($`\ge 3/4`$ mass within $`2\sigma`$) + Cauchy–Schwarz lower bound give $`\sum P_i^2 \asymp k^{-1/2}`$ and degree-$`k`$ sum $`\asymp k^{-d-1/2}`$, summable iff $`d > 1/2`$. The axis cases $`u=0`$ ($`\mathrm{Re}\,\gamma < 0`$) and $`v=0`$ ($`d > \mathrm{Re}\,B`$) are exact.
   - **Right-hand side and combined boundary (lines 229–277):** Note that (2) transforms each single $`(m,n)`$ term without inter-term cancellation, so absolute convergence of (1) equals absolute convergence of (3). The uniform estimate $`\frac{|(S+k+n)_{k-n}|}{(k-n)!} \asymp \binom{2k}{k-n}`$ holds because the ratio is $`\prod_{j=k+n}^{2k-1}|1+(S-1)/(j+1)|`$ and $`\sum_{j=k}^{2k-1} |S-1|/(j+1) < |S-1|`$. For $`|t|<1`$, $`\sum_{n \ge 0} \frac{|(S-\gamma)_n|}{n!}|t|^n < \infty`$ settles $`R<1`$, $`R>1`$, and $`R=1`$ ($`d>1/2`$).
   - **Corner case $`|t|=1, L \le 1`$ (lines 255–277):** Equality in $`1=|t| \le |u|+|v| \le L^2 \le 1`$ forces $`uv=0`$. If $`u=0`$ then $`R=0<1`$. If $`v=0`$ then $`|u|=1`$ and $`R = 4/|1-t|^2 \ge 1`$ with equality iff $`t=-1`$ ($`x=-1`$). There, both $`L=1, v=0`$ ($`d > \mathrm{Re}\,B`$) and $`R=1`$ ($`d > 1/2`$) apply. The bound $`\binom{2k}{k-n}/\binom{2k}{k} = \prod_{j=1}^n (1-\frac{2j-1}{k+j}) \le \exp(-n^2/(2k))`$ is exact (since $`k+j \le 2k`$ and $`\sum_{j=1}^n(2j-1)=n^2`$). Summing $`(n+1)^{\delta-1}\exp(-n^2/(2k))`$ with $`\delta = d+\mathrm{Re}\,B - 1`$ over $`n \le \sqrt{k}`$ and dyadic bands $`2^j\sqrt{k} < n \le 2^{j+1}\sqrt{k}`$ yields $`O(k^{-d-1/2+\delta/2}) = O(k^{-1-(d-\mathrm{Re}\,B)/2})`$ for $`\delta>0`$ (summable since $`d>\mathrm{Re}\,B`$) and $`O(k^{-d-1/2}\log(k+2))`$ for $`\delta \le 0`$ (summable since $`d>1/2`$). Necessity of both $`d>\mathrm{Re}\,B`$ (from LHS) and $`d>1/2`$ (from RHS $`n=0`$) holds, as confirmed empirically in our boundary diagnostics.
   - **Note on the RHS-alone domain when $`|t|>1`$:** If one ignored the LHS condition $`L \le 1`$, the RHS series with $`|t|>1`$ would have exponential growth rate $`((|t|+1)^2/|t|)|w|`$ rather than $`4|w|`$ (attained near $`n/k = (|t|-1)/(|t|+1)`$). Because `proof.md` only characterizes the *common* domain where $`L \le 1`$ already forces $`|t| \le 1`$, this does not affect any claim in the report.

5. **Alternative derivation via DLMF 16.16.5_5 and Horn's $`H_3`$ (`proof.md` lines 343–381):**
   - Verified algebraically from `16.16.E5_5.tex` and `16.16.E8.tex`: $`X=(1-t)\zeta`$, $`Y=X+t`$, $`\zeta(1-\zeta)=w`$ give $`X(1-Y)=u`$, $`Y(1-X)=v`$, and reduce $`F_4(\gamma,B;A,B;u,v)`$ to $`(1-Y)^{-\gamma} F_1(\gamma; S-\gamma, \gamma+1-A; A; X, X/(1-Y))`$.
   - The subsequent $`F_1 \to H_3`$ reduction follows directly from applying the quadratic $`{}_2F_1`$ transformation (DLMF 15.8.13) to the inner $`{}_2F_1(\gamma+q, \gamma+1-A; A+q; b)`$ series, yielding $`(1-t)^{-\gamma} H_3(\gamma, S-\gamma; A; w, (1-t)w)`$, and the binomial convolution in lines 373–376 recovers (3) exactly. We also verified the $`H_3`$ equality directly as an exact rational $`t`$-series through $`t^{16}`$. *(Caveat: Ernst & Karlsson's 2025 Authorea preprint could not be opened, so its equation numbers (22), (42)–(43) were not directly checked against that PDF.)*

---

### 1.3 Independent Recomputation & Planted-Error Tests

Two standalone pure-Python-3 standard-library scripts were written from scratch directly from `cohl-f4.jpg` (without reading or reusing the proof's own check programs):

- **Exact rational verification:** `checks/verify_gemini_exact.py` (output in `checks/gemini-exact-output.txt`), using exact `fractions.Fraction` arithmetic on the literal uncancelled formula in (1):
  - **11 parameter sets** checked through $`t^{20}`$ (and $`t^{18}`$ for the near-pole case), including generic positive rationals, awkward negative non-integer $`(\alpha,\beta,\gamma)`$ such as $`(-7/3,-5/2,-11/5)`$ and $`(-17/5,-13/3,-9/2)`$, special values $`x \in \{-1, 0, 1, 5, 13/4\}`$, terminating $`\gamma = -5`$ (also verifying all coefficients for $`r>5`$ vanish identically), specialization $`\gamma=\alpha+\beta+1`$, and near-pole parameters $`\alpha=-1999/1000, \beta=-2990/997`$. **All 11 passed with zero residual.**
  - **Exact polynomial-in-$`x`$ check:** verified all 153 exact coefficients $`[t^r x^d]`$ ($`0 \le d \le r \le 16`$) for $`\alpha=-7/3, \beta=-5/2, \gamma=-11/5`$. **Zero residual.**
  - **Horn $`H_3`$ identity check:** verified `proof.md` lines 360–363 through $`t^{16}`$ for $`\alpha=-7/3, \beta=-5/2, \gamma=-11/5, x=-2/3`$. **Zero residual.**
  - **4 planted deliberate errors** (tested at $`\alpha=-7/3, \beta=-5/2, \gamma=-11/5, x=-2/3`$) — all 4 caught immediately:
    1. `sign_flip` (drop minus sign in second base: $`-2t^2(x-1) \to +2t^2(x-1)`$): caught at $`t^2`$, exact residual $`\text{LHS}-\text{RHS} = -539/120`$.
    2. `denom_m_plus_n` (subscript typo $`(\alpha+\beta+1)_{m+2n} \to (\alpha+\beta+1)_{m+n}`$): caught at $`t^2`$, exact residual $`-12397/1440`$.
    3. `shift_param` (shift typo $`(\alpha+\beta-\gamma+1)_n \to (\alpha+\beta-\gamma)_n`$): caught at $`t^2`$, exact residual $`11/8`$.
    4. `half_param_shift` (half-parameter typo $`(\frac{1}{2}(\alpha+\beta+2))_{m+n} \to (\frac{1}{2}(\alpha+\beta+3))_{m+n}`$): caught at $`t^1`$, exact residual $`11/8`$.

- **High-precision numerical & boundary verification:** `checks/verify_gemini_numeric.py` (output in `checks/gemini-numeric-output.txt`), using 110-digit `decimal.Decimal` complex arithmetic (self-tested against the binomial series for $`(1-t)^{-\gamma}`$ to $`2.98 \times 10^{-107}`$):
  - **7 interior test points** (positive real, negative non-integer real, complex parameters/arguments, closer-to-boundary complex with $`L=0.7845, R=0.6562`$, and $`x = 0, -1, +1`$) evaluated at two total-degree cutoffs $`K_1 \to K_2`$: **all 7 passed** with maximum scaled discrepancy $`|\text{LHS}-\text{RHS}|/\max(1,|\text{LHS}|,|\text{RHS}|) = 1.101 \times 10^{-85}`$ (down to $`3.155 \times 10^{-108}`$).
  - **Removable $`S \in \mathbb{Z}_{\le 0}`$ pole check:** at $`\alpha=-2.3, \beta=-1.7`$ ($`S=-3`$), the literal RHS raises `InvalidOperation` ($`0/0`$), while the cancelled form (3) matches LHS to $`2.467 \times 10^{-87}`$.
  - **2 planted errors in high-precision numerics** (complex interior point): `sign_flip` produces discrepancy $`8.003 \times 10^{-2}`$; `denom_m_plus_n` produces discrepancy $`1.000`$. Both caught cleanly.
  - **Boundary & divergence diagnostics:**
    - $`L=R=1, uv \ne 0`$ ($`x=5/4, t=1/2`$) with $`d=5.5>1/2`$: LHS and RHS converge together (discrepancy $`1.39 \times 10^{-5}`$ at $`K=120`$, degree-$`k`$ terms decaying as $`k^{-6}`$).
    - Corner $`x=-1, t=-1`$ ($`L=R=1, v=0`$) with $`d=5.5 > \max(\mathrm{Re}\,B=1.5, 1/2)`$: both converge to Gauss's exact value $`\frac{8192}{2205\pi} \approx 1.18258258`$.
    - Corner $`x=-1, t=-1`$ with $`A=0.75, B=-0.5, \gamma=0.5`$ ($`d=0.25`$: satisfies $`d>\mathrm{Re}\,B`$ but violates $`d>1/2`$): LHS degree-$`k`$ ratio $`k=120/k=60`$ is $`0.2964 \approx 2^{-1.75}=0.2973`$ (converges), while RHS ratio is $`0.5972 \approx 2^{-0.75}=0.5946`$ (diverges), confirming both table rows are needed.
    - Outside common domain $`x=-1, t=0.5`$ ($`L=0.7071<1, R=8>1`$): LHS degree-60 term is $`6.30 \times 10^{-20}`$ while RHS degree-60 term blows up to $`3.35 \times 10^{51}`$.
    - **Conditional-convergence examples** (relevant to the reply text below):
      - $`x=1, t=-1, \gamma=1/2`$ ($`L=1, u=0, \mathrm{Re}\,\gamma=1/2 \not< 0`$): LHS absolute sum diverges ($`\sum k^{-1/2}`$), yet the alternating LHS series converges conditionally to $`1/\sqrt{2} = \text{RHS}`$.
      - $`x=-1, t=3-2\sqrt{2}, \alpha=0, \beta=0, \gamma=1/2`$ ($`L=0.4142<1, R=1, d=1/2`$): LHS converges absolutely to $`\cos(\pi/8) \approx 0.9238795325`$, while RHS absolute sum diverges logarithmically ($`\sim c/k`$), yet its total-degree partial sums alternate around and converge conditionally to LHS ($`K=150`$: $`0.92495`$, $`K=151`$: $`0.92281`$, average $`0.923883`$).

---

## 2. Does the Plan Make Sense?

**Replying makes sense, but the reply text should be changed (`SENSE: CHANGE TEXT`).**

### 2.1 Is the Identity Already in Print?

- **Component pieces are classical:**
  - The LHS is a classical Jacobi-polynomial generating function, $`F_4(\gamma,\beta+1;\alpha+1,\beta+1; t(x-1)/2, t(x+1)/2) = \sum_{j\ge 0} \frac{(\gamma)_j}{(\alpha+1)_j} P_j^{(\alpha,\beta)}(x) t^j`$ (e.g., Chan–Chen–Chyan–Srivastava 2004, eq. (71), and earlier Jacobi generating-function literature).
  - The specialization $`\gamma = \alpha+\beta+1`$ (where the $`n`$-sum collapses to $`n=0`$) is classical: DLMF 18.12.3 (after $`(x,\alpha,\beta,t) \mapsto (-x,\beta,\alpha,-t)`$) and Cohl & MacKenzie (2013, `arXiv:1210.0039`, eq. (7)). Cohl's 2023 post is clearly a general-$`\gamma`$ extension of that formula.
  - The repeated-parameter reduction $`F_4(a,b;c,b; \dots) \to F_1`$ is classical (Bailey 1933; Schlosser survey eq. (43); Shpot 2007 eq. (21); added to DLMF as `16.16.5_5` in V1.1.12 on 2023-12-15, three months *after* Cohl's 2023-09-08 post), and the further reduction to Horn's $`H_3`$ is also classical (via DLMF 15.8.13; Erdélyi 1948; Srivastava & Panda 1974; Srivastava & Manocha 1984 p. 57 as cited by Shpot 2007; Ernst & Karlsson 2025 preprint).
- **The exact posted double series was not found in accessible print:** Searching DLMF (§16.16, §18.12), all 13 of Cohl's 2023–2026 arXiv preprints (full LaTeX source), Cohl's personal bibliography (`hcohl.sdf.org/bibliography.html`), NIST publication list, Costas-Santos's bibliography, arXiv API, zbMATH Open API, Schlosser's survey, Shpot (2007), Srivastava (1973), Ernst & Karlsson's open-access 2024–2025 papers, and MathOverflow/Math.SE turned up **no prior display of Cohl's exact Srivastava–Daoust RHS** with $`(\alpha+\beta+1)_{m+2n}`$ in the denominator and arguments $`\frac{2t(x-1)}{(t-1)^2}, -\frac{2t^2(x-1)}{(t-1)^2}`$.
- **Honest caveat:** Several relevant books and articles could not be opened in full text (listed in §3 below, notably Srivastava–Manocha 1984, Srivastava–Karlsson 1985, Srivastava–Panda 1974, Chan et al. 2004, and Ernst–Karlsson's 2025 Authorea preprint), so prior appearance in a non-open-access source cannot be completely ruled out. The reply should not claim novelty — and neither the draft nor the revised text does.

### 2.2 Is It Useful to Cohl, or Too Trivial?

**Useful, not silly.** Cohl explicitly asked *"Here's a conjecture. Can anyone prove it?"* on 2023-09-08, and public post metadata (`api.fxtwitter.com`) still shows **0 replies and 0 quotes** as of 2026-11-21. Although the formal proof is short and elementary (~half a page), finding the cleanest direct path (Chu–Vandermonde + Euler's transformation), nailing down the exact absolute-convergence boundary including the subtle corner $`x=t=-1`$, and tying the series to DLMF 16.16.5_5 (added after his tweet) and Horn's $`H_3`$ is genuinely helpful to a NIST DLMF editor. It would only look silly if hyped.

### 2.3 Problems with the Planned Reply Text

Original draft:
> It's true: as a formal power series in t, and wherever both sides converge. One Vandermonde sum, then Euler's transformation (DLMF 15.8.1). Proof and code you can run: <link to this repository>
> Written by GPT-6.2 Sol, checked independently by Claude Opus 5.5 and GLM 5.3.

Three issues require changing it:

1. **Mathematical overclaim (`"wherever both sides converge"`):** The report only proves equality **wherever both sides converge absolutely** (with $`t \ne 1`$). For double series, bare convergence depends on summation order and can hold conditionally outside the proven domain — e.g. at $`(x,t)=(1,-1), \gamma=1/2`$ and at $`(x,t)=(-1, 3-2\sqrt{2}), \alpha=\beta=0, \gamma=1/2`$, both sides converge (the latter by total-degree partial sums) while one side fails absolute convergence, so `proof.md`'s dominated-convergence proof does not cover them. Adding `"for t ≠ 1 wherever both sides converge absolutely"` makes the claim exact.
2. **Brittle / promotional model roll call:** `"Written by GPT-6.2 Sol, checked independently by Claude Opus 5.5 and GLM 5.3."` reads like AI-benchmark marketing under a specialist's math post and hardcodes a specific reviewer list that does not reflect all actual checks (or future edits). A plain disclosure (`"AI-drafted proof and checks"`) is transparent without hype; full provenance belongs in the linked repository README.
3. **Raw length and missed DLMF/Horn pointer:** The original draft is **290 raw characters** (269 X t.co-weighted). Trimming the model roll call frees space to stay strictly under 280 raw characters while adding a high-signal note for Cohl (the link to DLMF 16.16.5_5 and Horn's $`H_3`$).

### 2.4 Proposed Replacement Text

> It holds as a formal power series in t, and for t ≠ 1 wherever both sides converge absolutely: one Vandermonde sum, then Euler's transformation (DLMF 15.8.1); also via DLMF 16.16.5_5 and Horn's H3. AI-drafted proof and checks: <link to this repository>

Length: **271 raw characters** (250 X t.co-weighted characters). Plain, mathematically exact, transparent about AI authorship, and no hype.

---

## 3. Sources Consulted

### Opened and Read

- Local image the repository's `cohl-f4.jpg` and original image `https://pbs.twimg.com/media/F5iLOskbAAA0o8A.jpg?name=orig` — **opened** (byte-identical SHA-256).
- Public X post metadata `https://api.fxtwitter.com/hcohl001/status/1700261948808700126` — **opened** (post text, date 2023-09-08, 0 replies, 0 quotes).
- NIST DLMF equations & pages:
  - `https://dlmf.nist.gov/16.13` and `https://dlmf.nist.gov/16.13.E4.tex` — **opened**.
  - `https://dlmf.nist.gov/15.8.E1.tex` — **opened**.
  - `https://dlmf.nist.gov/5.11.E12.tex` — **opened**.
  - `https://dlmf.nist.gov/16.16`, `https://dlmf.nist.gov/16.16.E1.tex`–`16.16.E10.tex` (incl. `16.16.E5_5.tex`, `16.16.E8.tex`) — **opened**.
  - `https://dlmf.nist.gov/18.12`, `https://dlmf.nist.gov/18.12.E1.tex`–`18.12.E5.tex` (incl. `18.12.E2_5.tex`, `18.12.E3_5.tex`) — **opened**.
  - `https://dlmf.nist.gov/errata/` (V1.1.12, 2023-12-15) — **opened**.
- Cohl publication lists & 2023–2026 papers:
  - `https://arxiv.org/a/cohl_h_1.html` — **opened**.
  - `https://www.nist.gov/publications/search_by_author/1159096?page=0` — **opened**.
  - `http://hcohl.sdf.org/` and `http://hcohl.sdf.org/bibliography.html` — **opened**.
  - `https://rscosan.com/papers` — **opened**.
  - Full LaTeX sources (`https://arxiv.org/e-print/<id>`) of Cohl's 13 arXiv preprints from 2023–2026: `2303.02235`, `2303.13680`, `2306.03035`, `2307.04884`, `2308.13652`, `2410.00246`, `2411.03571`, `2411.03574`, `2411.09040`, `2504.21782`, `2508.08085`, `2508.08162`, `2602.15824` — **opened**.
  - H. S. Cohl & C. MacKenzie (2013), `https://arxiv.org/e-print/1210.0039` and `https://arxiv.org/pdf/1210.0039` — **opened**.
- Related hypergeometric literature:
  - M. A. Shpot (2007), `https://arxiv.org/e-print/0711.2742` and `https://arxiv.org/pdf/0711.2742` — **opened**.
  - M. J. Schlosser, *Multiple hypergeometric series — Appell series and beyond*, `https://www.mat.univie.ac.at/~schlosse/pdf/appell-survey.pdf` — **opened**.
  - H. M. Srivastava (1973), *On the reducibility of Appell's function F₄*, `https://doi.org/10.4153/CMB-1973-049-9` (Cambridge Core OA PDF) — **opened**.
  - T. Ernst & P. W. Karlsson (2025), *Reduction and summation formulas for Appell functions in the spirit of Burchnall–Chaundy*, `https://doi.org/10.7153/jca-2025-25-04` / `http://files.ele-math.com/articles/jca-25-04.pdf` — **opened**.
  - P. W. Karlsson & T. Ernst (2024), *Applications of quadratic and cubic hypergeometric transformations*, `https://doi.org/10.17951/a.2024.78.1.37-73` / `https://journals.umcs.pl/a/article/download/18000/11666` — **opened**.
  - arXiv LaTeX sources (`https://arxiv.org/e-print/<id>`) for `2010.16093`, `0804.0655`, `2005.07170` — **opened**.
- Search APIs:
  - arXiv API `https://export.arxiv.org/api/query` — **opened**.
  - zbMATH Open API `https://api.zbmath.org/v1/document/` (records `3421169`, `3849611`, `3879293`, `2085965`, `3047921`, `7540065`, `7905942`, `7960861`, `7904194`, `8022932`, and Cohl 2023–2026 records) — **opened** (bibliographic metadata/open summaries).
  - StackExchange API `https://api.stackexchange.com/2.3/search/advanced` (MathOverflow and Math.SE) — **opened**.

### Not Read (with Reason)

- `https://x.com/hcohl001/status/1700261948808700126` and `https://x.com/hcohl001/status/1700261948808700126/quotes` — **not read** (X requires login/JS rendering; used `api.fxtwitter.com` snapshot instead, so live quote-tweets cannot be independently confirmed).
- T. Ernst & W. P. Karlsson (2025-02-17 preprint), `https://doi.org/10.22541/au.173980497.70111458/v1` (`https://www.authorea.com/doi/full/10.22541/au.173980497.70111458/v1`) and `https://d197for5662m48.cloudfront.net/documents/publicationstatus/245830/preprint_pdf/793746cb6d0de25bbcca2c3cbce48018.pdf` — **not read** (Authorea returned HTTP 403; CloudFront host failed DNS resolution, curl exit code 6).
- H. M. Srivastava & H. L. Manocha (1984), *A Treatise on Generating Functions*, `https://zbmath.org/3849611` — **not read** (copyrighted monograph, no open-access full text).
- H. M. Srivastava & P. W. Karlsson (1985), *Multiple Gaussian Hypergeometric Series*, `https://zbmath.org/3879293` — **not read** (copyrighted monograph, no open-access full text).
- H. M. Srivastava & R. Panda (1974), *Some hypergeometric transformations involving Horn's function H₃*, `https://zbmath.org/3421169` — **not read** (full text unavailable online; only zbMATH bibliographic record opened).
- W.-C. C. Chan, K.-Y. Chen, C.-J. Chyan & H. M. Srivastava (2004), `https://doi.org/10.1016/j.jmaa.2004.01.015` — **not read** (paywalled on ScienceDirect; zbMATH summary restricted).
- A. Erdélyi (1948), `https://zbmath.org/3047921`, and A. Erdélyi et al. (1953), *Higher Transcendental Functions*, Vol. 1 (`https://authors.library.caltech.edu/records/cnd32-h9x80/files/Volume_1.pdf`) — **not read** (1948 article paywalled; Caltech S3 URL returned HTTP 403 on HEAD and full scanned volume was not downloaded).
- M. I. Qureshi et al. (2023), `https://www.mdpi.com/2504-3110/7/10/700`; H. Liu & W. Wang (2014), `https://www.sciencedirect.com/science/article/pii/S0022247X13006264`; Y. Vyas & K. Fatawat (2024), `https://link.springer.com/chapter/10.1007/978-3-031-56307-2_15` — **not read** (MDPI returned HTTP 403; ScienceDirect and Springer paywalled).

---

## 4. Files Produced

- `checks/verify_gemini_exact.py` — independent exact rational (`fractions.Fraction`) verifier + symbolic polynomial-in-$`x`$ check + Horn $`H_3`$ check + 4 planted-error checks.
- `checks/gemini-exact-output.txt` — recorded output of `verify_gemini_exact.py`.
- `checks/verify_gemini_numeric.py` — independent 110-digit `decimal.Decimal` complex numeric verifier + removable $`S`$-pole check + boundary/divergence/conditional-convergence diagnostics + 2 planted-error checks.
- `checks/gemini-numeric-output.txt` — recorded output of `verify_gemini_numeric.py`.
