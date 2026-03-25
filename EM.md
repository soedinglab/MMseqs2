# EM Formulation in `reclassify_taxonomy.cpp`

## Overview

`reclassify_taxonomy.cpp` reorders MMseqs2 alignment hits by estimating posterior support for each target using:

- alignment score
- target abundance
- entropy-based coverage penalty

The method is an EM-like iterative update accelerated with SQUAREM.

## Notation

Let:

- \( q \): a query
- \( H_q \): the set of hits for query \(q\)
- \( h \in H_q \): one hit
- \( t(h) \): the target of hit \(h\)
- \( Q \): total number of queries
- \( s_h \): alignment score of hit \(h\)
- \( A_t \): abundance of target \(t\)
- \( E_t \): entropy penalty of target \(t\)
- \( \lambda, \alpha, \gamma \): model parameters

## 1. Initial Abundance

For each query \(q\), compute the sum of scores over all hits:

\[
S_q = \sum_{h \in H_q} s_h
\]

Each hit contributes a normalized score fraction to its target:

\[
c_h = \frac{s_h}{S_q}
\]

Initial abundance of target \(t\) is:

\[
A_t^{(0)} = \frac{1}{Q} \sum_q \sum_{h \in H_q,\ t(h)=t} \frac{s_h}{S_q}
\]

## 2. Coverage Weight for Entropy

For each hit \(h\), define target-covered length:

\[
\ell_h = dbEndPos_h - dbStartPos_h + 1
\]

The code assigns each hit a coverage weight:

\[
m_h = \frac{e^{\lambda s_h}}{\ell_h}
\]

For each covered target position \(p\), accumulate:

\[
C_t(p) = \sum_{h \text{ covering } p} m_h
\]

## 3. Position Probability on a Target

Let total coverage mass on target \(t\) be:

\[
Z_t = \sum_p C_t(p)
\]

Then the normalized positional probability is:

\[
p_t(p) = \frac{C_t(p)}{Z_t}
\]

## 4. Target Entropy

Coverage entropy for target \(t\) is:

\[
H_t = - \sum_p p_t(p)\log_2 p_t(p)
\]

If \(Z_t = 0\), then:

\[
H_t = 0
\]

## 5. Entropy Penalty

Let the total entropy over all targets be:

\[
H_{\mathrm{sum}} = \sum_t H_t
\]

Then the penalty used in scoring is:

\[
E_t =
\begin{cases}
1 - \frac{H_t}{H_{\mathrm{sum}}}, & H_{\mathrm{sum}} > 0 \\
0, & H_{\mathrm{sum}} = 0
\end{cases}
\]

## 6. Score Term for One Hit

For each query \(q\), define:

\[
S_{\max,q} = \max_{h \in H_q} s_h
\]

In the current implementation, this is the maximum observed bit score among the query's candidate hits.

With

\[
\varepsilon = 10^{-12}
\]

the abundance and entropy penalty are clipped from below:

\[
A_h = \max(A_{t(h)}, \varepsilon)
\]

\[
E_h = \max(E_{t(h)}, \varepsilon)
\]

The score term is:

\[
f_h =
\exp\left(\lambda \frac{s_h}{S_{\max,q}}\right)
\cdot A_h^{\alpha}
\cdot E_h^{\gamma}
\]

## 7. E-step: Posterior Probability

For each query \(q\), normalize the score terms over all its hits:

\[
Z_q = \sum_{h \in H_q} f_h
\]

Then the posterior probability of hit \(h\) is:

\[
P(h \mid q) =
\begin{cases}
\frac{f_h}{Z_q}, & Z_q > 0 \\
0, & Z_q = 0
\end{cases}
\]

## 8. M-step: Abundance Update

The updated abundance of target \(t\) is the average posterior mass assigned to it:

\[
A_t^{\mathrm{new}} =
\frac{1}{Q}
\sum_q
\sum_{h \in H_q,\ t(h)=t}
P(h \mid q)
\]

This is the EM abundance update used in the code.

## 9. Log-Likelihood

For each query \(q\), the code evaluates:

\[
\ell_q =
\sum_{h \in H_q} P(h \mid q)\log f_h
-
\log\left(\sum_{h \in H_q} f_h\right)
\]

The average log-likelihood is:

\[
LL = \frac{1}{Q} \sum_q \ell_q
\]

This value is used to monitor whether the accelerated update is acceptable.

## 10. Simplex Projection

After acceleration, abundance values may become negative or fail to sum to 1.

The code projects back to the simplex by:

1. clipping negatives:

\[
x_i' = \max(x_i, 0)
\]

2. renormalizing if the total is positive:

\[
x_i'' = \frac{x_i'}{\sum_j x_j'}
\]

## 11. SQUAREM Acceleration

Let:

\[
x_1 = EM(x_0), \qquad x_2 = EM(x_1)
\]

Define:

\[
r = x_1 - x_0
\]

\[
v = x_2 - x_1 - r
\]

Their Euclidean norms are:

\[
\|r\| = \sqrt{\sum_i r_i^2}, \qquad
\|v\| = \sqrt{\sum_i v_i^2}
\]

The acceleration factor is:

\[
a =
\begin{cases}
-1, & \|v\| = 0 \\
-\frac{\|r\|}{\|v\|}, & \|v\| > 0
\end{cases}
\]

Then it is clipped to:

\[
a \in [-1, 1]
\]

The accelerated update is:

\[
x_{\mathrm{new}} = x_0 - 2ar + a^2 v
\]

Finally, `projectSimplex(x_new)` is applied.

## 12. Likelihood Safeguard

If the accelerated update lowers the log-likelihood, the code falls back to the second EM step:

\[
x_{\mathrm{new}} \leftarrow x_2
\]

This keeps the iteration stable.

## 13. Convergence Criterion

The stopping criterion is the maximum absolute parameter change:

\[
\Delta = \max_i |x_{\mathrm{new}, i} - x_{0,i}|
\]

Convergence is declared when:

\[
\Delta < \mathrm{tol}
\]

after at least 6 iterations.

## 14. Final Ranking

After convergence, hits for each query are sorted by posterior probability:

\[
h_i \prec h_j \quad \text{if} \quad P(h_i \mid q) > P(h_j \mid q)
\]

If two hits have equal posterior, MMseqs2's default hit comparison is used as a tie-breaker.

## Summary

The model implemented in `reclassify_taxonomy.cpp` is:

\[
f_h =
\exp\left(\lambda \frac{s_h}{S_{\max,q}}\right)
\cdot
A_{t(h)}^{\alpha}
\cdot
E_{t(h)}^{\gamma}
\]

with EM updates:

\[
P(h \mid q) = \frac{f_h}{\sum_{h' \in H_q} f_{h'}}
\]

\[
A_t^{\mathrm{new}} =
\frac{1}{Q}
\sum_q
\sum_{h \in H_q,\ t(h)=t}
P(h \mid q)
\]

and SQUAREM acceleration applied to the abundance vector.
