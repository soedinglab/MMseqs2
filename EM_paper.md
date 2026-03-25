# A Paper-Style Description of the EM Reclassification Model

## Abstract

The reclassification procedure implemented in `reclassify_taxonomy.cpp` refines target ranking for each query by combining local alignment evidence with global target-level support. The model assigns each hit a latent posterior probability and iteratively updates target abundances using an expectation-maximization framework. To improve convergence speed, the abundance update is accelerated using SQUAREM.

## Model Setup

Let \(q \in \{1, \dots, Q\}\) index queries, and let \(H_q\) denote the set of candidate hits for query \(q\). For each hit \(h \in H_q\), define:

- \(t(h)\): target assigned to hit \(h\)
- \(s_h\): alignment score

For each target \(t\), the model estimates:

- abundance \(A_t\)
- entropy-derived penalty \(E_t\)

The free model parameters are:

\[
\lambda, \alpha, \gamma
\]

which control the relative effect of alignment score, abundance, and entropy penalty.

## Hit Likelihood Term

For each query \(q\), define the query-specific normalization term:

\[
S_{\max,q} = \max_{h \in H_q} s_h
\]

In the current implementation, this is the maximum observed bit score among the candidate hits of query \(q\).

For each hit \(h\), abundance and entropy penalty are clipped from below:

\[
\tilde{A}_h = \max(A_{t(h)}, \varepsilon)
\]

\[
\tilde{E}_h = \max(E_{t(h)}, \varepsilon)
\]

with

\[
\varepsilon = 10^{-12}
\]

The unnormalized compatibility score of hit \(h\) is:

\[
f_h =
\exp\left(\lambda \frac{s_h}{S_{\max,q}}\right)
\cdot \tilde{A}_h^{\alpha}
\cdot \tilde{E}_h^{\gamma}
\]

This quantity defines the relative support for assigning query \(q\) to target \(t(h)\).

## Entropy-Based Penalty

For a hit \(h\), let its target-covered length be:

\[
\ell_h = dbEndPos_h - dbStartPos_h + 1
\]

The hit contributes position-wise target coverage mass:

\[
m_h = \frac{\exp(\lambda s_h)}{\ell_h}
\]

For target \(t\), the accumulated coverage at target position \(p\) is:

\[
C_t(p) = \sum_{h \text{ covering } p,\ t(h)=t} m_h
\]

The normalized positional distribution is:

\[
p_t(p) = \frac{C_t(p)}{\sum_{p'} C_t(p')}
\]

The target entropy is then:

\[
H_t = - \sum_p p_t(p)\log_2 p_t(p)
\]

Let

\[
H_{\mathrm{sum}} = \sum_t H_t
\]

The penalty assigned to target \(t\) is:

\[
E_t =
\begin{cases}
1 - \frac{H_t}{H_{\mathrm{sum}}}, & H_{\mathrm{sum}} > 0 \\
0, & H_{\mathrm{sum}} = 0
\end{cases}
\]

This term modulates target support according to the spatial distribution of its matched regions.

## Initialization

The abundance is initialized by normalized query-local alignment support. For each query \(q\), define:

\[
S_q = \sum_{h \in H_q} s_h
\]

Then the initial abundance of target \(t\) is:

\[
A_t^{(0)} =
\frac{1}{Q}
\sum_q
\sum_{h \in H_q,\ t(h)=t}
\frac{s_h}{S_q}
\]

This initialization gives each query unit mass, distributed proportionally to hit score.

## E-step

Given the current abundances \(A_t\), posterior probabilities are assigned to hits within each query:

\[
P(h \mid q) =
\frac{f_h}{\sum_{h' \in H_q} f_{h'}}
\]

whenever the denominator is positive; otherwise the posterior is set to zero.

Thus, the E-step computes a soft assignment of each query to its candidate targets.

## M-step

The target abundance is updated by averaging posterior mass across all queries:

\[
A_t^{\mathrm{new}} =
\frac{1}{Q}
\sum_q
\sum_{h \in H_q,\ t(h)=t}
P(h \mid q)
\]

The entropy penalty is kept fixed during the EM iteration after its initial computation.

## Objective Function

The code evaluates an average query-normalized log-likelihood-like objective:

\[
\ell_q =
\sum_{h \in H_q} P(h \mid q)\log f_h
-
\log\left(\sum_{h \in H_q} f_h\right)
\]

and

\[
LL = \frac{1}{Q}\sum_q \ell_q
\]

This objective is used to reject unstable accelerated updates.

## SQUAREM Acceleration

Let \(EM(x)\) denote one EM abundance update applied to abundance vector \(x\). Starting from \(x_0\), compute:

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

with Euclidean norms:

\[
\|r\| = \sqrt{\sum_i r_i^2}, \qquad
\|v\| = \sqrt{\sum_i v_i^2}
\]

The acceleration parameter is:

\[
a =
\begin{cases}
-1, & \|v\| = 0 \\
-\frac{\|r\|}{\|v\|}, & \|v\| > 0
\end{cases}
\]

and is clipped into the interval:

\[
a \in [-1,1]
\]

The accelerated proposal is:

\[
x_{\mathrm{new}} = x_0 - 2ar + a^2 v
\]

Since \(x_{\mathrm{new}}\) must remain a valid abundance vector, it is projected back onto the simplex.

## Simplex Projection

Given a proposed abundance vector \(x\), projection is implemented as:

1. clip negative entries:

\[
x_i' = \max(x_i, 0)
\]

2. normalize:

\[
x_i'' = \frac{x_i'}{\sum_j x_j'}
\]

if the denominator is positive.

This ensures abundances are nonnegative and sum to 1.

## Safeguard Step

If the accelerated proposal decreases the objective, the algorithm falls back to the second ordinary EM iterate:

\[
x_{\mathrm{new}} \leftarrow x_2
\]

This provides a monotonicity safeguard against unstable extrapolation.

## Convergence Criterion

The iteration stops when the maximum coordinate-wise parameter change becomes sufficiently small:

\[
\Delta = \max_i |x_{\mathrm{new},i} - x_{0,i}|
\]

Convergence is declared when:

\[
\Delta < \mathrm{tol}
\]

after an initial burn-in of at least six iterations.

## Final Ranking

After convergence, hits are sorted by posterior probability:

\[
h_i \prec h_j
\quad \Longleftrightarrow \quad
P(h_i \mid q) > P(h_j \mid q)
\]

with MMseqs2's default hit comparison used to break ties.

## Final Model Summary

The reclassification model can be summarized as:

\[
f_h =
\exp\left(\lambda \frac{s_h}{S_{\max,q}}\right)
\cdot
A_{t(h)}^{\alpha}
\cdot
E_{t(h)}^{\gamma}
\]

with posterior assignment:

\[
P(h \mid q) =
\frac{f_h}{\sum_{h' \in H_q} f_{h'}}
\]

and abundance update:

\[
A_t^{\mathrm{new}} =
\frac{1}{Q}
\sum_q
\sum_{h \in H_q,\ t(h)=t}
P(h \mid q)
\]

This framework combines local alignment quality with global target prevalence to produce a query-specific posterior ranking of targets.
