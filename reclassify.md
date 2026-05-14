# Reclassify EM Notes

This document summarizes the equations and outputs implemented in `src/util/EM_reclassify.cpp`.

What reclassify adds:
- It does not treat each query in isolation.
- It estimates target abundance over all matched queries.
- It uses that abundance to boost hits to globally supported targets and downweight isolated or ambiguous hits.
- It penalizes high-entropy or diffuse mappings.
- It iterates this with an EM-style procedure, accelerated with SQUAREM, until convergence.

You can use `mmseqs reclassify` after running `mmseqs search` and (optional) `mmseqs align -a`.

## Example usage (search -> align -> reclassify)
- mmseqs createdb query.fasta queryDB
- mmseqs search queryDB targetDB alignDB tmp -a
- mmseqs reclassify queryDB targetDB alignDB newDB
- mmseqs convertalis queryDB targetDB newDB reclassify_result.m8
- mmseqs abundance queryDB targetDB newDB abundance.report (--taxonomy 1)

Notes:
- `newDB` is an alignment DB in the same format as `alignDB`, but ordered by posterior probability per query.
- The posterior probability replaces the `seqId` field in `newDB` to keep the column count unchanged.

## Notation

- Query index: $q$
- Target/protein index: $t$
- Hit list for query $q$: $H(q)$
- Bit score: $s_{qt}$
- Query max score: $S_q = \max_{t \in H(q)} s_{qt}$
- Target abundance: $A_t$
- Coverage confidence: $C_t$
- Posterior for hit $(q,t)$: $R_{qt}$
- Parameters: $\lambda, \alpha, \gamma$
- Small constant: $\varepsilon = 10^{-12}$

## Abundance initialization

For each query $q$, compute the sum of scores over its hit list:

$$
S^{\text{sum}}_q = \sum_{t \in H(q)} s_{qt}
$$

If $S^{\text{sum}}_q > 0$, each target gets a fractional count:

$$
c_{qt} = \frac{s_{qt}}{S^{\text{sum}}_q}
$$

Aggregate counts across all queries and normalize by $|Q|$:

$$
A_t = \frac{1}{|Q|} \sum_q c_{qt}
$$

If $|Q| = 0$, all $A_t$ are treated as $0$.

## Coverage confidence

For each target $t$, find its global alignment bounds:

$$
	ext{minPos}_t = \min_q \text{dbStart}_{qt}, \quad
	ext{maxPos}_t = \max_q \text{dbEnd}_{qt}
$$

Coverage is accumulated over the interval $[\text{minPos}_t,\text{maxPos}_t]$.
For each hit $(q,t)$:

$$
\ell_{qt} = \text{dbEnd}_{qt} - \text{dbStart}_{qt} + 1
$$

$$
m_{qt} = \frac{\exp(\lambda s_{qt})}{\ell_{qt}}
$$

For each position $p$ covered by the hit:

$$
	ext{cov}_t(p) \mathrel{+}= m_{qt}
$$

The per-target coverage confidence $C_t$ is the fraction of target positions where the accumulated per-hit weights reach at least 1.0. Values are clipped to $[0,1]$.

## Score term

Normalized score for a hit:

$$
\text{ns}_{qt} = \frac{s_{qt}}{S_q}
$$

Score term used for posteriors:

$$
T_{qt} = \exp(\lambda \text{ns}_{qt}) \cdot \max(A_t,\varepsilon)^{\alpha} \cdot \max(C_t,\varepsilon)^{\gamma}
$$

If $S_q \le 0$, then $T_{qt} = 0$ for all $t \in H(q)$.

## Posterior for each query

For each query $q$:

$$
Z_q = \sum_{t \in H(q)} T_{qt}
$$

$$
R_{qt} = \frac{T_{qt}}{Z_q}
$$

If $Z_q = 0$, then $R_{qt} = 0$.

## EM abundance update

Given posteriors, the next abundance is:

$$
\hat{A}_t = \frac{1}{|Q|} \sum_q R_{qt}
$$

## Log-likelihood (per query)

$$
\mathcal{L} = \frac{1}{|Q|} \sum_q \left( \sum_t R_{qt} \log(T_{qt}) - \log(Z_q) \right)
$$

## SQUAREM step (acceleration)

Two EM steps produce $x_1$ and $x_2$ from $x_0$:

$$
\mathbf{r} = x_1 - x_0, \quad \mathbf{v} = x_2 - x_1 - \mathbf{r}
$$

Acceleration scalar (clipped to $[-1, 1]$):

$$
a = -\frac{\|\mathbf{r}\|}{\|\mathbf{v}\|}
$$

Proposed update and simplex projection:

$$
\mathbf{x}_{\text{new}} = x_0 - 2a\mathbf{r} + a^2\mathbf{v}, \quad \mathbf{x}_{\text{new}} \in \Delta
$$

## Simplex projection

Clamp negatives and renormalize:

$$
x_i \leftarrow \max(0, x_i), \quad x_i \leftarrow \frac{x_i}{\sum_j x_j}
$$
