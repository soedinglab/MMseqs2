# Reclassify 0406 Notes

This document summarizes the logic and formulas in `src/util/reclassify_taxonomy.cpp`.

## Overview

- Input: an alignment-result DB (from `mmseqs search -a` or `mmseqs align`) with per-hit fields.
- Output:
  - `new_alignment_result.m8` (re-ranked by posterior)
  - `protein_abundance.tsv` (per-target stats)
  - `taxonomy_abundance.tsv` (optional taxonomy aggregation)

The pipeline is:
1) Read alignments into per-query hit lists.
2) Initialize abundance and coverage confidence.
3) Run EM with SQUAREM acceleration.
4) Drop low-abundance targets by tail cutoff.
5) Write outputs.

## Ex usage order
- mmseqs createdb query.fasta queryDB
- mmseqs search queryDB targetDB alignDB tmp -a
- mmseqs reclassify queryDB targetDB alignmentDB newDB
- mmseqs convertalis queryDB targetDB newnDB(or alignDB) reclassify_result.m8
- mmseqs abundance queryDB targetDB newDB abundance.report (--taxonomy 1)
* Plain custom FASTA target DB: createdb + createtaxdb
* Prebuilt taxonomy-ready MMseqs database: createtaxdb usually not needed

## Data structures

- Query key: the DB key for a query sequence.
- Target key: the DB key for a target sequence.
- Hit: a `Matcher::result_t` record (alignment fields, scores, positions).
- Per-hit state:
  - abundance $A_t$
  - posterior $R_{qt}$
  - coverage confidence $C_t$

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

## Step 1: Load alignment DB

- Each DB entry is parsed into per-query hit lists.
- Targets are tracked in a global set.

## Step 2: Abundance initialization

For each query $q$, sum scores over hits:

$$
S^{\text{sum}}_q = \sum_{t \in H(q)} s_{qt}
$$

If $S^{\text{sum}}_q > 0$, assign fractional counts:

$$
c_{qt} = \frac{s_{qt}}{S^{\text{sum}}_q}
$$

Aggregate over all queries and normalize by $|Q|$:

$$
A_t = \frac{1}{|Q|} \sum_q c_{qt}
$$

## Step 3: Coverage confidence (per target)

Threads usage:
- The coverage-confidence calculation is parallelized across targets with OpenMP.
- The loop over target IDs runs with `#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)`.
- Other steps in this pipeline are single-threaded in this file.

For each target $t$, compute global bounds:

$$
\text{minPos}_t = \min_q \text{dbStart}_{qt}, \quad
\text{maxPos}_t = \max_q \text{dbEnd}_{qt}
$$

For each hit $(q,t)$, define alignment length:

$$
\ell_{qt} = \text{dbEnd}_{qt} - \text{dbStart}_{qt} + 1
$$

Define per-hit weight:

$$
\text{expScore}_{qt} = \exp(\lambda s_{qt}), \quad
w_{qt} = \frac{\text{expScore}_{qt}}{\sum_{t' \in H(q)} \text{expScore}_{qt'}}
$$

Coverage accumulation over the target interval:

$$
\text{cov}_t(p) \mathrel{+}= \frac{\text{expScore}_{qt}}{\ell_{qt}}
$$

Coverage confidence accumulation:

$$
\text{conf}_t(p) \mathrel{+}= w_{qt}
$$

Coverage confidence (per target):

$$
C_t = \mathrm{clamp}_{[0,1]}\left(\frac{\sum_p \min(1, \text{conf}_t(p))}{\text{targetLen}_t}\right)
$$

## Step 4: Score term and posterior

Normalized score:

$$
\text{ns}_{qt} = \frac{s_{qt}}{S_q}
$$

Score term:

$$
T_{qt} = \exp(\lambda \text{ns}_{qt}) \cdot \max(A_t, \varepsilon)^{\alpha} \cdot \max(C_t, \varepsilon)^{\gamma}
$$

Posterior for query $q$:

$$
Z_q = \sum_{t \in H(q)} T_{qt}, \quad
R_{qt} = \frac{T_{qt}}{Z_q}
$$

## Step 5: EM update

Abundance update:

$$
\hat{A}_t = \frac{1}{|Q|} \sum_q R_{qt}
$$

## Step 6: Log-likelihood (per query)

$$
\mathcal{L} = \frac{1}{|Q|} \sum_q \left( \sum_t R_{qt} \log(T_{qt}) - \log(Z_q) \right)
$$

## Step 7: SQUAREM acceleration

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

Simplex projection:

$$
x_i \leftarrow \max(0, x_i), \quad x_i \leftarrow \frac{x_i}{\sum_j x_j}
$$

If the accelerated step decreases log-likelihood, the code falls back to $x_2$.

## Step 8: Target filtering (low-abundance drop)

- Compute a cutoff on target abundances using a largest-gap rule or tail quantile rule.
- Drop up to `maxDropPercentage` of targets (with a minimum tail size).
- Remove dropped targets from the mapping table and recompute output stats.

## Outputs

- `new_alignment_result.m8`
  - Hits re-sorted by posterior.
  - Alignment summary derived from backtrace if present.
- `protein_abundance.tsv`
  - Per-target: abundance %, coverage confidence, drop flag, mapped interval(s).
- `taxonomy_abundance.tsv` (if taxonomy enabled)
  - Aggregated by taxonomy ID.

## Notes

- `C_t` (coverage confidence) is used directly in the score term.
- Parallelism is used in coverage confidence initialization across targets.
- No other steps in this file use OpenMP; EM and output phases are serial here.
