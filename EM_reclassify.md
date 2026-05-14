# EM_reclassify.cpp Summary (04/12)

## Overview
- This module performs EM-based reclassification of alignment hits.
- It updates target abundances and posterior probabilities, optionally drops low-abundance targets, and writes a new alignment DB where hits are sorted by posterior.
- Primary input: alignment result DB (`par.db3`). Primary output: reclassified alignment DB (`par.db4`).

## Example Usage Order
- `mmseqs createdb query.fasta queryDB`
- `mmseqs search queryDB targetDB alignDB tmp -a`
- `mmseqs reclassify queryDB targetDB alignDB newDB`
- `mmseqs convertalis queryDB targetDB newDB reclassify_result.m8`
- `mmseqs abundance queryDB targetDB newDB abundance.tsv --taxonomy 1`

Notes:
- Plain custom FASTA target DB: use `createdb` + `createtaxdb` as needed.
- Prebuilt taxonomy-ready MMseqs DB: `createtaxdb` is usually not needed.

## Key Data Structures
- `ReclassTaxEntry`: per-hit working record with `abundance`, `posterior`, `coverageConfidence`.
- `MappingTable`: `queryKey -> vector<ReclassTaxEntry>`.
- `ReclassTaxContext`: full working set (mapping table, query order, target set, query count, output-format flags).
- `TargetStats`: per-target aggregate stats used for filtering/reporting.

## Mathematical Definitions
Let `q` be a query and `t` a target hit for `q`.

- Normalized score:

$$
s_{q,t} = \frac{\text{score}_{q,t}}{\max_{t'} \text{score}_{q,t'}}
$$

- Score term used in posterior update:

$$
\text{scoreTerm}_{q,t} = \exp(\lambda \cdot s_{q,t}) \cdot a_t^{\alpha} \cdot c_t^{\gamma}
$$

where `a_t` is abundance, `c_t` is coverage confidence.

- Posterior:

$$
p_{q,t} = \frac{\text{scoreTerm}_{q,t}}{\sum_{t'} \text{scoreTerm}_{q,t'}}
$$

- Abundance update:

$$
a_t^{\text{new}} = \frac{1}{|Q|}\sum_{q \in Q} p_{q,t}
$$

- Log-likelihood (average per query):

$$
\mathcal{L} = \frac{1}{|Q|}\sum_{q \in Q}\left(\sum_t p_{q,t}\log(\text{scoreTerm}_{q,t}) - \log\sum_t \text{scoreTerm}_{q,t}\right)
$$

## Coverage Confidence
Coverage confidence is computed once at initialization and then kept fixed during EM.

### 1) Query-level hit weight (score normalization)
For each query `q`, define hit weight by plain bit-score normalization (no `exp`, no `lambda`):

$$
w_{q,h} = \frac{\text{score}_{q,h}}{\sum_{h'} \text{score}_{q,h'}}
$$

### 2) Position-wise accumulation on each target
For each target position `p`, accumulate weighted support:

$$
\text{covConf}_t(p) = \sum_{h: p \in [\text{start}_h, \text{end}_h]} w_{q,h}
$$

Clip each position to avoid over-crediting stacked/repeat mappings at one locus:

$$
\tilde{c}_t(p) = \min(1, \text{covConf}_t(p))
$$

### 3) Base coverage over observed span
Observed span length:

$$
L_t = \max(1, \text{endPos}_t - \text{startPos}_t + 1)
$$

Base coverage fraction:

$$
f_t = \frac{1}{L_t}\sum_{p=1}^{L_t}\tilde{c}_t(p)
$$

### 4) Concentration penalty (HHI-based)
Define concentration (Herfindahl-Hirschman Index over clipped position mass):

$$
\text{HHI}_t = \frac{\sum_{p=1}^{L_t}\tilde{c}_t(p)^2}{\left(\sum_{p=1}^{L_t}\tilde{c}_t(p)\right)^2}
$$

with the edge case:

$$
\text{HHI}_t = 1 \quad \text{if } \sum_{p}\tilde{c}_t(p)=0
$$

Convert to dispersion reward (low when concentrated, high when spread):

$$
\text{penalty}_t = 1 - \text{HHI}_t
$$

### 5) Final coverage confidence

$$
c_t = \operatorname{clamp}_{[0,1]}\left(f_t \cdot \text{penalty}_t\right)
$$

Important details:
- Position contribution is clipped by `min(1, ...)`.
- Normalization length is **observed span** (`maxEnd - minStart + 1`), not full `dbLen`.
- Final value includes HHI-based concentration penalty and is clamped to `[0, 1]`.
- `lambda` still affects posterior via `scoreTerm`, but no longer affects coverage-confidence initialization.
- Behavior intent:
  - Many hits stacked in one local region: high concentration -> high HHI -> small `(1-HHI)` -> lower `c_t`.
  - Hits distributed across many positions: lower concentration -> lower HHI -> larger `(1-HHI)` -> higher `c_t`.

## Processing Pipeline
1. Load alignment DB
- `loadAlignmentDb()` parses alignment rows into `mappingTable`.
- Keeps backtrace/ORF layout flags for output compatibility.

2. Initialize
- `initAbundance()` initializes abundance from query-normalized hit scores.

### Abundance Initialization (초기값 계산)

각 쿼리 $q$에 대해, 해당 쿼리의 모든 target hit $t$에 대해 bit score를 정규화합니다:

$$
w_{q,t} = \frac{\text{score}_{q,t}}{\sum_{t'} \text{score}_{q,t'}}
$$

각 target $t$의 초기 abundance는, 해당 target이 등장한 모든 쿼리에서의 $w_{q,t}$ 값을 평균낸 값입니다:

$$
a_t^{(0)} = \frac{1}{|Q|} \sum_{q \in Q} w_{q,t}
$$

즉, 각 target의 abundance는 쿼리별로 score-normalized weight를 합산한 뒤 전체 쿼리 수로 나누어 초기화합니다.
- `initCoverageConfidence()` computes per-target `coverageConfidence` using `score/sum(score)` weights.

3. EM + SQUAREM
- `computePosterior()` computes per-query posterior.
- `emUpdate()` updates abundance; coverage confidence remains fixed.
- `squarem()` accelerates EM with extrapolation and simplex projection.
- If LL decreases, fallback to conservative step (`x2`).

4. Optional target filtering
- `collectTargetStats()` builds target-level abundance/coverage/interval stats.
- `selectDroppedTargets()` chooses low-tail drops under mass cap.
- `applyDroppedTargets()` removes dropped targets from mapping/query/target sets.

5. Write output DB
- `writeReclassifiedDb()` writes hits sorted by posterior.
- Posterior is written into output `seqId` field.

## Parameters (Current Defaults)
- `--lambda` (`reclassifyLambda`): `2`
- `--alpha` (`reclassifyAlpha`): `1.0`
- `--gamma` (`reclassifyGamma`): `1.0`
- `--max-iter` (`reclassifyMaxIterations`): `100`
- `--tol` (`reclassifyTolerance`): `1e-5`
- `--drop-percentage` (`reclassifyMaxDropPercentage`): `10.0`

## Output/Storage Clarification
- `reclassify` output DB stores posterior in `seqId`.
- `coverageConfidence` is not persisted in alignment DB columns.
- `abundance` recomputes `coverageConfidence` from input alignments and writes it to `abundance.tsv`.

## Consistency Across Modules
The updated coverage-confidence logic (`score/sum(score)` + observed-span normalization + HHI penalty) is applied consistently in:
- `src/util/EM_reclassify.cpp`
- `src/util/EM_abundnace.cpp`

## Logging
- Coverage-confidence initialization complete message.
- Iteration LL and parameter delta.
- Number of dropped targets and applied abundance cutoff.
