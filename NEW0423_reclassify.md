# EM_reclassify.cpp Equations (NEW0423, KaTeX)

본 문서는 `src/util/EM_reclassify.cpp`의 현재 구현 수식을 KaTeX 친화 문법으로 정리한 것이다.

## 0. Constants
- $\varepsilon = 10^{-12}$ (`EPS`)
- $r_{\min}=-60,\ r_{\max}=60$ (`LOG_COMPATIBILITY_MIN`, `LOG_COMPATIBILITY_MAX`)
- $\tau_{\alpha}=3$ (`ABUNDANCE_EXP_TAU`)
- $\epsilon_\pi=10^{-8}$ (`ABUNDANCE_SMOOTH_EPS`)
- SQUAREM step bound: $[-1,1]$ (`STEP_MIN`, `STEP_MAX`)
- Drop filter activation threshold: $N_{\min}=20$ (`MIN_FILTER_TARGETS`)
- Minimum tail size for cutoff: $k_{\min}=2$ (`MIN_TAIL_TARGETS`)

---

## 1. Abundance Initialization (`initAbundance`)
쿼리 $q$의 hit 집합을 $H(q)$라 하자.

쿼리 내 score 합:
$$
S_q = \sum_{t\in H(q)} \max(\mathrm{score}_{q,t}, 0)
$$

초기 가중치:
$$
w_{q,t}^{(0)} =
\begin{cases}
\dfrac{\max(\mathrm{score}_{q,t},0)}{S_q}, & S_q > 0 \\
0, & S_q \le 0
\end{cases}
$$

타깃별 누적:
$$
\tilde{\pi}_t^{(0)} = \sum_q w_{q,t}^{(0)}
$$

쿼리 수 평균:
$$
\hat{\pi}_t^{(0)} =
\begin{cases}
\dfrac{\tilde{\pi}_t^{(0)}}{|Q|}, & |Q| > 0 \\
\tilde{\pi}_t^{(0)}, & |Q| = 0
\end{cases}
$$

simplex 정규화:
$$
\pi_t^{(0)} = \frac{\hat{\pi}_t^{(0)}}{\sum_{t'} \hat{\pi}_{t'}^{(0)}}
$$
(분모가 0이 아니면 수행)

---

## 2. Coverage Confidence (`initCoverageConfidence`)
타깃 $t$의 관측 span 길이:
$$
L_t = \max\left(1,\ \mathrm{maxEnd}_t - \mathrm{minStart}_t + 1\right)
$$

쿼리 내 hit 정규화 weight:
$$
\omega_{q,h} =
\begin{cases}
\dfrac{\mathrm{score}_{q,h}}{\sum_{h'\in H(q)} \mathrm{score}_{q,h'}}, & \sum \mathrm{score} > 0 \\
0, & \text{otherwise}
\end{cases}
$$

위치별 누적:
$$
\mathrm{covConf}_t(p) = \sum_{h:\ p\in[start_h,end_h]} \omega_{q(h),h}
$$

클리핑:
$$
\tilde{c}_t(p) = \min\left(1,\mathrm{covConf}_t(p)\right)
$$

coverage fraction:
$$
f_t = \frac{1}{L_t} \sum_{p=1}^{L_t} \tilde{c}_t(p)
$$

HHI:
$$
\mathrm{HHI}_t =
\begin{cases}
\dfrac{\sum_p \tilde{c}_t(p)^2}{\left(\sum_p \tilde{c}_t(p)\right)^2}, & \sum_p \tilde{c}_t(p) > 0 \\
1, & \sum_p \tilde{c}_t(p) = 0
\end{cases}
$$

penalty:
$$
\mathrm{penalty}_t = 1 - \mathrm{HHI}_t
$$

최종 confidence:
$$
c_t = \operatorname{clamp}_{[0,1]}\left(f_t\cdot\mathrm{penalty}_t\right)
$$

여기서
$$
\operatorname{clamp}_{[0,1]}(x) = \max(0, \min(1, x))
$$

---

## 3. Query-Target Compatibility (`compatibilityLogTerm`)
best bit score:
$$
b_q^{\max} = \max_{t\in H(q)} \mathrm{score}_{q,t}
$$

coverage feature:
$$
\mathrm{cov}_{q,t} = \operatorname{clamp}_{[0,1]}\left(\min(qcov_{q,t}, dbcov_{q,t})\right)
$$

bit score 차이:
$$
\Delta b_{q,t} = \mathrm{score}_{q,t} - b_q^{\max}
$$

compatibility log-term:
$$
r_{q,t} = \beta_{bit}\Delta b_{q,t} + \beta_{id}\,id_{q,t} + \beta_{cov}\,\mathrm{cov}_{q,t}
$$

clip:
$$
r_{q,t} \leftarrow \min\left(r_{\max},\max\left(r_{\min}, r_{q,t}\right)\right)
$$

compatibility:
$$
\phi_{q,t} = \exp(r_{q,t})
$$

기본 호출값:
- $\beta_{id}=1.0$
- $\beta_{cov}=0.25$
- $\beta_{bit}=\texttt{par.reclassifyLambda}$

---

## 4. Annealed Exponent (`squarem`)
iteration $m$ (`iter+1`)에서:
$$
\alpha^{(m)} = \alpha_{\max}\left(1 - e^{-m/\tau_\alpha}\right)
$$
- $\alpha_{\max}=\texttt{par.reclassifyAlpha}$
- $\tau_\alpha=3$

---

## 5. E-step Posterior (`computePosterior`)
수치 안정화 abundance:
$$
\bar{\pi}_t = \max(\pi_t, \epsilon_\pi)
$$

가중치:
$$
N_{q,t} = \phi_{q,t}\cdot\bar{\pi}_t^{\alpha^{(m)}}
$$

posterior:
$$
p_{q,t} = \frac{N_{q,t}}{\sum_{t'\in H(q)} N_{q,t'}}
$$

---

## 6. M-step Abundance Update (`emUpdate`)
responsibility 합:
$$
R_t = \sum_q p_{q,t}
$$

coverage prior smoothing:
$$
\tilde{\pi}_t = R_t + \kappa c_t + \epsilon_\pi
$$
- $\kappa=\texttt{coveragePriorWeight}$ (코드에서는 `par.reclassifyGamma` 전달)

정규화:
$$
\pi_t^{new} = \frac{\tilde{\pi}_t}{\sum_{t'} \tilde{\pi}_{t'}}
$$

---

## 7. Log-likelihood (`logLikelihood`)
쿼리별 mixture:
$$
M_q = \sum_{t\in H(q)} \phi_{q,t}\,\bar{\pi}_t^{\alpha^{(m)}}
$$

평균 로그우도:
$$
\mathcal{L} = \frac{1}{|Q|}\sum_q \log\left(\max(M_q, 10^{-300})\right)
$$

---

## 8. Simplex Projection (`projectSimplex`)
비음수화:
$$
x_i^+ = \max(x_i,0)
$$

합 $S=\sum_i x_i^+$가 양수면:
$$
\Pi_i = \frac{x_i^+}{S}
$$

---

## 9. SQUAREM Extrapolation (`squarem`)
두 번의 EM 결과로:
$$
r = x_1 - x_0,\qquad v = x_2 - x_1 - r
$$

노름:
$$
\|r\|_2 = \sqrt{\sum_i r_i^2},\qquad \|v\|_2 = \sqrt{\sum_i v_i^2}
$$

가속계수:
$$
a =
\begin{cases}
-1, & \|v\|_2 = 0 \\
-\dfrac{\|r\|_2}{\|v\|_2}, & \text{otherwise}
\end{cases}
$$

clip:
$$
a \leftarrow \min(1,\max(-1,a))
$$

외삽:
$$
x_{new} = x_0 - 2ar + a^2 v
$$

LL 감소 fallback:
$$
\text{if } \mathcal{L}(x_{new}) < \mathcal{L}_{prev} - 10^{-9},\ \ x_{new} \leftarrow x_2
$$

수렴량:
$$
\Delta = \max_i |x_{new,i} - x_{0,i}|
$$

조건:
$$
\Delta < \mathrm{tol}\ \text{and}\ \mathrm{iter} > 5
$$

---

## 10. Target Drop Cutoff (`selectDroppedTargets`)
$$
\mathrm{maxTailFraction} = \operatorname{clamp}_{[0,1]}\left(\frac{\mathrm{maxDropPercentage}}{100}\right)
$$

활성 조건:
$$
n_{\text{targets}} \ge 20,\quad T=\sum_i a_i > \varepsilon,\quad \mathrm{maxTailFraction}>0
$$

전체 abundance 질량:
$$
T = \sum_i a_i,\qquad T_{maxTail} = \mathrm{maxTailFraction}\cdot T
$$

tail-quantile cutoff:
- low-tail: $\mathrm{cutoff}=a_{(k)}$
- high-tail: $\mathrm{cutoff}=a_{(n-k+1)}$
- 여기서 $k$는 누적 질량이 $T_{maxTail}$를 넘지 않는 최대 tail 크기이며, $2 \le k < n$ 필요

(현재 reclassify 경로에서는 low-tail abundance drop 사용)

최종 drop 집합은 abundance 오름차순으로 tail을 다시 선택하며, 선택 질량이 $T_{maxTail}$를 넘기기 직전까지 포함한다(최소 2개는 허용). 전부 drop되는 경우에는 drop을 취소한다.

---

## 11. Aux
제거 비율 로그:
$$
\mathrm{removedPct} = 100\cdot\frac{|\mathrm{dropped}|}{\mathrm{totalTargets}}
$$
(분모 0이면 0 처리)

출력 정렬 우선순위:
$$
\mathrm{posterior\ desc} \rightarrow \mathrm{bitScore\ desc}
$$

출력 시 `seqId`에 posterior 저장:
$$
\mathrm{seqId}_{out} = p_{q,t}
$$

---

## 12. Update Notes (Code-aligned)
- `squarem`/`emUpdate`/`logLikelihood`에서 compatibility 가중치 호출값은 $\beta_{cov}=0.25$로 고정되어 있다.
- 타깃 drop은 coverage가 아니라 abundance low-tail 기준으로만 수행된다.
- drop cutoff는 tail-quantile 방식만 사용한다.
