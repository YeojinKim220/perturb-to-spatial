# Phase 3 Review: 독립 정의와 evidence funnel 재점검

## One-line judgement
이번 추가 결과로 기존의 큰 오류들은 상당히 교정됐지만, clean target-specific spatial story는 아직 완성되지 않았다. 따라서 추가 target-specific spatial map은 만들지 않았다.

## Figure 1: independent state x environment context
Figure 1은 combined feature-niche clustering을 본문 근거로 쓰지 않고, intrinsic state와 environment context를 독립적으로 정의했다.

- Intrinsic state: naive-annotated CD4 T cell의 expression score만 사용했다.
- Environment context: tumor distance/contact와 neighbor cell-type composition만 사용했다.
- Myeloid/fibroblast/B-cell context는 section별 positive-abundance top-quantile rule로 다시 정의했고, context count QC를 `F1_context_qc.csv`에 저장했다.
- 총 15개 state-context association 중 8개가 R1-R3에서 같은 방향이었다.

관련 산출물: `fig1_state_context.png`, `F1_state_context_assoc.csv`, `F1_context_qc.csv`.

## Figure 2: activation-dependent perturbation with reliability/viability
기존 activation-context 결과는 유지하되, 후보 해석은 FDR 방식과 donor/viability 필터를 함께 보아야 한다.

- 기존 pooled-null cytotoxic FDR<0.1 후보는 23개였다.
- per-target add-one permutation 기준 cytotoxic FDR<0.1 후보는 127개였다.
- per-target FDR<0.1이면서 donor reliable, no viability flag인 후보는 3개였다.
- 기존 pooled-null 23개 중 donor reliable, no viability flag를 동시에 만족한 후보는 0개였다.

ZMAT1 FDR 0.024는 pooled-null result로만 표기해야 하며, 확정 hit로 사용하지 않는다. `n_downstream=0`은 measured response가 없다는 뜻이 아니라 statistically significant downstream DE genes가 0개라는 뜻이다.

관련 산출물: `FUN_fdr_hit_audit.csv`, `FUN_evidence_funnel.csv`, `FUN_clean_shortlist.csv`.

## Figure 3: corrected spatial contrast x perturbation integration
Contact-vs-far signature는 section별 cutoff, state 내 matching, `n_counts`, state score, local density, tumor-neighbor fraction 보정 후 다시 만들었다.

- Stress/stromal flag와 2/3 section direction agreement를 유지했다.
- Raw vs corrected target ranking은 Spearman 0.974-0.995로 높았지만, 이것은 correction이 모든 문제를 해결했다는 뜻이 아니라 leading rank structure가 크게 바뀌지 않았다는 뜻이다.
- Leave-one-section-out target-rank Spearman 평균은 0.129, top-20 recurrence 평균은 0.11로 낮았다.
- Evidence funnel 결과, clean shortlist는 0개였다.

관련 산출물: `fig3_corrected_integration.png`, `S2c_target_x_corrected_context.csv`, `S2c_bootstrap_rank_ci.csv`, `REP_leave_one_section.csv`.

## Moran's I re-check
Moran's I는 FOV-safe max-distance graph와 `n_counts` residualized program score로 재검정했다.

- BH 보정 후 유의한 test는 0/15개였다.
- 평균 absolute Moran's I는 0.031였다.
- 따라서 현재 해석은 "공간 구조가 약하거나 graph/residualization 설정에 민감하다"가 가장 안전하다.

관련 산출물: `fig_moran_recheck.png`, `MI_morans_recheck.csv`.

## Conclusion and limitations
현재 결과는 더 많은 target을 주장하기보다, 기존 주장의 취약한 부분을 정리하는 데 의미가 있다. R1-R3은 one-patient sections이므로 biological replicates로 해석하면 안 된다. 최종 메시지는 "spatial context can generate testable CRISPRi knockdown hypotheses, but current clean target-specific support is insufficient"로 두는 것이 맞다.
