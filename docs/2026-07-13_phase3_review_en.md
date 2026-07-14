# Phase 3 Review: Independent Definitions and Evidence Filtering

## One-line judgement
The additional analyses fixed several major issues, but the clean target-specific spatial story is still not complete. For that reason, no additional target-specific spatial map was generated.

## Figure 1: independent state x environment context
Figure 1 no longer uses combined feature-niche clustering as primary evidence. Intrinsic state and environment context are defined independently.

- Intrinsic state uses expression scores within naive-annotated CD4 T cells.
- Environment context uses tumor distance/contact and neighbor cell-type composition.
- Myeloid/fibroblast/B-cell contexts use a section-wise positive-abundance top-quantile rule, with context count QC saved in `F1_context_qc.csv`.
- 8 of 15 state-context associations had the same direction across R1-R3.

Outputs: `fig1_state_context.png`, `F1_state_context_assoc.csv`, `F1_context_qc.csv`.

## Figure 2: activation-dependent perturbation with reliability/viability
The activation-context result is retained, but target interpretation must combine FDR definition, donor reliability, and viability flags.

- The old pooled-null cytotoxic FDR<0.1 set contained 23 targets.
- The per-target add-one permutation cytotoxic FDR<0.1 set contained 127 targets.
- Per-target FDR<0.1 plus donor reliable plus no viability flag left 3 targets.
- Among the old pooled-null 23 targets, donor reliable plus no viability flag left 0 targets.

ZMAT1 FDR 0.024 should be described only as a pooled-null result, not as a confirmed hit. `n_downstream=0` means zero statistically significant downstream DE genes, not absence of measured expression response.

Outputs: `FUN_fdr_hit_audit.csv`, `FUN_evidence_funnel.csv`, `FUN_clean_shortlist.csv`.

## Figure 3: corrected spatial contrast x perturbation integration
The contact-vs-far signature was rebuilt using section-specific cutoffs, within-state matching, and adjustment by `n_counts`, state score, local density, and tumor-neighbor fraction.

- Stress/stromal flags and >=2/3 section direction agreement are retained.
- Raw vs corrected target ranking Spearman values were 0.974-0.995; this means the rank structure was stable, not that the confound is fully solved.
- Leave-one-section-out target-rank Spearman averaged 0.129, and top-20 recurrence averaged 0.11.
- After the full evidence funnel, the clean shortlist was empty.

Outputs: `fig3_corrected_integration.png`, `S2c_target_x_corrected_context.csv`, `S2c_bootstrap_rank_ci.csv`, `REP_leave_one_section.csv`.

## Moran's I re-check
Moran's I was re-tested with a FOV-safe max-distance graph and `n_counts`-residualized program scores.

- BH-significant tests: 0/15.
- Mean absolute Moran's I: 0.031.
- The safest interpretation is that spatial structure is weak or sensitive to graph/residualization choices.

Outputs: `fig_moran_recheck.png`, `MI_morans_recheck.csv`.

## Conclusion and limitations
The value of the current result is not a larger target list, but a more defensible audit of the central claim. R1-R3 are one-patient sections and should not be interpreted as biological replicates. The final message should be: spatial context can generate testable CRISPRi knockdown hypotheses, but current clean target-specific support is insufficient.
