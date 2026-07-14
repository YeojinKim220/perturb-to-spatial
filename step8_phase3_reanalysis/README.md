# Step 8 — Phase-3 re-analysis (independent definitions, FOV-safe spatial, evidence funnel)

The final, most defensible pass. It removes the circularity that let step 6 look
positive, re-tests the spatial signal with a stricter graph, and runs every
candidate through an evidence funnel. The result is a **calibrated negative**:
the pipeline works end-to-end, but no endpoint yields a clean, reproducible,
target-specific spatial shortlist.

## Scripts

| Script | What it does |
|---|---|
| [`phase3_core.py`](phase3_core.py) | Runs analyses **F1, MI, S2c, REP, FUN** (below). Uses `Path(__file__)` to locate the project base, `seed=0`. |
| [`phase3_report.py`](phase3_report.py) | Generates the phase-3 figures and the English/Korean Markdown reports from the core outputs. |

### Analyses in `phase3_core.py`

- **F1 — independent state × context.** Intrinsic state (expression scores within
  naive-annotated CD4⁺ cells) and environment context (tumor distance/contact,
  neighbor composition) are defined **independently**, instead of being read off
  one combined clustering — so a state–context association can no longer be an
  artifact of the clustering that defined both.
- **MI — Moran's I re-check.** Rebuilt with a FOV-safe maximum-distance graph
  (edges beyond ~40 µm and across FOV boundaries removed) and each program score
  residualized against total counts, then BH-corrected.
- **S2c — stress-corrected contact-vs-far signatures.** The niche contrast is
  recomputed with section-specific cutoffs, within-state matching, and adjustment
  for total counts, state score, local density, and tumor-neighbor fraction.
- **REP — leave-one-section-out reproducibility.** Define the contrast on two
  sections, transfer to the third, and measure the perturbation-rank Spearman —
  the honest cross-section estimate (R1–R3 are three sections of **one patient**,
  not biological replicates).
- **FUN — evidence funnel.** Pass every candidate through FDR + reliability +
  viability + spatial evidence and count what survives.

## How to run

Submit as a CPU job (needs `numpy`, `pandas`, `h5py`, `scipy`, `scanpy`,
`statsmodels`, `matplotlib`; see [`../SLURM_TEMPLATE.md`](../SLURM_TEMPLATE.md)):

```bash
python phase3_core.py       # F1, MI, S2c, REP, FUN -> 2-Results/phase3/
python phase3_report.py     # figures + EN/KO markdown
```

## Results — the calibrated negative

- **Moran's I:** 0 of 15 section×program tests reach significance after BH
  correction (mean |I| = 0.031). The earlier "significant but weak"
  autocorrelation was partly an artifact of long-range / cross-FOV edges.
- **Cross-section reproducibility:** leave-one-section-out perturbation-rank
  Spearman averages ~0.13 (top-20 recurrence ~0.11) — far below the ~0.9 obtained
  when the niche was defined jointly on all three sections; the independent value
  is the honest one.
- **State × context:** only the tumor-distance axis carries reproducible signal
  (Th1-high cells enriched at tumor contact, cytotoxic-high cells depleted, both
  sign-consistent across R1–R3); myeloid/fibroblast/B-cell contexts are
  degenerate for naive CD4⁺ cells (a sampling limitation, stated as such).
- **Evidence funnel:** after FDR + reliability + viability + spatial evidence,
  the clean target-specific shortlist is **empty**.

## Final message

The value of this round is not a longer target list but a more defensible audit
of the central claim. Spatial context can *generate* testable CRISPRi-knockdown
hypotheses, but current clean, target-specific spatial support is insufficient —
a genuinely useful, calibrated negative result. R1–R3 are one-patient sections
and must not be read as biological replicates. This remains a hypothesis-
generation instrument, not a source of causal or predictive claims; every
projected map is labeled exploratory.

## Outputs

`fig1_state_context.png`, `F1_state_context_assoc.csv`, `F1_context_qc.csv`,
`fig_moran_recheck.png`, `MI_morans_recheck.csv`, `fig3_corrected_integration.png`,
`S2c_target_x_corrected_context.csv`, `REP_leave_one_section.csv`,
`FUN_evidence_funnel.csv`, `FUN_clean_shortlist.csv`, and the phase-3 EN/KO
reports.
