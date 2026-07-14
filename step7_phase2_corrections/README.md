# Step 7 — Phase-2 corrections (proper nulls, FDR, reliability & viability filters)

A methodology review of step 6 flagged that the exploratory relevance ranking
lacked a proper statistical null and multiple-testing control. This step rebuilds
the validation with permutation nulls and BH-FDR, and adds reliability and
viability filters — and the corrected result is substantially more conservative.

## Scripts

| Script | What it does |
|---|---|
| [`phase2_analysis.py`](phase2_analysis.py) | Runs stages **S1–S5** (below), each writing its own outputs so a late failure keeps earlier ones. GRCh38 symbols, `seed=0`, overlap reported per score. |
| [`s1_fdr_audit.py`](s1_fdr_audit.py) | The decisive audit: recomputes the held-out validation with a **per-target add-one permutation null** (1,000 permutations) and compares it against the old pooled-null definition, applying BH-FDR per endpoint. |
| [`phase2_figures.py`](phase2_figures.py) | Builds the compact summary JSON and the niche-contrast figure from the stage outputs. |

### Stages in `phase2_analysis.py`

- **S1** — held-out, direction-corrected validation on the 471 markerless genes,
  with a permutation null and BH-FDR (reports `rank_pct` and FDR, not a raw
  percentile).
- **S2** — within-state niche-contrast signature (tumor-contact vs tumor-far
  cells) and the resulting perturbation rank-shift.
- **S3** — Perturb-seq reliability (donor concordance, downstream-DE count,
  on-target significance) and a viability / general-disruption filter.
- **S4** — spatial quantification: Moran's I per program per section,
  neighborhood enrichment with a permutation null, minimum niche-size filter, and
  R1–R3 rank reproducibility.
- **S5** — recomputation of the Rest-vs-Stim48 flip screen from a single frozen
  run, so all downstream flip counts come from one internally consistent source.

## How to run

Submit as a CPU job (needs `numpy`, `pandas`, `h5py`, `scipy`, `statsmodels`,
`scanpy`; see [`../SLURM_TEMPLATE.md`](../SLURM_TEMPLATE.md)):

```bash
python phase2_analysis.py     # S1-S5 -> 2-Results/phase2/
python s1_fdr_audit.py        # per-target permutation FDR audit
python phase2_figures.py      # summary JSON + figures
```

## Results — the correction bites

- **Per-target permutation FDR (S1 / audit):** under the corrected null, only
  `cytotoxic_enhancement` retains any targets at FDR < 0.1 (23 under the old
  pooled null; 127 under the per-target null). Every other endpoint's minimum FDR
  is ≫ 0.1 (activation 0.73, Th1 0.86, Treg reduction 0.99, exhaustion reversal
  0.99). **The earlier PTPN2 / exhaustion / Treg confirmations do not reproduce
  and are withdrawn.**
- **Reliability (S3):** only ~10% of perturbations (1,144 / 11,287) are donor-
  reliable with ≥ 5 downstream DE genes and a significant on-target effect; ~17%
  carry a viability / broad-disruption flag.
- **Niche contrast (S2):** the contact-vs-far signature is dominated by
  heat-shock / stress / stromal genes — a tumor-proximity stress axis rather than
  a clean niche-resolution signal; batch/stress correction is the needed next
  step (done in step 8).
- **Spatial (S4):** Moran's I is statistically detectable but very weak
  (I ≈ 0.01–0.02) under the 6-nearest-neighbor graph used here — a signal that
  step 8 re-tests with a stricter graph.

## Takeaway

The proper null converts the exploratory step-6 shortlist into a mostly
**negative, calibrated** result: FDR alone does not select clean targets here.
Step 8 pushes this further with independent state/context definitions and a
FOV-safe spatial graph.

## Outputs

`S1_heldout_direction_corrected_fdr.csv`, `S2_*`, `S3_reliability_viability.csv`,
`S4_morans_I_per_section.csv`, `S4_section_rank_reproducibility.json`,
`S5_flip_qc.json`, `phase2_report_summary.json`, `FUN_fdr_hit_audit.csv`, and the
phase-2 figures.
