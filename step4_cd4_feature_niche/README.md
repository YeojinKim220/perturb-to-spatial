# Step 4 — CD4⁺ feature-niches (4-block, interpretable)

Refine the niche definition into an interpretable feature vector on the CD4⁺
substrate, and validate that intrinsic state scores are non-degenerate.

## Script

[`cd4_feature_niche.py`](cd4_feature_niche.py). Substrate: all CD4⁺ T cells
(naive + memory + Treg), LUAD-5 replicates R1–R3 (25,493 cells). Each cell gets a
**4-block feature vector**:

1. **Intrinsic state scores** — expression-based scores for 7 programs (naive,
   activated, cytotoxic, exhausted, Treg, Tfh, Th1) using panel-present markers.
2. **Neighbor cell-type composition** — the Gaussian-weighted mix of neighboring
   cell types (macrophage, fibroblast, B cell, endothelial, tumor, …) within
   150 µm, i.e. what kinds of cells surround this T cell.
3. **Neighbor suppressive/ligand programs** — expression of suppressive and
   ligand programs in the neighborhood.
4. **Distance to tumor** — log distance from each T cell to its nearest tumor
   cell.

Each block is standardized and then down-weighted by the square root of its
feature count, so a block with many features (e.g. neighbor composition) does not
drown out a single-column block (tumor distance). Then PCA → kNN →
Leiden(res=1.0).

## Why all-CD4 (not naive-only)

Restricting to naive cells collapses the exhaustion / Treg / Th1 endpoints to
≈ 0 (those states are absent in naive cells), which would make any relevance
scoring on them pure noise. Expanding to all CD4⁺ restores non-degenerate program
scores — this decision is validated as a testable claim inside the script.

## How to run

```bash
python cd4_feature_niche.py
```

Needs `scanpy`, `anndata`, `scipy`, `leidenalg`, `matplotlib`. Submit as a CPU
job (see [`../SLURM_TEMPLATE.md`](../SLURM_TEMPLATE.md)).

## Results

- **19 feature-niches.**
- **Non-degeneracy validated:** fraction of cells scoring positive — activated
  0.65, exhausted 0.53, cytotoxic 0.35, Treg 0.35, Th1 0.27 (all ≈ 0 in a
  naive-only substrate). Treg cells carry the highest Treg (1.60) and exhausted
  (0.59) scores, as expected.
- **Mature states separate into dedicated niches** — one niche is 86.8% Treg,
  another 80.2% memory; others differentiate along program axes (a Tfh-high
  niche, Th1-high niches) and by distance to the tumor boundary.

## Outputs

`cd4_allcd4_feature_niche_{assignment,profile}.csv`,
`cd4_allcd4_niche_celltype_composition.csv`,
`cd4_allcd4_program_{positivity,scores_by_subtype}.csv`,
`cd4_allcd4_feature_niche_{spatial,umap}.png`,
`cd4_allcd4_program_scores_by_subtype.png`. The per-cell niche assignment feeds
the corrected analyses in [`step7`](../step7_phase2_corrections/) and
[`step8`](../step8_phase3_reanalysis/).
