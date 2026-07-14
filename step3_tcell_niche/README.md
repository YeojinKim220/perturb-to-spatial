# Step 3 — T-cell niches from neighborhood expression

Build T-cell niches bottom-up rather than reusing the atlas's provided niche
labels: cluster each T cell by the expression of its spatial neighborhood, then
compare the result against the reference labels.

## Script

[`tcell_niche.py`](tcell_niche.py). For each T cell *i*, compute a
Gaussian-kernel-weighted neighborhood expression vector
`NE_i(g) = Σ_j w_ij · x_j(g)`, with `w_ij = exp(-d_ij² / 2σ²)`, over neighbors
found by a KD-tree radius query **within the same sample**. Then
scale → PCA(30) → kNN(15) → Leiden(res=1.0), and compare to the reference niche
labels by adjusted Rand index (ARI) and normalized mutual information (NMI).

**Data-driven bandwidth:** σ = 50 µm and radius = 3σ = 150 µm are set from the
*measured* median nearest-neighbor distance (46 µm), not a guessed constant.
Neighbors include all cell types (so a niche is the microenvironment, not just
other T cells).

## How to run

```bash
python tcell_niche.py
```

Needs `scanpy`, `anndata`, `scipy`, `scikit-learn`, `leidenalg`, `matplotlib`.
Submit as a CPU job (~5 min; see [`../SLURM_TEMPLATE.md`](../SLURM_TEMPLATE.md)).

## Results

- **20 T-cell niches** on 66,825 T cells.
- **ARI = 0.208, NMI = 0.338** vs the reference niches — moderate agreement; the
  new niches *subdivide* the reference (the TLS and immune niches each split into
  several sub-niches).

## Caveat surfaced here

The per-sample spatial maps show a strong batch effect — several samples collapse
to one or two niches (e.g. LUAD-12 is 94.6% one niche). Cause: clustering on
**non-batch-corrected** neighborhood expression, so sample-level expression
differences dominate. This motivates the tighter substrate and the interpretable
feature blocks used in [`step4`](../step4_cd4_feature_niche/), and batch
correction remains an open next step.

## Outputs

`tcell_niche_vs_cellcharter_heatmap.png`, `tcell_niche_umap.png`,
`tcell_niche_spatial.png`, `tcell_niche_annotation.csv`,
`tcell_niche_vs_cellcharter_{counts,rownorm}.csv`, `tcell_niche_assignment.csv`.
