# Step 2 — Spatial atlas visualization & T-cell-rich niche flagging

Visualize the CosMx NSCLC atlas by cell type and quantify which niches are
enriched for T cells, to confirm the planned positive control (the TLS-like
"lymphoid structure" niche) before any perturbation work.

## Script

[`spatial_viz.py`](spatial_viz.py) — reads only `.obs` + `.obsm['spatial']`
(backed mode, no expression matrix loaded) and produces:

1. an 8-panel per-sample spatial scatter colored by the 22 cell types (shared
   legend);
2. the T-cell fraction per CellCharter niche, with enrichment vs the global
   T-cell fraction;
3. a highlight map of the flagged T-cell-rich niches against a greyed background.

## How to run

```bash
python spatial_viz.py
```

Needs `anndata`, `numpy`, `pandas`, `matplotlib`. On a cluster, submit as a
short CPU job (see [`../SLURM_TEMPLATE.md`](../SLURM_TEMPLATE.md)).

## Results

- **Global T-cell fraction = 8.73%.** A niche is flagged as T-cell-rich at ≥ 1.5×
  the global fraction.
- **Flagged:** *lymphoid structure* (T fraction 29.7%, 3.41×) and *immune*
  (27.5%, 3.15×). The lymphoid-structure (TLS-like) niche ranks #1 — confirming
  it as the positive control for organized T-cell programs.
- The two flagged niches differ in subtype composition (lymphoid structure is
  CD4-naive-dominated; the immune niche carries more CD4 memory and Treg),
  suggesting they carry different activation contexts.

## Note on provenance

The `niche` column is **pre-provided** in the data (CellCharter clustering from
the source). The *T-cell-rich flag* is what this step computes (per-niche
fraction vs global, 1.5× threshold; the same two niches hold at a 2× threshold).

## Outputs

`cosmx_nsclc_celltype_spatial.png`, `cosmx_nsclc_tcell_niche_highlight.png`,
`tcell_fraction_per_niche.csv`.
