# Step 5 — Cross-platform state matching (which Perturb-seq condition do tumor CD4⁺ cells resemble?)

Establish which Perturb-seq condition (Rest / Stim8hr / Stim48hr) the tumor CD4⁺
T cells are closest to, so the downstream projection knows which condition to
lean on.

## Script

[`state_match.py`](state_match.py). Both platforms are brought back to a common
starting point (raw counts → log-normalized CPM) on the **912 shared genes**
(CosMx ∩ pseudobulk). Each tumor CD4-naive cell's expression profile is
correlated (Pearson) against the pseudobulk profile of each condition, and the
best-matching condition is taken per cell and per niche.

## How to run

```bash
python state_match.py
```

Needs `scanpy`, `anndata`, `scipy`, `matplotlib`.

## Result

Tumor CD4-naive cells are closest to the Perturb-seq **Rest** state: the
best-matching condition is Rest for **81.7%** of cells and for all niches (mean
Pearson: Rest 0.194 > Stim48hr 0.175 > Stim8hr 0.170). This is the biologically
expected result (tumor-resident naive cells are not acutely TCR-stimulated).

## Consequence for the rest of the project

Because per-cell matching collapses so strongly onto Rest, a per-cell "soft"
projection would be almost identical to a plain Rest effect and carries little
extra information. The activation-context question is therefore anchored not on a
soft projection but on the **Rest-vs-Stim direction comparison** — regulators
whose effect *reverses* between resting and stimulated states — which is the
subject of [`step7`](../step7_phase2_corrections/) and
[`step8`](../step8_phase3_reanalysis/). Rest is used as the primary condition
throughout; Stim48hr is always reported alongside with a direction-flip flag.

## Outputs

`cd4_state_similarity_percell.csv`, `cd4_state_similarity_summary.csv`,
`cd4_state_similarity_violin.png`, `cd4_state_umap_overlay.png`.
