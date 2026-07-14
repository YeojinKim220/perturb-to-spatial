import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

H5 = "/storage/scratch1/9/ykim3030/hackathon/0-Data/ST_CosMx_NSCLC/cosmx_human_nsclc_clustered.h5ad"

# ---- load only obs + spatial (no expression needed for these figures) ----
a = ad.read_h5ad(H5, backed="r")
obs = a.obs[["sample", "cell_type", "niche"]].copy()
sp = np.asarray(a.obsm["spatial"])
obs["x"] = sp[:, 0]
obs["y"] = sp[:, 1]
del a

# sample order (fixed, readable)
sample_order = ["LUAD-5 R1", "LUAD-5 R2", "LUAD-5 R3", "LUSC-6",
                "LUAD-9 R1", "LUAD-9 R2", "LUAD-12", "LUAD-13"]
sample_order = [s for s in sample_order if s in set(obs["sample"])]

T_TYPES = ["T CD4 naive", "T CD4 memory", "T CD8 naive", "T CD8 memory", "Treg"]

# ---- deterministic cell-type color map (22 types) ----
cell_types = sorted(obs["cell_type"].unique().tolist())
# put T cells on vivid reds/oranges, tumors on greys, rest on tab20
import matplotlib.cm as cm
from matplotlib.colors import to_hex
base = [to_hex(x) for x in (list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors))]
ct_color = {}
# T cells: distinct warm palette
t_palette = {"T CD4 naive": "#e41a1c", "T CD4 memory": "#ff7f00",
             "T CD8 naive": "#984ea3", "T CD8 memory": "#f781bf", "Treg": "#a65628"}
tumor_grey = ["#8c8c8c", "#a6a6a6", "#bdbdbd", "#d9d9d9"]
ti = 0; bi = 0
for ct in cell_types:
    if ct in t_palette:
        ct_color[ct] = t_palette[ct]
    elif ct.lower().startswith("tumor"):
        ct_color[ct] = tumor_grey[ti % len(tumor_grey)]; ti += 1
    else:
        ct_color[ct] = base[bi % len(base)]; bi += 1

# =========================================================
# FIGURE 1: multi-panel spatial map colored by cell type
# =========================================================
ncol = 4; nrow = 2
fig, axes = plt.subplots(nrow, ncol, figsize=(20, 11))
axes = axes.ravel()
for i, s in enumerate(sample_order):
    ax = axes[i]
    d = obs[obs["sample"] == s]
    cols = d["cell_type"].astype(str).map(ct_color).values
    ax.scatter(d["x"], d["y"], s=0.8, c=cols, linewidths=0, rasterized=True)
    ax.set_title(f"{s}  (n={len(d):,})", fontsize=11)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.invert_yaxis()
for j in range(len(sample_order), len(axes)):
    axes[j].axis("off")
handles = [Line2D([0],[0], marker='o', linestyle='', markersize=7,
                  markerfacecolor=ct_color[ct], markeredgewidth=0, label=ct)
           for ct in cell_types]
fig.legend(handles=handles, loc="center right", fontsize=9, frameon=False,
           title="Cell type", title_fontsize=11, ncol=1,
           bbox_to_anchor=(1.0, 0.5))
fig.suptitle("CosMx NSCLC — spatial map by cell type (8 samples)", fontsize=15, y=0.995)
fig.tight_layout(rect=[0, 0, 0.86, 0.98])
fig.savefig("cosmx_nsclc_celltype_spatial.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("FIG1 done")

# =========================================================
# TABLE: T-cell fraction per niche
# =========================================================
obs["is_T"] = obs["cell_type"].isin(T_TYPES)
g = obs.groupby("niche", observed=True)
tab = pd.DataFrame({
    "n_cells": g.size(),
    "n_Tcells": g["is_T"].sum(),
})
tab["Tcell_fraction"] = tab["n_Tcells"] / tab["n_cells"]
global_frac = obs["is_T"].mean()
tab["global_Tcell_fraction"] = global_frac
tab["enrichment_vs_global"] = tab["Tcell_fraction"] / global_frac
# per-subtype fraction within niche
for tt in T_TYPES:
    frac = obs.assign(f=obs["cell_type"].eq(tt)).groupby("niche", observed=True)["f"].mean()
    tab[f"frac_{tt}"] = frac
tab = tab.sort_values("Tcell_fraction", ascending=False)
# flag T-cell-rich: enrichment >= 1.5x global
tab["Tcell_rich_flag"] = tab["enrichment_vs_global"] >= 1.5
tab = tab.reset_index()
tab.to_csv("tcell_fraction_per_niche.csv", index=False)
flagged = tab.loc[tab["Tcell_rich_flag"], "niche"].tolist()
print("GLOBAL T-cell fraction: %.4f" % global_frac)
print(tab[["niche","n_cells","n_Tcells","Tcell_fraction","enrichment_vs_global","Tcell_rich_flag"]].to_string(index=False))
print("FLAGGED_NICHES:", flagged)

# =========================================================
# FIGURE 2: highlight flagged T-cell-rich niches per sample
# =========================================================
niche_all = sorted(obs["niche"].unique().tolist())
nfl = len(flagged)
niche_palette = [to_hex(x) for x in (list(plt.get_cmap("Set1").colors) + list(plt.get_cmap("Set2").colors))]
niche_color = {n: niche_palette[i % len(niche_palette)] for i, n in enumerate(flagged)}

fig, axes = plt.subplots(nrow, ncol, figsize=(20, 11))
axes = axes.ravel()
for i, s in enumerate(sample_order):
    ax = axes[i]
    d = obs[obs["sample"] == s]
    infl = d["niche"].isin(flagged)
    # background greyed
    ax.scatter(d.loc[~infl,"x"], d.loc[~infl,"y"], s=0.5, c="#e8e8e8",
               linewidths=0, rasterized=True)
    # highlighted niche cells
    dh = d[infl]
    if len(dh):
        cols = dh["niche"].astype(str).map(niche_color).values
        ax.scatter(dh["x"], dh["y"], s=1.2, c=cols, linewidths=0, rasterized=True)
    ax.set_title(f"{s}  (T-niche n={int(infl.sum()):,})", fontsize=11)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.invert_yaxis()
for j in range(len(sample_order), len(axes)):
    axes[j].axis("off")
handles = [Line2D([0],[0], marker='o', linestyle='', markersize=8,
                  markerfacecolor=niche_color[n], markeredgewidth=0, label=n)
           for n in flagged]
handles.append(Line2D([0],[0], marker='o', linestyle='', markersize=8,
                       markerfacecolor="#e8e8e8", markeredgewidth=0, label="other niches"))
fig.legend(handles=handles, loc="center right", fontsize=10, frameon=False,
           title="T-cell-rich niche", title_fontsize=11,
           bbox_to_anchor=(1.0, 0.5))
fig.suptitle("CosMx NSCLC — flagged T-cell-rich niches (enrichment \u2265 1.5\u00d7 global)",
             fontsize=15, y=0.995)
fig.tight_layout(rect=[0, 0, 0.85, 0.98])
fig.savefig("cosmx_nsclc_tcell_niche_highlight.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("FIG2 done")
print("ALL DONE")
