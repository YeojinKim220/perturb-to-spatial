import scanpy as sc, anndata as ad, numpy as np, pandas as pd
import scipy.sparse as sp_sparse
from scipy.spatial import cKDTree
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

H5 = "/storage/scratch1/9/ykim3030/hackathon/0-Data/ST_CosMx_NSCLC/cosmx_human_nsclc_clustered.h5ad"
SIGMA = 50.0          # Gaussian bandwidth (um) ~ median 1-NN dist 46
RADIUS = 3 * SIGMA    # 150 um cutoff
T_TYPES = ["T CD4 naive", "T CD4 memory", "T CD8 naive", "T CD8 memory", "Treg"]

print("loading...", flush=True)
a = sc.read_h5ad(H5)
X = a.X.tocsr() if sp_sparse.issparse(a.X) else sp_sparse.csr_matrix(a.X)
genes = a.var_names.to_numpy()
sp = np.asarray(a.obsm["spatial"])
samples = a.obs["sample"].astype(str).to_numpy()
celltypes = a.obs["cell_type"].astype(str).to_numpy()
niche_ref = a.obs["niche"].astype(str).to_numpy()
is_T = np.isin(celltypes, T_TYPES)
print("total cells", a.n_obs, "T cells", int(is_T.sum()), flush=True)

# ---- neighborhood expression for each T cell, from ALL neighbors, Gaussian weighted ----
NE_rows = []      # dense neighborhood expression per T cell
T_index = []      # global index of each T cell
for s in pd.unique(samples):
    smask = samples == s
    idx_all = np.where(smask)[0]
    coords_all = sp[smask]
    Xs = X[smask]                        # cells_in_sample x genes (normalized)
    tmask_local = is_T[smask]
    t_local = np.where(tmask_local)[0]
    if len(t_local) == 0:
        continue
    tree = cKDTree(coords_all)
    qpts = coords_all[t_local]
    neigh = tree.query_ball_point(qpts, r=RADIUS, workers=-1)
    for k, nb in enumerate(neigh):
        nb = np.asarray(nb)
        d = np.linalg.norm(coords_all[nb] - qpts[k], axis=1)
        w = np.exp(-(d**2) / (2 * SIGMA**2))
        w = w / w.sum()
        ne = (Xs[nb].multiply(w[:, None])).sum(axis=0)   # 1 x genes
        NE_rows.append(np.asarray(ne).ravel())
        T_index.append(idx_all[t_local[k]])
    print(f"[{s}] T={len(t_local)} done", flush=True)

NE = np.vstack(NE_rows).astype(np.float32)
T_index = np.asarray(T_index)
print("NE matrix", NE.shape, flush=True)

# ---- build AnnData of T-cell neighborhoods, cluster with Leiden ----
tn = ad.AnnData(X=NE, obs=pd.DataFrame({
    "sample": samples[T_index],
    "cell_type": celltypes[T_index],
    "niche_ref": niche_ref[T_index],
}, index=[f"T{i}" for i in T_index]))
tn.var_names = genes
sc.pp.scale(tn, max_value=10)
sc.tl.pca(tn, n_comps=30)
sc.pp.neighbors(tn, n_neighbors=15, use_rep="X_pca")
sc.tl.leiden(tn, resolution=1.0, key_added="tcell_niche", flavor="igraph", n_iterations=2, directed=False)
tn.obs["tcell_niche"] = "TN" + tn.obs["tcell_niche"].astype(str)
n_niche = tn.obs["tcell_niche"].nunique()
print("N T-cell niches (leiden res=1.0):", n_niche, flush=True)
print(tn.obs["tcell_niche"].value_counts().to_string(), flush=True)

# UMAP for viz
sc.tl.umap(tn)

# ---- compare with reference CellCharter niche ----
ari = adjusted_rand_score(tn.obs["niche_ref"], tn.obs["tcell_niche"])
nmi = normalized_mutual_info_score(tn.obs["niche_ref"], tn.obs["tcell_niche"])
print(f"ARI vs CellCharter niche = {ari:.4f}", flush=True)
print(f"NMI vs CellCharter niche = {nmi:.4f}", flush=True)

# crosstab (row-normalized): each new niche's composition over reference niches
ct = pd.crosstab(tn.obs["tcell_niche"], tn.obs["niche_ref"])
ct_norm = ct.div(ct.sum(1), axis=0)
ct.to_csv("tcell_niche_vs_cellcharter_counts.csv")
ct_norm.round(4).to_csv("tcell_niche_vs_cellcharter_rownorm.csv")

# annotate each new niche by dominant reference niche + dominant T subtype
anno = []
for tnid in sorted(tn.obs["tcell_niche"].unique(), key=lambda x:int(x[2:])):
    sub = tn.obs[tn.obs["tcell_niche"] == tnid]
    dom_ref = sub["niche_ref"].value_counts().idxmax()
    dom_ref_frac = sub["niche_ref"].value_counts(normalize=True).max()
    dom_t = sub["cell_type"].value_counts().idxmax()
    dom_t_frac = sub["cell_type"].value_counts(normalize=True).max()
    anno.append({"tcell_niche": tnid, "n": len(sub),
                 "dominant_ref_niche": dom_ref, "ref_frac": round(dom_ref_frac,3),
                 "dominant_T_subtype": dom_t, "T_frac": round(dom_t_frac,3)})
anno_df = pd.DataFrame(anno)
anno_df.to_csv("tcell_niche_annotation.csv", index=False)
print(anno_df.to_string(index=False), flush=True)

# save the T-niche assignment back (global index + label)
assign = pd.DataFrame({"global_cell_index": T_index,
                       "sample": samples[T_index],
                       "cell_type": celltypes[T_index],
                       "niche_ref": niche_ref[T_index],
                       "tcell_niche": tn.obs["tcell_niche"].to_numpy()})
assign.to_csv("tcell_niche_assignment.csv", index=False)

# ---- FIGURE A: heatmap crosstab (row-normalized) ----
fig, ax = plt.subplots(figsize=(10, max(4, 0.5*n_niche)))
import numpy as _np
mat = ct_norm.values
im = ax.imshow(mat, aspect="auto", cmap="viridis", vmin=0, vmax=1)
ax.set_xticks(range(ct_norm.shape[1])); ax.set_xticklabels(ct_norm.columns, rotation=45, ha="right", fontsize=9)
ax.set_yticks(range(ct_norm.shape[0])); ax.set_yticklabels(ct_norm.index, fontsize=9)
ax.set_xlabel("CellCharter niche (reference)"); ax.set_ylabel("T-cell niche (new, Leiden)")
ax.set_title(f"T-cell niche vs CellCharter niche (row-normalized)\nARI={ari:.3f}  NMI={nmi:.3f}")
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if mat[i,j] >= 0.15:
            ax.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center",
                    color="white" if mat[i,j]<0.6 else "black", fontsize=7)
fig.colorbar(im, ax=ax, label="fraction of T-cell-niche cells")
fig.tight_layout()
fig.savefig("tcell_niche_vs_cellcharter_heatmap.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# ---- FIGURE B: UMAP colored by new niche and by reference niche ----
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
for ax, key, ttl in [(axes[0],"tcell_niche","T-cell niche (Leiden)"),
                     (axes[1],"niche_ref","CellCharter niche (reference)")]:
    cats = sorted(tn.obs[key].unique())
    pal = [to_hex(c) for c in (list(plt.get_cmap("tab20").colors)+list(plt.get_cmap("tab20b").colors))]
    cmap = {c: pal[i % len(pal)] for i, c in enumerate(cats)}
    U = tn.obsm["X_umap"]
    cols = pd.Series(tn.obs[key].to_numpy()).map(cmap).to_numpy()
    ax.scatter(U[:,0], U[:,1], s=1.5, c=cols, linewidths=0, rasterized=True)
    from matplotlib.lines import Line2D
    h=[Line2D([0],[0],marker='o',linestyle='',markersize=6,markerfacecolor=cmap[c],markeredgewidth=0,label=c) for c in cats]
    ax.legend(handles=h, fontsize=7, loc="center left", bbox_to_anchor=(1.0,0.5), frameon=False)
    ax.set_title(ttl); ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("T-cell neighborhood-expression UMAP", fontsize=14)
fig.tight_layout()
fig.savefig("tcell_niche_umap.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# ---- FIGURE C: spatial map of new T-cell niches, per sample ----
sample_order = ["LUAD-5 R1","LUAD-5 R2","LUAD-5 R3","LUSC-6","LUAD-9 R1","LUAD-9 R2","LUAD-12","LUAD-13"]
sample_order = [s for s in sample_order if s in set(assign["sample"])]
cats = sorted(tn.obs["tcell_niche"].unique(), key=lambda x:int(x[2:]))
pal = [to_hex(c) for c in (list(plt.get_cmap("tab20").colors)+list(plt.get_cmap("tab20b").colors))]
cmap = {c: pal[i % len(pal)] for i,c in enumerate(cats)}
spT = sp[T_index]
fig, axes = plt.subplots(2,4, figsize=(20,11)); axes=axes.ravel()
for i,s in enumerate(sample_order):
    ax=axes[i]; m = assign["sample"].to_numpy()==s
    cols = pd.Series(assign["tcell_niche"].to_numpy()[m]).map(cmap).to_numpy()
    ax.scatter(spT[m,0], spT[m,1], s=3, c=cols, linewidths=0, rasterized=True)
    ax.set_title(f"{s} (T={int(m.sum()):,})", fontsize=11)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.invert_yaxis()
for j in range(len(sample_order),len(axes)): axes[j].axis("off")
from matplotlib.lines import Line2D
h=[Line2D([0],[0],marker='o',linestyle='',markersize=7,markerfacecolor=cmap[c],markeredgewidth=0,label=c) for c in cats]
fig.legend(handles=h, fontsize=9, loc="center right", bbox_to_anchor=(1.0,0.5), frameon=False, title="T-cell niche")
fig.suptitle("T-cell niches (neighborhood-expression Leiden) — spatial map", fontsize=15, y=0.995)
fig.tight_layout(rect=[0,0,0.9,0.98])
fig.savefig("tcell_niche_spatial.png", dpi=200, bbox_inches="tight")
plt.close(fig)

print("ALL DONE", flush=True)
