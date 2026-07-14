import scanpy as sc, anndata as ad, numpy as np, pandas as pd
import scipy.sparse as sp_sparse
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

COSMX="/storage/scratch1/9/ykim3030/hackathon/0-Data/ST_CosMx_NSCLC/cosmx_human_nsclc_clustered.h5ad"
PB="/storage/scratch1/9/ykim3030/hackathon/0-Data/perturbseq_marson2025/GWCD4i.pseudobulk_merged.h5ad"
SAMPLES=["LUAD-5 R1","LUAD-5 R2","LUAD-5 R3"]; ANCHOR="T CD4 naive"
CONDS=["Rest","Stim8hr","Stim48hr"]

def lognorm_cpm(M):  # M: cells x genes sparse counts -> dense log1p CPM
    M=M.tocsr().astype(np.float64)
    lib=np.asarray(M.sum(1)).ravel(); lib[lib==0]=1
    M=M.multiply((1e4/lib)[:,None]).tocsr()
    return np.log1p(M.toarray())

# ---- CosMx CD4 naive: use RAW counts (layers['counts']) so both platforms start from counts ----
print("load cosmx", flush=True)
ac=sc.read_h5ad(COSMX)
sm=ac.obs["sample"].astype(str).isin(SAMPLES).values
ac=ac[sm].copy()
amask=(ac.obs["cell_type"].astype(str)==ANCHOR).values
anchor_subset_idx=np.where(amask)[0]   # indices into the 3-sample subset
counts = ac.layers["counts"] if "counts" in ac.layers else ac.X
cos_genes=ac.var_names.to_numpy()
cd4_counts = counts[anchor_subset_idx]
cd4_counts = cd4_counts.tocsr() if sp_sparse.issparse(cd4_counts) else sp_sparse.csr_matrix(cd4_counts)
# align niche labels by subset index (assignment may be grouped by sample)
_asg = pd.read_csv("cd4_niche_assignment.csv").set_index("global_cell_index_in_subset")
cd4_niche = _asg.reindex(anchor_subset_idx).reset_index()
assert cd4_niche["cd4_niche"].notna().all(), "niche label alignment failed"
print("cd4 cells", cd4_counts.shape[0], "niche rows", len(cd4_niche), flush=True)

# ---- Perturb-seq NTC per-condition reference (sum counts across NTC pseudobulks -> normalize) ----
print("load pseudobulk", flush=True)
ap=ad.read_h5ad(PB, backed="r")
# PB var_names are Ensembl IDs; use var['gene_name'] (symbols) to match CosMx symbols
pb_genes=ap.var["gene_name"].astype(str).to_numpy()
ntc = (ap.obs["guide_type"].astype(str)=="non-targeting").values
cond = ap.obs["culture_condition"].astype(str).values
ref_counts={}
for cnd in CONDS:
    idx=np.where(ntc & (cond==cnd))[0]
    sub=ap[idx].X
    sub=sub.tocsr() if sp_sparse.issparse(sub) else sp_sparse.csr_matrix(sub)
    ref_counts[cnd]=np.asarray(sub.sum(0)).ravel()  # summed counts across NTC pseudobulks (genes,)
    print(f"  NTC {cnd}: n_pseudobulk={len(idx)}", flush=True)

# ---- shared genes ----
shared=np.array(sorted(set(cos_genes)&set(pb_genes)))
print("shared genes", len(shared), flush=True)
cos_pos={g:i for i,g in enumerate(cos_genes)}
pb_pos={}
for i,g in enumerate(pb_genes):
    if g not in pb_pos: pb_pos[g]=i  # first occurrence
ci=[cos_pos[g] for g in shared]; pi=[pb_pos[g] for g in shared]

# CD4 naive expression on shared genes (log CPM)
cd4_shared = cd4_counts[:, ci]
cd4_log = lognorm_cpm(cd4_shared)   # cells x shared
# reference log CPM on shared genes
ref_log={}
for cnd in CONDS:
    v=ref_counts[cnd][pi].astype(np.float64)
    v=v*(1e4/max(v.sum(),1)); ref_log[cnd]=np.log1p(v)
R=np.vstack([ref_log[c] for c in CONDS])  # 3 x shared

# ---- similarity: Pearson & cosine of each CD4 cell vs each state ref ----
def zscore_rows(A):
    A=A-A.mean(1,keepdims=True); s=A.std(1,keepdims=True); s[s==0]=1; return A/s
Xz=zscore_rows(cd4_log); Rz=zscore_rows(R)
pear = Xz @ Rz.T / cd4_log.shape[1]      # cells x 3 (Pearson r)
def l2n(A): n=np.linalg.norm(A,axis=1,keepdims=True); n[n==0]=1; return A/n
cos = l2n(cd4_log) @ l2n(R).T            # cells x 3 cosine

sim=pd.DataFrame(pear, columns=[f"pearson_{c}" for c in CONDS])
for j,c in enumerate(CONDS): sim[f"cosine_{c}"]=cos[:,j]
sim["cd4_niche"]=cd4_niche["cd4_niche"].to_numpy()
sim["sample"]=cd4_niche["sample"].to_numpy()
sim["argmax_pearson"]=[CONDS[k] for k in pear.argmax(1)]
sim["argmax_cosine"]=[CONDS[k] for k in cos.argmax(1)]
sim.to_csv("cd4_state_similarity_percell.csv", index=False)

print("=== argmax(Pearson) distribution ===", flush=True)
print(sim["argmax_pearson"].value_counts().to_string(), flush=True)
print("=== argmax(cosine) distribution ===", flush=True)
print(sim["argmax_cosine"].value_counts().to_string(), flush=True)
print("=== mean similarity per state ===", flush=True)
summ=pd.DataFrame({
  "mean_pearson":[pear[:,j].mean() for j in range(3)],
  "mean_cosine":[cos[:,j].mean() for j in range(3)],
  "argmax_pearson_frac":[ (np.array(sim['argmax_pearson'])==c).mean() for c in CONDS],
}, index=CONDS)
summ.to_csv("cd4_state_similarity_summary.csv")
print(summ.to_string(), flush=True)

# per-niche mean pearson to each state
niche_state = sim.groupby("cd4_niche")[[f"pearson_{c}" for c in CONDS]].mean()
niche_state["argmax"]=niche_state.idxmax(1).str.replace("pearson_","")
niche_state.round(4).to_csv("cd4_niche_state_similarity.csv")

# ---- FIG A: violin/box of pearson per state ----
fig, ax = plt.subplots(figsize=(7,5))
data=[pear[:,j] for j in range(3)]
parts=ax.violinplot(data, showmeans=True, showextrema=False)
ax.set_xticks([1,2,3]); ax.set_xticklabels(CONDS)
ax.set_ylabel("Pearson r (CD4 naive cell vs state NTC ref, 501 shared genes)")
ax.set_title("CD4 naive similarity to Perturb-seq states\n(LUAD-5 R1-3)")
fig.tight_layout(); fig.savefig("cd4_state_similarity_violin.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# ---- FIG B: UMAP overlay of CD4 naive + 3 state refs in shared-gene space ----
allX=np.vstack([cd4_log, R])
lab=np.array(["CD4 naive"]*cd4_log.shape[0]+[f"NTC {c}" for c in CONDS])
emb=ad.AnnData(X=allX.astype(np.float32))
sc.pp.scale(emb, max_value=10); sc.tl.pca(emb, n_comps=max(2, min(30, allX.shape[0]-1, allX.shape[1]-1)))
sc.pp.neighbors(emb, n_neighbors=15); sc.tl.umap(emb)
U=emb.obsm["X_umap"]
fig, ax=plt.subplots(figsize=(8.5,7))
m=lab=="CD4 naive"
# color CD4 cells by argmax state
cmap_state={"Rest":"#2ca02c","Stim8hr":"#ff7f0e","Stim48hr":"#d62728"}
cellcol=pd.Series(sim["argmax_pearson"].to_numpy()).map(cmap_state).to_numpy()
ax.scatter(U[m,0],U[m,1],s=4,c=cellcol,linewidths=0,rasterized=True,alpha=0.5,label=None)
# ref points big
for j,c in enumerate(CONDS):
    rm=lab==f"NTC {c}"
    ax.scatter(U[rm,0],U[rm,1],s=420,marker="*",c=cmap_state[c],edgecolors="black",linewidths=1.5,zorder=5,label=f"NTC {c}")
from matplotlib.lines import Line2D
h=[Line2D([0],[0],marker='o',linestyle='',markersize=6,markerfacecolor=cmap_state[c],markeredgewidth=0,label=f"CD4 argmax={c}") for c in CONDS]
h+=[Line2D([0],[0],marker='*',linestyle='',markersize=13,markerfacecolor=cmap_state[c],markeredgecolor='black',label=f"NTC {c} ref") for c in CONDS]
ax.legend(handles=h, fontsize=8, loc="best", frameon=True)
ax.set_title("CD4 naive + Perturb-seq NTC state refs\n(shared 501-gene space, log-CPM)"); ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout(); fig.savefig("cd4_state_umap_overlay.png", dpi=180, bbox_inches="tight"); plt.close(fig)
print("ALL DONE", flush=True)
