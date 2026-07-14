import scanpy as sc, anndata as ad, numpy as np, pandas as pd
import scipy.sparse as sp_sparse
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from matplotlib.lines import Line2D

H5="/storage/scratch1/9/ykim3030/hackathon/0-Data/ST_CosMx_NSCLC/cosmx_human_nsclc_clustered.h5ad"
SIGMA=50.0; RADIUS=150.0; ANCHOR="T CD4 naive"; SAMPLES=["LUAD-5 R1","LUAD-5 R2","LUAD-5 R3"]
TUMOR=["tumor 5","tumor 6","tumor 9","tumor 12","tumor 13","epithelial"]

# state signatures (only panel-present genes, from diagnostic)
STATE_SIGS={
 "cytotoxic":["GZMB","GZMK","GZMA","GZMH","PRF1","NKG7","IFNG"],
 "exhausted":["PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","TOX","ENTPD1","VSIR"],
 "activated":["CD69","IL2RA","ICOS","TNFRSF9","CD28","TNFRSF4","MKI67","HLA-DRA"],
 "naive":["CCR7","SELL","IL7R"],
 "Treg":["FOXP3","IL2RA","CTLA4"],
 "Tfh":["CXCR5","PDCD1","ICOS"],
 "Th1":["TBX21","IFNG","CXCR3","STAT1"],
}
# suppressive/ligand program aggregated over NEIGHBORS
SUPP_GENES=["CD274","PDCD1LG2","TGFB1","IL10","IDO1","ARG1","CD163","LGALS9","VSIR",
            "CXCL9","CXCL10","CXCL12","CCL2","CCL5","CCL19","CCL21"]

print("load", flush=True)
a=sc.read_h5ad(H5)
a=a[a.obs["sample"].astype(str).isin(SAMPLES)].copy()
sp=np.asarray(a.obsm["spatial"]); samples=a.obs["sample"].astype(str).to_numpy()
ct=a.obs["cell_type"].astype(str)
X=a.X.tocsr() if sp_sparse.issparse(a.X) else sp_sparse.csr_matrix(a.X)  # log-norm
genes=a.var_names.to_numpy(); gpos={g:i for i,g in enumerate(genes)}
is_anchor=(ct.to_numpy()==ANCHOR); is_tumor=ct.isin(TUMOR).to_numpy()

# ---- Block 1: T-cell state signature scores (self) via sc.tl.score_genes ----
for name,gl in STATE_SIGS.items():
    gl=[g for g in gl if g in gpos]
    sc.tl.score_genes(a, gl, score_name=f"sig_{name}", use_raw=False)
B1_cols=[f"sig_{k}" for k in STATE_SIGS]
B1 = a.obs.loc[is_anchor, B1_cols].to_numpy()
anchor_idx=np.where(is_anchor)[0]

# ---- Blocks 2,3,4 per anchor via within-sample KDTree ----
supp_pos=[gpos[g] for g in SUPP_GENES if g in gpos]
ct_cats=sorted(ct.unique()); ct_idx={c:i for i,c in enumerate(ct_cats)}
ct_code=ct.map(ct_idx).to_numpy()
B2=np.zeros((len(anchor_idx),len(ct_cats)))            # neighbor composition
B3=np.zeros((len(anchor_idx),len(supp_pos)))           # neighbor suppressive program
B4=np.zeros(len(anchor_idx))                           # dist to nearest tumor
pos_in_anchor={g:k for k,g in enumerate(anchor_idx)}
row=0
for s in SAMPLES:
    smask=samples==s; idx_all=np.where(smask)[0]; coords=smask_coords=sp[smask]
    codes=ct_code[smask]; Xs=X[smask]; tum_local=is_tumor[smask]
    a_local=np.where(is_anchor[smask])[0]
    tree=cKDTree(coords)
    # neighbor query for composition + suppressive
    nb=tree.query_ball_point(coords[a_local], r=RADIUS, workers=-1)
    # tumor tree for distance
    if tum_local.sum()>0:
        ttree=cKDTree(coords[tum_local])
        dd,_=ttree.query(coords[a_local], k=1, workers=-1)
    else:
        dd=np.full(len(a_local), np.nan)
    Xs_supp=Xs[:,supp_pos].toarray()
    for k,ind in enumerate(nb):
        ind=np.asarray(ind); d=np.linalg.norm(coords[ind]-coords[a_local[k]],axis=1)
        w=np.exp(-(d**2)/(2*SIGMA**2)); w/=w.sum()
        comp=np.zeros(len(ct_cats)); np.add.at(comp, codes[ind], w)
        B2[row]=comp
        B3[row]=(Xs_supp[ind]*w[:,None]).sum(0)
        B4[row]=dd[k]
        row+=1
    print(f"[{s}] anchor={len(a_local)} tumor={int(tum_local.sum())} done", flush=True)

B2_cols=[f"nbr__{c}" for c in ct_cats]
B3_cols=[f"supp__{genes[p]}" for p in supp_pos]
# log-transform distance (skewed), fill nan with max
B4=np.where(np.isnan(B4), np.nanmax(B4), B4); B4log=np.log1p(B4)

# ---- assemble: z-score each block, concat ----
def z(A):
    A=np.asarray(A,float); m=A.mean(0); s=A.std(0); s[s==0]=1; return (A-m)/s
# block-count normalization: each z-scored block divided by sqrt(n_features)
# so all 4 blocks contribute equally to Euclidean distance (distance block, n=1, gets full weight)
def zb(A):
    A=z(A); n=A.shape[1] if A.ndim>1 else 1
    return A/np.sqrt(n)
F=np.hstack([zb(B1), zb(B2), zb(B3), zb(B4log[:,None])])
Fcols=B1_cols+B2_cols+B3_cols+["dist_tumor_log"]
print("feature matrix", F.shape, "block-count normalized (each block sqrt(n)-scaled)", flush=True)

# ---- Leiden ----
fn=ad.AnnData(X=F.astype(np.float32), obs=pd.DataFrame({
    "sample":samples[anchor_idx],"niche_ref":a.obs["niche"].astype(str).to_numpy()[anchor_idx],
}, index=[f"C{i}" for i in anchor_idx]))
fn.var_names=Fcols
sc.pp.pca(fn, n_comps=min(30,F.shape[1]-1))
sc.pp.neighbors(fn, n_neighbors=15)
sc.tl.leiden(fn, resolution=1.0, key_added="fniche", flavor="igraph", n_iterations=2, directed=False)
fn.obs["fniche"]="FN"+fn.obs["fniche"].astype(str)
sc.tl.umap(fn)
nN=fn.obs["fniche"].nunique(); print("N feature-niches:", nN, flush=True)
print(fn.obs["fniche"].value_counts().to_string(), flush=True)

# ---- per-niche feature means (for naming) ----
featdf=pd.DataFrame(F, columns=Fcols); featdf["fniche"]=fn.obs["fniche"].to_numpy()
featdf["dist_tumor_um"]=B4  # raw distance
niche_prof=featdf.groupby("fniche").mean(numeric_only=True)
# keep interpretable columns: states, top neighbor types, suppressive mean, distance
tum_cols=[c for c in B2_cols if any(t in c for t in ["tumor","epithelial"])]
key_cols=B1_cols+["nbr__macrophage","nbr__monocyte","nbr__mDC","nbr__fibroblast",
                  "nbr__endothelial","nbr__B-cell","nbr__T CD4 naive"]+["dist_tumor_um"]
key_cols=[c for c in key_cols if c in niche_prof.columns]
niche_prof["supp_mean"]=featdf.groupby("fniche")[B3_cols].mean().mean(1)
niche_prof["tumor_nbr_frac"]=featdf.groupby("fniche")[tum_cols].mean().sum(1)
niche_prof.round(4).to_csv("cd4_feature_niche_profile.csv")

# auto-name by dominant axis
def name_niche(r):
    tags=[]
    st={k:r[f"sig_{k}"] for k in STATE_SIGS}
    top_state=max(st,key=st.get)
    tags.append(top_state)
    if r["tumor_nbr_frac"]>niche_prof["tumor_nbr_frac"].median(): tags.append("tumor-contact")
    if r["dist_tumor_um"]>niche_prof["dist_tumor_um"].median(): tags.append("tumor-far")
    if r["nbr__macrophage"]+r.get("nbr__monocyte",0)>0.15: tags.append("myeloid-rich")
    if r["nbr__fibroblast"]>0.15: tags.append("fibroblast-rich")
    if r.get("nbr__B-cell",0)>0.15: tags.append("Bcell/TLS")
    if r["supp_mean"]>niche_prof["supp_mean"].median(): tags.append("suppressive-hi")
    return "; ".join(tags)
niche_prof["auto_label"]=niche_prof.apply(name_niche,axis=1)
lab_df=niche_prof[["auto_label","tumor_nbr_frac","dist_tumor_um","supp_mean"]+B1_cols].round(3)
lab_df.to_csv("cd4_feature_niche_labels.csv")
print(lab_df["auto_label"].to_string(), flush=True)

# assignment
assign=pd.DataFrame({"global_cell_index_in_subset":anchor_idx,"sample":samples[anchor_idx],
                     "niche_ref":a.obs["niche"].astype(str).to_numpy()[anchor_idx],
                     "fniche":fn.obs["fniche"].to_numpy(),"dist_tumor_um":B4})
for cc in B1_cols: assign[cc]=B1[:, B1_cols.index(cc)]
assign.to_csv("cd4_feature_niche_assignment.csv", index=False)

# ---- FIG A: heatmap of per-niche z-scored feature profile (states + key nbrs + dist + supp) ----
show_cols=B1_cols+["nbr__macrophage","nbr__fibroblast","nbr__endothelial","nbr__B-cell",
                   "tumor_nbr_frac","supp_mean","dist_tumor_um"]
show_cols=[c for c in show_cols if c in niche_prof.columns]
Z=niche_prof[show_cols].copy()
Zz=(Z-Z.mean())/Z.std().replace(0,1)
fig,ax=plt.subplots(figsize=(max(8,0.7*len(show_cols)), max(5,0.4*nN)))
im=ax.imshow(Zz.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
ax.set_xticks(range(len(show_cols))); ax.set_xticklabels([c.replace("sig_","").replace("nbr__","nbr:") for c in show_cols], rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(nN)); ax.set_yticklabels([f"{i} {niche_prof.loc[i,'auto_label'][:32]}" for i in niche_prof.index], fontsize=7)
fig.colorbar(im,ax=ax,label="z-score across niches"); ax.set_title("CD4 naive feature-niche profiles (LUAD-5 R1-3)")
fig.tight_layout(); fig.savefig("cd4_feature_niche_heatmap.png",dpi=200,bbox_inches="tight"); plt.close(fig)

# ---- FIG B: UMAP colored by fniche + by distance ----
fig,axes=plt.subplots(1,2,figsize=(16,7)); U=fn.obsm["X_umap"]
cats=sorted(fn.obs["fniche"].unique(),key=lambda x:int(x[2:]))
pal=[to_hex(c) for c in (list(plt.get_cmap("tab20").colors)+list(plt.get_cmap("tab20b").colors))]
cmap={c:pal[i%len(pal)] for i,c in enumerate(cats)}
axes[0].scatter(U[:,0],U[:,1],s=3,c=pd.Series(fn.obs["fniche"].to_numpy()).map(cmap).to_numpy(),linewidths=0,rasterized=True)
h=[Line2D([0],[0],marker='o',linestyle='',markersize=6,markerfacecolor=cmap[c],markeredgewidth=0,label=c) for c in cats]
axes[0].legend(handles=h,fontsize=7,loc="center left",bbox_to_anchor=(1,0.5),frameon=False); axes[0].set_title("feature-niche"); axes[0].set_xticks([]);axes[0].set_yticks([])
sc2=axes[1].scatter(U[:,0],U[:,1],s=3,c=np.clip(B4,0,np.percentile(B4,99)),cmap="viridis",linewidths=0,rasterized=True)
fig.colorbar(sc2,ax=axes[1],label="dist to tumor (um)"); axes[1].set_title("distance to tumor"); axes[1].set_xticks([]);axes[1].set_yticks([])
fig.suptitle("CD4 naive feature-niche UMAP (LUAD-5 R1-3)",fontsize=14)
fig.tight_layout(); fig.savefig("cd4_feature_niche_umap.png",dpi=180,bbox_inches="tight"); plt.close(fig)

# ---- FIG C: spatial per sample ----
spA=sp[anchor_idx]
fig,axes=plt.subplots(1,3,figsize=(18,6.5)); axes=axes.ravel()
for i,s in enumerate(SAMPLES):
    ax=axes[i]; m=assign["sample"].to_numpy()==s
    ax.scatter(spA[m,0],spA[m,1],s=6,c=pd.Series(assign["fniche"].to_numpy()[m]).map(cmap).to_numpy(),linewidths=0,rasterized=True)
    ax.set_title(f"{s} (n={int(m.sum()):,})",fontsize=11); ax.set_aspect("equal"); ax.set_xticks([]);ax.set_yticks([]); ax.invert_yaxis()
h=[Line2D([0],[0],marker='o',linestyle='',markersize=7,markerfacecolor=cmap[c],markeredgewidth=0,label=c) for c in cats]
fig.legend(handles=h,fontsize=8,loc="center right",bbox_to_anchor=(1,0.5),frameon=False,title="feature-niche")
fig.suptitle("CD4 naive feature-niches — spatial (LUAD-5 R1-3)",fontsize=15,y=1.0)
fig.tight_layout(rect=[0,0,0.9,0.97]); fig.savefig("cd4_feature_niche_spatial.png",dpi=200,bbox_inches="tight"); plt.close(fig)
print("ALL DONE", flush=True)
