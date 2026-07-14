
import h5py, numpy as np
base="/storage/home/hcoda1/9/ykim3030/scratch/hackathon/0-Data"
def decode(a): return [x.decode() if isinstance(x,bytes) else str(x) for x in a]
def cats(f, col):
    # anndata categorical: group with 'categories' and 'codes'
    g=f["obs"][col]
    if isinstance(g,h5py.Group):
        c=decode(g["categories"][:]); return c
    return decode(np.unique(g[:]))
# DE genes
with h5py.File(base+"/perturbseq_marson2025/GWCD4i.DE_stats.h5ad","r") as f:
    de_names=set(decode(f["var"]["gene_name"][:]))
    de_conds=cats(f,"culture_condition")
    n_pert=f["obs"]["target_contrast_gene_name"]
    npert = len(decode(n_pert["categories"][:])) if isinstance(n_pert,h5py.Group) else "?"
print("DE genes measured:", len(de_names))
print("DE culture conditions:", de_conds)
print("DE n unique perturbed genes:", npert)
# CosMx genes + celltypes + niches
with h5py.File(base+"/ST_CosMx_NSCLC/cosmx_human_nsclc_clustered.h5ad","r") as f:
    cx_genes=set(decode(f["var"]["_index"][:]))
    ct=cats(f,"cell_type"); ni=cats(f,"niche"); samp=cats(f,"sample")
    ncells=f["obs"]["_index"].shape[0]
print("CosMx panel genes:", len(cx_genes))
print("CosMx n cells:", ncells)
print("CosMx cell_types:", ct)
print("CosMx niches:", ni)
print("CosMx samples:", samp)
# negative controls in CosMx panel (SystemControl/NegPrb)
negs=[g for g in cx_genes if g.lower().startswith(("negprb","systemcontrol","neg_","control"))]
print("CosMx control/neg probes:", len(negs), negs[:8])
ov = de_names & cx_genes
print("OVERLAP DE x CosMx:", len(ov), "of", len(cx_genes), "panel genes")
