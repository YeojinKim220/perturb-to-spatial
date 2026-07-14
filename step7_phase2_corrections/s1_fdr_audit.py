#!/usr/bin/env python3
"""FDR AUDIT: recompute held-out validation with a PROPER PER-TARGET permutation null.
Compares three p-value definitions per (endpoint,target):
  p_pertarget : add-one per-target permutation p = (1 + #{null_k >= obs_k}) / (NPERM+1)   [PRIMARY]
  p_pooled    : fraction of the pooled (all-target) permutation null >= obs                [old, sensitivity]
BH-FDR applied per endpoint across targets, separately for each p definition.
Also dumps the nominal p-value distribution of the previous '23 cytotoxic hits'.
Seed fixed. Writes audit CSV + JSON summary."""
import os, json, numpy as np, pandas as pd, h5py, scipy.sparse as sp
from statsmodels.stats.multitest import multipletests
import scanpy as sc
np.random.seed(0)
BASE="/storage/home/hcoda1/9/ykim3030/scratch/hackathon"
DEF=f"{BASE}/0-Data/perturbseq_marson2025/GWCD4i.DE_stats.h5ad"
H5=f"{BASE}/0-Data/ST_CosMx_NSCLC/cosmx_human_nsclc_clustered.h5ad"
OUT=f"{BASE}/2-Results/phase2"; os.makedirs(OUT,exist_ok=True)
SAMPLES=["LUAD-5 R1","LUAD-5 R2","LUAD-5 R3"]
NPERM=1000
def log(*a): print(*a,flush=True)
STATE_SIGS={
 "cytotoxic":["GZMB","GZMK","GZMA","GZMH","PRF1","NKG7","IFNG"],
 "exhausted":["PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","TOX","ENTPD1","VSIR"],
 "activated":["CD69","IL2RA","ICOS","TNFRSF9","CD28","TNFRSF4","MKI67","HLA-DRA"],
 "Treg":["FOXP3","IL2RA","CTLA4"],"Th1":["TBX21","IFNG","CXCR3","STAT1"]}
ENDPOINTS={"activation_enhancement":("activated",+1),"Th1_enhancement":("Th1",+1),
 "cytotoxic_enhancement":("cytotoxic",+1),"exhaustion_reversal":("exhausted",-1),
 "Treg_reduction":("Treg",-1)}

log("== load DE ==")
f=h5py.File(DEF,"r")
degenes=np.array([g.decode() if isinstance(g,bytes) else g for g in f["var"]["gene_name"][:]])
dgpos={g:i for i,g in enumerate(degenes)}
cc=f["obs"]["culture_condition"]; ccat=[c.decode() for c in cc["categories"][:]]
cond=np.array([ccat[c] for c in cc["codes"][:]])
tg=f["obs"]["target_contrast_gene_name"]; tcat=[c.decode() for c in tg["categories"][:]]
target=np.array([tcat[c] if c>=0 else "NA" for c in tg["codes"][:]])
LFC=f["layers"]["log_fc"]
rest_rows=np.where(cond=="Rest")[0]; rest_t=target[rest_rows]
uniq_t=[t for t in pd.unique(rest_t) if t!="NA"]
rowmap={t:rest_rows[np.where(rest_t==t)[0][0]] for t in uniq_t}

log("== load CosMx ==")
A=sc.read_h5ad(H5); A=A[A.obs["sample"].astype(str).isin(SAMPLES)].copy()
naive=(A.obs["cell_type"].astype(str).to_numpy()=="T CD4 naive")
cg=A.var_names.to_numpy(); cgpos={g:i for i,g in enumerate(cg)}
X=A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(A.X)
for name,gl in STATE_SIGS.items():
    sc.tl.score_genes(A,[g for g in gl if g in cgpos],score_name=f"sig_{name}",use_raw=False)
shared=set(l.strip() for l in open(f"{BASE}/shared_genes_genelevel.txt"))
markers=set(g for gl in STATE_SIGS.values() for g in gl)
heldout=[g for g in shared if g in dgpos and g in cgpos and g not in markers]
ho_de=np.array([dgpos[g] for g in heldout]); ho_cx=np.array([cgpos[g] for g in heldout])
log("heldout genes:",len(heldout))

Xn=X[naive]; Sp={}
for p in STATE_SIGS:
    s=A.obs.loc[naive,f"sig_{p}"].to_numpy(); s=s-s.mean()
    col=np.asarray(Xn[:,ho_cx].todense()); col=col-col.mean(0)
    den=(np.linalg.norm(col,axis=0)*np.linalg.norm(s)); den[den==0]=1
    Sp[p]=(col.T@s)/den
Dmat=np.vstack([np.nan_to_num(LFC[rowmap[t],:][ho_de]) for t in uniq_t])
Dn=Dmat/(np.linalg.norm(Dmat,axis=1,keepdims=True)+1e-12)
rng=np.random.default_rng(0)
T=len(uniq_t)
recs=[]
for ep,(prog,d) in ENDPOINTS.items():
    sig=Sp[prog]; obs=d*(Dn@(sig/(np.linalg.norm(sig)+1e-12)))     # (T,)
    # build (NPERM, T) null matrix ONCE
    NM=np.empty((NPERM,T),dtype=np.float32)
    for j in range(NPERM):
        sp_=sig[rng.permutation(len(sig))]; sp_=sp_/(np.linalg.norm(sp_)+1e-12)
        NM[j]=d*(Dn@sp_)
    # PER-TARGET add-one p
    ge_pt=(NM>=obs[None,:]).sum(0)                  # per-target count
    p_pt=(1.0+ge_pt)/(NPERM+1.0)
    # POOLED p (old method)
    pool=np.sort(NM.ravel())
    idx=np.searchsorted(pool,obs,side="left"); p_pool=(len(pool)-idx)/len(pool)
    p_pool=np.clip(p_pool,1.0/len(pool),1.0)
    for k,t in enumerate(uniq_t):
        recs.append(dict(endpoint=ep,program=prog,direction=d,target=t,
                         obs_favorable=float(obs[k]),perm_ge_count=int(ge_pt[k]),
                         p_pertarget=float(p_pt[k]),p_pooled=float(p_pool[k])))
S=pd.DataFrame(recs)
for col,fcol in [("p_pertarget","fdr_pertarget"),("p_pooled","fdr_pooled")]:
    S[fcol]=np.nan
    for ep in S.endpoint.unique():
        m=S.endpoint==ep; S.loc[m,fcol]=multipletests(S.loc[m,col],method="fdr_bh")[1]
S=S.sort_values(["endpoint","p_pertarget"])
S.to_csv(f"{OUT}/S1_fdr_audit.csv",index=False)

# summary
summ={"NPERM":NPERM,"n_targets":T,"n_heldout":len(heldout),"min_p_pertarget_possible":1.0/(NPERM+1)}
for ep,g in S.groupby("endpoint"):
    summ[ep]={
      "pertarget_fdr_lt0.1":int((g.fdr_pertarget<0.1).sum()),
      "pertarget_fdr_lt0.25":int((g.fdr_pertarget<0.25).sum()),
      "pooled_fdr_lt0.1":int((g.fdr_pooled<0.1).sum()),
      "min_p_pertarget":float(g.p_pertarget.min()),
      "min_fdr_pertarget":float(g.fdr_pertarget.min()),
      "min_fdr_pooled":float(g.fdr_pooled.min()),
      "top5":[{"target":t,"p_pt":round(pp,4),"fdr_pt":round(fp,3),"p_pool":round(ppl,6),"fdr_pool":round(fpl,4)}
              for t,pp,fp,ppl,fpl in zip(g.target.head(5),g.p_pertarget.head(5),g.fdr_pertarget.head(5),
                                         g.p_pooled.head(5),g.fdr_pooled.head(5))]}
json.dump(summ,open(f"{OUT}/S1_fdr_audit_summary.json","w"),indent=1)
print(json.dumps(summ,indent=1))
log("== AUDIT DONE ==")
