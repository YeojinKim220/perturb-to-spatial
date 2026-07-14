#!/usr/bin/env python3
"""Phase-2 additional analyses (methodology-review corrections).
Stages (each writes its own outputs; late failure keeps earlier ones):
 S1 held-out direction-corrected validation + BH-FDR (rank_pct, not p)
 S2 within-state niche-contrast signature (tumor-contact vs tumor-far) + perturbation rank-shift
 S3 perturb-seq reliability + viability/general-disruption filter annotation
 S4 spatial quantification: Moran's I per program/section, neighborhood enrichment + permutation null,
    min-niche-size filter, section (R1-R3) rank reproducibility
All: GRCh38 symbols, seed=0, shared-gene overlap reported per score.
"""
import os, json, numpy as np, pandas as pd, h5py, scipy.sparse as sp
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests
np.random.seed(0)
BASE="/storage/home/hcoda1/9/ykim3030/scratch/hackathon"
DEF=f"{BASE}/0-Data/perturbseq_marson2025/GWCD4i.DE_stats.h5ad"
H5=f"{BASE}/0-Data/ST_CosMx_NSCLC/cosmx_human_nsclc_clustered.h5ad"
ASSIGN=f"{BASE}/2-Results/cd4_feature_niche/cd4_feature_niche_assignment.csv"
OUT=f"{BASE}/2-Results/phase2"; os.makedirs(OUT,exist_ok=True)
SAMPLES=["LUAD-5 R1","LUAD-5 R2","LUAD-5 R3"]
def log(*a): print(*a,flush=True)

STATE_SIGS={
 "cytotoxic":["GZMB","GZMK","GZMA","GZMH","PRF1","NKG7","IFNG"],
 "exhausted":["PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","TOX","ENTPD1","VSIR"],
 "activated":["CD69","IL2RA","ICOS","TNFRSF9","CD28","TNFRSF4","MKI67","HLA-DRA"],
 "Treg":["FOXP3","IL2RA","CTLA4"],
 "Th1":["TBX21","IFNG","CXCR3","STAT1"],
}
ENDPOINTS={  # endpoint -> (program, favorable direction: +1 enhance, -1 reverse/reduce)
 "activation_enhancement":("activated",+1),"Th1_enhancement":("Th1",+1),
 "cytotoxic_enhancement":("cytotoxic",+1),"exhaustion_reversal":("exhausted",-1),
 "Treg_reduction":("Treg",-1)}
def cos(a,b):
    na=np.linalg.norm(a); nb=np.linalg.norm(b)
    return float(a@b/(na*nb)) if na>0 and nb>0 else 0.0

# ---------- load DE_stats (Rest) ----------
log("== load DE_stats ==")
f=h5py.File(DEF,"r")
degenes=np.array([g.decode() if isinstance(g,bytes) else g for g in f["var"]["gene_name"][:]])  # symbols, not Ensembl _index
dgpos={g:i for i,g in enumerate(degenes)}
cc=f["obs"]["culture_condition"]; ccat=[c.decode() for c in cc["categories"][:]]
cond=np.array([ccat[c] for c in cc["codes"][:]])
tg=f["obs"]["target_contrast_gene_name"]; tcat=[c.decode() for c in tg["categories"][:]]
target=np.array([tcat[c] if c>=0 else "NA" for c in tg["codes"][:]])
LFC=f["layers"]["log_fc"]  # lazy
rest_rows=np.where(cond=="Rest")[0]; rest_t=target[rest_rows]
uniq_t=[t for t in pd.unique(rest_t) if t!="NA"]
rowmap={t:rest_rows[np.where(rest_t==t)[0][0]] for t in uniq_t}
log("targets in Rest:",len(uniq_t))

# ---------- load CosMx subset ----------
log("== load CosMx subset ==")
import anndata as ad, scanpy as sc
A=sc.read_h5ad(H5)
A=A[A.obs["sample"].astype(str).isin(SAMPLES)].copy()
ct=A.obs["cell_type"].astype(str)
naive=(ct.to_numpy()=="T CD4 naive")
cg=A.var_names.to_numpy(); cgpos={g:i for i,g in enumerate(cg)}
X=A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(A.X)  # log-norm
# program scores on all cells
for name,gl in STATE_SIGS.items():
    gl=[g for g in gl if g in cgpos]
    sc.tl.score_genes(A,gl,score_name=f"sig_{name}",use_raw=False)
prog_cols=[f"sig_{k}" for k in STATE_SIGS]
shared_gl=set(l.strip() for l in open(f"{BASE}/shared_genes_genelevel.txt"))
markers=set(g for gl in STATE_SIGS.values() for g in gl)
heldout=[g for g in shared_gl if g in dgpos and g in cgpos and g not in markers]
log("held-out genes:",len(heldout),"| shared_genelevel:",len(shared_gl))

# ================= S1 held-out direction-corrected + FDR =================
log("== S1 held-out validation ==")
ho_de=np.array([dgpos[g] for g in heldout]); ho_cx=np.array([cgpos[g] for g in heldout])
# program co-expression signature on held-out genes, from naive CosMx cells (non-circular: markers excluded)
Xn=X[naive]; 
Sp={}
for p in STATE_SIGS:
    s=A.obs.loc[naive,f"sig_{p}"].to_numpy()
    s=s-s.mean()
    col=np.asarray(Xn[:,ho_cx].todense()); col=col-col.mean(0)
    denom=(np.linalg.norm(col,axis=0)*np.linalg.norm(s)); denom[denom==0]=1
    Sp[p]=(col.T@s)/denom  # per-heldout-gene corr with program score
# precompute candidate response matrix on held-out genes (targets x genes)
Dmat=np.vstack([np.nan_to_num(LFC[rowmap[t],:][ho_de]) for t in uniq_t])  # (T, G)
Dn=Dmat/ (np.linalg.norm(Dmat,axis=1,keepdims=True)+1e-12)
NPERM=1000; rng=np.random.default_rng(0)
recs=[]
for ep,(prog,d) in ENDPOINTS.items():
    sig=Sp[prog]; sn=sig/(np.linalg.norm(sig)+1e-12)
    obs=Dn@sn                      # cosine of each target vs signature
    fav=d*obs
    order=fav.argsort(); pct={uniq_t[i]:(r+1)/len(fav) for r,i in enumerate(order)}
    # permutation null: shuffle signature gene labels -> null favorable-cosine distribution (pooled)
    nullvals=[]
    for _ in range(NPERM):
        sp_=sig[rng.permutation(len(sig))]; sp_=sp_/(np.linalg.norm(sp_)+1e-12)
        nullvals.append(d*(Dn@sp_))
    nullpool=np.concatenate(nullvals)                    # (NPERM*T,)
    nullpool.sort()
    # one-sided p per target: fraction of null favorable-cos >= observed
    import numpy as _np
    idx=_np.searchsorted(nullpool,fav,side="left")
    p=(len(nullpool)-idx)/len(nullpool)
    p=_np.clip(p,1.0/len(nullpool),1.0)
    for k,t in enumerate(uniq_t):
        recs.append(dict(endpoint=ep,program=prog,direction=d,target=t,cosine=float(obs[k]),
                         favorable=float(fav[k]),rank_pct=pct[t],p_perm=float(p[k]),n_heldout=len(heldout)))
S1=pd.DataFrame(recs); S1["fdr"]=np.nan
for ep in S1.endpoint.unique():
    m=S1.endpoint==ep; S1.loc[m,"fdr"]=multipletests(S1.loc[m,"p_perm"],method="fdr_bh")[1]
S1=S1.sort_values(["endpoint","rank_pct"],ascending=[True,False])
S1.to_csv(f"{OUT}/S1_heldout_direction_corrected_fdr.csv",index=False)
S1.groupby("endpoint").head(10).to_csv(f"{OUT}/S1_heldout_top10.csv",index=False)
log(S1.groupby("endpoint").head(3)[["endpoint","target","cosine","rank_pct","fdr"]].to_string())

# ================= S2 within-state niche-contrast + rank-shift =================
log("== S2 niche-contrast ==")
asg=pd.read_csv(ASSIGN)
# restrict to naive cells (assignment is already naive-anchored); align by order
asg=asg.reset_index(drop=True)
# stratify naive cells by dominant intrinsic program among 3 tumor-relevant states
strata={"Th1":"sig_Th1","exhausted":"sig_exhausted","cytotoxic":"sig_cytotoxic"}
sp_pos=[cgpos[g] for g in heldout]  # contrast on held-out+shared genes for cleanliness
shared_cx=np.array([cgpos[g] for g in shared_gl if g in cgpos and g in dgpos])
shared_de=np.array([dgpos[g] for g in shared_gl if g in cgpos and g in dgpos])
shared_syms=[g for g in shared_gl if g in cgpos and g in dgpos]
dist=asg["dist_tumor_um"].to_numpy(); samp=asg["sample"].to_numpy()
# map assignment rows to CosMx X rows: assignment has global_cell_index_in_subset
gidx=asg["global_cell_index_in_subset"].to_numpy()
contrast_sigs={}; contrast_meta=[]
for stname,scol in strata.items():
    sv=asg[scol].to_numpy()
    hi=sv>=np.quantile(sv,0.66)  # cells expressing this program
    per_sec=[]
    for s in SAMPLES:
        m=hi&(samp==s)
        if m.sum()<40: continue
        dloc=dist[m]; g_=gidx[m]
        contact=g_[dloc<=np.quantile(dloc,0.33)]; far=g_[dloc>=np.quantile(dloc,0.67)]
        if len(contact)<15 or len(far)<15: continue
        mc=np.asarray(X[contact][:,shared_cx].mean(0)).ravel()
        mf=np.asarray(X[far][:,shared_cx].mean(0)).ravel()
        per_sec.append(mc-mf)
    if not per_sec: continue
    S_r=np.mean(per_sec,axis=0)  # section-averaged contact-vs-far signature
    contrast_sigs[stname]=S_r
    contrast_meta.append(dict(stratum=stname,n_sections=len(per_sec),n_genes=len(shared_syms)))
    pd.DataFrame({"gene":shared_syms,"contact_minus_far":S_r}).sort_values("contact_minus_far").to_csv(
        f"{OUT}/S2_signature_{stname}.csv",index=False)
# align every target to each contrast signature (cosine on shared genes, Rest)
align=[]
for st,S_r in contrast_sigs.items():
    for t in uniq_t:
        dg=np.nan_to_num(LFC[rowmap[t],:][shared_de])
        align.append(dict(stratum=st,target=t,cos=cos(dg,S_r)))
S2=pd.DataFrame(align)
S2p=S2.pivot(index="target",columns="stratum",values="cos")
S2p.to_csv(f"{OUT}/S2_target_x_contrast_cosine.csv")
# rank-shift: rank of each target within each contrast; show top movers
ranks=S2p.rank(ascending=False)
ranks.to_csv(f"{OUT}/S2_target_ranks.csv")
if ranks.shape[1]>=2:
    cols=list(ranks.columns)
    rs=ranks.copy(); rs["rank_spread"]=ranks.max(1)-ranks.min(1)
    rs.sort_values("rank_spread",ascending=False).head(30).to_csv(f"{OUT}/S2_top_rank_shift.csv")
    log("contrast strata:",cols,"| rank spread top:",
        rs.sort_values("rank_spread",ascending=False).head(5).index.tolist())
json.dump(contrast_meta,open(f"{OUT}/S2_contrast_meta.json","w"),indent=1)

# ================= S3 reliability + viability =================
log("== S3 reliability + viability ==")
def oarr(k):
    g=f["obs"][k]
    if isinstance(g,h5py.Group):
        cat=[c.decode() if isinstance(c,bytes) else c for c in g["categories"][:]]
        return np.array([cat[c] if c>=0 else None for c in g["codes"][:]])
    return g[:]
rel_fields=["donor_correlation_hits_mean","donor_correlation_all_mean","guide_correlation_signif",
            "guide_n_signif_ontarget","ontarget_significant","ontarget_effect_size","low_target_gex",
            "n_downstream","n_total_de_genes","n_cells_target","distal_offtarget_flag"]
d={"target":target[rest_rows],"condition":"Rest"}
for k in rel_fields:
    try: d[k]=oarr(k)[rest_rows]
    except Exception as e: log("skip",k,e)
S3=pd.DataFrame(d)
# reliability score: donor+guide consistency, on-target, min trans-effect
for c in ["donor_correlation_hits_mean","ontarget_effect_size","n_downstream","n_total_de_genes","n_cells_target"]:
    S3[c]=pd.to_numeric(S3[c],errors="coerce")
S3["reliable"]=((S3["donor_correlation_hits_mean"]>0.2)&(S3["n_downstream"]>=5)&
                (S3["ontarget_significant"].astype(str).isin(["True","1","1.0"]) if "ontarget_significant" in S3 else True))
# viability / general-disruption penalty: very high global DE (top decile) or very low cell count (depletion)
hi_de=S3["n_total_de_genes"].quantile(0.90)
lo_cells=S3["n_cells_target"].quantile(0.10)
BROAD={"CDK9","ATP2A2","CHD4","STK11","PRKAR1A","POLR2A","RPL","EIF"}
S3["viability_flag"]=((S3["n_total_de_genes"]>=hi_de)|(S3["n_cells_target"]<=lo_cells)|
                      (S3["target"].isin(BROAD)))
S3.to_csv(f"{OUT}/S3_reliability_viability.csv",index=False)
log("reliable targets:",int(S3["reliable"].sum()),"/",len(S3),"| viability-flagged:",int(S3["viability_flag"].sum()))

# ================= S4 spatial quantification =================
log("== S4 spatial stats ==")
def morans_I_nn(vals,nn):
    n=len(vals); k=nn.shape[1]; z=vals-vals.mean(); denom=(z**2).sum()
    if denom==0: return 0.0
    num=(z[:,None]*z[nn]).sum(); W=n*k
    return float((n/W)*(num/denom))
rows=[]
for s in SAMPLES:
    m=(A.obs["sample"].astype(str)==s).to_numpy()&naive
    coords=np.asarray(A.obsm["spatial"])[m]
    if m.sum()<50: continue
    _,nn=cKDTree(coords).query(coords,k=7); nn=nn[:,1:]
    for p in STATE_SIGS:
        v=A.obs.loc[m,f"sig_{p}"].to_numpy()
        obs=morans_I_nn(v,nn)
        null=np.array([morans_I_nn(np.random.permutation(v),nn) for _ in range(199)])
        pval=float((np.sum(null>=obs)+1)/(len(null)+1))
        rows.append(dict(sample=s,program=p,morans_I=obs,perm_p=pval,n_cells=int(m.sum())))
S4=pd.DataFrame(rows); S4.to_csv(f"{OUT}/S4_morans_I_per_section.csv",index=False)
log(S4.to_string())
# niche size filter + section presence
nsz=asg.groupby("fniche").size().rename("n_cells")
sec=asg.groupby(["fniche","sample"]).size().unstack(fill_value=0)
nfilt=pd.concat([nsz,sec],axis=1); nfilt["keep_min50"]=nsz>=50
nfilt.to_csv(f"{OUT}/S4_niche_size_section.csv")
# section rank reproducibility of program burden per niche
prof=pd.read_csv(f"{BASE}/2-Results/cd4_feature_niche/cd4_feature_niche_profile.csv").set_index("fniche")
sec_rankcorr={}
for pa in [("LUAD-5 R1","LUAD-5 R2"),("LUAD-5 R1","LUAD-5 R3"),("LUAD-5 R2","LUAD-5 R3")]:
    # mean program score per niche per section
    def nichemean(s):
        sub=asg[asg["sample"]==s]; return sub.groupby("fniche")[prog_cols].mean()
    a1=nichemean(pa[0]); a2=nichemean(pa[1]); common=a1.index.intersection(a2.index)
    from scipy.stats import spearmanr
    r=spearmanr(a1.loc[common].values.ravel(),a2.loc[common].values.ravel()).correlation
    sec_rankcorr["_vs_".join(pa)]=float(r)
json.dump(sec_rankcorr,open(f"{OUT}/S4_section_rank_reproducibility.json","w"),indent=1)
log("section program-burden rank corr:",sec_rankcorr)

# ================= S5 flip-screen recomputation (Rest vs Stim48) + QC breakdown =================
log("== S5 flip recompute ==")
stim_rows=np.where(cond=="Stim48hr")[0]; stim_t=target[stim_rows]
srowmap={t:stim_rows[np.where(stim_t==t)[0][0]] for t in pd.unique(stim_t) if t!="NA"}
nd=None
try: nd_rest=oarr("n_downstream")[rest_rows]
except Exception: nd_rest=None
# program effect per (target, program) = mean log_fc over program markers present in DE panel
prog_de={p:[dgpos[g] for g in gl if g in dgpos] for p,gl in STATE_SIGS.items()}
common_t=[t for t in uniq_t if t in srowmap]
# n_downstream per condition
ndr={}; nds={}
try:
    a_nd=oarr("n_downstream")
    for t in common_t:
        ndr[t]=float(a_nd[rowmap[t]]); nds[t]=float(a_nd[srowmap[t]])
except Exception as e: log("n_downstream unavailable",e)
frecs=[]
for t in common_t:
    r_vec=np.nan_to_num(LFC[rowmap[t],:]); s_vec=np.nan_to_num(LFC[srowmap[t],:])
    for p,idx in prog_de.items():
        if not idx: continue
        er=float(r_vec[idx].mean()); es=float(s_vec[idx].mean())
        flip=(abs(er)>=0.10 and abs(es)>=0.10 and np.sign(er)!=np.sign(es))
        incomplete=(ndr.get(t,1)==0) or (nds.get(t,1)==0)
        frecs.append(dict(target=t,program=p,effect_rest=er,effect_stim48=es,
                          flip=flip,nd_rest=ndr.get(t,np.nan),nd_stim=nds.get(t,np.nan),
                          incomplete=incomplete))
S5=pd.DataFrame(frecs); S5.to_csv(f"{OUT}/S5_flip_recompute.csv",index=False)
n_flip=int(S5["flip"].sum())
n_flip_incomplete=int((S5["flip"]&S5["incomplete"]).sum())
n_both_zero=int(((S5["nd_rest"]==0)&(S5["nd_stim"]==0)).sum())
n_any_zero=int(((S5["nd_rest"]==0)|(S5["nd_stim"]==0)).sum())
stable=S5[~S5["flip"]]; n_stable=len(stable)
n_stable_incomplete=int(stable["incomplete"].sum())
qc=dict(n_pairs=len(S5),n_flip=n_flip,n_flip_incomplete=n_flip_incomplete,
        n_any_nd0=n_any_zero,n_both_nd0=n_both_zero,
        n_stable=n_stable,n_stable_incomplete=n_stable_incomplete,
        cutoff=0.10,n_targets=len(common_t))
json.dump(qc,open(f"{OUT}/S5_flip_qc.json","w"),indent=1)
log("flip QC:",qc)

log("== DONE ==",OUT)
