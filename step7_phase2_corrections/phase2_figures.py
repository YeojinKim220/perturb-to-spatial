#!/usr/bin/env python3
"""Build report-ready compact summaries + Figure 3 (niche-contrast core) from phase2 outputs."""
import json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
BASE="/storage/home/hcoda1/9/ykim3030/scratch/hackathon"; P=f"{BASE}/2-Results/phase2"

# ---- compact summary json (all values the report needs) ----
s1=pd.read_csv(f"{P}/S1_heldout_direction_corrected_fdr.csv")
summ={"heldout":{}}
for ep,g in s1.groupby("endpoint"):
    g=g.sort_values("rank_pct",ascending=False)
    summ["heldout"][ep]={
        "n_fdr_lt_0.1":int((g.fdr<0.1).sum()),
        "n_fdr_lt_0.25":int((g.fdr<0.25).sum()),
        "min_fdr":float(g.fdr.min()),
        "top5":[{"target":t,"cosine":round(float(c),3),"rank_pct":round(float(r),4),"fdr":round(float(f),3)}
                 for t,c,r,f in zip(g.target.head(5),g.cosine.head(5),g.rank_pct.head(5),g.fdr.head(5))]}
qc=json.load(open(f"{P}/S5_flip_qc.json")); summ["flip_qc"]=qc
mor=pd.read_csv(f"{P}/S4_morans_I_per_section.csv")
summ["morans"]={"per_program_mean":{p:round(float(mor[mor.program==p].morans_I.mean()),4) for p in mor.program.unique()},
                "max_perm_p":round(float(mor.perm_p.max()),4),"n_sig_0.05":int((mor.perm_p<0.05).sum()),"n_tests":len(mor)}
summ["section_reproducibility"]=json.load(open(f"{P}/S4_section_rank_reproducibility.json"))
cm=json.load(open(f"{P}/S2_contrast_meta.json")); summ["niche_contrast_meta"]=cm
s3=pd.read_csv(f"{P}/S3_reliability_viability.csv")
summ["reliability"]={"n_targets":len(s3),"n_reliable":int(s3.reliable.sum()),"n_viability_flag":int(s3.viability_flag.sum())}
# annotate top cytotoxic candidates with reliability/viability
rel=s3.set_index("target")
cyto=summ["heldout"]["cytotoxic_enhancement"]["top5"]
for c in cyto:
    if c["target"] in rel.index:
        c["reliable"]=bool(rel.loc[c["target"],"reliable"]); c["viability_flag"]=bool(rel.loc[c["target"],"viability_flag"])
json.dump(summ,open(f"{P}/phase2_report_summary.json","w"),indent=1)
print(json.dumps(summ,indent=1)[:2500])

# ---- Figure 3: within-state niche-contrast (core) ----
# panel A: contact-vs-far signature top genes per stratum; panel B: target rank-shift across strata
sigs={s:pd.read_csv(f"{P}/S2_signature_{s}.csv") for s in ["Th1","exhausted","cytotoxic"]}
rs=pd.read_csv(f"{P}/S2_top_rank_shift.csv",index_col=0)
fig,axes=plt.subplots(1,4,figsize=(18,5.2))
for ax,(st,df) in zip(axes[:3],sigs.items()):
    df=df.sort_values("contact_minus_far")
    top=pd.concat([df.head(8),df.tail(8)])
    colors=["#c0392b" if v>0 else "#2c7fb8" for v in top.contact_minus_far]
    ax.barh(range(len(top)),top.contact_minus_far,color=colors)
    ax.set_yticks(range(len(top))); ax.set_yticklabels(top.gene,fontsize=7)
    ax.axvline(0,color="k",lw=0.6)
    ax.set_title(f"{st}-high cells\ntumor-contact vs tumor-far",fontsize=10)
    ax.set_xlabel("Δ mean expr (contact − far)",fontsize=8)
# panel D: rank spread of top movers
strata=[c for c in rs.columns if c!="rank_spread"]
mv=rs.sort_values("rank_spread",ascending=False).head(10)
ax=axes[3]
x=np.arange(len(mv))
for i,st in enumerate(strata):
    ax.scatter(mv[st],x,label=st,s=40)
ax.set_yticks(x); ax.set_yticklabels(mv.index,fontsize=7); ax.invert_yaxis()
ax.set_xlabel("perturbation rank (1=best aligned)",fontsize=8)
ax.set_title("Perturbation rank-shift\nacross within-state niches",fontsize=10)
ax.legend(fontsize=7)
fig.suptitle("Figure 3 — Within-state niche contrast (tumor-contact vs tumor-far) and perturbation rank-shift  ·  exploratory projection — not experimentally validated",fontsize=11)
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig(f"{P}/fig3_niche_contrast.png",dpi=150)
print("SAVED fig3_niche_contrast.png")

# ---- Figure: held-out direction-corrected + FDR ----
fig2,ax=plt.subplots(figsize=(9,5.5))
order=["cytotoxic_enhancement","activation_enhancement","Th1_enhancement","exhaustion_reversal","Treg_reduction"]
for i,ep in enumerate(order):
    g=s1[s1.endpoint==ep].sort_values("rank_pct",ascending=False).head(15)
    y=[i]*len(g)
    sig=g.fdr<0.1
    ax.scatter(g.rank_pct[~sig],[i]*(~sig).sum(),c="#bbb",s=25)
    ax.scatter(g.rank_pct[sig],[i]*sig.sum(),c="#c0392b",s=45,zorder=3)
    # label top target
    tt=g.iloc[0]; ax.annotate(f"{tt.target} (fdr={tt.fdr:.2f})",(tt.rank_pct,i),fontsize=7,xytext=(3,3),textcoords="offset points")
ax.set_yticks(range(len(order))); ax.set_yticklabels(order)
ax.set_xlabel("genome-wide favorable rank percentile (held-out genes)")
ax.set_title("Held-out validation, direction-corrected + permutation FDR\nred = BH-FDR < 0.1  ·  only cytotoxic-enhancement survives correction",fontsize=10)
ax.set_xlim(0.9,1.001)
fig2.tight_layout(); fig2.savefig(f"{P}/fig_heldout_corrected.png",dpi=150)
print("SAVED fig_heldout_corrected.png")
