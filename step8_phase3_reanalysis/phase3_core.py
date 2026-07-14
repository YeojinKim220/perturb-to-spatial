#!/usr/bin/env python3
"""Phase-3 core re-analyses per methodology review.

Outputs are written under hackathon/2-Results/phase3.

Analyses:
  F1   independent intrinsic-state x environment-only-context association
  MI   Moran's I re-check with max-distance, FOV-safe graph, n_counts residuals, BH
  S2c  stress-corrected contact-vs-far signatures and target rankings
  REP  leave-one-section-out ranking reproducibility
  FUN  FDR/reliability/viability/spatial evidence funnel
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


np.random.seed(0)
rng = np.random.default_rng(0)

BASE = Path(__file__).resolve().parents[1]
H5 = BASE / "0-Data/ST_CosMx_NSCLC/cosmx_human_nsclc_clustered.h5ad"
DEF = BASE / "0-Data/perturbseq_marson2025/GWCD4i.DE_stats.h5ad"
PHASE2 = BASE / "2-Results/phase2"
OUT = BASE / "2-Results/phase3"
OUT.mkdir(parents=True, exist_ok=True)

SAMPLES = ["LUAD-5 R1", "LUAD-5 R2", "LUAD-5 R3"]
STATE_SIGS = {
    "cytotoxic": ["GZMB", "GZMK", "GZMA", "GZMH", "PRF1", "NKG7", "IFNG"],
    "exhausted": ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "ENTPD1", "VSIR"],
    "activated": ["CD69", "IL2RA", "ICOS", "TNFRSF9", "CD28", "TNFRSF4", "MKI67", "HLA-DRA"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4"],
    "Th1": ["TBX21", "IFNG", "CXCR3", "STAT1"],
}
STATES_FOR_CONTRAST = ["Th1", "exhausted", "cytotoxic"]
TUMOR_TYPES = {"tumor 5", "tumor 6", "tumor 9", "tumor 12", "tumor 13", "epithelial"}
NBR_TYPES = {
    "macrophage": ["macrophage", "monocyte", "mDC", "pDC", "neutrophil"],
    "fibroblast": ["fibroblast"],
    "Bcell": ["B-cell", "plasmablast"],
    "endothelial": ["endothelial"],
    "tumor": list(TUMOR_TYPES),
}
RADIUS_UM = 50.0
CONTACT_UM = 30.0
BOOT_N = 200

TECH_GENES = {
    "HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1", "HSPA8", "HSPB1", "DNAJB1",
    "FOS", "JUN", "JUNB", "EGR1", "HSPA6", "BAG3", "DNAJA1", "VIM", "FN1",
    "COL1A1", "COL1A2", "COL3A1", "S100A4", "S100A6", "ACTB", "B2M", "MALAT1",
}


def log(*parts: object) -> None:
    print(*parts, flush=True)


def is_tech_gene(gene: str) -> bool:
    if gene in TECH_GENES:
        return True
    return gene.startswith(("RPL", "RPS", "MT-", "MRPL", "MRPS"))


def normalize_vec(x: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(x)
    return x / (n + 1e-12)


def ranks_desc(values: np.ndarray) -> np.ndarray:
    order = np.argsort(-values)
    ranks = np.empty(len(values), dtype=int)
    ranks[order] = np.arange(1, len(values) + 1)
    return ranks


def top_tertile_mask(values: np.ndarray) -> np.ndarray:
    return values >= np.nanquantile(values, 2.0 / 3.0)


def positive_top_mask(values: np.ndarray, q: float = 2.0 / 3.0) -> tuple[np.ndarray, float, int]:
    positive = np.isfinite(values) & (values > 0)
    n_positive = int(positive.sum())
    if n_positive == 0:
        return np.zeros(len(values), dtype=bool), np.nan, 0
    threshold = float(np.nanquantile(values[positive], q)) if n_positive >= 10 else float(values[positive].min())
    return positive & (values >= threshold), threshold, n_positive


def context_mask(frame: pd.DataFrame, context: str) -> tuple[np.ndarray, float, int]:
    if context == "tumor_contact":
        values = frame["dist_tumor"].to_numpy(float)
        mask = np.isfinite(values) & (values <= CONTACT_UM)
        return mask, CONTACT_UM, int(np.isfinite(values).sum())
    if context == "tumor_far":
        values = frame["dist_tumor"].to_numpy(float)
        threshold = float(np.nanquantile(values, 2.0 / 3.0))
        return values >= threshold, threshold, int(np.isfinite(values).sum())
    col = {
        "myeloid_rich": "nbr_macrophage",
        "fibroblast_rich": "nbr_fibroblast",
        "Bcell_rich": "nbr_Bcell",
    }[context]
    return positive_top_mask(frame[col].to_numpy(float))


def load_cosmx() -> tuple[sc.AnnData, pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    log("== load CosMx ==")
    adata = sc.read_h5ad(H5)
    adata = adata[adata.obs["sample"].astype(str).isin(SAMPLES)].copy()

    genes = adata.var_names.to_numpy()
    gene_pos = {g: i for i, g in enumerate(genes)}
    for name, gene_list in STATE_SIGS.items():
        present = [g for g in gene_list if g in gene_pos]
        sc.tl.score_genes(adata, present, score_name=f"sig_{name}", use_raw=False)

    obs = adata.obs
    xy = np.asarray(adata.obsm["spatial"])
    cell_type = obs["cell_type"].astype(str).to_numpy()
    sample = obs["sample"].astype(str).to_numpy()
    fov = obs["fov"].astype(str).to_numpy()
    n_counts = obs["n_counts"].to_numpy().astype(float)
    naive = cell_type == "T CD4 naive"
    log("naive-annotated CD4 T cells:", int(naive.sum()), "total cells:", adata.n_obs)
    return adata, obs, xy, cell_type, sample, fov, n_counts


def build_context(
    adata: sc.AnnData,
    obs: pd.DataFrame,
    xy: np.ndarray,
    cell_type: np.ndarray,
    sample: np.ndarray,
    fov: np.ndarray,
    n_counts: np.ndarray,
) -> pd.DataFrame:
    log("== build environment-only context ==")
    dist_tumor = np.full(adata.n_obs, np.nan)
    local_density = np.zeros(adata.n_obs, dtype=float)
    comp = {k: np.zeros(adata.n_obs, dtype=float) for k in NBR_TYPES}

    for section in SAMPLES:
        section_idx = np.where(sample == section)[0]
        coords = xy[section_idx]
        tree = cKDTree(coords)

        tumor_local = np.isin(cell_type[section_idx], list(TUMOR_TYPES))
        if tumor_local.sum() > 0:
            tumor_tree = cKDTree(coords[tumor_local])
            dist_tumor[section_idx], _ = tumor_tree.query(coords, k=1)

        for local_i, global_i in enumerate(section_idx):
            nbr_local = tree.query_ball_point(coords[local_i], RADIUS_UM)
            nbr = [section_idx[j] for j in nbr_local if section_idx[j] != global_i]
            local_density[global_i] = len(nbr)
            if not nbr:
                continue
            nbr_types = cell_type[nbr]
            for key, members in NBR_TYPES.items():
                comp[key][global_i] = np.isin(nbr_types, members).mean()
        log(f"  {section}: cells={len(section_idx):,}")

    ctx = pd.DataFrame(index=np.arange(adata.n_obs))
    for key in NBR_TYPES:
        ctx[f"nbr_{key}"] = comp[key]
    ctx["dist_tumor"] = dist_tumor
    ctx["log_dist_tumor"] = np.log1p(dist_tumor)
    ctx["contact"] = (dist_tumor <= CONTACT_UM).astype(int)
    ctx["local_density"] = local_density
    for state in STATE_SIGS:
        ctx[f"sig_{state}"] = obs[f"sig_{state}"].to_numpy()
    ctx["sample"] = sample
    ctx["naive"] = cell_type == "T CD4 naive"
    ctx["n_counts"] = n_counts
    ctx["fov"] = fov

    ctx_naive = ctx[ctx.naive].copy()
    ctx_naive.to_csv(OUT / "naive_state_context.csv.gz", index=False, compression="gzip")
    log("context rows:", len(ctx_naive))
    return ctx_naive


def run_f1(ctx_naive: pd.DataFrame) -> pd.DataFrame:
    log("== F1 state x context association ==")
    contexts = ["tumor_contact", "tumor_far", "myeloid_rich", "fibroblast_rich", "Bcell_rich"]
    qc_rows = []
    rows = []

    for context in contexts:
        for section in SAMPLES:
            d = ctx_naive[ctx_naive["sample"] == section]
            inc, threshold, n_positive = context_mask(d, context)
            qc_rows.append({
                "context": context,
                "sample": section,
                "n_total": int(len(d)),
                "n_context": int(inc.sum()),
                "fraction_context": float(inc.mean()) if len(d) else np.nan,
                "threshold": threshold,
                "n_positive_or_finite": int(n_positive),
            })
    pd.DataFrame(qc_rows).to_csv(OUT / "F1_context_qc.csv", index=False)

    for state in STATES_FOR_CONTRAST:
        for context in contexts:
            effects = []
            counts = {}
            for section in SAMPLES:
                d = ctx_naive[ctx_naive["sample"] == section]
                high = top_tertile_mask(d[f"sig_{state}"].to_numpy(float))
                inc, _, _ = context_mask(d, context)

                a = int((high & inc).sum())
                b = int((high & ~inc).sum())
                c = int((~high & inc).sum())
                dd = int((~high & ~inc).sum())
                odds_ratio = ((a + 0.5) * (dd + 0.5)) / ((b + 0.5) * (c + 0.5))
                effects.append(float(np.log(odds_ratio)))
                tag = section.split()[-1]
                counts[f"n_context_{tag}"] = int(inc.sum())
                counts[f"n_state_high_in_context_{tag}"] = a

            signs = np.sign(effects)
            sign_consistent = int((signs[0] == signs[1] == signs[2]) and signs[0] != 0)
            rows.append({
                "state": state,
                "context": context,
                "logOR_R1": effects[0],
                "logOR_R2": effects[1],
                "logOR_R3": effects[2],
                "logOR_mean": float(np.mean(effects)),
                "sign_consistent": sign_consistent,
                **counts,
            })

    f1 = pd.DataFrame(rows)
    f1.to_csv(OUT / "F1_state_context_assoc.csv", index=False)
    log("F1 associations:", len(f1), "sign-consistent:", int(f1.sign_consistent.sum()))
    return f1


def morans_i(vals: np.ndarray, coords: np.ndarray, fov: np.ndarray, maxd: float = 40.0, k: int = 6) -> float:
    tree = cKDTree(coords)
    dist, idx = tree.query(coords, k=k + 1)
    z = vals - vals.mean()
    den = float((z ** 2).sum())
    if den == 0:
        return np.nan
    num = 0.0
    weight = 0
    for i in range(len(vals)):
        for jj in range(1, k + 1):
            j = idx[i, jj]
            if dist[i, jj] > maxd:
                continue
            if fov[i] != fov[j]:
                continue
            num += z[i] * z[j]
            weight += 1
    if weight == 0:
        return np.nan
    return float((len(vals) / weight) * (num / den))


def run_moran(ctx_naive: pd.DataFrame, xy: np.ndarray, fov: np.ndarray) -> pd.DataFrame:
    log("== MI Moran re-check ==")
    rows = []
    for section in SAMPLES:
        d = ctx_naive[ctx_naive["sample"] == section]
        global_idx = d.index.to_numpy()
        coords = xy[global_idx]
        fov_section = fov[global_idx]
        n_counts = d["n_counts"].to_numpy(float)
        design = np.c_[np.ones_like(n_counts), n_counts]

        for state in STATE_SIGS:
            values = d[f"sig_{state}"].to_numpy(float)
            beta = np.linalg.lstsq(design, values, rcond=None)[0]
            resid = values - design @ beta
            observed = morans_i(resid, coords, fov_section)
            if np.isnan(observed):
                pval = 1.0
            else:
                null = np.array([
                    morans_i(resid[rng.permutation(len(resid))], coords, fov_section)
                    for _ in range(199)
                ])
                pval = float((1 + np.sum(null >= observed)) / 200.0)
            rows.append({
                "sample": section,
                "program": state,
                "morans_I": float(observed),
                "perm_p": pval,
                "n": int(len(d)),
            })

    mi = pd.DataFrame(rows)
    mi["fdr_bh"] = multipletests(mi["perm_p"].fillna(1.0), method="fdr_bh")[1]
    mi.to_csv(OUT / "MI_morans_recheck.csv", index=False)
    log("MI BH<0.05:", int((mi.fdr_bh < 0.05).sum()), "/", len(mi))
    return mi


def read_de() -> tuple[h5py.File, np.ndarray, list[str], dict[str, int], list[str], np.ndarray, np.ndarray, np.ndarray]:
    log("== load Perturb-seq DE_stats ==")
    h5 = h5py.File(DEF, "r")
    de_genes = np.array([g.decode() if isinstance(g, bytes) else g for g in h5["var"]["gene_name"][:]])
    de_pos = {g: i for i, g in enumerate(de_genes)}

    cond_group = h5["obs"]["culture_condition"]
    cond_cat = [x.decode() for x in cond_group["categories"][:]]
    cond = np.array([cond_cat[x] for x in cond_group["codes"][:]])

    target_group = h5["obs"]["target_contrast_gene_name"]
    target_cat = [x.decode() for x in target_group["categories"][:]]
    target = np.array([target_cat[x] if x >= 0 else "NA" for x in target_group["codes"][:]])

    rest_rows = np.where(cond == "Rest")[0]
    rest_target = target[rest_rows]
    targets = [t for t in pd.unique(rest_target) if t != "NA"]
    rowmap = {t: rest_rows[np.where(rest_target == t)[0][0]] for t in targets}
    return h5, de_genes, targets, rowmap, cond.tolist(), target, rest_rows, h5["layers"]["log_fc"]


def section_contrast_vector(
    state: str,
    section: str,
    ctx_naive: pd.DataFrame,
    xlog: sp.csr_matrix,
    shared_cx: np.ndarray,
) -> tuple[np.ndarray | None, dict[str, int]]:
    d = ctx_naive[ctx_naive["sample"] == section]
    high = d[f"sig_{state}"] >= np.nanquantile(d[f"sig_{state}"], 2.0 / 3.0)
    dd = d[high].copy()
    if len(dd) < 60:
        return None, {"n_state_high": int(len(dd))}

    dist = dd["dist_tumor"].to_numpy(float)
    contact = dd.index[dist <= np.nanquantile(dist, 1.0 / 3.0)].to_numpy()
    far = dd.index[dist >= np.nanquantile(dist, 2.0 / 3.0)].to_numpy()
    if len(contact) < 15 or len(far) < 15:
        return None, {"n_state_high": int(len(dd)), "n_contact": int(len(contact)), "n_far": int(len(far))}

    cov_cols = ["n_counts", f"sig_{state}", "local_density", "nbr_tumor"]
    cov = ctx_naive.loc[np.r_[contact, far], cov_cols].to_numpy(float)
    cov[:, 0] = np.log1p(cov[:, 0])
    mu = np.nanmean(cov, axis=0)
    sd = np.nanstd(cov, axis=0)
    sd[sd == 0] = 1.0
    cov = (cov - mu) / sd
    cov_contact = cov[: len(contact)]
    cov_far = cov[len(contact):]

    tree = cKDTree(cov_far)
    _, nearest = tree.query(cov_contact, k=1)
    matched_far = far[nearest]

    contact_mean = np.asarray(xlog[contact][:, shared_cx].mean(0)).ravel()
    far_mean = np.asarray(xlog[matched_far][:, shared_cx].mean(0)).ravel()
    stats = {
        "n_state_high": int(len(dd)),
        "n_contact": int(len(contact)),
        "n_far": int(len(far)),
        "n_matched_far_unique": int(len(np.unique(matched_far))),
    }
    return contact_mean - far_mean, stats


def signature_from_section_matrix(lfc_by_section: np.ndarray, tech_mask: np.ndarray, apply_filter: bool) -> np.ndarray:
    mean = lfc_by_section.mean(axis=0)
    if lfc_by_section.shape[0] >= 2:
        score = mean / (lfc_by_section.std(axis=0) + 1e-9)
        agree_needed = 2 if lfc_by_section.shape[0] == 3 else lfc_by_section.shape[0]
        agree = (np.sign(lfc_by_section) == np.sign(mean)).sum(axis=0)
        if apply_filter:
            score[agree < agree_needed] = 0.0
    else:
        score = mean.copy()
    score[tech_mask] = 0.0
    return normalize_vec(score)


def run_s2c(
    adata: sc.AnnData,
    ctx_naive: pd.DataFrame,
    de_pos: dict[str, int],
    targets: list[str],
    rowmap: dict[str, int],
    lfc_layer,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, np.ndarray, list[str], np.ndarray]:
    log("== S2c corrected contact-vs-far signatures ==")
    genes = adata.var_names.to_numpy()
    cosmx_pos = {g: i for i, g in enumerate(genes)}
    shared = [g.strip() for g in open(BASE / "shared_genes_genelevel.txt")]
    shared = [g for g in shared if g in de_pos and g in cosmx_pos]
    shared_de = np.array([de_pos[g] for g in shared], dtype=int)
    shared_cx = np.array([cosmx_pos[g] for g in shared], dtype=int)
    tech_mask = np.array([is_tech_gene(g) for g in shared], dtype=bool)

    xlog = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    dmat = np.vstack([np.nan_to_num(lfc_layer[rowmap[t], :][shared_de]) for t in targets])
    dnorm = dmat / (np.linalg.norm(dmat, axis=1, keepdims=True) + 1e-12)

    section_mats = {}
    sig_raw = {}
    sig_corr = {}
    meta_rows = []
    match_rows = []

    for state in STATES_FOR_CONTRAST:
        vectors = []
        for section in SAMPLES:
            vec, stats = section_contrast_vector(state, section, ctx_naive, xlog, shared_cx)
            stats.update({"state": state, "sample": section})
            match_rows.append(stats)
            if vec is not None:
                vectors.append(vec)
        if len(vectors) != 3:
            log("  skip", state, "sections with valid contrast:", len(vectors))
            continue

        mat = np.vstack(vectors)
        section_mats[state] = mat
        mean = mat.mean(axis=0)
        z = mean / (mat.std(axis=0) + 1e-9)
        agree = (np.sign(mat) == np.sign(mean)).sum(axis=0)
        keep = agree >= 2

        sig = pd.DataFrame({
            "gene": shared,
            "mean_lfc": mean,
            "z": z,
            "agree": agree,
            "keep": keep,
            "tech": tech_mask,
            "state": state,
        }).sort_values("z", ascending=False)
        sig.to_csv(OUT / f"S2c_signature_{state}.csv", index=False)

        raw = normalize_vec(z.copy())
        corr_score = z.copy()
        corr_score[tech_mask] = 0.0
        corr_score[~keep] = 0.0
        corrected = normalize_vec(corr_score)
        sig_raw[state] = raw
        sig_corr[state] = corrected

        meta_rows.append({
            "state": state,
            "n_keep": int(keep.sum()),
            "n_tech_in_top20": int(sig.head(20)["tech"].sum()),
            "top_genes": ",".join(sig.head(8)["gene"]),
        })

    meta = pd.DataFrame(meta_rows)
    meta.to_csv(OUT / "S2c_meta.csv", index=False)
    pd.DataFrame(match_rows).to_csv(OUT / "S2c_matching_qc.csv", index=False)

    target_rows = []
    target_array = np.array(targets)
    for state in STATES_FOR_CONTRAST:
        if state not in sig_corr:
            continue
        raw_align = dnorm @ sig_raw[state]
        corr_align = dnorm @ sig_corr[state]
        raw_rank = ranks_desc(raw_align)
        corr_rank = ranks_desc(corr_align)

        section_top20_count = np.zeros(len(targets), dtype=int)
        for i, section in enumerate(SAMPLES):
            single_sig = signature_from_section_matrix(section_mats[state][i : i + 1], tech_mask, apply_filter=False)
            single_rank = ranks_desc(dnorm @ single_sig)
            section_top20_count += single_rank <= 20

        for i, target in enumerate(targets):
            target_rows.append({
                "state": state,
                "target": target,
                "raw_alignment": float(raw_align[i]),
                "corrected_alignment": float(corr_align[i]),
                "raw_rank": int(raw_rank[i]),
                "corrected_rank": int(corr_rank[i]),
                "rank_shift_raw_minus_corrected": int(raw_rank[i] - corr_rank[i]),
                "section_top20_recurrence": int(section_top20_count[i]),
                "stress_correction_stable": bool((corr_rank[i] <= 20) and (raw_rank[i] <= 50)),
                "section_consistent": bool(section_top20_count[i] >= 2),
            })

    target_x = pd.DataFrame(target_rows)
    target_x.to_csv(OUT / "S2c_target_x_corrected_context.csv", index=False)

    sens_rows = []
    for state in STATES_FOR_CONTRAST:
        if state not in sig_corr:
            continue
        state_rows = target_x[target_x["state"] == state]
        top_raw = set(state_rows.nsmallest(20, "raw_rank")["target"])
        top_corr = set(state_rows.nsmallest(20, "corrected_rank")["target"])
        sens_rows.append({
            "state": state,
            "spearman_raw_vs_corr": float(spearmanr(state_rows["raw_alignment"], state_rows["corrected_alignment"]).correlation),
            "top20_overlap": int(len(top_raw & top_corr)),
            "top_corr": ",".join(state_rows.nsmallest(10, "corrected_rank")["target"]),
        })
    sensitivity = pd.DataFrame(sens_rows)
    sensitivity.to_csv(OUT / "S2c_sensitivity.csv", index=False)

    log("S2c states:", meta.to_dict("records"))
    return target_x, meta, sensitivity, dnorm, targets, tech_mask


def run_rep(
    target_x: pd.DataFrame,
    dnorm: np.ndarray,
    targets: list[str],
    tech_mask: np.ndarray,
    ctx_naive: pd.DataFrame,
    adata: sc.AnnData,
) -> pd.DataFrame:
    log("== REP leave-one-section-out ==")
    genes = adata.var_names.to_numpy()
    cosmx_pos = {g: i for i, g in enumerate(genes)}
    shared = [g.strip() for g in open(BASE / "shared_genes_genelevel.txt")]
    shared = [g for g in shared if g in cosmx_pos]
    shared_cx = np.array([cosmx_pos[g] for g in shared], dtype=int)
    xlog = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    target_array = np.array(targets)

    def sig_for_sections(state: str, sections: list[str], apply_filter: bool) -> np.ndarray | None:
        vecs = []
        for section in sections:
            vec, _ = section_contrast_vector(state, section, ctx_naive, xlog, shared_cx)
            if vec is None:
                return None
            vecs.append(vec)
        return signature_from_section_matrix(np.vstack(vecs), tech_mask, apply_filter=apply_filter)

    rows = []
    top20_records = []
    for state in STATES_FOR_CONTRAST:
        for held_out in SAMPLES:
            train_sections = [s for s in SAMPLES if s != held_out]
            train_sig = sig_for_sections(state, train_sections, apply_filter=True)
            test_sig = sig_for_sections(state, [held_out], apply_filter=False)
            if train_sig is None or test_sig is None:
                continue
            train_align = dnorm @ train_sig
            test_align = dnorm @ test_sig
            rho = spearmanr(train_align, test_align).correlation
            train_top20 = set(target_array[np.argsort(-train_align)[:20]])
            test_top20 = set(target_array[np.argsort(-test_align)[:20]])
            rows.append({
                "state": state,
                "held_out": held_out,
                "spearman": float(rho),
                "top20_recurrence": int(len(train_top20 & test_top20)),
            })
            for target in sorted(train_top20 | test_top20):
                top20_records.append({
                    "state": state,
                    "held_out": held_out,
                    "target": target,
                    "in_train_top20": target in train_top20,
                    "in_test_top20": target in test_top20,
                })

    rep = pd.DataFrame(rows)
    rep.to_csv(OUT / "REP_leave_one_section.csv", index=False)
    pd.DataFrame(top20_records).to_csv(OUT / "REP_top20_targets_by_split.csv", index=False)
    log("REP mean Spearman:", float(rep.spearman.mean()) if len(rep) else np.nan)
    return rep


def run_bootstrap(
    target_x: pd.DataFrame,
    dnorm: np.ndarray,
    targets: list[str],
    tech_mask: np.ndarray,
    adata: sc.AnnData,
    ctx_naive: pd.DataFrame,
) -> pd.DataFrame:
    log("== bootstrap corrected ranks ==")
    genes = adata.var_names.to_numpy()
    cosmx_pos = {g: i for i, g in enumerate(genes)}
    shared = [g.strip() for g in open(BASE / "shared_genes_genelevel.txt")]
    shared = [g for g in shared if g in cosmx_pos]
    shared_cx = np.array([cosmx_pos[g] for g in shared], dtype=int)
    xlog = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    target_array = np.array(targets)
    target_to_idx = {t: i for i, t in enumerate(targets)}

    audit = pd.read_csv(PHASE2 / "S1_fdr_audit.csv")
    rel = pd.read_csv(PHASE2 / "S3_reliability_viability.csv")
    cyto = audit[audit.endpoint == "cytotoxic_enhancement"].merge(rel[["target", "reliable", "viability_flag"]], on="target", how="left")
    evidence_targets = set(cyto[(cyto.fdr_pertarget < 0.1) & (cyto.reliable) & (~cyto.viability_flag)]["target"])

    rows = []
    for state in STATES_FOR_CONTRAST:
        state_rows = target_x[target_x["state"] == state]
        focus = set(state_rows.nsmallest(200, "corrected_rank")["target"]) | evidence_targets
        focus = [t for t in focus if t in target_to_idx]
        if not focus:
            continue

        section_vecs = []
        for section in SAMPLES:
            vec, _ = section_contrast_vector(state, section, ctx_naive, xlog, shared_cx)
            if vec is not None:
                section_vecs.append(vec)
        if len(section_vecs) != 3:
            continue
        mat = np.vstack(section_vecs)

        rank_samples = {t: [] for t in focus}
        for _ in range(BOOT_N):
            boot_idx = rng.integers(0, mat.shape[0], size=mat.shape[0])
            sig = signature_from_section_matrix(mat[boot_idx], tech_mask, apply_filter=True)
            ranks = ranks_desc(dnorm @ sig)
            for target in focus:
                rank_samples[target].append(int(ranks[target_to_idx[target]]))

        for target, vals in rank_samples.items():
            arr = np.array(vals)
            rows.append({
                "state": state,
                "target": target,
                "rank_median": float(np.median(arr)),
                "rank_q10": float(np.quantile(arr, 0.10)),
                "rank_q90": float(np.quantile(arr, 0.90)),
                "n_bootstrap": BOOT_N,
            })

    boot = pd.DataFrame(rows)
    boot.to_csv(OUT / "S2c_bootstrap_rank_ci.csv", index=False)
    log("bootstrap rows:", len(boot))
    return boot


def run_fun(target_x: pd.DataFrame, boot: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    log("== FUN evidence funnel ==")
    old = pd.read_csv(PHASE2 / "S1_heldout_direction_corrected_fdr.csv")
    audit = pd.read_csv(PHASE2 / "S1_fdr_audit.csv")
    rel = pd.read_csv(PHASE2 / "S3_reliability_viability.csv")

    old_cyto = old[old.endpoint == "cytotoxic_enhancement"][["target", "cosine", "p_perm", "fdr"]].rename(
        columns={"cosine": "old_pooled_cosine", "p_perm": "old_pooled_p", "fdr": "old_pooled_fdr"}
    )
    audit_cyto = audit[audit.endpoint == "cytotoxic_enhancement"].copy()
    rel = rel[["target", "reliable", "viability_flag", "n_downstream", "n_total_de_genes", "n_cells_target"]]

    spatial = target_x[target_x.state == "cytotoxic"].copy()
    spatial = spatial[[
        "target", "raw_rank", "corrected_rank", "section_top20_recurrence",
        "stress_correction_stable", "section_consistent",
    ]]
    boot_cyto = boot[boot.state == "cytotoxic"][["target", "rank_median", "rank_q10", "rank_q90"]]

    fun = audit_cyto.merge(old_cyto, on="target", how="left")
    fun = fun.merge(rel, on="target", how="left")
    fun = fun.merge(spatial, on="target", how="left")
    fun = fun.merge(boot_cyto, on="target", how="left")
    fun["old_pooled_hit"] = fun["old_pooled_fdr"] < 0.1
    fun["per_target_hit"] = fun["fdr_pertarget"] < 0.1
    fun["donor_reliable"] = fun["reliable"].fillna(False).astype(bool)
    fun["no_viability_flag"] = ~fun["viability_flag"].fillna(True).astype(bool)
    fun["bootstrap_stable"] = fun["rank_q90"].fillna(np.inf) <= 100
    fun["clean_candidate"] = (
        fun["per_target_hit"]
        & fun["donor_reliable"]
        & fun["no_viability_flag"]
        & fun["stress_correction_stable"].fillna(False).astype(bool)
        & fun["section_consistent"].fillna(False).astype(bool)
        & fun["bootstrap_stable"]
    )

    audit_out = fun[fun["old_pooled_hit"] | fun["per_target_hit"]].sort_values(
        ["old_pooled_hit", "per_target_hit", "fdr_pertarget"], ascending=[False, False, True]
    )
    audit_out.to_csv(OUT / "FUN_fdr_hit_audit.csv", index=False)
    clean = fun[fun["clean_candidate"]].sort_values("fdr_pertarget")
    clean.to_csv(OUT / "FUN_clean_shortlist.csv", index=False)

    steps = [
        ("old pooled-null cytotoxic FDR < 0.1", fun["old_pooled_hit"]),
        ("old pooled-null hit + donor reliable + no viability flag", fun["old_pooled_hit"] & fun["donor_reliable"] & fun["no_viability_flag"]),
        ("per-target cytotoxic FDR < 0.1", fun["per_target_hit"]),
        ("per-target hit + donor reliable", fun["per_target_hit"] & fun["donor_reliable"]),
        ("per-target hit + donor reliable + no viability flag", fun["per_target_hit"] & fun["donor_reliable"] & fun["no_viability_flag"]),
        ("+ stable after stress correction", fun["per_target_hit"] & fun["donor_reliable"] & fun["no_viability_flag"] & fun["stress_correction_stable"].fillna(False).astype(bool)),
        ("+ section-consistent corrected rank", fun["per_target_hit"] & fun["donor_reliable"] & fun["no_viability_flag"] & fun["stress_correction_stable"].fillna(False).astype(bool) & fun["section_consistent"].fillna(False).astype(bool)),
        ("+ bootstrap rank q90 <= 100", fun["clean_candidate"]),
    ]
    funnel = pd.DataFrame({
        "step": [s for s, _ in steps],
        "n_targets": [int(mask.sum()) for _, mask in steps],
    })
    funnel.to_csv(OUT / "FUN_evidence_funnel.csv", index=False)

    summary = {
        "old_pooled_cytotoxic_fdr_lt_0_1": int(fun["old_pooled_hit"].sum()),
        "per_target_cytotoxic_fdr_lt_0_1": int(fun["per_target_hit"].sum()),
        "per_target_reliable_no_viability": int((fun["per_target_hit"] & fun["donor_reliable"] & fun["no_viability_flag"]).sum()),
        "old_pooled_reliable_no_viability": int((fun["old_pooled_hit"] & fun["donor_reliable"] & fun["no_viability_flag"]).sum()),
        "clean_shortlist": int(fun["clean_candidate"].sum()),
        "clean_shortlist_targets": clean["target"].tolist(),
    }
    with open(OUT / "phase3_core_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    log("funnel:", funnel.to_dict("records"))
    return fun, funnel, clean


def main() -> None:
    adata, obs, xy, cell_type, sample, fov, n_counts = load_cosmx()
    ctx_naive = build_context(adata, obs, xy, cell_type, sample, fov, n_counts)
    run_f1(ctx_naive)
    run_moran(ctx_naive, xy, fov)

    h5, _, targets, rowmap, _, _, _, lfc_layer = read_de()
    de_genes = np.array([g.decode() if isinstance(g, bytes) else g for g in h5["var"]["gene_name"][:]])
    de_pos = {g: i for i, g in enumerate(de_genes)}

    target_x, _, _, dnorm, targets, tech_mask = run_s2c(adata, ctx_naive, de_pos, targets, rowmap, lfc_layer)
    run_rep(target_x, dnorm, targets, tech_mask, ctx_naive, adata)
    boot = run_bootstrap(target_x, dnorm, targets, tech_mask, adata, ctx_naive)
    run_fun(target_x, boot)
    log("== ALL PHASE3 CORE STAGES COMPLETE ==")


if __name__ == "__main__":
    main()
