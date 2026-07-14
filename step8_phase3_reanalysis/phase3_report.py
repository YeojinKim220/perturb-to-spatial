#!/usr/bin/env python3
"""Generate phase-3 figures and Korean/English Markdown reports."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm


BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "2-Results/phase3"
DOC = BASE / "1-Doc"
DOC.mkdir(parents=True, exist_ok=True)


def read_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(OUT / name)


def save_fig1() -> None:
    f1 = read_csv("F1_state_context_assoc.csv")
    mat = f1.pivot(index="state", columns="context", values="logOR_mean")
    cons = f1.pivot(index="state", columns="context", values="sign_consistent")
    contexts = ["tumor_contact", "tumor_far", "myeloid_rich", "fibroblast_rich", "Bcell_rich"]
    states = ["Th1", "exhausted", "cytotoxic"]
    mat = mat.loc[states, contexts]
    cons = cons.loc[states, contexts]

    fig, ax = plt.subplots(figsize=(9.5, 4.5))
    vmax = max(0.1, np.nanmax(np.abs(mat.values)))
    im = ax.imshow(mat.values, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax))
    ax.set_xticks(range(len(contexts)))
    ax.set_xticklabels([c.replace("_", "\n") for c in contexts], fontsize=9)
    ax.set_yticks(range(len(states)))
    ax.set_yticklabels(states, fontsize=10)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            label = f"{mat.values[i, j]:.2f}"
            if cons.values[i, j] == 1:
                label += "*"
            ax.text(j, i, label, ha="center", va="center", fontsize=9, color="black")
    ax.set_title("Figure 1 - Independent intrinsic state x environment context")
    fig.colorbar(im, ax=ax, label="mean section log odds ratio")
    fig.tight_layout()
    fig.savefig(OUT / "fig1_state_context.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_moran() -> None:
    mi = read_csv("MI_morans_recheck.csv")
    programs = ["cytotoxic", "exhausted", "activated", "Treg", "Th1"]
    samples = ["LUAD-5 R1", "LUAD-5 R2", "LUAD-5 R3"]
    colors = {"LUAD-5 R1": "#4C78A8", "LUAD-5 R2": "#F58518", "LUAD-5 R3": "#54A24B"}
    x = np.arange(len(programs))

    fig, ax = plt.subplots(figsize=(9.5, 4.3))
    width = 0.22
    for i, sample in enumerate(samples):
        d = mi[mi["sample"] == sample].set_index("program").loc[programs]
        ax.bar(x + (i - 1) * width, d["morans_I"], width=width, color=colors[sample], label=sample)
    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(programs, rotation=20, ha="right")
    ax.set_ylabel("Moran's I after n_counts residualization")
    n_sig = int((mi["fdr_bh"] < 0.05).sum())
    ax.set_title(f"Moran re-check with FOV-safe max-distance graph: BH-significant {n_sig}/{len(mi)}")
    ax.legend(frameon=False, ncol=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_moran_recheck.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_fig3() -> None:
    target_x = read_csv("S2c_target_x_corrected_context.csv")
    rep = read_csv("REP_leave_one_section.csv")
    funnel = read_csv("FUN_evidence_funnel.csv")
    boot = read_csv("S2c_bootstrap_rank_ci.csv")
    hit_audit = read_csv("FUN_fdr_hit_audit.csv")

    selected = []
    for state in ["Th1", "exhausted", "cytotoxic"]:
        d = target_x[target_x["state"] == state].nsmallest(5, "corrected_rank")
        selected.extend(d["target"].tolist())
    evidence = hit_audit[
        (hit_audit["per_target_hit"])
        & (hit_audit["donor_reliable"])
        & (hit_audit["no_viability_flag"])
    ]["target"].tolist()
    selected = list(dict.fromkeys(selected + evidence))[:18]

    heat = target_x[target_x["target"].isin(selected)].pivot(
        index="target", columns="state", values="corrected_alignment"
    ).fillna(0)
    heat = heat[["Th1", "exhausted", "cytotoxic"]]

    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 0.85], height_ratios=[1, 1], wspace=0.35, hspace=0.35)

    ax0 = fig.add_subplot(gs[0, 0])
    vmax = max(0.01, np.nanmax(np.abs(heat.values)))
    im = ax0.imshow(heat.values, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), aspect="auto")
    ax0.set_xticks(range(heat.shape[1]))
    ax0.set_xticklabels(heat.columns)
    ax0.set_yticks(range(heat.shape[0]))
    ax0.set_yticklabels(heat.index, fontsize=8)
    ax0.set_title("Corrected target alignment by state-context contrast")
    fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.03, label="cosine alignment")

    ax1 = fig.add_subplot(gs[0, 1])
    rep_plot = rep.copy()
    rep_plot["label"] = rep_plot["state"] + "\nheld " + rep_plot["held_out"].str.split().str[-1]
    ax1.barh(np.arange(len(rep_plot)), rep_plot["spearman"], color="#6C8EBF")
    ax1.set_yticks(np.arange(len(rep_plot)))
    ax1.set_yticklabels(rep_plot["label"], fontsize=8)
    ax1.axvline(0, color="black", linewidth=0.7)
    ax1.set_xlabel("rank Spearman")
    ax1.set_title("Leave-one-section-out target-rank reproducibility")

    ax2 = fig.add_subplot(gs[1, 0])
    if len(boot):
        boot_focus = boot[boot["target"].isin(selected)].copy()
        boot_focus = boot_focus.sort_values(["state", "rank_median"]).head(20)
        y = np.arange(len(boot_focus))
        ax2.errorbar(
            boot_focus["rank_median"],
            y,
            xerr=[
                boot_focus["rank_median"] - boot_focus["rank_q10"],
                boot_focus["rank_q90"] - boot_focus["rank_median"],
            ],
            fmt="o",
            color="#333333",
            ecolor="#999999",
            markersize=4,
        )
        ax2.set_yticks(y)
        ax2.set_yticklabels(boot_focus["target"] + " / " + boot_focus["state"], fontsize=8)
        ax2.invert_yaxis()
        ax2.set_xlabel("bootstrap rank interval")
    ax2.set_title("Bootstrap rank uncertainty for displayed targets")

    ax3 = fig.add_subplot(gs[1, 1])
    funnel_plot = funnel.iloc[::-1].copy()
    ax3.barh(np.arange(len(funnel_plot)), funnel_plot["n_targets"], color="#C44E52")
    ax3.set_yticks(np.arange(len(funnel_plot)))
    ax3.set_yticklabels(funnel_plot["step"], fontsize=8)
    ax3.set_xlabel("target count")
    ax3.set_title("Evidence funnel")
    for i, value in enumerate(funnel_plot["n_targets"]):
        ax3.text(value + max(funnel["n_targets"]) * 0.01, i, str(int(value)), va="center", fontsize=8)

    fig.suptitle("Figure 3 - Corrected spatial contrast and perturbation evidence filtering", fontsize=14)
    fig.savefig(OUT / "fig3_corrected_integration.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def report_numbers() -> dict:
    core = json.load(open(OUT / "phase3_core_summary.json"))
    mi = read_csv("MI_morans_recheck.csv")
    f1 = read_csv("F1_state_context_assoc.csv")
    s2meta = read_csv("S2c_meta.csv")
    sens = read_csv("S2c_sensitivity.csv")
    rep = read_csv("REP_leave_one_section.csv")
    funnel = read_csv("FUN_evidence_funnel.csv")
    clean = read_csv("FUN_clean_shortlist.csv")
    return {
        **core,
        "f1_consistent": int(f1["sign_consistent"].sum()),
        "f1_total": int(len(f1)),
        "mi_sig": int((mi["fdr_bh"] < 0.05).sum()),
        "mi_total": int(len(mi)),
        "mi_mean_abs": float(mi["morans_I"].abs().mean()),
        "s2meta": s2meta.to_dict("records"),
        "sensitivity": sens.to_dict("records"),
        "rep_mean": float(rep["spearman"].mean()) if len(rep) else float("nan"),
        "rep_top20_mean": float(rep["top20_recurrence"].mean()) if len(rep) else float("nan"),
        "funnel": funnel.to_dict("records"),
        "clean_n": int(len(clean)),
    }


def write_reports() -> None:
    n = report_numbers()
    clean_msg_ko = "clean shortlist는 0개였다" if n["clean_n"] == 0 else f"clean shortlist는 {n['clean_n']}개였다"
    clean_msg_en = "the clean shortlist was empty" if n["clean_n"] == 0 else f"the clean shortlist contained {n['clean_n']} targets"

    ko = f"""# Phase 3 Review: 독립 정의와 evidence funnel 재점검

## One-line judgement
이번 추가 결과로 기존의 큰 오류들은 상당히 교정됐지만, clean target-specific spatial story는 아직 완성되지 않았다. 따라서 추가 target-specific spatial map은 만들지 않았다.

## Figure 1: independent state x environment context
Figure 1은 combined feature-niche clustering을 본문 근거로 쓰지 않고, intrinsic state와 environment context를 독립적으로 정의했다.

- Intrinsic state: naive-annotated CD4 T cell의 expression score만 사용했다.
- Environment context: tumor distance/contact와 neighbor cell-type composition만 사용했다.
- Myeloid/fibroblast/B-cell context는 section별 positive-abundance top-quantile rule로 다시 정의했고, context count QC를 `F1_context_qc.csv`에 저장했다.
- 총 {n['f1_total']}개 state-context association 중 {n['f1_consistent']}개가 R1-R3에서 같은 방향이었다.

관련 산출물: `fig1_state_context.png`, `F1_state_context_assoc.csv`, `F1_context_qc.csv`.

## Figure 2: activation-dependent perturbation with reliability/viability
기존 activation-context 결과는 유지하되, 후보 해석은 FDR 방식과 donor/viability 필터를 함께 보아야 한다.

- 기존 pooled-null cytotoxic FDR<0.1 후보는 {n['old_pooled_cytotoxic_fdr_lt_0_1']}개였다.
- per-target add-one permutation 기준 cytotoxic FDR<0.1 후보는 {n['per_target_cytotoxic_fdr_lt_0_1']}개였다.
- per-target FDR<0.1이면서 donor reliable, no viability flag인 후보는 {n['per_target_reliable_no_viability']}개였다.
- 기존 pooled-null 23개 중 donor reliable, no viability flag를 동시에 만족한 후보는 {n['old_pooled_reliable_no_viability']}개였다.

ZMAT1 FDR 0.024는 pooled-null result로만 표기해야 하며, 확정 hit로 사용하지 않는다. `n_downstream=0`은 measured response가 없다는 뜻이 아니라 statistically significant downstream DE genes가 0개라는 뜻이다.

관련 산출물: `FUN_fdr_hit_audit.csv`, `FUN_evidence_funnel.csv`, `FUN_clean_shortlist.csv`.

## Figure 3: corrected spatial contrast x perturbation integration
Contact-vs-far signature는 section별 cutoff, state 내 matching, `n_counts`, state score, local density, tumor-neighbor fraction 보정 후 다시 만들었다.

- Stress/stromal flag와 2/3 section direction agreement를 유지했다.
- Raw vs corrected target ranking은 Spearman 0.974-0.995로 높았지만, 이것은 correction이 모든 문제를 해결했다는 뜻이 아니라 leading rank structure가 크게 바뀌지 않았다는 뜻이다.
- Leave-one-section-out target-rank Spearman 평균은 {n['rep_mean']:.3f}, top-20 recurrence 평균은 {n['rep_top20_mean']:.2f}로 낮았다.
- Evidence funnel 결과, {clean_msg_ko}.

관련 산출물: `fig3_corrected_integration.png`, `S2c_target_x_corrected_context.csv`, `S2c_bootstrap_rank_ci.csv`, `REP_leave_one_section.csv`.

## Moran's I re-check
Moran's I는 FOV-safe max-distance graph와 `n_counts` residualized program score로 재검정했다.

- BH 보정 후 유의한 test는 {n['mi_sig']}/{n['mi_total']}개였다.
- 평균 absolute Moran's I는 {n['mi_mean_abs']:.3f}였다.
- 따라서 현재 해석은 "공간 구조가 약하거나 graph/residualization 설정에 민감하다"가 가장 안전하다.

관련 산출물: `fig_moran_recheck.png`, `MI_morans_recheck.csv`.

## Conclusion and limitations
현재 결과는 더 많은 target을 주장하기보다, 기존 주장의 취약한 부분을 정리하는 데 의미가 있다. R1-R3은 one-patient sections이므로 biological replicates로 해석하면 안 된다. 최종 메시지는 "spatial context can generate testable CRISPRi knockdown hypotheses, but current clean target-specific support is insufficient"로 두는 것이 맞다.
"""

    en = f"""# Phase 3 Review: Independent Definitions and Evidence Filtering

## One-line judgement
The additional analyses fixed several major issues, but the clean target-specific spatial story is still not complete. For that reason, no additional target-specific spatial map was generated.

## Figure 1: independent state x environment context
Figure 1 no longer uses combined feature-niche clustering as primary evidence. Intrinsic state and environment context are defined independently.

- Intrinsic state uses expression scores within naive-annotated CD4 T cells.
- Environment context uses tumor distance/contact and neighbor cell-type composition.
- Myeloid/fibroblast/B-cell contexts use a section-wise positive-abundance top-quantile rule, with context count QC saved in `F1_context_qc.csv`.
- {n['f1_consistent']} of {n['f1_total']} state-context associations had the same direction across R1-R3.

Outputs: `fig1_state_context.png`, `F1_state_context_assoc.csv`, `F1_context_qc.csv`.

## Figure 2: activation-dependent perturbation with reliability/viability
The activation-context result is retained, but target interpretation must combine FDR definition, donor reliability, and viability flags.

- The old pooled-null cytotoxic FDR<0.1 set contained {n['old_pooled_cytotoxic_fdr_lt_0_1']} targets.
- The per-target add-one permutation cytotoxic FDR<0.1 set contained {n['per_target_cytotoxic_fdr_lt_0_1']} targets.
- Per-target FDR<0.1 plus donor reliable plus no viability flag left {n['per_target_reliable_no_viability']} targets.
- Among the old pooled-null 23 targets, donor reliable plus no viability flag left {n['old_pooled_reliable_no_viability']} targets.

ZMAT1 FDR 0.024 should be described only as a pooled-null result, not as a confirmed hit. `n_downstream=0` means zero statistically significant downstream DE genes, not absence of measured expression response.

Outputs: `FUN_fdr_hit_audit.csv`, `FUN_evidence_funnel.csv`, `FUN_clean_shortlist.csv`.

## Figure 3: corrected spatial contrast x perturbation integration
The contact-vs-far signature was rebuilt using section-specific cutoffs, within-state matching, and adjustment by `n_counts`, state score, local density, and tumor-neighbor fraction.

- Stress/stromal flags and >=2/3 section direction agreement are retained.
- Raw vs corrected target ranking Spearman values were 0.974-0.995; this means the rank structure was stable, not that the confound is fully solved.
- Leave-one-section-out target-rank Spearman averaged {n['rep_mean']:.3f}, and top-20 recurrence averaged {n['rep_top20_mean']:.2f}.
- After the full evidence funnel, {clean_msg_en}.

Outputs: `fig3_corrected_integration.png`, `S2c_target_x_corrected_context.csv`, `S2c_bootstrap_rank_ci.csv`, `REP_leave_one_section.csv`.

## Moran's I re-check
Moran's I was re-tested with a FOV-safe max-distance graph and `n_counts`-residualized program scores.

- BH-significant tests: {n['mi_sig']}/{n['mi_total']}.
- Mean absolute Moran's I: {n['mi_mean_abs']:.3f}.
- The safest interpretation is that spatial structure is weak or sensitive to graph/residualization choices.

Outputs: `fig_moran_recheck.png`, `MI_morans_recheck.csv`.

## Conclusion and limitations
The value of the current result is not a larger target list, but a more defensible audit of the central claim. R1-R3 are one-patient sections and should not be interpreted as biological replicates. The final message should be: spatial context can generate testable CRISPRi knockdown hypotheses, but current clean target-specific support is insufficient.
"""

    (DOC / "2026-07-13_phase3_review_ko.md").write_text(ko, encoding="utf-8")
    (DOC / "2026-07-13_phase3_review_en.md").write_text(en, encoding="utf-8")


def main() -> None:
    save_fig1()
    save_moran()
    save_fig3()
    write_reports()
    print("PHASE3_REPORT_DONE", flush=True)


if __name__ == "__main__":
    main()
