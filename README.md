# Perturb-to-Spatial

⭐ **Research Report:** [report/hackathon_report_en.md](report/hackathon_report_en.md)

## 🎬 Video walkthrough (≈2 min)

<!-- INLINE PLAYER: drag `video/out/report_web.mp4` (8 MB) into a GitHub issue or
     PR comment box, wait for it to upload, copy the
     https://github.com/YeojinKim220/perturb-to-spatial/assets/... URL it inserts,
     and paste that URL on the line below as a bare URL (no [](...) wrapper).
     GitHub then renders an inline video player here. -->

_Inline player: paste the uploaded `assets/...mp4` URL on the line above this note._

Direct files: [report_web.mp4](https://github.com/user-attachments/assets/c7580294-0315-4806-8870-bd3c29b1861f) (8 MB, compressed) · [report.mp4](video/out/report.mp4) (18 MB, full quality) — a figure-by-figure walkthrough of the report.

**Where inside a tumor would a T-cell perturbation matter — and does it depend on
whether the T cell is resting or activated?**

This project connects a genome-scale CD4⁺ T-cell Perturb-seq atlas (what
happens inside a T cell when each gene is knocked down, measured in resting and
stimulated states) to tumor spatial transcriptomics (where T-cell states sit
relative to tumor, myeloid, and stromal neighbors), to produce a prioritized,
activation-context-annotated shortlist of candidate regulators per immune-program
endpoint per spatial niche.

> **What this is.** A **hypothesis-generation instrument for target discovery** —
> not a source of causal or predictive claims. Every projected map is labeled
> *exploratory — not experimentally validated*. After proper permutation nulls,
> multiple-testing correction, and reliability/viability filtering, the most
> honest headline of the current round is a **calibrated negative result**: FDR
> alone does not select a clean target-specific spatial shortlist here. That
> negative — and the audit that produced it — is the deliverable.

---

## Pipeline

Each step is a self-contained folder with its own README. Run them in order.

| Step | Folder | What it does |
|---|---|---|
| 0 | [`step0_data`](step0_data/) | Download the two public datasets (perturbation + spatial) and reconcile gene identifiers |
| 1 | [`step1_data_inspection`](step1_data_inspection/) | Inspect both AnnData schemas; freeze the 912 / 501 shared-gene sets |
| 2 | [`step2_spatial_visualization`](step2_spatial_visualization/) | Visualize the CosMx atlas by cell type; flag T-cell-rich niches (TLS-like positive control) |
| 3 | [`step3_tcell_niche`](step3_tcell_niche/) | Build T-cell niches from neighborhood expression (Leiden); compare to reference labels |
| 4 | [`step4_cd4_feature_niche`](step4_cd4_feature_niche/) | Interpretable 4-block CD4⁺ feature-niches; validate non-degeneracy of state scores |
| 5 | [`step5_state_matching`](step5_state_matching/) | Match tumor CD4⁺ cells to a Perturb-seq condition (→ Rest); anchor the projection |
| 6 | [`step6_perturbation_screen_relevance`](step6_perturbation_screen_relevance/) | Direction-flip screen + first niche × target relevance map *(exploratory; superseded)* |
| 7 | [`step7_phase2_corrections`](step7_phase2_corrections/) | Permutation nulls, BH-FDR, reliability & viability filters — the correction that bites |
| 8 | [`step8_phase3_reanalysis`](step8_phase3_reanalysis/) | Independent state/context definitions, FOV-safe spatial test, evidence funnel → calibrated negative |

Steps 2–8 are the *analysis order*; steps 6→7→8 are also a *chronological
refinement* — an exploratory first pass (6), a statistical correction (7), and a
final independent-definition re-analysis (8). Reading them in order shows how an
apparently positive shortlist became a defensible negative under proper controls.

## Organizing logic

1. **Classify each tumor T cell by what it *is*** — its intrinsic
   activation/differentiation state (7 programs).
2. **Use the spatial niche as a *context layer*** telling us *where* a given
   state sits — not as the unit of classification.
3. **Transfer the Perturb-seq response at the pathway/module level**, because the
   imaging panel shares only a few hundred genes with the whole-transcriptome
   perturbation axis.
4. **Score resting and stimulated states separately**, so regulators whose effect
   *reverses* with activation context become visible rather than averaged away.

## Data

Two public, credential-free datasets — see [`step0_data/README.md`](step0_data/README.md)
for full provenance, sizes, and download commands.

- **Perturbation:** Marson–Pritchard genome-scale CD4⁺ Perturb-seq (Zhu, Dann
  *et al.* 2025; bioRxiv doi:10.64898/2025.12.23.696273; portal on the CZI
  Virtual Cells Platform; GEO `GSE314342`). Driven by the `DE_stats` response
  vectors (~57 GB core pull).
- **Spatial:** CosMx human NSCLC (He *et al.* 2022), 765,771 cells × 960 genes,
  clustered AnnData with CellCharter niche labels (~2.6 GB).

## Conventions

- Human, **GRCh38** gene symbols throughout.
- **Random seed fixed** (`seed=0`) in every scored step.
- Resting vs stimulated states scored **separately**; direction-flipping
  regulators flagged.
- Every cross-platform score reports the **shared-gene overlap size** it used.
- Significance via **permutation / spatial-rotation nulls** with FDR per
  (program, state).
- **R1–R3 are three sections of one patient** — within-patient reproducibility,
  not biological replicates.

## Environment

Two conda environments were used on the HPC cluster (PACE-Phoenix, SLURM):

- **spatial stack** (steps 2–4): `scanpy`, `squidpy`, `anndata`, `scipy`,
  `scikit-learn`, `leidenalg`, `igraph`, `matplotlib`.
- **analysis stack** (steps 1, 5, 7, 8): the above plus `h5py`, `pandas`,
  `statsmodels`.

See [`requirements.txt`](requirements.txt) for a pip-installable superset. Any
step over ~1 minute was run via `sbatch` on the cluster, with outputs written to
scratch.

## Reproducing

```bash
# 0. get the data (use a large-capacity location, e.g. cluster scratch)
bash step0_data/download_data.sh /path/to/0-Data

# 1-8. run each step in order (submit as SLURM jobs on a cluster)
#      edit the H5 / BASE path at the top of each script to point at your 0-Data
python step1_data_inspection/overlap.py
python step2_spatial_visualization/spatial_viz.py
# ... etc.
```

> **Path note.** The scripts were written with absolute cluster paths at the top
> (`BASE = /storage/.../hackathon`, `H5 = .../cosmx_...h5ad`). Edit those to your
> local `0-Data` location before running. The path is always the first few lines
> of each script.
