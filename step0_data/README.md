# Step 0 — Data acquisition

This project joins two public datasets that were never designed to be linked.
Neither requires credentials. Run [`download_data.sh`](download_data.sh) to pull
them, or follow the manual instructions below.

```bash
bash download_data.sh /path/to/large/storage/0-Data   # use scratch, not $HOME
```

---

## Dataset 1 — Perturbation atlas (Marson–Pritchard CD4⁺ Perturb-seq)

A genome-scale CRISPRi Perturb-seq atlas in primary human CD4⁺ T cells. Every
expressed gene is knocked down one at a time and the transcriptome-wide
consequence is read out by single-cell RNA-seq, **in both a resting and a
stimulated state** — the feature that makes the activation-context question in
this project possible.

| | |
|---|---|
| **Paper** | Zhu R., Dann E., *et al.* (2025) *Genome-scale perturb-seq in primary human CD4+ T cells maps context-specific regulators of T cell programs and human immune traits.* bioRxiv, doi:10.64898/2025.12.23.696273 (CC-BY 4.0) |
| **Data portal** | https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq |
| **Authors' analysis code** | https://github.com/emdann/GWT_perturbseq_analysis_2025 |
| **Raw reads** | SRA `SRP643211` / GEO `GSE314342` |
| **S3 bucket** | `s3://genome-scale-tcell-perturb-seq/marson2025_data/` (public, `--no-sign-request`; ~1.7 TB, 32 objects) |

**Scale (from the paper):** 22M primary CD4⁺ T cells perturbed; 4 donors (D1–D4);
3 culture conditions (`Rest`, `Stim8hr`, `Stim48hr`); 12,748 genes targeted,
~2 guides/gene; DE computed per condition with DESeq2 (pseudobulk by
donor×condition×guide, regressing out donor and cell count).

**Files this project uses (~57 GB core):**

| File | Size | Role |
|---|---|---|
| `GWCD4i.DE_stats.h5ad` | 15.6 GB | perturbation × gene `log_fc` / `zscore` per condition — **the response vectors** (10,282 measured genes, 33,983 perturbation×condition rows) |
| `GWCD4i.pseudobulk_merged.h5ad` | 41.5 GB | pseudobulk profiles for QC and for cross-platform state matching (18,129 genes) |
| `suppl_tables/*.csv` | ~15 MB | sample metadata, guide-library annotation, flat DE summary |

The 1.6 TB cell-level `D{1-4}_{Rest,Stim8hr,Stim48hr}.assigned_guide.h5ad`
matrices are **not** needed here — this pipeline runs on the perturbation
response vectors, not raw cells.

---

## Dataset 2 — Spatial atlas (CosMx human NSCLC)

Single-cell-resolution spatial transcriptomics of non-small-cell lung cancer.

| | |
|---|---|
| **Source data** | He, S. *et al.* (2022) *High-plex imaging of RNA and proteins at subcellular resolution in fixed tissue by spatial molecular imaging.* Nat. Biotechnol. — NanoString CosMx NSCLC public release |
| **Distribution used** | `cosmx_human_nsclc_clustered.h5ad` (~2.6 GB), the clustered AnnData carrying the CellCharter niche labels used as a spatial-context reference here |
| **Scale** | 765,771 cells × 960-gene panel, 8 samples, 22 annotated cell types |

> **Note on provenance.** The clustered `.h5ad` with CellCharter niche labels was
> obtained from the CellCharter NSCLC distribution rather than assembled from raw
> flat files in this project. If your copy has a different filename, point the
> `H5` path at the top of each spatial step's script at it.

**Schema this project relies on** (verified in [`step1_data_inspection`](../step1_data_inspection/)):

- `.X` — log-normalized expression; raw UMI counts in `.layers['counts']`.
- `.obsm['spatial']` — 2-D cell coordinates (each sample occupies a disjoint
  coordinate region; neighbors are always computed within a sample).
- `.obs['cell_type']` — 22 types, including 5 T-cell subtypes: T CD4 naive
  (27,988), T CD4 memory (16,836), T CD8 naive (10,452), T CD8 memory (5,503),
  Treg (6,046).
- `.obs['niche']` — 9 CellCharter-provided niche labels (tumor interior, stroma,
  tumor–stroma boundary, myeloid-enriched stroma, **lymphoid structure** [TLS-like],
  immune, neutrophils, plasmablast-enriched stroma, macrophages).
- `.obs['sample']` — 8 samples; the primary analysis substrate is the three
  LUAD-5 replicate sections (R1–R3).

---

## Gene-identifier reconciliation (important)

The two platforms use different gene identifiers, which silently produced empty
intersections early in the project:

- CosMx `.var_names` are **HGNC symbols** (AATK, ABL1, …).
- Perturb-seq `.var_names` are **Ensembl IDs** (ENSG…); symbols live in
  `.var['gene_name']`.

All cross-platform steps map CosMx symbols through the Perturb-seq `gene_name`
column. Two overlaps are frozen and kept **separate** (see
[`step1_data_inspection`](../step1_data_inspection/)):

- **912 genes** — CosMx ∩ pseudobulk symbols (used for state matching).
- **501 genes** — CosMx ∩ `DE_stats` symbols, a proper subset of the 912 (used
  for gene-level / held-out confirmation).

Reporting these as one number would overstate the validation power, so every
score reports the overlap it used.
