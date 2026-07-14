#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# SpatialPerturb — data download
#
# Two public datasets are joined in this project. Neither requires credentials.
# Set DEST to a large-capacity location (scratch on an HPC cluster); the
# perturbation core pull is ~57 GB and the CosMx spatial file is ~2.6 GB.
# ---------------------------------------------------------------------------
set -euo pipefail

DEST="${1:-./0-Data}"
mkdir -p "$DEST/perturbseq_marson2025/suppl_tables" "$DEST/ST_CosMx_NSCLC"

# ===========================================================================
# 1. Perturbation atlas — Marson–Pritchard genome-scale CD4+ Perturb-seq
#    Paper : Zhu R., Dann E. et al. (2025) bioRxiv doi:10.64898/2025.12.23.696273 (CC-BY 4.0)
#    Portal: https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
#    Code  : https://github.com/emdann/GWT_perturbseq_analysis_2025
#    Raw   : SRA SRP643211 / GEO GSE314342
#    Bucket: s3://genome-scale-tcell-perturb-seq/marson2025_data/  (public, anonymous)
# ---------------------------------------------------------------------------
# This project is driven by the perturbation *response vectors*, not raw cells,
# so only the DE-stats + pseudobulk + supplementary tables are needed (~57 GB).
# Do NOT pull the 1.6 TB cell-level matrices unless single-cell resolution is
# required.
PB="s3://genome-scale-tcell-perturb-seq/marson2025_data"
PBD="$DEST/perturbseq_marson2025"

# Option A — AWS CLI (recommended). Always use --no-sign-request (anonymous).
aws s3 cp --no-sign-request "$PB/GWCD4i.DE_stats.h5ad"          "$PBD/"   # 15.6 GB — perturbation x gene log_fc/zscore per condition (the response vectors)
aws s3 cp --no-sign-request "$PB/GWCD4i.pseudobulk_merged.h5ad" "$PBD/"   # 41.5 GB — pseudobulk profiles for QC / state matching
aws s3 cp --no-sign-request --recursive "$PB/suppl_tables/"     "$PBD/suppl_tables/"   # ~15 MB — metadata, guide library, flat DE summary

# Option B — no AWS CLI, plain HTTPS (per-object direct URL):
#   curl -L -o "$PBD/GWCD4i.DE_stats.h5ad" \
#     https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.DE_stats.h5ad

# ===========================================================================
# 2. Spatial atlas — CosMx human NSCLC (He et al. 2022)
#    765,771 cells x 960-gene panel, 8 samples, 22 cell types, single-cell
#    resolution with provided CellCharter niche labels (incl. a TLS-like
#    "lymphoid structure" niche used as a positive control).
# ---------------------------------------------------------------------------
# The clustered AnnData (cosmx_human_nsclc_clustered.h5ad, ~2.6 GB) is the file
# every spatial step reads. It was obtained from the CellCharter NSCLC tutorial
# distribution. If your copy is under a different name, point CFG_H5 in each
# step at it. See step0_data/README.md for the exact provenance and the .obs /
# .obsm schema this project relies on.
echo "Place cosmx_human_nsclc_clustered.h5ad under $DEST/ST_CosMx_NSCLC/"
echo "  (see step0_data/README.md for the source link and expected schema)"

echo "Done. Core perturbation pull -> $PBD ; spatial file -> $DEST/ST_CosMx_NSCLC/"
