# Hackathon Progress — Report Draft

Running summary of **progress, results, and methods** only (failures and dead
ends are excluded — those live in the dated logs under `1-Doc/`). This file is
the working draft for the final report and is updated continuously.

## Objective
Test whether T cell-intrinsic regulators act on spatial immune programs in an
activation-context-specific manner, and whether spatial structure predicts
where each regulator acts. Hypothesis generation for target discovery — not
causal claims.

## Data
- **Perturbation:** Marson–Pritchard genome-scale CD4+ T cell Perturb-seq
  atlas (CRISPRi, 4 donors, conditions Rest/Stim8hr/Stim48hr). Author-provided
  `GWCD4i.DE_stats.h5ad` (per-perturbation × gene logFC/zscore) and
  `GWCD4i.pseudobulk_merged.h5ad` (NTC + targeting pseudobulks).
- **Spatial:** CosMx NSCLC (He 2022), 765,771 cells × 960-gene panel, 8 samples.
  Primary substrate: **all CD4 T cells** (naive + memory + Treg) in LUAD-5 R1-3.

## Design decisions
- **Cell-identity first, niche as context.** T cells are classified by intrinsic
  state (7 programs); spatial niche is the layer that says *where* a matched
  state sits, not the classification unit.
- **Scope = all CD4** (not naive-only): naive-only makes exhaustion/Treg/Th1
  endpoints degenerate (scores ≈ 0), so the substrate was expanded.
- **Primary result = Rest vs Stim direction-flip.** Per-cell soft context
  projection is exploratory only, because per-cell state-matching collapses to
  Rest (argmax=Rest 81.7%).
- **Two gene overlaps, kept separate:** state-matching uses 912 genes
  (CosMx ∩ pseudobulk); gene-level/held-out confirmation uses 501
  (CosMx ∩ DE_stats, ⊂ 912). Overlap size reported with every score.

## Methods
1. **Gene-ID harmonization.** Map CosMx HGNC symbols to Perturb-seq
   `var['gene_name']`; freeze 912 (state-matching) and 501 (confirmation) shared
   genes. Zero duplicate symbol→Ensembl mappings.
2. **Spatial niches (all-CD4).** Per cell, a 4-block feature vector — (1) 7
   intrinsic state signatures (`score_genes`), (2) Gaussian-weighted neighbor
   cell-type composition (150µm), (3) neighbor suppressive/ligand program,
   (4) log distance to nearest malignant cell — each block z-scored and
   sqrt(n)-normalized so all contribute equally; PCA → kNN → Leiden.
3. **Perturbation program effects.** From DE_stats logFC, E_{g,c,program} =
   mean logFC over the program's markers (target gene excluded), per condition
   c ∈ {Rest, Stim8hr, Stim48hr}, on QC-pass perturbations
   (on-target significant, no distal off-target, adequate target expression).
4. **Direction-flip classification.** For each (target, program), compare
   E_Rest vs E_Stim48hr: direction_flip (sign change, both |E|≥0.10),
   rest_specific, stim_specific, or stable.

## Results
- **Gene overlaps frozen:** 912 (state-matching), 501 (gene-level, ⊂ 912).
- **All-CD4 substrate validated (non-degeneracy):** across 25,493 CD4 cells,
  program scores are now real — fraction-positive exhausted 0.53, cytotoxic
  0.35, Treg 0.35, Th1 0.27 (≈0 under naive-only). Treg cells carry the highest
  Treg (1.60) and exhausted (0.59) scores.
- **19 all-CD4 feature-niches.** Mature states separate into dedicated niches:
  FN1 = Treg (86.8%), FN2 = memory (80.2%); niches also differentiate on
  program axes (Tfh-high FN10; Th1-high FN4/FN15) and tumor-boundary distance.
- **Direction-flip screen (PRIMARY):** 7,828 perturbations QC-passed; of
  54,796 target×program pairs, **5,125 are Rest↔Stim48hr direction-flips**,
  concentrated in Tfh (1266), cytotoxic (980), Treg (817), exhausted (719).
  Positive control: TCR-signaling KOs (ZAP70, CD3E, VAV1, PTPRC) flip the
  cytotoxic program from negative (Rest) to positive (Stim48hr), consistent
  with known activation biology.

- **Niche × target relevance map (tracks combined):** Relevance(niche, target,
  endpoint) = max(0, niche program burden) × max(0, direction·E_Rest). Target
  niche per endpoint = argmax-burden niche (exhaustion FN10, Treg FN1, Th1 FN4,
  activation FN10, cytotoxic FN4). 586,476 positive relevance records.
- **Held-out gene confirmation (circularity control):** recomputing niche↔
  perturbation cosine on 471 genes (501-shared minus 30 program markers), vs a
  matched null of all 6,734 QC-pass Rest perturbations, **confirms** the
  activation (null-pct 0.93), Th1 (0.89; PTPN2 0.98) and cytotoxic (0.86)
  endpoints on independent genes, but does **not** confirm exhaustion (0.14) or
  Treg (0.39) — those ranks were marker-driven and are down-weighted.
- **Literature support (OpenAlex-retrieved):** PTPN2 (ABBV-CLS-484 2023; primary
  human T-cell CRISPR screens 2018) and BHLHE40 (CD8 exhaustion CRISPR 2023).
  Most defensible single hit = **PTPN2** (literature + held-out).

## Status
- 2026-07-08: Documentation structure established.
- 2026-07-11: Niche pipeline built; scope set to all-CD4; 912/501 frozen;
  Track A (niche + validation) and Track B (direction-flip, primary) complete.
- 2026-07-11 (cont.): tracks combined into relevance map; held-out
  confirmation run (Th1/activation/cytotoxic hold, exhaustion/Treg do not);
  top candidates literature-checked (PTPN2 strongest).
- **Next:** direction-flip trajectory figures for confirmed candidates
  (PTPN2 etc.); revisit exhaustion/Treg endpoints with a marker-free niche
  definition; cross-patient validation remains open (single-patient LUAD-5).
