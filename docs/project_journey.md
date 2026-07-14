# SpatialPerturb — Project Journey & Process Log

> A narrative account of **how this project actually unfolded** — the path taken,
> the problems that forced detours, the decisions that were hard to make, and
> where the work is heading. This is the "story" companion to the failure-free
> `hackathon_progress.md` and the dated `1-Doc/` attempt logs; it exists to let
> anyone (including a future session) understand not just *what* was done but
> *why the project looks the way it does now*.
>
> EN master. Korean twin: `project_journey_ko.md`.
> Last updated: 2026-07-11.

---

## 0. What the project is (in one paragraph)

**SpatialPerturb** asks a single practical question: *where inside a tumor would
a T-cell perturbation actually matter?* We connect two datasets that were never
designed to be linked — the **Marson–Pritchard genome-scale CD4⁺ T cell
Perturb-seq atlas** (what a gene knockdown does to a CD4⁺ T cell, measured in
**Rest / Stim8hr / Stim48hr** states) and **tumor spatial transcriptomics**
(where T-cell states sit relative to tumor, myeloid, and stromal neighbors).
The goal is **hypothesis generation for target discovery**, not causal claims:
produce a *target × niche* priority map, scored **separately for resting and
stimulated states** so that regulators whose effect flips with activation
context surface explicitly. Every projection is labeled exploratory.

---

## 1. The arc so far — five phases

The project did not proceed in a straight line. It went through five recognizable
phases, and two of them were pivots forced by what we learned.

### Phase 1 — Planning & framing (Day 0, ~07-06/07)
A seven-day plan was written (`hackathon_plan.md`, ~28 KB), notable for putting
**validity and over-interpretation risk first** (§1) rather than as an
afterthought. The plan explicitly names the project's inferential ceiling:
valid as a *prioritization / hypothesis-generation* exercise, **not** as
evidence that a knockout causally remodels a niche. The single largest scientific
risk was identified up front — **cell-type confounding** (a "T-cell program"
spatial signal can be driven by T-cell *abundance*, not *state*) — and the
single largest methodological risk — **unstable transfer on a small shared gene
set**.

### Phase 2 — Data reconnaissance & dataset lock (07-08/09)
This phase resolved the four "Day-1 blocker" questions before any heavy compute:
1. Does the Perturb-seq distribution ship precomputed per-condition signatures?
   → **Yes** (`DE_stats.h5ad`, 10,282 measured genes) — so the 22M-cell raw data
   never had to be reprocessed. This is the decision that made a one-week project
   feasible at all.
2. Which spatial dataset? → **CosMx NSCLC locked as primary** (see §2, D1).
3. Deconvolution needed? → **No** — CosMx is single-cell with shipped `cell_type`.
4. Is the shared gene space large enough? → **501 shared genes** at that point,
   above the stability floor, but small enough that **pathway/module-level
   transfer became the default** and gene-level scores must always report the
   overlap size.

### Phase 3 — The reframe: "outcomes" → "where would it matter" (07-11)
The project's central question was **reframed by the user** from *predict the
outcome of a perturbation in tissue* to *prioritize which Perturb-seq targets
matter in which tumor niche*. This is more than wording. It changed the method
order: **niche definition became Step 1**, and Moran's I / spatial-autocorrelation
was **demoted from "the way programs are extracted" to "a validation that the
niches are genuinely spatially organized."** The scoring object became a product:
`Spatial Priority = Response Compatibility × Desired Program Change × Niche
Specificity`. This reframe is why the current work centers on niches, not on a
Moran's-I gene list.

### Phase 4 — Building niches, bottom-up (07-11, Sessions 1–2)
Rather than trust the atlas's shipped CellCharter niche labels, we built T-cell
niches from scratch: per-T-cell **Gaussian-kernel neighborhood expression**
(σ=50 µm, 3σ=150 µm cutoff, *set from the measured median nearest-neighbor
distance of 46 µm*, not a guessed constant) → PCA → kNN → Leiden. The result was
20 T-cell niches; compared to CellCharter, ARI=0.208 / NMI=0.338 — moderate
agreement, with the new niches **subdividing** the reference (TLS and immune
niches each split into several sub-niches). An **unprompted caveat** was raised:
the niches show a strong per-sample/batch effect because neighborhood expression
was not batch-corrected.

### Phase 5 — Narrowing to CD4 naive + state matching (07-11, Session 3)
Scope was **deliberately narrowed** to CD4 naive cells in LUAD-5 replicates
R1–R3 only (see §2, decision D3). On that anchor: 21 CD4 naive niches, and a
**state-matching** analysis showing tumor CD4 naive cells are closest to the
Perturb-seq **Rest** state — argmax Pearson = Rest for **81.7%** of cells and for
**all 21 niches** (mean Pearson Rest 0.194 > Stim48hr 0.175 > Stim8hr 0.170).
This is the biologically expected result (tumor-resident naive cells are not
acutely TCR-stimulated) and it anchors *which* Perturb-seq condition the
downstream projection should lean on.

---

## 2. The decisions that were genuinely hard

These are the forks where the answer was not obvious and real deliberation (often
with the user) was required.

**D1 — Which spatial dataset, and how many.**
The requirement itself *changed mid-project*: an initial Visium HD colorectal
choice was dropped once we realized Visium HD is bin-grid data with **no
single-cell segmentation and no cell-type annotation** (Space Ranger 3.0 ships no
cell boundaries). The requirement hardened to "single-cell resolution *with*
existing cell-type annotation." A CosMx **breast** option was investigated and
abandoned — no clean public download. The landing point: **CosMx NSCLC as
primary** (public, annotated, has CD4 subsets, ships a TLS-like 'lymphoid
structure' niche as a built-in positive control) **+ CosMx BCC as secondary**.
The difficulty here was balancing "clear immune architecture" against "actually
publicly downloadable and annotated" — several attractive datasets failed the
second test.

**D2 — Classify by niche, or by cell identity?** *(a real tension)*
A core conceptual decision: when connecting Perturb-seq to space, is the unit of
classification the **niche** (where the cell is) or the **cell identity/state**
(what the cell is)? The resolution — grounded in ProjecTILs, Yeh et al.
Nat Immunol 2024, and SpatialProp — was **cell-identity-first**: classify T cells
by state, and let the niche become the *spatial-context layer* (where the matched
states land), not the classification unit. This decision is still the live
organizing principle and is why Session 3 pivoted to "CD4 naive cells" as the
anchor rather than "the TLS niche" as the anchor.

**D3 — How far to narrow scope.**
Whether to keep all 8 samples and all 5 T subtypes, or narrow. Narrowed hard, and
for a defensible reason: **Perturb-seq is a CD4-axis experiment — CD8 was never
profiled**, so projecting perturbation responses onto CD8 cells is unjustified;
and Treg appears only as an in-vitro-arising state. So: **CD4 first, naive first,
LUAD-5 replicates only** (within-patient reproducibility before cross-patient
generalization). The cost is obvious (we are looking at one histology, one T
subtype); the benefit is that every downstream claim stays inside what the
Perturb-seq data can actually support.

**D4 — Pathway-level vs gene-level transfer.**
Because the CosMx panel shares only a few hundred genes with the whole-
transcriptome Perturb-seq axis, gene-level transfer is statistically shaky. The
decision: **default to pathway/module-level transfer** (activation, cytotoxicity,
exhaustion, IFN, proliferation, antigen presentation, metabolic stress), with
gene-level cosine/GSEA on the shared genes as a *secondary* check that always
reports the overlap size.

**D5 — Reference-projection method: ProjecTILs vs signature scoring.**
ProjecTILs offers ready human CD4⁺ TIL reference atlases, but it is designed for
full-transcriptome scRNA-seq; on a 960-gene imaging panel with sparse counts its
projection quality may degrade. The **recommended main pipeline** leans toward
**robust signature scoring** for the CD4 identity axis, with ProjecTILs as an
alternative — and SpatialProp/Celcomen-style GNN tissue-propagation noted as a
*more sophisticated follow-up deliberately deferred* as too heavy for the
hackathon timeframe.

---

## 3. Where problems bit us (and how they were resolved)

The dead ends and the fixes — the part usually invisible in a clean report.

**P1 — `shared genes = 0` (the silent killer).**
The first state-matching run produced a **blank** violin plot. Root cause, found
in the log: **CosMx `var_names` are HGNC gene symbols** (AATK, ABL1…) while
**pseudobulk `var_names` are Ensembl IDs** (ENSG…), with symbols hidden in
`var['gene_name']`. A naive intersection was empty. Fix: map CosMx symbols
through the pseudobulk `gene_name` column → **912 shared genes**. This also
reconciled a confusing "501 vs 912" discrepancy in the notes (501 was an earlier
filter against `DE_stats`; 912 is the symbol intersection with the pseudobulk).
Standing rule adopted: **always verify shared-gene count > 0 before running any
similarity**, because an empty gene set silently produces blank plots and later a
PCA `n_components=-1` crash.

**P2 — Batch effect in the niches.**
The first (8-sample) niche run showed several samples collapsing to a single
niche (LUAD-12 → 94.6% one niche; LUSC-6 → 88.8% one niche). Diagnosed correctly
as clustering on **non-batch-corrected** neighborhood expression, so sample-level
expression differences dominate. This was surfaced as a caveat, not hidden;
narrowing to LUAD-5 replicates (Session 3) reduced it but did not eliminate it
(some CN niches are still replicate-specific). Batch correction (harmony on the
neighborhood-expression PCA) remains a candidate next step.

**P3 — Auditor rejection: neighbor composition vs sample distribution.**
An automated review flagged that the niche annotation had used **per-sample
distribution as a stand-in for neighbor cell-type composition** — a different
quantity than the plan required. Corrected by computing the real Gaussian-weighted
neighbor composition per niche in a follow-up job; the fix validated itself
biologically (TLS-mapped niches show high B-cell neighbor fractions, consistent
with TLS B/T zones).

**P4 — Infrastructure friction.**
Several non-scientific problems cost real time and shaped how work is now done:
- **PACE-Phoenix login node intermittently unresponsive** (even `echo alive`
  times out). Mitigation: save deliverables as artifacts (server-independent) and
  defer server file writes until a connectivity probe succeeds.
- **HOME is 20 GB and ~77% full** → all conda envs and outputs go to the
  **scratch** prefix, never HOME.
- **matplotlib RGB-tuple colors** passed to a pandas categorical `.map()` raised
  `TypeError: Expected tuple, got str`; fix is `to_hex` + `.astype(str)`.
- **Large base64 blobs on the command line** hit "Broken pipe" (arg-length limit)
  → write files individually via per-file heredocs.
- **Korean matplotlib labels** rendered as tofu until AppleSDGothicNeo was
  registered.
- **Translation truncation**: long docs translated via LLM can silently truncate
  at `max_tokens` mid-document, dropping bibliography entries — now caught by an
  EN-vs-KO item-count diff at the bullet level, not just headers.

---

## 4. Current state (as of 2026-07-11)

**Done:**
- Plan, data-description, related-work, and background literature review docs
  (all bilingual EN master + KO twin).
- Data downloaded (~75 GB on scratch): Perturb-seq `DE_stats` + pseudobulk +
  supplementary tables; CosMx NSCLC clustered `.h5ad`; CosMx BCC `.rds`.
- `env_spatial` built and verified (scanpy 1.11.5, squidpy 1.8.2, anndata
  0.12.19).
- Perturb-seq EDA / sanity check passed (on-target KD ~62%; TCR genes are
  stim-specific; ~101 Stim-specific and ~152 Rest-specific regulators).
- CosMx atlas visualized; T-cell-rich niches flagged (lymphoid structure 3.4×,
  immune 3.1× enrichment).
- T-cell niches built two ways (8-sample all-T; CD4-naive LUAD-5).
- State matching: tumor CD4 naive ↔ Perturb-seq **Rest** established.

**Not yet done (the real remaining science):**
- **Step 2 proper** — rank Perturb-seq CD4 perturbations that strengthen/suppress
  each spatially organized program, **scored Rest vs Stim separately**, with
  context-flip flagging, permutation nulls, and FDR.
- **Step 3** — the `Spatial Priority` target × niche map itself (the primary
  deliverable), plus the framework figure and the labeled spatial atlas.
- **Moran's I validation** that the niches/programs are spatially organized (now
  a validation step, per the reframe).
- Optional: batch correction; BCC reproduction; `pert2state_model` / CONCERT
  comparison; Streamlit demo.

---

## 5. Where the project is heading

The organizing logic is now settled: **classify T cells by identity/state
(CD4 axis, naive first), use the niche as the spatial-context layer, transfer
Perturb-seq responses at the pathway level, and score Rest vs Stim separately so
context-flippers are the headline output.** The immediate next move is Step 2 —
turning the response-vector matrix into per-program, per-condition perturbation
rankings — followed by the Step 3 target × niche map that is the project's actual
deliverable. Everything downstream of a minimum viable submission (one T subtype,
one histology, ranked per-state perturbations with nulls, one labeled niche ×
target map, honest limitations) is strengthening, and the plan's cut-lines are
explicit if time runs short.

The honest framing has been preserved throughout: this is a **prioritization
instrument and hypothesis generator**, every map is an **exploratory projection**,
and the in-vitro→in-tissue gap plus the ~900-gene shared space are stated with
every score — not buried in a limitations paragraph.
