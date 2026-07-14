# Claude Use — SpatialPerturb (Researcher Track)

> Submission criterion (25%): *How creatively did this team use Claude Code? Did they go beyond a basic application? Did they surface capabilities that surprised even us?*
>
> This document is a running log of **how Claude Science (Claude Code) was used** to drive the project. Updated every work session. Newest entries on top.

## Summary of the mode of use

We did not use Claude as a chatbot that hands back snippets. We used **Claude Science as an autonomous research operator** that:
- **Drives a remote HPC cluster end to end** (PACE-Phoenix, Slurm). Claude reads `compute_details`, builds conda environments at a scratch prefix, writes Slurm job scripts (correct account/QOS/partition directives), submits via `sbatch`, waits asynchronously for completion notifications, harvests outputs back, and records what it learned about the host for the next session — with no manual SSH from us.
- **Self-corrects across the submit→fail→fix→resubmit loop.** When a job failed (e.g. matplotlib RGB-tuple colors breaking a pandas categorical `.map()`), Claude read the traceback, diagnosed the root cause, patched the script, and resubmitted — recording the gotcha to durable memory so it never repeats.
- **Persists project knowledge across sessions** via a memory system (data paths, env locations, niche-definition decisions, host quirks) and a bilingual (EN master + KO review) documentation convention, so each session resumes with full context instead of re-discovering the setup.
- **Runs a plan → approval → step-tracked execution workflow**, keeping the human in the loop at decision points (niche-definition method, neighbor scope, clustering algorithm) via structured multiple-choice questions rather than open-ended back-and-forth.

## Session log

### 2026-07-11 — Session 4: idea-doc review, scope decision, two-track HPC execution

**What we asked Claude to do:** review a revised soft-matching pipeline draft, tell us what to fix, then decide the verifications/decisions and start the analysis.

**How Claude used its capabilities (beyond basic):**
- **Reviewed a design against its own memory of the project, not just the text.** Claude pulled prior state-matching results from memory (argmax=Rest 81.7%, mean Pearson values) and *read the actual per-cell CSV back* to confirm the numbers before critiquing — catching that the draft's headline per-cell projection would collapse to a plain Rest effect, and recommending the direction-flip analysis as the defensible primary instead.
- **Turned a "why do we need this?" question into a one-shot cluster verification.** Rather than argue the 912-vs-501 gene-count contradiction in prose, Claude wrote an h5py var-group reader, ran it on the login node (seconds, under the 4 GiB cap), and *froze both numbers from the real files* — 912 and 501, zero duplicate mappings — then corrected the two contradicting memory notes.
- **Ran two independent analyses as parallel Slurm jobs** (all-CD4 niche redefinition; Rest/Stim direction-flip on the 16 GB DE_stats), submitted together, harvested asynchronously.
- **Reused prior work through lineage, not re-authoring.** Claude extracted the frozen 7-program marker set and the 4-block niche code from an existing artifact's lineage and adapted only the anchor (naive → all-CD4), keeping the method identical across tracks.
- **Built the scope decision as a testable hypothesis and then tested it.** The all-CD4 expansion was justified by a claim ("naive-only makes exhaustion/Treg scores degenerate"); Claude added a non-degeneracy validation to the job and confirmed it (exhausted positivity 0→0.53).
- **Sanity-checked its own primary result against known biology** — flagged that TCR-signaling KOs flip the cytotoxic program Rest→Stim exactly as expected, as evidence the pipeline isn't spurious.
- **Self-recovered a job failure** (matplotlib `labels`→`tick_labels` API rename) with a one-line fix and resubmit, and — when it claimed an extra defensive guard it hadn't written — accepted the auditor's correction rather than back-fill.


**Same session, continued (relevance map + validation):**
- **Combined the two tracks into a niche × target relevance map** — niche baseline program burden × condition-specific perturbation effect — without re-deriving either input, then ranked candidate regulators per endpoint.
- **Retrieved real literature to check its own candidates.** After surfacing PTPN2, Claude queried OpenAlex and returned verifiable references (ABBV-CLS-484 2023; primary human T-cell CRISPR screens 2018) rather than asserting a citation from memory.
- **Caught and corrected its own fabricated citation.** Claude had attributed PTPN2 to a specific "Yeh 2024" without a lookup; when the auditor flagged it, Claude retracted the attribution, ran the actual retrieval, and logged the correction.
- **Designed and ran a circularity control unprompted.** Recognizing that the same marker genes defined both the niches and the ranking, Claude built a held-out confirmation (471 genes = shared minus program markers) with a matched-null percentile, ran it on Slurm, and reported the honest split: three endpoints confirm on independent genes, two do not — down-weighting its own earlier candidates rather than defending them.

### 2026-07-11 — Session 3: CD4 naive scope, cell-identity pivot, 4-block feature niches
**What we asked Claude to do:** narrow to CD4 naive / LUAD-5 R1-3; reconsider whether to link Perturb-seq at niche vs cell-identity level; then define a tumor spatial niche combining T-cell state + neighbor composition + suppressive program + tumor-boundary distance.

**How Claude used its capabilities (beyond basic):**
- **Grounded a methods debate in a live literature review.** When the human questioned whether niche-level was the right linkage unit, Claude ran targeted web searches and synthesized three established paradigms (ProjecTILs reference-projection, Yeh 2024 HGSC spatial+Perturb-seq, SpatialProp GNN), each with its data resources and its fit/limits for our 960-gene panel — then recommended a concrete pipeline (cell-identity classification → Yeh-style Perturb-seq ranking → spatial mapping) rather than just answering yes/no.
- **Bug diagnosis by isolation.** A state-matching job returned 0 shared genes; instead of guessing, Claude submitted a one-shot diagnostic that compared var_names and var columns across both AnnDatas, pinpointing the mismatch (gene symbols vs Ensembl IDs) and the fix (pseudobulk `gene_name` column → 912 genes).
- **Diagnostic-before-commit, again.** Before the feature-niche pipeline, Claude submitted a job to enumerate cell_type labels (to define which are malignant for the distance block) and to measure signature-marker coverage in the 960-gene panel — then used only panel-present markers.
- **Interpretable feature engineering with a defensible weighting.** Claude built a 4-block feature vector, and when the human asked whether the block weighting was principled, Claude answered honestly (it was implicit/uneven due to differing feature counts), explained the consequence (distance buried among 46 features), and re-ran with block-count (sqrt-n) normalization the human chose — after which the tumor-distance axis spanned 67-1121µm and the expected niche types (tumor-contacting exhausted, myeloid-rich Treg, fibroblast-rich excluded, TLS-active) emerged cleanly.
- **Publication-grade figure refinement.** On request, Claude reordered the niche heatmap by hierarchical clustering of both axes and ran a geometric overlap self-check before saving.
- **Honest status-keeping.** When asked whether the daily docs were being kept current, Claude said plainly they were not yet updated for this session, rather than implying otherwise — then updated them.

**Result:** 21 CD4 naive neighborhood niches (Phase 0) + 16 interpretable 4-block feature-niches with an explicit, principled block weighting; a literature-grounded decision to pivot to cell-identity linkage; state-matching bug root-caused with a known fix.

### 2026-07-11 — Session 2: T-cell-specific niche construction (neighborhood-expression clustering)
**What we asked Claude to do:** build T-cell-specific niches from scratch (not rely on the paper's CellCharter labels): per-T-cell Gaussian-kernel neighborhood expression → cluster → annotate → compare to the existing niche labels.

**How Claude used its capabilities (beyond basic):**
- **Diagnostic-before-commit.** Before writing the pipeline, Claude submitted a small Slurm diagnostic job to measure what it needed to parameterize the method correctly: matrix normalization state (found X is log-normalized, raw in `layers['counts']`), spatial-coordinate structure (samples occupy disjoint coordinate regions; FOVs tile within a sample → neighbors computed within-sample), and the median nearest-neighbor distance (46 µm) — which it used to set the Gaussian bandwidth σ=50 µm and a 3σ=150 µm cutoff **from the data**, not a guessed constant.
- **Human-in-the-loop design choices** surfaced as structured questions: neighbor scope (all cell types vs T-only) and clustering algorithm (Leiden vs fixed-k KMeans) — the human picked, Claude implemented exactly that.
- **Custom spatial method, not an off-the-shelf call.** Gaussian-kernel-weighted neighborhood expression over a KD-tree radius query, per sample, vectorized with sparse matrices; then scanpy PCA→kNN→Leiden on the neighborhood matrix; then a quantitative comparison to the reference CellCharter niche via ARI/NMI + a row-normalized confusion heatmap.
- **Honest interpretation.** Claude flagged an unprompted finding: the spatial maps show a clear per-sample/batch effect (some samples collapse to a single niche), correctly attributing it to clustering on non-batch-corrected neighborhood expression — a caveat, surfaced without being asked, that sets up the next methodological step.

**Result:** 20 T-cell niches (Leiden res=1.0); TLS-like (lymphoid structure) and immune niches cleanly resolved into multiple sub-niches; ARI=0.208 / NMI=0.338 vs CellCharter (moderate agreement — new niches subdivide the reference). Artifacts: heatmap, UMAP, per-sample spatial map, annotation + crosstab + per-cell assignment tables.

### 2026-07-11 — Session 1: env + atlas visualization + T-cell-niche flagging
**What we asked Claude to do:** create the spatial analysis environment on the cluster, visualize the CosMx NSCLC atlas by cell type, and define the T-cell niche.

**How Claude used its capabilities (beyond basic):**
- **Stood up a fresh conda environment on PACE via Slurm** (`env_spatial`: scanpy/squidpy/anndata) at a scratch prefix — reasoning on its own that HOME was ~77% full and redirecting the install to scratch — then verified imports and recorded the activation recipe to `compute_details` for reuse.
- **Inspected a 2.6 GB, 765k-cell AnnData without loading it into local memory** (backed read on the cluster) to recover the exact schema before writing any figure code.
- **Rendered an 8-panel single-cell spatial atlas** (765,771 cells, 22 cell types, shared legend) and computed the **T-cell fraction per niche** with enrichment vs global, flagging T-cell-rich niches (lymphoid structure 3.4×, immune 3.1×) — confirming the planned TLS positive control quantitatively.
- **Answered the human's conceptual questions** (what is a niche, what does resting vs stimulated mean, why Moran's I) grounded in the actual data schema, distinguishing what was *provided* in the data (CellCharter niche labels) from what Claude *computed* (the T-cell-rich flag).

**Result:** `env_spatial` ready; cell-type atlas + T-cell-fraction table + highlight map delivered as artifacts and mirrored to the server.

## Capabilities worth highlighting for judges
- Autonomous **Slurm orchestration** with async job harvesting and per-host learning.
- **Data-driven parameterization** (measuring NN distance to set a kernel bandwidth) rather than hardcoded constants.
- **Cross-session memory + bilingual documentation** that makes a multi-day project resumable.
- **Self-diagnosed failure recovery** with the fix written to memory.
- **Unprompted scientific caveats** (batch effect) that shape the next step, not just green checkmarks.
