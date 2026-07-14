# Step 1 — Data inspection & gene-space reconciliation

Before any analysis, confirm the schema of both AnnData files and freeze the
shared gene sets. This step exists because a naive gene intersection between the
two platforms is **empty** (CosMx uses HGNC symbols, Perturb-seq uses Ensembl
IDs) — a mismatch that otherwise surfaces only as blank downstream plots.

## Scripts

| Script | What it does |
|---|---|
| [`peek_var.py`](peek_var.py) | Prints the root keys, `.var`/`.obs` group structure, and AnnData encoding attributes of both files. Use it to confirm where gene symbols live before writing any loader. |
| [`overlap.py`](overlap.py) | Reads both files with `h5py`, decodes categoricals, and reports: DE measured-gene count, DE culture conditions, number of unique perturbed genes, CosMx panel size, cell count, cell types, niches, and samples. |

## How to run

Both are lightweight `h5py` reads (no full matrix load) — safe on a login node:

```bash
python peek_var.py
python overlap.py
```

## Key facts frozen here

- **Perturb-seq `DE_stats`:** 10,282 measured genes; conditions `Rest` /
  `Stim8hr` / `Stim48hr`; symbols in `.var['gene_name']`, Ensembl IDs in
  `.var['gene_ids']` / `_index`.
- **CosMx panel:** 960 genes, 765,771 cells, 22 cell types, 9 niches, 8 samples.
- **Shared-gene overlaps** (kept separate throughout the project):
  - **912** genes = CosMx ∩ pseudobulk symbols → state matching.
  - **501** genes = CosMx ∩ `DE_stats` symbols (⊂ 912) → gene-level / held-out
    confirmation.
  - Zero duplicate symbol→Ensembl mappings, so no gene-merging is needed.

## Output

Console output only, plus the frozen shared-gene lists used downstream
(`shared_gene_map.csv` with a `in_genelevel_501` flag; `shared_genes_statematch.txt`,
912 entries; `shared_genes_genelevel.txt`, 501 entries).
