# Step 6 — Direction-flip screen & first niche × target relevance map (exploratory)

This step produced the project's first integrated deliverable: a per-endpoint
target ranking combining niche program burden (step 4) with a state-resolved
perturbation effect. **It has since been superseded** by the statistically
corrected re-analysis in [`step7`](../step7_phase2_corrections/) and
[`step8`](../step8_phase3_reanalysis/), which apply proper permutation nulls and
multiple-testing correction and reach a more conservative conclusion. It is kept
here because it is where the direction-flip idea and the relevance-map structure
were defined.

> **Code note.** This step ran as inline job scripts in an earlier session rather
> than as committed standalone `.py` files, so this folder holds the method
> description and the exact recomputable definitions rather than a script. The
> recomputed, corrected versions of these quantities live in step 7
> (`s1_fdr_audit.py`, `phase2_analysis.py` stage S5) and are the ones to run.

## What it computed

**Perturbation program effect.** From the `DE_stats` `log_fc` layer, for each
QC-passing perturbation `g`, condition `c ∈ {Rest, Stim8hr, Stim48hr}`, and
program `p`:

```
E(g, c, p) = mean log_fc over program p's marker genes,  excluding the target gene g itself
```

QC = on-target knockdown significant, no distal off-target flag, adequate target
expression.

**Direction-flip classification.** For each (target, program), compare `E` at
Rest vs Stim48hr and label: `direction_flip` (sign reversal, both |E| ≥ 0.10),
`rest_specific`, `stim_specific`, or `stable`.

**Niche × target relevance.**
`relevance(niche r, target g, endpoint) = max(0, niche program burden) ×
max(0, favorable-direction · E(g, Rest, p))`, ranked per endpoint, with each
endpoint anchored to its highest-burden niche.

**Held-out confirmation (circularity control).** Because the same program markers
define both the niches and the ranking, top candidates were re-scored on the 471
held-out genes (501 shared minus the 30 program markers) against a matched null,
and the null-percentile reported.

## Original (exploratory) findings — read step 7/8 for the corrected result

- ~5,100 (target × program) pairs flagged as Rest↔Stim48 direction-flips,
  concentrated in Tfh > cytotoxic > Treg > exhausted.
- Positive control held: TCR-signaling knockouts (ZAP70, CD3E, VAV1, PTPRC) flip
  the cytotoxic program from negative at rest to positive when stimulated.
- The relevance map ranked candidates across five endpoints (exhaustion reversal,
  Treg reduction, Th1 enhancement, activation enhancement, cytotoxic enhancement).

**Why superseded:** the held-out check used a percentile against a pooled null
rather than a per-target permutation null with FDR. When step 7 replaced that
with a proper per-target permutation null and BH-FDR, most of the apparent hits
did not survive — see step 7.

## Outputs (archived)

`condition_direction_flips.csv`, `perturbation_program_effects_long.csv`,
`spatialperturb_tracks_summary.png`, `niche_target_relevance_map.png`,
`niche_target_relevance_top10.csv`, `niche_target_relevance_long.csv`,
`heldout_confirmation.png`, `heldout_gene_confirmation.csv`,
`candidate_literature_evidence.csv`.
