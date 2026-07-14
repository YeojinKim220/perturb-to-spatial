# Running the steps on SLURM

Every step that takes more than ~1 minute was submitted as a SLURM batch job on
an HPC cluster (login-node use was limited to quick `ls` / light inspection).
Heavy outputs were written to **scratch**, never to a size-limited home
directory.

## CPU job template

Copy this, edit the env name / script path / resources, and `sbatch` it:

```bash
#!/usr/bin/env bash
#SBATCH --job-name=spatialperturb
#SBATCH --account=<your-account>
#SBATCH --partition=<cpu-partition>
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH -o log/%j.out

mkdir -p log
module load anaconda3        # adjust to your cluster's module
conda activate <env_name>    # e.g. a scratch-prefix env with scanpy + statsmodels
cd "$SLURM_SUBMIT_DIR"

python <step_folder>/<script>.py
```

## Notes

- **Write outputs to scratch.** Each script's `BASE` / `OUT` path points at a
  scratch location; keep it there — the AnnData files and intermediate CSVs are
  large.
- **Memory.** The spatial steps load a 765k-cell AnnData; 64 GB is comfortable.
  The `DE_stats` steps read via `h5py` slice-wise and need less.
- **One job per script.** The steps were run as individual jobs, not an array.
- **Reproducibility.** Every scored step sets `seed=0` (and permutation nulls use
  a fixed generator), so reruns match.
