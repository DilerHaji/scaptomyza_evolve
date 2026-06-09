#!/bin/bash
#SBATCH --job-name=bp_wild
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --qos=rosalind_htc4_normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER_EMAIL}

export PYTHONUNBUFFERED=1

. ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate PATH/TO/CONDA_ENVS/snakemake4

cd PATH/TO/DATA/v5_dedup_trim2

snakemake --unlock || true

mkdir -p slurm logs benchmarks

snakemake all_baypass_wild \
  --cluster "sbatch {resources.resources}" \
  --keep-going \
  --rerun-incomplete \
  --use-conda \
  --conda-frontend conda \
  --nolock \
  --latency-wait 60 \
  --rerun-triggers mtime \
  -j 100 \
  --scheduler greedy
