#!/bin/bash
#SBATCH --job-name=founder_pi
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --qos=rosalind_htc4_normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER_EMAIL}

export PYTHONUNBUFFERED=1
. ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate PATH/TO/CONDA_ENVS/snakemake4

cd PATH/TO/DATA/v5_dedup_trim2

snakemake all_founder_diversity \
  --cluster "sbatch {resources.resources}" \
  --keep-going \
  --rerun-incomplete \
  --use-conda \
  --conda-frontend conda \
  -j 20
