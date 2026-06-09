#!/bin/bash
#SBATCH --job-name=traj_div_driver
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --qos=savio_lowprio
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER_EMAIL}
#SBATCH --output=slurm/%j.out
#SBATCH --error=slurm/%j.err


cd PATH/TO/DATA/v5_dedup_trim2

. ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate PATH/TO/CONDA_ENVS/snakemake4

set -eo pipefail

mkdir -p slurm

rm -rf grenfst/diversity_bins/trajectory/390000 \
       grenfst/diversity_bins/trajectory/200000 \
       grenfst/diversity/trajectory_pi_390000.csv \
       grenfst/diversity/trajectory_pi_200000.csv

snakemake -s Snakefile_diversity \
    --cluster "sbatch {resources.resources}" \
    --use-conda --conda-frontend conda \
    --rerun-incomplete --rerun-triggers mtime \
    --keep-going --latency-wait 60 \
    --nolock \
    -j 8 \
    all_trajectory_diversity

ls -lh grenfst/diversity/trajectory_pi_200000.csv
