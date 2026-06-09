#!/bin/bash
#SBATCH --job-name=pool_var_af
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --qos=savio_lowprio
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER_EMAIL}

. ~/.bashrc
eval "$(conda shell.bash hook)"
conda config --set channel_priority strict
conda activate PATH/TO/CONDA_ENVS/snakemake8_new
CONDA_EXE=PATH/TO/CONDA_ENVS/snakemake8_new/bin/conda

snakemake -s Snakefile_pool_comparison \
    --configfile config_pool_comparison.yml \
    -j 200 \
    --executor slurm \
    --profile profiles/slurm \
    --use-conda --conda-frontend conda \
    --rerun-incomplete
