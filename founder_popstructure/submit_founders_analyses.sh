#!/bin/bash
#SBATCH --job-name=founders_analyses
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --qos=savio_normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER_EMAIL}
#SBATCH --output=slurm/%j_founders_analyses.out
#SBATCH --error=slurm/%j_founders_analyses.err

find slurm/ -name "*.out" -o -name "*.err" | sort | head -n -5 | xargs -r rm -f
find .snakemake/slurm_logs/ -name "*.log" -mtime +3 -delete 2>/dev/null || true

. ~/.bashrc
eval "$(conda shell.bash hook)"
conda config --set channel_priority strict
conda activate PATH/TO/CONDA_ENVS/snakemake8
snakemake --unlock -s Snakefile_founders_analyses
CONDA_EXE=PATH/TO/CONDA_ENVS/snakemake8/bin/conda

snakemake \
    -s Snakefile_founders_analyses \
    -j 200 \
    --executor slurm \
    --profile PATH/TO/DATA/profiles/slurm \
    --use-conda \
    --conda-frontend conda \
    --rerun-incomplete
