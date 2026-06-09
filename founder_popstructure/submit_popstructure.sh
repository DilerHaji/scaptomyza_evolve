#!/bin/bash
#SBATCH --job-name=popstructure
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --qos=savio_normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER_EMAIL}
#SBATCH --output=slurm/%j_popstructure.out
#SBATCH --error=slurm/%j_popstructure.err

find slurm/ -name "*.out" -o -name "*.err" | sort | head -n -5 | xargs -r rm -f
find .snakemake/slurm_logs/ -name "*.log" -mtime +3 -delete 2>/dev/null || true

. ~/.bashrc
eval "$(conda shell.bash hook)"
conda config --set channel_priority strict
conda activate PATH/TO/CONDA_ENVS/snakemake8


snakemake --unlock -s Snakefile_popstructure --configfile config/popstructure.yaml

CONDA_EXE=PATH/TO/CONDA_ENVS/snakemake8/bin/conda

snakemake \
    -s Snakefile_popstructure \
    --configfile config/popstructure.yaml \
    -j 50 \
    --executor slurm \
    --profile PATH/TO/DATA/profiles/slurm \
    --use-conda \
    --conda-frontend conda \
    --rerun-incomplete
