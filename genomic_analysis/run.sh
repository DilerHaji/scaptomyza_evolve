#!/bin/bash
#SBATCH --job-name=v5
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=50:00:00
#SBATCH --qos=savio_normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER_EMAIL}
#SBATCH --exclude=n0032.savio4

export PYTHONUNBUFFERED=1

. ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate PATH/TO/CONDA_ENVS/snakemake4
snakemake --unlock || true

rm logs -r
rm benchmarks -r
rm slurm -r
mkdir slurm

#snakemake -j 1000 --touch

snakemake all_timeslice_multiplots \
--cluster "sbatch {resources.resources}" \
--keep-going \
--rerun-incomplete \
--use-conda \
--conda-frontend conda \
--nolock \
--latency-wait 60 \
-j 2000 \
--rerun-triggers mtime

#--use-singularity \
#--singularity-args "--nv --disable-cache" \
#--conda-prefix PATH/TO/DATA/popgen/.snakemake/conda


#snakemake --cluster "sbatch {resources.resources}" -j 600 --use-conda --rerun-incomplete --use-singularity --singularity-args "--nv --disable-cache" --latency-wait 10
#snakemake --cluster "sbatch {resources.resources}" -j 1000 --use-conda --rerun-incomplete --rerun-triggers mtime --use-singularity --singularity-args "--nv --disable-cache" --latency-wait 60
#snakemake --cluster "sbatch {resources.resources}" -j 1000 --keep-going --use-conda --rerun-incomplete --rerun-triggers mtime --use-singularity --singularity-args "--nv --disable-cache" --latency-wait 60
#snakemake --cluster "sbatch {resources.resources}" -j 11 --use-conda --rerun-incomplete --use-singularity --singularity-args "--nv --disable-cache"
