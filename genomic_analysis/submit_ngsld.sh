#!/bin/bash
#SBATCH --job-name=ngsld_coord
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --qos=rosalind_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_ngsld_coord_%j.out
#SBATCH --error=slurm_ngsld_coord_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER_EMAIL}


cd PATH/TO/DATA/IND

. ~/.bashrc
conda activate snakemake8

snakemake -s Snakefile_ngsld \
    --configfile config/ngsld_config.yaml \
    --use-conda \
    --rerun-incomplete \
    -j 20 \
    --executor slurm \
    --default-resources slurm_account=${SLURM_ACCOUNT} slurm_partition=savio4_htc \
                        mem_mb=32000 runtime=480 \
    --keep-going
