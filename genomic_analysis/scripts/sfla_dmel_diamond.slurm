#!/bin/bash
#SBATCH --job-name=sfla_dmel_diamond
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/sfla_dmel_diamond_%j.out
#SBATCH --error=logs/sfla_dmel_diamond_%j.err

set -euo pipefail

WORKDIR=PATH/TO/DATA/synteny_xenrichment
cd "$WORKDIR"
mkdir -p logs

QUERY=sfla_v2_proteome.faa
DB=dmel.dmnd
OUT=sfla_v2_dmel_blastp.tsv

module load gcc/11.4.0 || true

if ! command -v diamond >/dev/null 2>&1; then
    source PATH/TO/DATA/.conda/etc/profile.d/conda.sh
    conda activate snakemake8
fi

diamond version

diamond blastp \
    --query    "$QUERY" \
    --db       "$DB" \
    --out      "$OUT" \
    --outfmt   6 qseqid sseqid stitle pident length evalue bitscore \
    --max-target-seqs 1 \
    --evalue   1e-5 \
    --threads  ${SLURM_CPUS_PER_TASK:-8} \
    --more-sensitive

wc -l "$OUT"
echo "done"
