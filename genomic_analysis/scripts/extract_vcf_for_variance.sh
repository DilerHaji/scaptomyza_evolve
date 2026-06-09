#!/bin/bash
#SBATCH --job-name=extract_vcf
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1


set -eo pipefail

export PS1="${PS1:-}"
eval "$(conda shell.bash hook)"
conda activate PATH/TO/CONDA_ENVS/bcftools

VCF="fvariants/fvOG_e10fe9w.fixed_no_neff.vcf.gz"
OUTDIR="variance_analysis"
mkdir -p "$OUTDIR"

bcftools query -l "$VCF" > "$OUTDIR/sample_list.txt"

bcftools query -f '%CHROM\t%POS\n' "$VCF" > "$OUTDIR/sites.tsv"

awk -v OFS='\t' '{print $1, $2-1, $2}' "$OUTDIR/sites.tsv" > "$OUTDIR/sites.bed"

bcftools query -f '%CHROM\t%POS[\t%DP]\n' "$VCF" > "$OUTDIR/merged_depth.tsv"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "$VCF" > "$OUTDIR/merged_ad.tsv"

ls -lh "$OUTDIR"
