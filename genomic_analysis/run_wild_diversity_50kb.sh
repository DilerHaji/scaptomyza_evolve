#!/bin/bash
#SBATCH --job-name=wild_pi_50kb
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --qos=savio_lowprio
#SBATCH --output=logs/wild_pi_50kb_%j.out
#SBATCH --error=logs/wild_pi_50kb_%j.err

# -----------------------------------------------------------------------------
# Wild-pool diversity at 50 kb non-overlapping windows for the chr_439
# balancing-selection scan. Same filter and window-policy settings as the
# existing 390-kb attrition_pi run -- only window size differs. Reuses the
# existing wild pileups (symlinked from baypass_wild/pileups).
#
# Output:  grenfst/diversity_attrition/wild_pi_50000.csv
# -----------------------------------------------------------------------------

cd PATH/TO/DATA/v5_dedup_trim2

# Source bashrc/conda BEFORE enabling strict mode -- /etc/bashrc references
# PS1 which is unset in non-interactive shells and trips `set -u`.
. ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate PATH/TO/CONDA_ENVS/grenedalf

set -eo pipefail
mkdir -p logs

WILD_PILEUP_DIR=PATH/TO/DATA/baypass_wild/pileups
OUTDIR=grenfst/diversity_attrition

# Pool-size file -- 6 wild pools, n=50 each
cat > "$OUTDIR/wild_pool_sizes.txt" << 'EOF'
AVB.1,100
AVT.1,100
PSB.1,100
PST.1,100
RMB.1,100
RMT.1,100
EOF

echo "=== wild diversity 50kb start: $(date) ==="
grenedalf diversity \
  --pileup-path "$WILD_PILEUP_DIR/AVB.mpileup" \
  --pileup-path "$WILD_PILEUP_DIR/AVT.mpileup" \
  --pileup-path "$WILD_PILEUP_DIR/PSB.mpileup" \
  --pileup-path "$WILD_PILEUP_DIR/PST.mpileup" \
  --pileup-path "$WILD_PILEUP_DIR/RMB.mpileup" \
  --pileup-path "$WILD_PILEUP_DIR/RMT.mpileup" \
  --pool-sizes "$OUTDIR/wild_pool_sizes.txt" \
  --filter-sample-min-count 2 \
  --filter-sample-min-read-depth 4 \
  --filter-sample-max-read-depth 500 \
  --window-type interval \
  --window-interval-width 50000 \
  --window-interval-stride 50000 \
  --window-average-policy valid-snps \
  --allow-file-overwriting \
  --out-dir "$OUTDIR" \
  --file-prefix wild_pi_50000_n100_ \
  --verbose

echo ""
echo "=== 50kb done: $(date) -- starting single-SNP scan ==="

# Per-SNP scan: --window-type single produces one row per variant site.
# π is well-defined per site; θW per site reduces to 1/a_{n-1} per segregating
# site (i.e., a constant), and Tajima's D per site is undefined / NaN. We
# keep all three columns for completeness and only use π in downstream code.
grenedalf diversity \
  --pileup-path "$WILD_PILEUP_DIR/AVB.mpileup" \
  --pileup-path "$WILD_PILEUP_DIR/AVT.mpileup" \
  --pileup-path "$WILD_PILEUP_DIR/PSB.mpileup" \
  --pileup-path "$WILD_PILEUP_DIR/PST.mpileup" \
  --pileup-path "$WILD_PILEUP_DIR/RMB.mpileup" \
  --pileup-path "$WILD_PILEUP_DIR/RMT.mpileup" \
  --pool-sizes "$OUTDIR/wild_pool_sizes.txt" \
  --filter-sample-min-count 2 \
  --filter-sample-min-read-depth 4 \
  --filter-sample-max-read-depth 500 \
  --window-type single \
  --window-average-policy valid-snps \
  --allow-file-overwriting \
  --out-dir "$OUTDIR" \
  --file-prefix wild_pi_persnp_n100_ \
  --verbose

echo ""
echo "=== all done: $(date) ==="
wc -l "$OUTDIR"/wild_pi_50000_n100_*.csv "$OUTDIR"/wild_pi_persnp_n100_*.csv || true
