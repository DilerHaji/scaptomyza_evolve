#!/bin/bash
#SBATCH --job-name=attrit_pi
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --time=8:00:00
#SBATCH --qos=savio_lowprio
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER_EMAIL}

# -----------------------------------------------------------------------------
# Diversity-attrition grenedalf run.
# Produces a single CSV with pi / theta_watterson / tajimas_d for:
#   6 wild pools     (AVB AVT PSB PST RMB RMT, n=50 haploid)
#   4 founder pools  (F1G00..F4G00, n=80)
#   12 G10 pools     (B1..B4 / T1..T4 / M1..M4, n=160)
# Same position bed, same 390 kb non-overlapping windows, same filters as
# the existing combined_pi_390000 run → rows line up window-for-window.
#
# Reuses existing wild + founder pileups via symlink; only G10 pileups are
# generated fresh.
# -----------------------------------------------------------------------------

cd PATH/TO/DATA/v5_dedup_trim2
. ~/.bashrc
eval "$(conda shell.bash hook)"

REF=sfla_v2.fa
POSITIONS=baypass_wild/variant_positions.bed
BAM_DIR=../v5_dedup_trim/mapping/filtered
OUTDIR=grenfst/diversity_attrition
PILEUP_DIR=$OUTDIR/pileups
mkdir -p "$PILEUP_DIR"

echo "=== attrition diversity ==="
echo "start: $(date)"

# ── Step 1: G10 pileups (parallel) ────────────────────────────────────────────
conda activate PATH/TO/CONDA_ENVS/samtools
set -euo pipefail

for SAMPLE in B1G10 B2G10 B3G10 B4G10 \
              T1G10 T2G10 T3G10 T4G10 \
              M1G10 M2G10 M3G10 M4G10; do
    BAM="$BAM_DIR/${SAMPLE}.bam"
    OUT="$PILEUP_DIR/${SAMPLE}.mpileup"
    if [ ! -s "$OUT" ]; then
        if [ ! -f "$BAM" ]; then
            echo "  WARNING: missing BAM for $SAMPLE ($BAM)"
            continue
        fi
        echo "  starting pileup: $SAMPLE"
        samtools mpileup -q 25 -Q 25 -l "$POSITIONS" -f "$REF" "$BAM" > "$OUT" &
    else
        echo "  exists, skipping: $SAMPLE ($(wc -l < $OUT) lines)"
    fi
done
echo "  waiting for parallel pileups..."
wait
echo "  G10 pileups done: $(date)"

for f in "$PILEUP_DIR"/?*G10.mpileup; do
    echo "    $(basename $f): $(wc -l < $f) lines"
done

# ── Step 2: symlink existing wild + founder pileups ───────────────────────────
WILD_PILEUP_DIR=PATH/TO/DATA/baypass_wild/pileups
FOUNDER_PILEUP_DIR=PATH/TO/DATA/grenfst/diversity_combined/pileups

# ── Step 3: pool-size file (grenedalf filename.1 convention) ──────────────────
cat > "$OUTDIR/attrition_pool_sizes.txt" << 'EOF'
AVB.1,50
AVT.1,50
PSB.1,50
PST.1,50
RMB.1,50
RMT.1,50
F1G00.1,80
F2G00.1,80
F3G00.1,80
F4G00.1,80
B1G10.1,160
B2G10.1,160
B3G10.1,160
B4G10.1,160
T1G10.1,160
T2G10.1,160
T3G10.1,160
T4G10.1,160
M1G10.1,160
M2G10.1,160
M3G10.1,160
M4G10.1,160
EOF

# ── Step 4: grenedalf diversity on all 22 pools ───────────────────────────────
conda activate PATH/TO/CONDA_ENVS/grenedalf

echo ""
echo "=== running grenedalf diversity on 22 pools ==="

grenedalf diversity \
  --pileup-path "$WILD_PILEUP_DIR/AVB.mpileup" \
  --pileup-path "$WILD_PILEUP_DIR/AVT.mpileup" \
  --pileup-path "$WILD_PILEUP_DIR/PSB.mpileup" \
  --pileup-path "$WILD_PILEUP_DIR/PST.mpileup" \
  --pileup-path "$WILD_PILEUP_DIR/RMB.mpileup" \
  --pileup-path "$WILD_PILEUP_DIR/RMT.mpileup" \
  --pileup-path "$FOUNDER_PILEUP_DIR/F1G00.mpileup" \
  --pileup-path "$FOUNDER_PILEUP_DIR/F2G00.mpileup" \
  --pileup-path "$FOUNDER_PILEUP_DIR/F3G00.mpileup" \
  --pileup-path "$FOUNDER_PILEUP_DIR/F4G00.mpileup" \
  --pileup-path "$PILEUP_DIR/B1G10.mpileup" \
  --pileup-path "$PILEUP_DIR/B2G10.mpileup" \
  --pileup-path "$PILEUP_DIR/B3G10.mpileup" \
  --pileup-path "$PILEUP_DIR/B4G10.mpileup" \
  --pileup-path "$PILEUP_DIR/T1G10.mpileup" \
  --pileup-path "$PILEUP_DIR/T2G10.mpileup" \
  --pileup-path "$PILEUP_DIR/T3G10.mpileup" \
  --pileup-path "$PILEUP_DIR/T4G10.mpileup" \
  --pileup-path "$PILEUP_DIR/M1G10.mpileup" \
  --pileup-path "$PILEUP_DIR/M2G10.mpileup" \
  --pileup-path "$PILEUP_DIR/M3G10.mpileup" \
  --pileup-path "$PILEUP_DIR/M4G10.mpileup" \
  --pool-sizes "$OUTDIR/attrition_pool_sizes.txt" \
  --filter-sample-min-count 2 \
  --filter-sample-min-read-depth 4 \
  --filter-sample-max-read-depth 500 \
  --window-type interval \
  --window-interval-width 390000 \
  --window-interval-stride 390000 \
  --window-average-policy valid-snps \
  --allow-file-overwriting \
  --out-dir "$OUTDIR" \
  --file-prefix attrition_pi_390000 \
  --verbose

echo ""
echo "=== done: $(date) ==="
wc -l "$OUTDIR"/attrition_pi_390000*.csv
head -2 "$OUTDIR"/attrition_pi_390000*.csv
