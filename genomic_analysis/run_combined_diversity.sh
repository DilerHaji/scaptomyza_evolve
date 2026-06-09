#!/bin/bash
#SBATCH --job-name=combined_pi
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --time=6:00:00
#SBATCH --qos=savio_lowprio
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER_EMAIL}

cd PATH/TO/DATA/v5_dedup_trim2

. ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate PATH/TO/CONDA_ENVS/samtools

set -euo pipefail

REF=sfla_v2.fa
POSITIONS=baypass_wild/variant_positions.bed
OUTDIR=grenfst/diversity_combined
PILEUP_DIR=$OUTDIR/pileups

mkdir -p $PILEUP_DIR

# Step 1: Generate all 4 founder pileups in parallel
echo "=== Generating founder pileups (parallel) ==="
for SAMPLE in F1G00 F2G00 F3G00 F4G00; do
    BAM=../v5_dedup_trim/mapping/filtered/${SAMPLE}.bam
    OUT=$PILEUP_DIR/${SAMPLE}.mpileup
    if [ ! -s "$OUT" ]; then
        echo "  Starting $SAMPLE..."
        samtools mpileup -q 25 -Q 25 -l $POSITIONS -f $REF $BAM > $OUT &
    else
        echo "  $SAMPLE exists ($(wc -l < $OUT) lines), skipping"
    fi
done

echo "  Waiting for all pileups to finish..."
wait
echo "  All pileups done: $(date)"

for f in $PILEUP_DIR/F*G00.mpileup; do
    echo "  $(basename $f): $(wc -l < $f) lines"
done

# Step 2: Pool sizes — use grenedalf's naming convention (filename.1)
cat > $OUTDIR/combined_pool_sizes.txt << 'EOF'
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
EOF

# Step 3: Run grenedalf
echo ""
echo "=== Running grenedalf diversity (10 samples, matched sites) ==="
conda activate PATH/TO/CONDA_ENVS/grenedalf || true

grenedalf diversity \
  --pileup-path baypass_wild/pileups/AVB.mpileup \
  --pileup-path baypass_wild/pileups/AVT.mpileup \
  --pileup-path baypass_wild/pileups/PSB.mpileup \
  --pileup-path baypass_wild/pileups/PST.mpileup \
  --pileup-path baypass_wild/pileups/RMB.mpileup \
  --pileup-path baypass_wild/pileups/RMT.mpileup \
  --pileup-path $PILEUP_DIR/F1G00.mpileup \
  --pileup-path $PILEUP_DIR/F2G00.mpileup \
  --pileup-path $PILEUP_DIR/F3G00.mpileup \
  --pileup-path $PILEUP_DIR/F4G00.mpileup \
  --pool-sizes $OUTDIR/combined_pool_sizes.txt \
  --filter-sample-min-count 2 \
  --filter-sample-min-read-depth 4 \
  --filter-sample-max-read-depth 500 \
  --window-type interval \
  --window-interval-width 390000 \
  --window-interval-stride 390000 \
  --window-average-policy valid-snps \
  --allow-file-overwriting \
  --out-dir $OUTDIR \
  --file-prefix combined_pi_390000 \
  --verbose

echo ""
echo "=== Done: $(date) ==="
wc -l $OUTDIR/combined_pi_390000*.csv
head -2 $OUTDIR/combined_pi_390000*.csv
