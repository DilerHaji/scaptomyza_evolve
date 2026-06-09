#!/bin/bash
#SBATCH --job-name=block_r2
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=ngsld/block_validation/slurm-%j.out
#SBATCH --error=ngsld/block_validation/slurm-%j.err

set -euo pipefail

OUT_DIR="ngsld/block_validation"
mkdir -p "$OUT_DIR"

BLOCKS_BED="$OUT_DIR/blocks.bed"
if [ ! -f "$BLOCKS_BED" ]; then
    echo "ERROR: $BLOCKS_BED not found."
    exit 1
fi
echo "Blocks: $(wc -l < "$BLOCKS_BED")"

SCAFFOLDS=(
  "chr_ScDA7r2_110_HRSCAF_295"
  "chr_ScDA7r2_126_HRSCAF_325"
  "chr_ScDA7r2_439_HRSCAF_779"
  "chr_ScDA7r2_597_HRSCAF_953"
)

process_one() {
    local cohort=$1
    local scaff=$2
    local ld_file="ngsld/ld/${cohort}/${scaff}.ld.gz"
    local shared_file="ngsld/shared_snp_analysis/${scaff}.shared.pos"
    local out_tsv="$OUT_DIR/blockr2_${cohort}_${scaff}.tsv"

    if [ ! -f "$ld_file" ] || [ ! -f "$shared_file" ]; then
        echo "[skip] $cohort $scaff — missing inputs"
        return
    fi

    echo "[start] $cohort $scaff — $(date +%T)"

    zcat "$ld_file" | awk -v LABEL="$cohort" -v SCAFF="$scaff" \
                          -v SHARED_FILE="$shared_file" \
                          -v BLOCKS_BED="$BLOCKS_BED" \
        'BEGIN {
            while ((getline line < SHARED_FILE) > 0) shared[line] = 1
            close(SHARED_FILE)

            n_blocks = 0
            while ((getline line < BLOCKS_BED) > 0) {
                split(line, a, "\t")
                if (a[1] != SCAFF) continue
                n_blocks++
                b_start[n_blocks] = a[2] + 0
                b_end[n_blocks]   = a[3] + 0
                b_tag[n_blocks]   = a[4]
                span = a[3] - a[2]
                f_start[n_blocks] = a[3] + 1
                f_end[n_blocks]   = a[3] + span
            }
            close(BLOCKS_BED)
        }
        NR == 1 { next }
        {
            split($1, a, ":"); p1 = a[2] + 0
            split($2, b, ":"); p2 = b[2] + 0
            if (!(p1 in shared) || !(p2 in shared)) next

            dist = $3 + 0
            r2 = $4 + 0

            if (dist < 1000)        bin = "1_<1kb"
            else if (dist < 10000)  bin = "2_1-10kb"
            else if (dist < 100000) bin = "3_10-100kb"
            else                    next

            for (i = 1; i <= n_blocks; i++) {
                if (p1 >= b_start[i] && p1 <= b_end[i] && \
                    p2 >= b_start[i] && p2 <= b_end[i]) {
                    key = LABEL "\t" b_tag[i] "\tinside\t" bin
                    n[key]++; sumr[key] += r2
                    break
                }
                if (p1 >= f_start[i] && p1 <= f_end[i] && \
                    p2 >= f_start[i] && p2 <= f_end[i]) {
                    key = LABEL "\t" b_tag[i] "\tflanking\t" bin
                    n[key]++; sumr[key] += r2
                    break
                }
            }
        }
        END {
            for (k in n) printf "%s\t%d\t%.6f\n", k, n[k], sumr[k] / n[k]
        }' > "$out_tsv"
}

for SCAFF in "${SCAFFOLDS[@]}"; do
    process_one founders "$SCAFF" &
    process_one T2G07    "$SCAFF" &
done
wait

SUMMARY="$OUT_DIR/block_r2_summary.tsv"
echo -e "cohort\tblock_tag\tregion\tdist_bin\tn_pairs\tmean_r2" > "$SUMMARY"
cat "$OUT_DIR"/blockr2_*.tsv >> "$SUMMARY"
wc -l "$SUMMARY"
head -5 "$SUMMARY"
