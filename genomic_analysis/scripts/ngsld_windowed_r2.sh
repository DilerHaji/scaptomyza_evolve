#!/bin/bash
#SBATCH --job-name=winr2
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=ngsld/block_validation/slurm-winr2-%j.out
#SBATCH --error=ngsld/block_validation/slurm-winr2-%j.err

set -euo pipefail

SCAFFOLDS=(
  "chr_ScDA7r2_110_HRSCAF_295"
  "chr_ScDA7r2_126_HRSCAF_325"
  "chr_ScDA7r2_439_HRSCAF_779"
  "chr_ScDA7r2_597_HRSCAF_953"
)

OUT_DIR="ngsld/block_validation"
WINSIZE=500000
SUMMARY="$OUT_DIR/windowed_r2.tsv"
echo -e "cohort\tscaffold\twin_start\twin_end\tn_pairs\tmean_r2" > "$SUMMARY"

process_one() {
    local cohort=$1
    local scaff=$2
    local ld_file="ngsld/ld/${cohort}/${scaff}.ld.gz"
    local shared_file="ngsld/shared_snp_analysis/${scaff}.shared.pos"
    local out_tsv="$OUT_DIR/winr2_${cohort}_${scaff}.tsv"

    if [ ! -f "$ld_file" ] || [ ! -f "$shared_file" ]; then
        echo "[skip] $cohort $scaff"
        return
    fi
    echo "[start] $cohort $scaff — $(date +%T)"

    zcat "$ld_file" | awk -v LABEL="$cohort" -v SCAFF="$scaff" \
                          -v WIN="$WINSIZE" \
                          -v SHARED_FILE="$shared_file" \
        'BEGIN {
            while ((getline line < SHARED_FILE) > 0) shared[line] = 1
            close(SHARED_FILE)
        }
        NR == 1 { next }
        {
            split($1, a, ":"); p1 = a[2] + 0
            split($2, b, ":"); p2 = b[2] + 0
            if (!(p1 in shared) || !(p2 in shared)) next
            dist = $3 + 0
            if (dist < 10000 || dist >= 100000) next
            r2 = $4 + 0
            # Window based on midpoint of the pair
            mid = int((p1 + p2) / 2)
            w = int(mid / WIN) * WIN
            key = LABEL "\t" SCAFF "\t" w "\t" (w + WIN)
            n[key]++
            sumr[key] += r2
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

cat "$OUT_DIR"/winr2_*.tsv >> "$SUMMARY"
wc -l "$SUMMARY"
head -5 "$SUMMARY"
