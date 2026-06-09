#!/bin/bash
#SBATCH --job-name=ngsld_shared
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --qos=rosalind_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=ngsld/shared_snp_analysis/slurm-%j.out
#SBATCH --error=ngsld/shared_snp_analysis/slurm-%j.err

set -euo pipefail

SCAFFOLDS=(
  "chr_ScDA7r2_110_HRSCAF_295"
  "chr_ScDA7r2_126_HRSCAF_325"
  "chr_ScDA7r2_439_HRSCAF_779"
  "chr_ScDA7r2_597_HRSCAF_953"
)

HOTSPOT_CHR="chr_ScDA7r2_439_HRSCAF_779"
HOTSPOT_START=2000000
HOTSPOT_END=12000000

OUT_DIR="ngsld/shared_snp_analysis"
mkdir -p "$OUT_DIR"

for SCAFF in "${SCAFFOLDS[@]}"; do
    POS_FND="ngsld/angsd/founders/${SCAFF}.pos"
    POS_T2="ngsld/angsd/T2G07/${SCAFF}.pos"
    SHARED="$OUT_DIR/${SCAFF}.shared.pos"

    if [ ! -f "$POS_FND" ] || [ ! -f "$POS_T2" ]; then
        echo "[WARN] POS missing for $SCAFF — skipping"
        continue
    fi

    if [ ! -f "$SHARED" ]; then
        echo "[intersect] $SCAFF"
        awk '$2 ~ /^[0-9]+$/ {print $2}' "$POS_FND" | sort > "$OUT_DIR/tmp_fnd_${SCAFF}.txt"
        awk '$2 ~ /^[0-9]+$/ {print $2}' "$POS_T2"  | sort > "$OUT_DIR/tmp_t2_${SCAFF}.txt"
        comm -12 "$OUT_DIR/tmp_fnd_${SCAFF}.txt" "$OUT_DIR/tmp_t2_${SCAFF}.txt" > "$SHARED"
        rm -f "$OUT_DIR/tmp_fnd_${SCAFF}.txt" "$OUT_DIR/tmp_t2_${SCAFF}.txt"
        wc -l "$SHARED"
    else
        echo "[intersect] $SCAFF already done ($(wc -l < "$SHARED") SNPs)"
    fi
done

process_one() {
    local cohort=$1
    local scaff=$2
    local ld_file="ngsld/ld/${cohort}/${scaff}.ld.gz"
    local shared_file="$OUT_DIR/${scaff}.shared.pos"
    local out_tsv="$OUT_DIR/partial_${cohort}_${scaff}.tsv"

    if [ ! -f "$ld_file" ] || [ ! -f "$shared_file" ]; then
        echo "[skip] $cohort $scaff — missing inputs"
        return
    fi

    echo "[start] $cohort $scaff — $(date +%T)"
    zcat "$ld_file" | awk -v LABEL="$cohort" -v SCAFF="$scaff" \
                          -v HSTART="$HOTSPOT_START" -v HEND="$HOTSPOT_END" \
                          -v HCHR="$HOTSPOT_CHR" \
                          -v SHARED_FILE="$shared_file" \
        'BEGIN {
            while ((getline line < SHARED_FILE) > 0) shared[line] = 1
            close(SHARED_FILE)
        }
        NR == 1 { next }
        {
            split($1, a, ":"); pos1 = a[2] + 0
            split($2, b, ":"); pos2 = b[2] + 0
            if (!(pos1 in shared) || !(pos2 in shared)) next

            dist = $3 + 0
            r2 = $4 + 0

            if (SCAFF == HCHR) {
                in_hot = (pos1 >= HSTART && pos1 <= HEND && \
                          pos2 >= HSTART && pos2 <= HEND) ? 1 : 0
            } else {
                in_hot = 0
            }

            if (dist < 1000)        bin = "1_<1kb"
            else if (dist < 10000)  bin = "2_1-10kb"
            else if (dist < 100000) bin = "3_10-100kb"
            else                    next

            key = LABEL "\t" SCAFF "\t" in_hot "\t" bin
            n[key]++
            sumr[key] += r2
        }
        END {
            for (k in n) printf "%s\t%d\t%.6f\n", k, n[k], sumr[k] / n[k]
        }' > "$out_tsv"
    echo "[done ] $cohort $scaff — $(date +%T) — $(wc -l < "$out_tsv") lines"
}

for SCAFF in "${SCAFFOLDS[@]}"; do
    process_one founders "$SCAFF" &
    process_one T2G07    "$SCAFF" &
done
wait

SUMMARY="$OUT_DIR/summary_shared_snps.tsv"
echo -e "cohort\tscaffold\tin_hotspot\tdist_bin\tn_pairs\tmean_r2" > "$SUMMARY"
cat "$OUT_DIR"/partial_*.tsv >> "$SUMMARY"

echo ""
column -t -s $'\t' "$SUMMARY"

echo ""
awk -F'\t' 'NR > 1 && $2 == "'"$HOTSPOT_CHR"'" && $3 == "1"' "$SUMMARY" \
  | awk -F'\t' '
    {
        r[$4, $1] = $6
        n[$4, $1] = $5
        bins[$4] = 1
    }
    END {
        printf "%-12s  %12s  %12s  %10s  %10s\n", "dist_bin", "founders", "T2G07", "delta", "pct"
        for (b in bins) {
            f = r[b, "founders"] + 0
            t = r[b, "T2G07"]    + 0
            d = t - f
            pct = (f != 0) ? 100 * d / f : 0
            printf "%-12s  %12.4f  %12.4f  %+10.4f  %+9.1f%%\n", b, f, t, d, pct
        }
    }' | sort

echo ""
awk -F'\t' 'NR > 1 && $2 == "'"$HOTSPOT_CHR"'"' "$SUMMARY" \
  | awk -F'\t' '
    {
        r[$1, $3, $4] = $6
    }
    END {
        printf "%-10s  %-10s  %12s  %12s  %10s\n",
               "cohort", "dist_bin", "hotspot", "flanking", "ratio"
        for (k in r) {
            split(k, a, SUBSEP)
            if (a[2] == "1") {
                h = r[a[1], "1", a[3]] + 0
                f = r[a[1], "0", a[3]] + 0
                ratio = (f != 0) ? h / f : 0
                printf "%-10s  %-10s  %12.4f  %12.4f  %10.3f\n",
                       a[1], a[3], h, f, ratio
            }
        }
    }' | sort
