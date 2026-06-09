#!/bin/bash

set -e

WORK_DIR="PATH/TO/DATA/v5_dedup_trim2"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

cd "$WORK_DIR"

# Create backup
mkdir -p "backups/optimization_${TIMESTAMP}"
cp rules/lrt_genome_wide.smk "backups/optimization_${TIMESTAMP}/"

sed -i "s/^LRT_NSIMS = 1000/LRT_NSIMS = 200  # OPTIMIZED/" rules/lrt_genome_wide.smk


sed -i 's/threads: 10$/threads: 20  # OPTIMIZED/' rules/lrt_genome_wide.smk


sed -i 's/min_coverage = 10,.*$/min_coverage = 15,      # OPTIMIZED: Stricter/' rules/lrt_genome_wide.smk
sed -i 's/min_samples = 4,.*$/min_samples = 6,         # OPTIMIZED: More samples required/' rules/lrt_genome_wide.smk
sed -i 's/min_distinct_freqs = 2.*$/min_distinct_freqs = 3   # OPTIMIZED: More variation/' rules/lrt_genome_wide.smk


sed -i 's/--time=10:00:00/--time=06:00:00/' rules/lrt_genome_wide.smk


# Submit jobs
snakemake \
  glm/lrt_gw/glmV1full/res_multi_001.csv \
  glm/lrt_gw/glmV1full/res_multi_002.csv \
  glm/lrt_gw/glmV1full/res_multi_003.csv \
  --cores 30 \
  --keep-going


WAIT_COUNT=0
MAX_WAIT=60  # Wait up to 2 hours (60 * 2 min)

while [ $WAIT_COUNT -lt $MAX_WAIT ]; do
    COMPLETE=0
    RUNNING=0
    PENDING=0

    for i in 001 002 003; do
        if [ -f "glm/lrt_gw/glmV1full/res_multi_$i.csv" ] && [ -s "glm/lrt_gw/glmV1full/res_multi_$i.csv" ]; then
            ((COMPLETE++))
        elif [ -f "logs/glm_lrt_gw/glmV1full/multi_$i.log" ]; then
            ((RUNNING++))
        else
            ((PENDING++))
        fi
    done

    echo "[$(date +%H:%M:%S)] Status: $COMPLETE completed, $RUNNING running, $PENDING pending"

    if [ $COMPLETE -eq 3 ]; then
        echo ""
        echo "✓ All test bins completed!"
        break
    fi

    sleep 120  # Wait 2 minutes
    ((WAIT_COUNT++))
done



TOTAL_INITIAL=0
TOTAL_REMOVED=0
TOTAL_REMAINING=0

for i in 001 002 003; do
    LOG="logs/glm_lrt_gw/glmV1full/multi_$i.log"

    if [ -f "$LOG" ]; then
        FILTER_LINE=$(grep -A 4 "Pre-filtering sites" "$LOG" 2>/dev/null | grep "Filtered out" || echo "")

        if [ -n "$FILTER_LINE" ]; then
            REMOVED=$(echo "$FILTER_LINE" | grep -oP '\d+(?=/\d+ sites)' || echo "0")
            INITIAL=$(echo "$FILTER_LINE" | grep -oP '(?<=/)\d+(?= sites)' || echo "0")
            TOTAL_INITIAL=$((TOTAL_INITIAL + INITIAL))
            TOTAL_REMOVED=$((TOTAL_REMOVED + REMOVED))
        fi

        RESULT_FILE="glm/lrt_gw/glmV1full/res_multi_$i.csv"
        if [ -f "$RESULT_FILE" ]; then
            SITES=$(tail -n +2 "$RESULT_FILE" | wc -l)
            CONVERGED=$(tail -n +2 "$RESULT_FILE" | cut -d',' -f6 | grep -c "true" || echo "0")
            TOTAL_REMAINING=$((TOTAL_REMAINING + SITES))
        fi

        START=$(grep "Starting R Conversion" "$LOG" 2>/dev/null | tail -1 | grep -oP '\w+ \d+ \d+:\d+:\d+' || echo "")
        FINISH=$(grep "Starting Julia Bootstrap" "$LOG" 2>/dev/null | tail -1 | grep -oP '\w+ \d+ \d+:\d+:\d+' || echo "")
        if [ -n "$START" ]; then
            echo "  Started: $START"
        fi
    fi

done


if [ $TOTAL_INITIAL -gt 0 ]; then
    REMOVAL_PCT=$(( (TOTAL_REMOVED * 100) / TOTAL_INITIAL ))
    if [ $REMOVAL_PCT -ge 70 ]; then
        :
    elif [ $REMOVAL_PCT -le 20 ]; then
        :
    else
        read -p "Launch full genome analysis now? (y/n) " -n 1 -r
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            snakemake glm_lrt_gw_final/glmV1full.csv --cores 40 --keep-going
        else
            :
        fi
    fi
fi