#!/bin/bash

set -euo pipefail

SCRIPTDIR="scripts"
SYNCDIR="hv_sync"
OUTBASE="hv_results"
SCAFFOLDS="chr_ScDA7r2_597_HRSCAF_953,chr_ScDA7r2_110_HRSCAF_295,chr_ScDA7r2_439_HRSCAF_779,chr_ScDA7r2_126_HRSCAF_325"

SINGLE_TRTS=("B" "T" "M")

PAIR_TRTS=("BT" "BM" "TM")

GENS="0,1,2,6,7,8,9"
NE=250

mkdir -p ${OUTBASE} ${SYNCDIR} slurm_logs


build_acer_samples() {
    local trt=$1
    local n_repl=$2
    local result=""

    if [ ${#trt} -eq 1 ]; then
        for r in $(seq 1 $n_repl); do
            [ -n "$result" ] && result="${result},"
            result="${result}F${r}G00"
            for g in 01 02 06 07 08 09; do
                result="${result},${trt}${r}G${g}"
            done
        done

    elif [ ${#trt} -eq 2 ]; then
        local t1=${trt:0:1}
        local t2=${trt:1:1}
        for r in $(seq 1 4); do
            [ -n "$result" ] && result="${result},"
            result="${result}F${r}G00"
            for g in 01 02 06 07 08 09; do
                result="${result},${t1}${r}G${g}"
            done
        done
        for r in $(seq 1 4); do
            result="${result},F${r}G00"
            for g in 01 02 06 07 08 09; do
                result="${result},${t2}${r}G${g}"
            done
        done
    fi
    echo "$result"
}


submit_hv_job() {
    local trt=$1
    local poolsize=$2
    local n_repl=$3
    local label=$4  # A_nominal or B_effective
    local sync_file="${SYNCDIR}/hv_${trt}.sync"
    local outdir="${OUTBASE}/hv_${trt}_${label}"
    local acer_samples=$(build_acer_samples $trt $n_repl)

    local jobname="hv_${trt}_${label}"

    sbatch --job-name=${jobname} \
           --account=${SLURM_ACCOUNT} \
           --partition=savio4_htc \
           --time=04:00:00 \
           --mem=32G \
           --cpus-per-task=1 \
           --output=slurm_logs/${jobname}_%j.out \
           --error=slurm_logs/${jobname}_%j.err \
           --wrap="
module load r
Rscript ${SCRIPTDIR}/hv_run_cluster.R \
  --sync ${sync_file} \
  --ad variance_analysis/merged_ad.tsv \
  --depth variance_analysis/merged_depth.tsv \
  --sample_list variance_analysis/sample_list.txt \
  --treatment ${trt} \
  --ne ${NE} \
  --poolsize ${poolsize} \
  --n_repl ${n_repl} \
  --outdir ${outdir} \
  --scaffolds ${SCAFFOLDS} \
  --max_cands_per_scaffold 50000 \
  --acer_samples ${acer_samples}
"
}



for trt in "${SINGLE_TRTS[@]}"; do
    if [ ! -f "${SYNCDIR}/hv_${trt}.sync" ]; then
        echo "  Building sync for ${trt}..."
        python3 ${SCRIPTDIR}/build_hv_sync.py \
            --treatment ${trt} \
            --scaffold ${SCAFFOLDS} \
            --ad variance_analysis/merged_ad.tsv \
            --sample_list variance_analysis/sample_list.txt \
            --out ${SYNCDIR}/hv_${trt}.sync
    else
        echo "  Sync exists: ${SYNCDIR}/hv_${trt}.sync"
    fi
done

for trt in "${PAIR_TRTS[@]}"; do
    if [ ! -f "${SYNCDIR}/hv_${trt}.sync" ]; then
        echo "  Building sync for ${trt}..."
        python3 ${SCRIPTDIR}/build_hv_sync.py \
            --treatment ${trt} \
            --scaffold ${SCAFFOLDS} \
            --ad variance_analysis/merged_ad.tsv \
            --sample_list variance_analysis/sample_list.txt \
            --out ${SYNCDIR}/hv_${trt}.sync
    fi
done



for trt in "${SINGLE_TRTS[@]}"; do
    echo "Treatment ${trt}:"
    submit_hv_job ${trt} 160 4 "A_nominal"    
    submit_hv_job ${trt} 58  4 "B_effective"   
done



for trt in "${PAIR_TRTS[@]}"; do
    if [ -f "${SYNCDIR}/hv_${trt}.sync" ]; then
        echo "Treatment ${trt}:"
        submit_hv_job ${trt} 160 8 "A_nominal"
        submit_hv_job ${trt} 58  8 "B_effective"
    fi
done
