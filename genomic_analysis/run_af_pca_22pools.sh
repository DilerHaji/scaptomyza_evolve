#!/bin/bash
#SBATCH --job-name=af_pca_22
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --qos=savio_lowprio
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER_EMAIL}

cd PATH/TO/DATA/v5_dedup_trim2

PY=/global/software/sl-7.x86_64/modules/langs/python/3.7/bin/python
[ -x "$PY" ] || { echo "FATAL: python not executable at $PY"; exit 1; }

"$PY" -c "import numpy, pandas, sklearn, matplotlib" || {
    echo "FATAL: python missing one of numpy/pandas/sklearn/matplotlib"
    exit 1
}

VCF="fvariants/fvOG_e10fe9w.fixed_no_neff.vcf.gz"
WILD_DIR="baypass_wild/pileups"
FOUNDER_DIR="grenfst/diversity_combined/pileups"
G10_DIR="grenfst/diversity_attrition/pileups"
OUTDIR="final_plots/wild"

mkdir -p "$OUTDIR"

EXPECTED=(
    AVB AVT PSB PST RMB RMT
    F1G00 F2G00 F3G00 F4G00
    B1G10 B2G10 B3G10 B4G10
    T1G10 T2G10 T3G10 T4G10
    M1G10 M2G10 M3G10 M4G10
)
MISSING=0
for S in "${EXPECTED[@]}"; do
    for D in "$WILD_DIR" "$FOUNDER_DIR" "$G10_DIR"; do
        if [ -f "$D/${S}.mpileup" ]; then
            continue 2
        fi
    done
    MISSING=$((MISSING+1))
done
if [ "$MISSING" -gt 0 ]; then
    exit 1
fi

export PYTHONUNBUFFERED=1
"$PY" -u scripts/af_pca_22pools.py \
    --vcf          "$VCF" \
    --wild-dir     "$WILD_DIR" \
    --founder-dir  "$FOUNDER_DIR" \
    --g10-dir      "$G10_DIR" \
    --outdir       "$OUTDIR" \
    --min-cov      10 \
    --max-cov      500 \
    --af-cache     "$OUTDIR/af_matrix_22pools.csv"

ls -lh "$OUTDIR"/af_pca_22pools* "$OUTDIR"/af_matrix_22pools.csv 2>/dev/null
