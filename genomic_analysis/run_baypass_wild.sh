#!/bin/bash
#SBATCH --job-name=bp_wild
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --qos=rosalind_htc4_normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER_EMAIL}

export PYTHONUNBUFFERED=1
. ~/.bashrc
eval "$(conda shell.bash hook)"

# ── paths ─────────────────────────────────────────────────────────────────────
PIPELINE_DIR="PATH/TO/DATA/v5_dedup_trim2"
BAM_DIR="PATH/TO/DATA/v5_dedup_trim/mapping/filtered"
REF="${PIPELINE_DIR}/sfla_v2.fa"

# VCF that defines the experimental variant positions + REF/ALT alleles.
# Same file used by grenedalf, PBS, and the experimental BayPass runs.
VCF="${PIPELINE_DIR}/fvariants/fvOG_e10fe9w.fixed_no_neff.vcf.gz"

POLARS_PYTHON="${PIPELINE_DIR}/.snakemake/conda/a20f2f036ca1245bd5eba6e9fba0a5ce_/bin/python"
BAYPASS_BIN="${PIPELINE_DIR}/baypass/sources/g_baypass"
NTHREADS=8

OUTDIR="${PIPELINE_DIR}/baypass_wild"
PILEUP_DIR="${OUTDIR}/pileups"
PREFIX="wild"

# Set to 1 to include S.flavaMA and S.flavaAZ as unspecialized populations.
# Their BAMs were not found in the standard BAM_DIR — set to 0 until located.
INCLUDE_UNSPEC=0

mkdir -p "${PILEUP_DIR}"

echo "=== BayPass Wild Analysis ==="
echo "Start: $(date)"
echo "Include unspecialized populations: ${INCLUDE_UNSPEC}"

# ── Step 1: generate positions BED from VCF ───────────────────────────────────
echo ""
echo "--- Step 1: Extract variant positions from VCF ---"
POS_BED="${OUTDIR}/variant_positions.bed"
if [ ! -f "${POS_BED}" ]; then
    bcftools view -H "${VCF}" \
      | awk '{print $1"\t"($2-1)"\t"$2}' \
      > "${POS_BED}"
    echo "  Written: ${POS_BED} ($(wc -l < ${POS_BED}) sites)"
else
    echo "  Exists, skipping: ${POS_BED}"
fi

# ── Step 2: generate per-sample mpileups ──────────────────────────────────────
echo ""
echo "--- Step 2: Per-sample mpileup at variant positions ---"

CORE_SAMPLES=("AVB" "AVT" "PSB" "PST" "RMB" "RMT")
if [ "${INCLUDE_UNSPEC}" -eq 1 ]; then
    ALL_SAMPLES=("AVB" "AVT" "PSB" "PST" "RMB" "RMT" "S.flavaMA" "S.flavaAZ")
else
    ALL_SAMPLES=("AVB" "AVT" "PSB" "PST" "RMB" "RMT")
fi

for SAMPLE in "${ALL_SAMPLES[@]}"; do
    MP="${PILEUP_DIR}/${SAMPLE}.mpileup"
    BAM="${BAM_DIR}/${SAMPLE}.bam"
    if [ ! -f "${MP}" ]; then
        if [ ! -f "${BAM}" ]; then
            echo "  WARNING: BAM not found for ${SAMPLE}: ${BAM}"
            continue
        fi
        echo "  Pileup: ${SAMPLE}"
        samtools mpileup \
            -q 25 -Q 25 \
            -B --ignore-RG \
            -l "${POS_BED}" \
            -f "${REF}" \
            "${BAM}" \
            > "${MP}"
        echo "    Done: $(wc -l < ${MP}) lines"
    else
        echo "  Exists, skipping: ${SAMPLE}"
    fi
done

# ── Step 3: prepare BayPass input ─────────────────────────────────────────────
echo ""
echo "--- Step 3: Prepare BayPass input ---"

UNSPEC_FLAG=""
if [ "${INCLUDE_UNSPEC}" -eq 1 ]; then
    UNSPEC_FLAG="--include-unspec"
fi

if [ ! -f "${OUTDIR}/${PREFIX}_pooldata.geno" ]; then
    ${POLARS_PYTHON} "${PIPELINE_DIR}/scripts/prepare_baypass_wild.py" \
        --vcf             "${VCF}" \
        --mpileup-dir     "${PILEUP_DIR}" \
        --output-dir      "${OUTDIR}" \
        --prefix          "${PREFIX}" \
        --min-cov         5 \
        --min-pop-cov     1.0 \
        --thin-step       16 \
        ${UNSPEC_FLAG}
else
    echo "  Input files exist, skipping preparation"
fi

cd "${OUTDIR}"

# ── Step 4: Omega estimation (thinned SNPs) ───────────────────────────────────
echo ""
echo "--- Step 4: Omega estimation (thinned SNPs) ---"
if [ ! -f "${PREFIX}_omega_mat_omega.out" ]; then
    "${BAYPASS_BIN}" \
        -pooldatafile ${PREFIX}_omega_pooldata.geno \
        -poolsizefile ${PREFIX}_poolsize.txt \
        -outprefix    ${PREFIX}_omega \
        -nthreads     ${NTHREADS}
    echo "  Omega done: $(date)"
else
    echo "  Omega exists, skipping"
fi

# ── Step 5: Covariate model (full SNPs, fixed Omega) ─────────────────────────
echo ""
echo "--- Step 5: Treatment covariate model (B/T, fixed Omega) ---"
if [ ! -f "${PREFIX}_trt_summary_betai_reg.out" ]; then
    "${BAYPASS_BIN}" \
        -pooldatafile ${PREFIX}_pooldata.geno \
        -poolsizefile ${PREFIX}_poolsize.txt \
        -omegafile    ${PREFIX}_omega_mat_omega.out \
        -efile        ${PREFIX}_treatment.cov \
        -outprefix    ${PREFIX}_trt \
        -nthreads     ${NTHREADS}
    echo "  Covariate model done: $(date)"
else
    echo "  Covariate model exists, skipping"
fi

# ── Step 6: Contrast model (full SNPs, fixed Omega) ──────────────────────────
echo ""
echo "--- Step 6: Contrast model (fixed Omega) ---"
if [ ! -f "${PREFIX}_contrast_summary_contrast.out" ]; then
    "${BAYPASS_BIN}" \
        -pooldatafile ${PREFIX}_pooldata.geno \
        -poolsizefile ${PREFIX}_poolsize.txt \
        -omegafile    ${PREFIX}_omega_mat_omega.out \
        -contrastfile ${PREFIX}_contrasts.con \
        -outprefix    ${PREFIX}_contrast \
        -nthreads     ${NTHREADS}
    echo "  Contrast model done: $(date)"
else
    echo "  Contrast model exists, skipping"
fi

cd "${PIPELINE_DIR}"

echo ""
echo "=== BayPass Wild complete: $(date) ==="
echo "Output: ${OUTDIR}/"
ls -lh "${OUTDIR}/${PREFIX}"_*.out 2>/dev/null || echo "(no .out files yet)"
