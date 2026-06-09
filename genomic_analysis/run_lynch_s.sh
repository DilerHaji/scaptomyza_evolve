#!/bin/bash

#SBATCH --job-name=lynch_s
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --qos=rosalind_htc4_normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${USER_EMAIL}

export PYTHONUNBUFFERED=1

. ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate .snakemake/conda/a20f2f036ca1245bd5eba6e9fba0a5ce_

python -c "import polars; import numpy; print(f'polars {polars.__version__}, numpy {numpy.__version__}')"

LYNCH_SCRIPT=scripts/lynch.py
SUMMARY_SCRIPT=scripts/summarize_lynch_s.py

mkdir -p lynch_s lynch_summary


declare -A SAMPLES
SAMPLES[dB1]="B1G01|B1G02|B1G03|B1G04|B1G05|B1G06|B1G07|B1G08|B1G09|B1G10"
SAMPLES[dB2]="B2G01|B2G02|B2G03|B2G04|B2G05|B2G06|B2G07|B2G08|B2G09|B2G10"
SAMPLES[dB3]="B3G01|B3G02|B3G06|B3G07|B3G08|B3G09|B3G10"
SAMPLES[dB4]="B4G01|B4G02|B4G06|B4G07|B4G08|B4G09|B4G10"
SAMPLES[dT1]="T1G01|T1G02|T1G03|T1G04|T1G05|T1G06|T1G07|T1G08|T1G09|T1G10"
SAMPLES[dT2]="T2G01|T2G02|T2G03|T2G04|T2G05|T2G06|T2G07|T2G08|T2G09|T2G10"
SAMPLES[dT3]="T3G01|T3G02|T3G06|T3G07|T3G08|T3G09|T3G10"
SAMPLES[dT4]="T4G01|T4G02|T4G06|T4G07|T4G08|T4G09|T4G10"
SAMPLES[dM1]="M1G01|M1G02|M1G03|M1G04|M1G05|M1G06|M1G07|M1G08|M1G09|M1G10"
SAMPLES[dM2]="M2G01|M2G02|M2G03|M2G04|M2G05|M2G06|M2G07|M2G08|M2G09|M2G10"
SAMPLES[dM3]="M3G01|M3G02|M3G06|M3G07|M3G08|M3G09|M3G10"
SAMPLES[dM4]="M4G01|M4G02|M4G06|M4G07|M4G08|M4G09|M4G10"

MEAN_S_FILES=""
LABELS=""

for delta in dB1 dB2 dB3 dB4 dT1 dT2 dT3 dT4 dM1 dM2 dM3 dM4; do
    TRAJ_OUT="lynch_s/${delta}_trajectory.csv"
    MEAN_OUT="lynch_s/${delta}_mean_s.csv"
    S_OUT="lynch_s/${delta}_per_interval_s.csv"

    if [ -s "$MEAN_OUT" ]; then
    else
        python $LYNCH_SCRIPT \
            --af_file "delta_tmp/${delta}/afmatsites.csv" \
            --neff_file "delta_tmp/${delta}/neffsites.csv" \
            --samples "${SAMPLES[$delta]}" \
            --output "$TRAJ_OUT" \
            --mean_output "$MEAN_OUT" \
            --s_output "$S_OUT"
    fi

    if [[ "$delta" == dB* ]] || [[ "$delta" == dT* ]]; then
        label=${delta#d}  # strip 'd' prefix
        MEAN_S_FILES="$MEAN_S_FILES $MEAN_OUT"
        LABELS="$LABELS $label"
    fi
done


python $SUMMARY_SCRIPT \
    --mean_s_files $MEAN_S_FILES \
    --labels $LABELS \
    --window_sizes 200000 500000 \
    --candidate_regions focused_mimicree/candidate_regions.tsv \
    --output_prefix lynch_summary/summary

