# ---------------------------------------------------------------------
# Per-window nucleotide diversity (θ_π, θ_W, Tajima's D) via grenedalf.
#
# Generalised over named sample sets so the same bin-and-gather machinery
# handles the founder-only run AND a full 112-pool wild → founder → G01–G10
# trajectory run. All sample sets share the same window geometry (200 kb,
# non-overlapping) and filter settings, chosen to match the PBS sweep grid
# and exceed the empirical LD-decay scale (see note on _DIV_WINDOW below).
#
# Outputs
#   grenfst/diversity/{sample_set}_pi_{window}.csv
#
# Targets
#   rule all_founder_diversity     — founder pools only (F1-F4)
#   rule all_trajectory_diversity  — 112 pools (wild + founder + G01-G10)
# ---------------------------------------------------------------------

import sys

# Reuse GREN_CONTIG_BINS and BIN_IDS from grenedalf_interval.smk.

# Chosen to (i) exceed the empirical LD-decay scale measured by ngsLD
# (background r^2 ≈ 0.02 already at ≥ 10 kb) by >20×, so adjacent windows
# are approximately independent for per-window diversity estimates, and
# (ii) align with the 200-kb grid point in the PBS sweep / timeslice
# analyses so window-based tables cross-reference exactly.
_DIV_WINDOW = 200000
_DIV_STRIDE = 1000

# -----------------------------------------------------------------------------
# Sample-set definitions. The value is the string grenedalf expects for
# --filter-samples-include: a comma-separated list of sample names (matching
# the VCF's column headers — no .1 suffix here).
# -----------------------------------------------------------------------------

_WILD_POOLS     = ["AVB", "AVT", "PSB", "PST", "RMB", "RMT"]
_FOUNDER_POOLS  = [f"F{i}G00" for i in range(1, 5)]
_FULL_REPS      = ["B1", "B2", "T1", "T2", "M1", "M2"]   # G01-G09 all present
_PARTIAL_REPS   = ["B3", "B4", "T3", "T4", "M3", "M4"]   # G01, G02, G06-G09
_FULL_GENS      = [f"{g:02d}" for g in range(1, 10)]
_PARTIAL_GENS   = ["01", "02", "06", "07", "08", "09"]
_G10_POOLS      = [f"{r}G10" for r in _FULL_REPS + _PARTIAL_REPS]

_TRAJECTORY_INTERMEDIATES = (
    [f"{r}G{g}" for r in _FULL_REPS    for g in _FULL_GENS] +
    [f"{r}G{g}" for r in _PARTIAL_REPS for g in _PARTIAL_GENS]
)

DIVERSITY_SAMPLE_SETS = {
    "founder":    ",".join(_FOUNDER_POOLS),
    # NOTE: wild pools (AVB, AVT, PSB, PST, RMB, RMT) are NOT in the main
    # experimental VCF — they live as per-sample pileups (see
    # run_attrition_diversity.sh). For the VCF-based trajectory run we include
    # only the 106 samples that exist in fvOG_e10fe9w.fixed_no_neff.vcf.gz:
    # 4 founders + 102 G01-G10 trajectory pools.
    "trajectory": ",".join(
        _FOUNDER_POOLS + _TRAJECTORY_INTERMEDIATES + _G10_POOLS
    ),
}


rule diversity_bin:
    input:
        vcf = "fvariants/fvOG_e10fe9w.fixed_no_neff.vcf.gz",
    output:
        csv = "grenfst/diversity_bins/{sample_set}/{window}/{bin_id}/diversity.csv",
    params:
        run_dir            = "grenfst/diversity_bins/{sample_set}/{window}/{bin_id}",
        window_width       = "{window}",
        window_stride      = _DIV_STRIDE,
        pool_sizes_file    = "pool_sizes.tsv",
        chrom_list_newline = lambda wc: GREN_CONTIG_BINS[wc.bin_id].replace(",", "\n"),
        samples            = lambda wc: DIVERSITY_SAMPLE_SETS[wc.sample_set],
    wildcard_constraints:
        sample_set = "|".join(DIVERSITY_SAMPLE_SETS.keys()),
        window     = r"[0-9]+",
        bin_id     = r"[0-9]+",
    resources:
        resources = config["default_resources"]
    conda:
        config["environments"]["grenedalf"]
    shell:
        """
        VCF_ABS=$(readlink -f {input.vcf})
        POOLS_ABS=$(readlink -f {params.pool_sizes_file})

        [[ -d {params.run_dir} ]] || mkdir -p {params.run_dir}
        cd {params.run_dir}

        echo -e "{params.chrom_list_newline}" > region_list.txt

        rm -f diversity.csv diversity_tmp.csv

        set +e
        grenedalf diversity \
            --vcf-path "$VCF_ABS" \
            --filter-samples-include "{params.samples}" \
            --pool-sizes "$POOLS_ABS" \
            --filter-sample-min-count 2 \
            --filter-sample-min-read-depth 4 \
            --filter-sample-max-read-depth 500 \
            --window-type interval \
            --window-interval-width {params.window_width} \
            --window-interval-stride {params.window_stride} \
            --window-average-policy valid-snps \
            --filter-region-list region_list.txt \
            --allow-file-overwriting \
            --file-suffix _tmp
        GREN_EXIT=$?
        set -e

        # Fail hard if grenedalf itself errored — regardless of whether a
        # partial _tmp file exists or not. An earlier version of this rule
        # silently created an empty diversity.csv when grenedalf aborted
        # before writing _tmp, which let snakemake report success on a
        # broken run.
        if [[ $GREN_EXIT -ne 0 ]]; then
            rm -f diversity_tmp.csv diversity.csv
            echo "ERROR: grenedalf diversity failed (exit $GREN_EXIT)" >&2
            exit 1
        fi

        if [[ -f diversity_tmp.csv ]]; then
            mv diversity_tmp.csv diversity.csv
        else
            # grenedalf exited 0 but produced no output (e.g. empty region).
            # Emit an empty file so gather still runs.
            touch diversity.csv
        fi
        """


rule diversity_gather:
    input:
        lambda wc: expand(
            "grenfst/diversity_bins/{sample_set}/{window}/{bin_id}/diversity.csv",
            sample_set=wc.sample_set,
            window=wc.window,
            bin_id=BIN_IDS,
        )
    output:
        "grenfst/diversity/{sample_set}_pi_{window}.csv"
    wildcard_constraints:
        sample_set = "|".join(DIVERSITY_SAMPLE_SETS.keys()),
        window     = r"[0-9]+",
    resources:
        resources = config["default_resources"]
    shell:
        """
        mkdir -p $(dirname {output})

        HEADER_FOUND=0
        for f in {input}; do
            if [[ -s "$f" ]]; then
                head -n 1 "$f" > {output}
                HEADER_FOUND=1
                break
            fi
        done

        if [[ $HEADER_FOUND -eq 0 ]]; then
            touch {output}
            exit 0
        fi

        for f in {input}; do
            if [[ -s "$f" ]]; then
                tail -n +2 "$f" >> {output}
            fi
        done
        """


# -----------------------------------------------------------------------------
# Target rules (named for backwards-compat with existing driver scripts).
# -----------------------------------------------------------------------------

rule all_founder_diversity:
    input:
        f"grenfst/diversity/founder_pi_{_DIV_WINDOW}.csv"


rule all_trajectory_diversity:
    input:
        f"grenfst/diversity/trajectory_pi_{_DIV_WINDOW}.csv"
