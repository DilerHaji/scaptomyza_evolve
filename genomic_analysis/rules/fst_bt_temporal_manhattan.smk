# ---------------------------------------------------------------------
# FST(B,T) TEMPORAL MANHATTAN PLOTS
#
# Generates stacked Manhattan plots of FST between B and T treatments
# at each generation, using the raw grenedalf FST queue CSV.
#
# Includes temporal stability analysis:
#   - CV of FST(B,T) across generations per window
#   - Mean FST vs CV scatter to identify balancing selection candidates
#
# Inputs:  grenfst/fst_queue/{analysis}.csv (already computed by grenedalf)
# Outputs: final_plots/fst_bt_temporal/{window}bp_fst_bt_temporal.png/svg
#
# No additional heavy computation — purely reads + plots.
# ---------------------------------------------------------------------

# Window sizes for FST Manhattan — use a subset that balances resolution and noise
_FST_MANHATTAN_WINDOWS = [200000, 500000]

import re as _re_fst

def _fst_queue_key_for_window(window):
    """Return the GRENFST_MULTIPLOT_DICT key for a given window size (full timecourse)."""
    for k in GRENFST_MULTIPLOT_DICT.keys():
        if "btwTB" in k and str(window) in k and not _re_fst.search(r'_\d{2}$', k):
            return k
    return None


def _fst_queue_csv_for_window(window):
    """Return the FST queue CSV path for a given window size."""
    k = _fst_queue_key_for_window(window)
    if k is None:
        return None
    return f"grenfst/fst_queue/{GRENFST_MULTIPLOT_DICT[k][0]}.csv"


def _fst_gens_for_window(window):
    """Return comma-separated generation list for a given window."""
    k = _fst_queue_key_for_window(window)
    if k is None:
        return "01,02,06,07,08,09,10"
    return ",".join(GRENFST_MULTIPLOT_DICT[k][3])


rule fst_bt_temporal_manhattan:
    input:
        csv = lambda w: _fst_queue_csv_for_window(int(w.window)),
    output:
        png = "final_plots/fst_bt_temporal/{window}bp_fst_bt_temporal.png",
        svg = "final_plots/fst_bt_temporal/{window}bp_fst_bt_temporal.svg",
        csv = "final_plots/fst_bt_temporal/{window}bp_fst_bt_per_gen.csv",
    log:
        "logs/fst_bt_temporal/{window}bp_fst_bt_temporal.log"
    benchmark:
        "benchmarks/fst_bt_temporal/{window}bp_fst_bt_temporal.log"
    params:
        script      = "scripts/fst_bt_temporal_manhattan.py",
        target_reps = "T1,T2,T3,T4",
        ref_reps    = "B1,B2,B3,B4",
        generations = lambda w: _fst_gens_for_window(int(w.window)),
    wildcard_constraints:
        window = "[0-9]+"
    conda:
        config["environments"]["polars"]
    resources:
        resources = config.get("default_resources", "")
    shell:
        """
        mkdir -p $(dirname {output.png})
        python {params.script} \
            --input           {input.csv} \
            --output          {output.png} \
            --output-svg      {output.svg} \
            --output-csv      {output.csv} \
            --target-reps     "{params.target_reps}" \
            --ref-reps        "{params.ref_reps}" \
            --generations     "{params.generations}" \
            --cv-panel \
        > {log} 2>&1
        """


rule all_fst_bt_temporal:
    input:
        expand("final_plots/fst_bt_temporal/{window}bp_fst_bt_temporal.png",
               window=_FST_MANHATTAN_WINDOWS)
