import math
import sys


def get_num_permutations(wildcards):
    if wildcards.grenfst_multiplot not in GRENFST_MULTIPLOT_DICT:
        return 1
        
    targets = GRENFST_MULTIPLOT_DICT[wildcards.grenfst_multiplot][1]
    refs    = GRENFST_MULTIPLOT_DICT[wildcards.grenfst_multiplot][2]
    
    if len(targets) != len(refs):
        return 1 
    
    return math.factorial(len(refs))


rule lmm_melt_perm:
    input:
        lambda w: (
            "grenfst/fst_queue/" + GRENFST_MULTIPLOT_DICT[w.grenfst_multiplot][0] + ".csv"
            if ("btwTM" in w.grenfst_multiplot or "btwBM" in w.grenfst_multiplot)
            else "grenfst/fst_queue/" + GRENFST_MULTIPLOT_DICT[w.grenfst_multiplot][0] + ".csv"
        )
#         lambda w: (
#             "grenfst/fst_queue-with-mixture/" + GRENFST_MULTIPLOT_DICT[w.grenfst_multiplot][0] + ".csv"
#             if ("btwTM" in w.grenfst_multiplot or "btwBM" in w.grenfst_multiplot)
#             else "grenfst/fst_queue/" + GRENFST_MULTIPLOT_DICT[w.grenfst_multiplot][0] + ".csv"
#         )
    output:
        "grenfst/lmm_perms/{grenfst_multiplot}/perm_{perm_id}.long.csv"
    params:
        target_reps = lambda w: ",".join(GRENFST_MULTIPLOT_DICT[w.grenfst_multiplot][1]), 
        ref_reps    = lambda w: ",".join(GRENFST_MULTIPLOT_DICT[w.grenfst_multiplot][2]), 
        generations = lambda w: ",".join(GRENFST_MULTIPLOT_DICT[w.grenfst_multiplot][3]),
        script      = "scripts/melt_permuted.py",
        perm_id     = "{perm_id}"
    conda:
        config["environments"]["polars"]
    resources:
        resources=config.get("default_resources_10cpus", "")
    shell:
        """
        python {params.script} \
        --input {input} \
        --output {output} \
        --perm-id {params.perm_id} \
        --target-reps "{params.target_reps}" \
        --ref-reps "{params.ref_reps}" \
        --generations "{params.generations}"
        """



rule lmm_run_perm_bin:
    input:
        "grenfst/lmm_perms/{grenfst_multiplot}/perm_{perm_id}.long.csv"
    output:
        "grenfst/lmm_perms/{grenfst_multiplot}/chunks/perm_{perm_id}_bin_{bin_id}.csv"
    params:
        script = "scripts/fst_lmm.py",
        chroms = get_chroms_for_bin
    conda:
        config["environments"]["statsmodels"]
    threads: 24 
    resources:
        resources=config.get("default_resources_24cpus", "")    
    shell:
        """
        python {params.script} \
        --input {input} \
        --output {output} \
        --chroms "{params.chroms}" \
        --threads {threads}  # CHANGE 3: Pass the Snakemake threads to Python
        """
        

rule lmm_gather_bins:
    input:
        lambda w: expand(
            "grenfst/lmm_perms/{{grenfst_multiplot}}/chunks/perm_{{perm_id}}_bin_{bin_id}.csv", 
            bin_id=BIN_IDS
        )
    output:
        "grenfst/lmm_perms/{grenfst_multiplot}/results/perm_{perm_id}_full.csv"
    conda:
        config["environments"]["polars"]
    resources:
        resources=config.get("default_resources", "")
    script:
        config["scripts"]["concat_bins"]
        
        

rule lmm_consensus:
    input:
        lambda w: expand(
            "grenfst/lmm_perms/{{grenfst_multiplot}}/results/perm_{perm_id}_full.csv", 
            perm_id=range(get_num_permutations(w))
        )
    output:
        "grenfst/lmm_results/{grenfst_multiplot}/consensus_lmm.csv"
    params:
        script = "scripts/aggregate_permutations.py"
    conda:
        config["environments"]["polars"]
    resources:
        resources=config.get("default_resources_24cpus", "")
    shell:
        """
        python {params.script} \
        --inputs {input} \
        --output {output}
        """



rule lmm_plot_consensus:
    input:
        "grenfst/lmm_results/{grenfst_multiplot}/consensus_lmm.csv"
    output:
        "grenfst/lmm_results/{grenfst_multiplot}/manhattan_consensus.png"
    params:
        script = config["scripts"]["plot_residuals"],
        width = 20, height = 6
    conda:
        config["environments"]["polars"]
    resources:
        resources=config["default_resources"]
    shell:
        """
        python {params.script} \
        --input {input} \
        --output {output} \
        --slope-col z_score_median \
        --r2-col z_score_sd \
        --ymin -5 --ymax 15 \
        --title "Replicated Divergence (Consensus Median Z-Score)"
        """
























def get_scale_space_inputs_lmm(wildcards):
    group = wildcards.scale_group
    matching_keys = [
        k for k in GRENFST_MULTIPLOT_DICT.keys()
        if k.startswith(f"{group}_")
    ]

    if not matching_keys:
        raise ValueError(f"No entries found in GRENFST_MULTIPLOT_DICT starting with '{group}_'")

    def extract_size(k):
        try:
            suffix = k.replace(f"{group}_", "")
            return int(suffix.split('_')[0])
        except:
            return 0 
    
    matching_keys.sort(key=extract_size)

    files = [f"grenfst/lmm_results/{k}/consensus_lmm.csv" for k in matching_keys]
        
    return files


rule grenfst_scale_space_plot_lmm:
    input:
        files = lambda w: get_scale_space_inputs_lmm(w)
    output:
        "grenfst/scale_space_plots_lmm/{scale_group}.png"
    params:
        script = config["scripts"]["plot_scale_space_filtered"],
        bin_size = 5000, 
        width = 15,
        height = 3,
        column = "z_score_median",
        # Plotting Top 10% (0.9 to 1.0)
        vmin = 0.90, 
        vmax = 1.0,   
        gap = 20
    conda:
        config["environments"]["polars"]
    resources:
        resources=config.get("default_resources_24cpus", "")
    shell:
        """
        python {params.script} \
        --inputs {input.files} \
        --output {output} \
        --bin-size {params.bin_size} \
        --width {params.width} \
        --height {params.height} \
        --column {params.column} \
        --top-n 10 \
        --vmin {params.vmin} \
        --vmax {params.vmax} \
        --gap {params.gap} \
        --percentile  # <--- NEW FLAG
        """


def get_all_lmm_inputs():
    return [
        f"grenfst/lmm_results/{k}/consensus_lmm.csv" 
        for k in GRENFST_MULTIPLOT_DICT.keys()
    ]

rule grenfst_global_multipanel_lmm:
    input:
        files = get_all_lmm_inputs()
    output:
        "grenfst/global_summary/combined_lmm_zscores.png"
    params:
        script = config["scripts"]["plot_scale_space_all"],
        bin_size = 1000, 
        width = 20,
        panel_height = 1.5,
        top_n = 10,
        column = "z_score_median",
        vmin = 5.0,
        vmax = 20.0, 
        noise_threshold = 1.5,
        smooth_sigma = 5.0, 
        gap = 100,   
        cmap = "mako" 
    conda:
        config["environments"]["polars"]
    resources:
        resources=config.get("default_resources_24cpus", "")
    shell:
        """
        python {params.script} \
        --inputs {input.files} \
        --output {output} \
        --bin-size {params.bin_size} \
        --width {params.width} \
        --panel-height {params.panel_height} \
        --column {params.column} \
        --top-n {params.top_n} \
        --vmin {params.vmin} \
        --vmax {params.vmax} \
        --gap {params.gap} \
        --common-only \
        --cmap {params.cmap} \
        --smooth-sigma {params.smooth_sigma} \
        --noise-threshold {params.noise_threshold} \
        --font-scale 1.5
        """