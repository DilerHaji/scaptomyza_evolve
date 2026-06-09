rule final_fst_glm_mixture_scan_Fig3:
    input:
        glm = "glm_lrt_gw_final/{glm}.csv",  # LRT with pre-filtering
        heatmaps = expand("grenfst/divergence_mixture_lmm/{analysis}/divergence_stats.csv", 
                          analysis=[k for k in GRENFST_MULTIPLOT_DICT.keys() if "btwTB" in k]),
        fasta = config["reference_genome"],
        repeats = "../popgen/repeat_content/masked_sfla_v2/sfla_v2.fa.out",
        chrom_map = "chrom_map.tsv",
        bed = config["palaez"]
    output:
        plot = "final_plots/fst_glm_mixture/final_aligned_landscape_{glm}_{stat}.png",
        candidates = "final_plots/fst_glm_mixture/final_aligned_landscape_{glm}_{stat}_candidate_sites_detailed.tsv"
    benchmark:
        "benchmarks/window_glm_mixture/final_aligned_landscape_{glm}_{stat}.log"
    log:
        "logs/window_glm_mixture/final_aligned_landscape_{glm}_{stat}.log"
    params:
        script = config["scripts"]["window_glm_aligned_v3"], 
        heatmap_col = lambda w: f"z_{w.stat}",          
        heatmap_col_ref = lambda w: f"z_{w.stat}_ref",  
        target_pval = "PB_p_val",  # LRT bootstrap p-value
        bin_size = 390000,
        ribbon_alpha = 0.15,
        min_alpha = 0.05,
        y_max = 1,
        y_min = 0.90,
        fst_top_pct = 3,
        glm_top_pct = 0.1,
        fig_width = 12,
        fig_height = 3,
        fig_ratios = "0.5 0.8 0.2 0.2"
    conda:
        config["environments"]["polars"]
    resources:
        resources = config.get("default_resources", "")
    shell:
        """
        python {params.script} \
        --heatmap-files {input.heatmaps} \
        --heatmap-col "{params.heatmap_col}" \
        --heatmap-col-ref "{params.heatmap_col_ref}" \
        --glm {input.glm} \
        --target-pval "{params.target_pval}" \
        --genome-fasta {input.fasta} \
        --repeat-masker {input.repeats} \
        --chrom-map {input.chrom_map} \
        --bed-file {input.bed} \
        --output {output.plot} \
        --percentile \
        --bin-size {params.bin_size} \
        --y-max {params.y_max} \
        --y-min {params.y_min} \
        --ribbon-alpha {params.ribbon_alpha} \
        --fst-top-pct {params.fst_top_pct} \
        --glm-top-pct {params.glm_top_pct} \
        --fig-width {params.fig_width} \
        --fig-height {params.fig_height} \
        --fig-ratios {params.fig_ratios} \
        --min-alpha {params.min_alpha} \
        --tracks heatmap glm overlap footer \
        --remove-repeats
        """
        
        

rule final_fst_glm_mixture_scan_interactive:
    input:
        glm = "glm_lrt_gw_final/{glm}.csv",  # LRT with pre-filtering
        heatmaps = expand("grenfst/divergence_mixture_lmm/{analysis}/divergence_stats.csv",
                          analysis=[k for k in GRENFST_MULTIPLOT_DICT.keys() if "btwTB" in k]),
        fasta = config["reference_genome"],
        repeats = "../popgen/repeat_content/masked_sfla_v2/sfla_v2.fa.out",
        chrom_map = "chrom_map.tsv",
        bed = config["palaez"]
    output:
        plot = "final_plots/fst_glm_mixture_interactive/{glm}/aligned.png",
        plot_clean = "final_plots/fst_glm_mixture_interactive/{glm}/aligned_clean.png",
        plot_svg = "final_plots/fst_glm_mixture_interactive/{glm}/aligned.svg",
        json = "final_plots/fst_glm_mixture_interactive/{glm}/aligned.json",
        candidates = "final_plots/fst_glm_mixture_interactive/{glm}/aligned_candidate_sites_detailed.tsv",
        candidates_summary = "final_plots/fst_glm_mixture_interactive/{glm}/aligned_candidates.tsv"
    benchmark:
        "benchmarks/window_glm_mixture/final_aligned_landscape_interactive_{glm}.log"
    log:
        "logs/window_glm_mixture/final_aligned_landscape_interactive_{glm}.log"
    params:
        script = config["scripts"]["window_glm_aligned_v3_interactive"],
        heatmap_col = "z_divergence",
        target_pval = "PB_p_val",  # LRT bootstrap p-value
        bin_size = 200000,
        ribbon_alpha = 0.3,
        min_alpha = 0.05,
        y_max = 0.75,
        y_min = 0.3,
        fst_top_pct = 2,
        glm_top_pct = 0.05,
        fig_width = 10,
        fig_height = 3,
        fig_ratios = "0.5 0.8 0.2 0.2"
    conda:
        config["environments"]["polars"]
    resources:
        resources = config.get("default_resources", "")
    shell:
        """
        mkdir -p $(dirname {output.plot})
        python {params.script} \
        --heatmap-files {input.heatmaps} \
        --heatmap-col "{params.heatmap_col}" \
        --glm {input.glm} \
        --target-pval "{params.target_pval}" \
        --genome-fasta {input.fasta} \
        --repeat-masker {input.repeats} \
        --chrom-map {input.chrom_map} \
        --bed-file {input.bed} \
        --output {output.plot} \
        --svg-output {output.plot_svg} \
        --json-output {output.json} \
        --percentile \
        --bin-size {params.bin_size} \
        --y-max {params.y_max} \
        --y-min {params.y_min} \
        --ribbon-alpha {params.ribbon_alpha} \
        --fst-top-pct {params.fst_top_pct} \
        --glm-top-pct {params.glm_top_pct} \
        --fig-width {params.fig_width} \
        --fig-height {params.fig_height} \
        --fig-ratios {params.fig_ratios} \
        --min-alpha {params.min_alpha} \
        --tracks heatmap glm overlap footer \
        --remove-repeats \
        --plot-bin-stats \
        --bin-stats-alpha 0.1 \
        --stats-lower-q 0.5 \
        --stats-upper-q 0.5 \
        > {log} 2>&1
        """

rule plot_divergence_sensitivity_mixture:
    input:
        candidates = "final_plots/fst_glm_mixture/final_aligned_landscape_{glm}_candidate_sites_detailed.tsv",
        af_data = "plotting_data/{glm}_raw_af.csv",
        metadata = "maps/sample_metadata.csv"
    output:
        out_dir = directory("final_plots/divergence_sensitivity_mixture_{glm}")
    log:
        "logs/divergence_sensitivity_mixture/{glm}.log"
    params:
        script = config["scripts"]["plot_divergence_sensitivity"],
        trt_1 = "B",
        trt_2 = "T",
        col_1 = "#2E86F0",
        col_2 = "#E8A60C",
        top_n = 100,
        ribbon_scale = 0.75,
        smoothness = 0.9,
        num_steps = 10,
        min_step = 100,
        white_fade_start = 0,
        white_fade_end = 1,
        max_dist = 1000000,
        divergence_ylim = 0.2,
        width = 5,
        height = 3,
        ratios = "1 1 1"
    conda:
        config["environments"]["polars"]
    resources:
        resources = config.get("default_resources", "")
    shell:
        """
        python {params.script} \
        --candidates {input.candidates} \
        --af-data {input.af_data} \
        --metadata {input.metadata} \
        --output-dir {output.out_dir} \
        --trt-1 {params.trt_1} \
        --trt-2 {params.trt_2} \
        --col-1 "{params.col_1}" \
        --col-2 "{params.col_2}" \
        --top-n {params.top_n} \
        --ribbon-scale {params.ribbon_scale} \
        --smoothness {params.smoothness} \
        --num-steps {params.num_steps} \
        --min-step {params.min_step} \
        --white-fade-start {params.white_fade_start} \
        --white-fade-end {params.white_fade_end} \
        --max-dist {params.max_dist} \
        --divergence-ylim {params.divergence_ylim} \
        --width {params.width} \
        --height {params.height} \
        --ratios {params.ratios} \
        > {log} 2>&1
        """

rule plot_bootstrap_zooms_mixture:
    input:
        lrt = "glm_bootstrap_final/{glm}.csv",
        heatmaps = expand("grenfst/divergence_mixture_lmm/{analysis}/divergence_stats.csv", 
                          analysis=[k for k in GRENFST_MULTIPLOT_DICT.keys() if "btwTB" in k]),
        bed = config["palaez"]
    output:
        out_dir = directory("final_plots/zoomed_plots_bootstrap_mixture/{glm}")
    log:
        "logs/bootstrap_zooms_mixture/{glm}.log"
    params:
        script = config["scripts"]["plot_lrt_zooms"],
        heatmap_col = "z_divergence",
        pval_col = "PB_p_val", 
        buffer = 2000000,
        smooth_type = "median",
        smooth_window = 10000,
        width = 3,
        height = 3,
        ratios = "1 1 0.1 0.1",
        point_size = 3,
        min_gene_width = 0.05
    conda:
        config["environments"]["polars"]
    resources:
        resources = config.get("default_resources", "")
    shell:
        """
        python3 {params.script} \
        --lrt-csv {input.lrt} \
        --pval-col {params.pval_col} \
        --heatmap-files {input.heatmaps} \
        --bed-file {input.bed} \
        --output-dir {output.out_dir} \
        --heatmap-col {params.heatmap_col} \
        --buffer {params.buffer} \
        --smooth-type {params.smooth_type} \
        --smooth-window {params.smooth_window} \
        --width {params.width} \
        --height {params.height} \
        --ratios {params.ratios} \
        --point-size {params.point_size} \
        --min-gene-width {params.min_gene_width} \
        > {log} 2>&1
        """