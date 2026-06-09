import os


JULIA_BIN = os.path.expanduser("~/software/julia-1.10.4/bin/julia")
LRT_NSIMS = 200  


rule pglm_lrt_gw_worker:
    input:
        rdata = "glm/haf/{glm}/gethaf_{entry}.Rdata",
        rds   = "glm/haf/{glm}/gethaf_{entry}_neff.RDS",
        sites = "glm/preprocessing/{glm}/pre_{entry}_sites.csv"
    output:
        csv = "glm/lrt_gw/{glm}/res_{entry}.csv"
    params:
        r_script = "scripts/convert_to_arrow.R",
        jl_script = "scripts/LRT_bootstrap_worker.jl",
        repName = config["glm"]["repName"],
        treatments = "B,T",
        nSims = LRT_NSIMS,
        julia = JULIA_BIN,
        min_coverage = 15,      # Stricter
        min_samples = 6,         #  More samples required
        min_distinct_freqs = 3   #  More variation
    conda:
        config["environments"]["r_arrow"]
    resources:
        resources = config["default_resources_10cpus"]
    threads: 20 
    benchmark:
        "benchmarks/glm_lrt_gw/{glm}/{entry}.log"
    log:
        "logs/glm_lrt_gw/{glm}/{entry}.log"
    shell:
        """
        if [ ! -s {input.rdata} ]; then
            echo "[$(date)] Input RData ({input.rdata}) is empty. Skipping bin {wildcards.entry}." > {log}
            touch {output.csv}
            exit 0
        fi
        TMP_ARROW=$(mktemp --tmpdir=. --suffix=.arrow)


        Rscript {params.r_script} \
            --haf {input.rdata} \
            --rd {input.rds} \
            --chunk {input.sites} \
            --repName {params.repName} \
            --selectTrts {params.treatments} \
            --min-coverage {params.min_coverage} \
            --min-samples {params.min_samples} \
            --min-distinct-freqs {params.min_distinct_freqs} \
            --out $TMP_ARROW >> {log} 2>&1


        export JULIA_NUM_THREADS={threads}

        {params.julia} {params.jl_script} \
            --arrow $TMP_ARROW \
            --nSims {params.nSims} \
            --out {output.csv} >> {log} 2>&1


        rm -f $TMP_ARROW
        """


def gather_lrt_gw(wildcards):
    entries = get_entries()
    return [f"glm/lrt_gw/{wildcards.glm}/res_{entry}.csv" for entry in entries]

rule pglm_lrt_gw_final:
    input:
        csv = gather_lrt_gw
    output:
        glm = "glm_lrt_gw_final/{glm}.csv"
    benchmark:
        "benchmarks/glm_lrt_gw_final/{glm}.log"
    log:
        "logs/glm_lrt_gw_final/{glm}.log"
    resources:
        resources = config["default_resources"]
    shell:
        """
        header_written=false

        for file in {input.csv}; do
            if [[ -s "$file" ]]; then
                head -n 1 "$file" > {output.glm}
                header_written=true
                break
            fi
        done

        for file in {input.csv}; do
            if [[ -s "$file" ]]; then
                tail -n +2 "$file" >> {output.glm}
            fi
        done
        """


rule final_fst_lrt_scan_genome_wide:
    input:
        glm = "glm_lrt_gw_final/{glm}.csv",
        heatmaps = expand("grenfst/lmm_results/{analysis}/consensus_lmm.csv",
                          analysis=[k for k in GRENFST_MULTIPLOT_DICT.keys() if "btwTB" in k]),
        fasta = config["reference_genome"],
        repeats = "../popgen/repeat_content/masked_sfla_v2/sfla_v2.fa.out",
        chrom_map = "chrom_map.tsv",
        bed = config["palaez"]
    output:
        plot = "final_plots/fst_lrt_gw/final_aligned_landscape_{glm}.png",
        candidates = "final_plots/fst_lrt_gw/final_aligned_landscape_{glm}_candidate_sites.tsv"
    benchmark:
        "benchmarks/window_lrt_gw/final_aligned_landscape_{glm}.log"
    log:
        "logs/window_lrt_gw/final_aligned_landscape_{glm}.log"
    params:
        script = config["scripts"]["window_glm_aligned_v3"],
        target_pval = "PB_p_val",
        bin_size = 390000,
        ribbon_alpha = 0.15,
        min_alpha = 0.05,
        y_max = 1.0,
        y_min = 0.0,
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




rule final_fst_heatmap_only:
    input:
        heatmaps = expand("grenfst/lmm_results/{analysis}/consensus_lmm.csv",
                          analysis=[k for k in GRENFST_MULTIPLOT_DICT.keys() if "btwTB" in k]),
        fasta = config["reference_genome"],
        repeats = "../popgen/repeat_content/masked_sfla_v2/sfla_v2.fa.out",
        chrom_map = "chrom_map.tsv",
        bed = config["palaez"]
    output:
        plot = "final_plots/fst_lrt_gw/final_aligned_landscape_heatmap_only.png"
    benchmark:
        "benchmarks/window_lrt_gw/final_aligned_landscape_heatmap_only.log"
    log:
        "logs/window_lrt_gw/final_aligned_landscape_heatmap_only.log"
    params:
        script = config["scripts"]["window_glm_aligned_v3"],
        bin_size = 390000,
        y_max = 1.0, 
        y_min = 0.0,
        fig_width = 12,
        fig_height = 3,
        fig_ratios = "0.7 0.15 0.1 0.05",
        tracks = "heatmap repeats genes footer"
    conda:
        config["environments"]["polars"]
    resources:
        resources = config.get("default_resources", "")
    shell:
        """
        python {params.script} \
        --heatmap-files {input.heatmaps} \
        --genome-fasta {input.fasta} \
        --repeat-masker {input.repeats} \
        --chrom-map {input.chrom_map} \
        --bed-file {input.bed} \
        --output {output.plot} \
        --percentile \
        --bin-size {params.bin_size} \
        --y-max {params.y_max} \
        --y-min {params.y_min} \
        --fig-width {params.fig_width} \
        --fig-height {params.fig_height} \
        --fig-ratios {params.fig_ratios} \
        --tracks {params.tracks} \
        --remove-repeats \
        # Only pass GLM args if they existed, here we intentionally omit --glm and --target-pval
        """
