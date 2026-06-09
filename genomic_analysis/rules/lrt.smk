import os

BOOTSTRAP_CHUNK_SIZE = 40
BOOTSTRAP_NSIMS = 1000   
JULIA_BIN = os.path.expanduser("~/software/julia-1.10.4/bin/julia")

CHUNK_DIR = "glm/bootstrap_chunks"


checkpoint pglm_split_candidates:
    input:
        candidates = "final_aligned_landscape_v4_candidate_sites_detailed.tsv"
    output:
        dir = directory(CHUNK_DIR)
    params:
        chunk_size = BOOTSTRAP_CHUNK_SIZE
    conda:
        config["environments"]["polars"]
    resources:
        resources=config.get("default_resources", "")
    log:
        "logs/glm_bootstrap/split.log"
    script:
        "../scripts/split_candidates_polars.py"


rule pglm_bootstrap_worker:
    input:
        rdata = "glm/haf/{glm}/merged_gethaf.Rdata",
        rds   = "glm/haf/{glm}/merged_gethaf_neff.RDS",
        chunk = os.path.join(CHUNK_DIR, "chunk_{chunk_id}.tsv")
    output:
        csv = "glm/bootstrap/{glm}/chunks/res_{chunk_id}.csv"
    params:
        r_script = "scripts/convert_to_arrow.R",
        jl_script = "scripts/LRT_bootstrap_worker.jl",
        repName = config["glm"]["repName"],
        treatments = "B,T",
        nSims = BOOTSTRAP_NSIMS,
        julia = JULIA_BIN
    conda:
        config["environments"]["r_arrow"]
    resources:
        resources=config["default_resources_10cpus"]
    threads: 10 
    log:
        "logs/glm_bootstrap/{glm}/{chunk_id}.log"
    shell:
        """
        TMP_ARROW=$(mktemp --tmpdir=. --suffix=.arrow)
        
        echo "[$(date)] Starting R Conversion..." > {log}
        
        Rscript {params.r_script} \
            --haf {input.rdata} \
            --rd {input.rds} \
            --chunk {input.chunk} \
            --repName {params.repName} \
            --selectTrts {params.treatments} \
            --out $TMP_ARROW >> {log} 2>&1
            
        
        export JULIA_NUM_THREADS={threads}
        
        {params.julia} {params.jl_script} \
            --arrow $TMP_ARROW \
            --nSims {params.nSims} \
            --out {output.csv} >> {log} 2>&1
            
        rm -f $TMP_ARROW
        """

def aggregate_bootstrap_inputs(wildcards):
    checkpoint_output = checkpoints.pglm_split_candidates.get(**wildcards).output.dir
    chunk_files = glob_wildcards(os.path.join(checkpoint_output, "chunk_{chunk_id}.tsv")).chunk_id
    return expand("glm/bootstrap/{glm}/chunks/res_{chunk_id}.csv", 
                  glm=wildcards.glm, 
                  chunk_id=chunk_files)



rule pglm_bootstrap_final:
    input:
        chunks = aggregate_bootstrap_inputs,
        candidates = "final_aligned_landscape_v4_candidate_sites_detailed.tsv"
    output:
        final_csv = "glm_bootstrap_final/{glm}.csv"
    conda:
        config["environments"]["polars"]
    resources:
        resources=config.get("default_resources", "")
    log:
        "logs/glm_bootstrap_final/{glm}.log"
    script:
        "../scripts/merge_bootstrap_polars.py"




rule plot_bootstrap_zooms_lrt:
    input:
        lrt = "glm_bootstrap_final/{glm}.csv",
        heatmaps = expand("grenfst/lmm_results/{analysis}/consensus_lmm.csv", 
                          analysis=[k for k in GRENFST_MULTIPLOT_DICT.keys() if "btwTB" in k]),
        bed = config["palaez"]
    output:
        out_dir = directory("final_plots/zoomed_plots_bootstrap_lrt/{glm}")
    log:
        "logs/plot_bootstrap_zooms_lrt/{glm}.log"
    params:
        script = config["scripts"]["plot_lrt_zooms"],
        heatmap_col = "z_score_median",
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