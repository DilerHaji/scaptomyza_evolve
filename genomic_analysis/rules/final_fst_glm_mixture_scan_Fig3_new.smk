_MANHATTAN_WIN = 200000

rule final_fst_glm_mixture_scan_Fig3_new:
    input:
        glm      = "glm_lrt_gw_final/{glm}.csv",
        heatmaps = lambda w: expand(
            "grenfst/divergence_mixture_lmm/{analysis}/"
            + ("divergence_stats_pbe2.csv" if w.stat == "pbe2" else "divergence_stats.csv"),
            analysis=[k for k in GRENFST_MULTIPLOT_DICT.keys() if "btwTB" in k]),
        gen_agg  = lambda w: next(
            f"grenfst/divergence_mixture_lmm/gen_agg/{k}.csv"
            for k in GRENFST_MULTIPLOT_DICT.keys()
            if "btwTB" in k and str(_MANHATTAN_WIN) in k),
        fasta     = config["reference_genome"],
        repeats   = "../popgen/repeat_content/masked_sfla_v2/sfla_v2.fa.out",
        chrom_map = "chrom_map.tsv",
        bed       = config["palaez"],
    output:
        plot       = "final_plots/fst_glm_mixture_new/final_aligned_landscape_{glm}_{stat}.png",
        svg        = "final_plots/fst_glm_mixture_new/final_aligned_landscape_{glm}_{stat}.svg",
        candidates = "final_plots/fst_glm_mixture_new/final_aligned_landscape_{glm}_{stat}_candidate_sites_detailed.tsv",
    benchmark:
        "benchmarks/window_glm_mixture_new/final_aligned_landscape_{glm}_{stat}.log"
    log:
        "logs/window_glm_mixture_new/final_aligned_landscape_{glm}_{stat}.log"
    params:
        script          = "scripts/window_glm_aligned_v3_pbs_manhat.py",
        heatmap_col     = lambda w: f"z_{w.stat}",
        heatmap_col_ref = lambda w: f"z_{w.stat}_ref",
        pbs_stat        = "fst_tb",
        target_pval     = "PB_p_val",
        bin_size        = 390000,
        manhattan_window = _MANHATTAN_WIN,
        ribbon_alpha    = 0.15,
        min_alpha       = 0.05,
        y_max           = 1.0,
        y_min           = 0.70,
        fst_top_pct     = 3,
        glm_top_pct     = 0.1,
        fig_width            = 12,
        fig_height           = 3.0,
        fig_ratios           = "0.5 0.2 0.8 0.2 0.15",
    conda:
        config["environments"]["polars"]
    resources:
        resources = config.get("default_resources", "")
    shell:
        """
        mkdir -p $(dirname {output.plot})

        python {params.script} \
            --heatmap-files   {input.heatmaps} \
            --heatmap-col     "{params.heatmap_col}" \
            --heatmap-col-ref "{params.heatmap_col_ref}" \
            --glm             {input.glm} \
            --target-pval     "{params.target_pval}" \
            --gen-agg-file    {input.gen_agg} \
            --pbs-stat        {params.pbs_stat} \
            --pbs-gen1        1 \
            --pbs-gen2        10 \
            --genome-fasta    {input.fasta} \
            --repeat-masker   {input.repeats} \
            --chrom-map       {input.chrom_map} \
            --bed-file        {input.bed} \
            --output          {output.plot} \
            --percentile \
            --bin-size        {params.bin_size} \
            --y-max           {params.y_max} \
            --y-min           {params.y_min} \
            --ribbon-alpha    {params.ribbon_alpha} \
            --fst-top-pct     {params.fst_top_pct} \
            --glm-top-pct     {params.glm_top_pct} \
            --fig-width       {params.fig_width} \
            --fig-height      {params.fig_height} \
            --fig-ratios      {params.fig_ratios} \
            --min-alpha       {params.min_alpha} \
            --tracks heatmap pbs_manhattan glm overlap footer \
            --remove-repeats \
        > {log} 2>&1

        python {params.script} \
            --heatmap-files   {input.heatmaps} \
            --heatmap-col     "{params.heatmap_col}" \
            --heatmap-col-ref "{params.heatmap_col_ref}" \
            --glm             {input.glm} \
            --target-pval     "{params.target_pval}" \
            --gen-agg-file    {input.gen_agg} \
            --pbs-stat        {params.pbs_stat} \
            --pbs-gen1        1 \
            --pbs-gen2        10 \
            --genome-fasta    {input.fasta} \
            --repeat-masker   {input.repeats} \
            --chrom-map       {input.chrom_map} \
            --bed-file        {input.bed} \
            --output          {output.svg} \
            --percentile \
            --bin-size        {params.bin_size} \
            --y-max           {params.y_max} \
            --y-min           {params.y_min} \
            --ribbon-alpha    {params.ribbon_alpha} \
            --fst-top-pct     {params.fst_top_pct} \
            --glm-top-pct     {params.glm_top_pct} \
            --fig-width       {params.fig_width} \
            --fig-height      {params.fig_height} \
            --fig-ratios      {params.fig_ratios} \
            --min-alpha       {params.min_alpha} \
            --tracks heatmap pbs_manhattan glm overlap footer \
            --remove-repeats \
        >> {log} 2>&1
        """
