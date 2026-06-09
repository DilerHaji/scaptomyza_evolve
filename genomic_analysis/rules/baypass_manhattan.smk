rule baypass_manhattan:
    input:
        snp_pos  = "baypass/g10_snp_positions.csv",
        xtx      = "baypass/g10_core_summary_pi_xtx.out",
        betai    = "baypass/g10_trt_summary_betai_reg.out",
        contrast = "baypass/g10_contrast_summary_contrast.out",
    output:
        png = "final_plots/baypass/baypass_manhattan.png",
        svg = "final_plots/baypass/baypass_manhattan.svg",
    log:
        "logs/baypass/baypass_manhattan.log"
    params:
        script = "scripts/baypass_manhattan.py",
        contrast_names = "B vs T,B vs M,T vs M",
    conda:
        config["environments"]["polars"]
    resources:
        resources = config.get("default_resources", "")
    shell:
        """
        mkdir -p $(dirname {output.png})
        python {params.script} \
            --snp-pos        {input.snp_pos} \
            --xtx            {input.xtx} \
            --betai          {input.betai} \
            --contrast        {input.contrast} \
            --contrast-names "{params.contrast_names}" \
            --output         {output.png} \
            --output-svg     {output.svg} \
        > {log} 2>&1
        """


rule baypass_heatmap_manhattan:
    input:
        snp_pos  = "baypass/g10_snp_positions.csv",
        xtx      = "baypass/g10_core_summary_pi_xtx.out",
        betai    = "baypass/g10_trt_summary_betai_reg.out",
        contrast = "baypass/g10_contrast_summary_contrast.out",
    output:
        png = "final_plots/baypass/baypass_heatmap_manhattan.png",
        svg = "final_plots/baypass/baypass_heatmap_manhattan.svg",
    log:
        "logs/baypass/baypass_heatmap_manhattan.log"
    params:
        script = "scripts/baypass_heatmap_manhattan.py",
        contrast_names = "B vs T,B vs M,T vs M",
        sweep_sizes = "50000,100000,200000,500000,1000000",
    conda:
        config["environments"]["polars"]
    resources:
        resources = config.get("default_resources", "")
    shell:
        """
        mkdir -p $(dirname {output.png})
        python {params.script} \
            --snp-pos        {input.snp_pos} \
            --xtx            {input.xtx} \
            --betai          {input.betai} \
            --contrast        {input.contrast} \
            --contrast-names "{params.contrast_names}" \
            --sweep-sizes    "{params.sweep_sizes}" \
            --output         {output.png} \
            --output-svg     {output.svg} \
        > {log} 2>&1
        """


_BAYPASS_PBS_WINDOWS = [50000, 110000, 200000, 500000, 1000000, 2000000]
_BAYPASS_PBS_DATADIR = "grenfst/divergence_mixture_lmm"

rule baypass_pbs_aligned:
    input:
        divergence = [f"{_BAYPASS_PBS_DATADIR}/e10Ffe9wGREN_btwTB_{w}_1000/divergence_stats.csv"
                      for w in _BAYPASS_PBS_WINDOWS],
        fasta      = config["reference_genome"],
        repeat_gff = "repeat_masker/sfla_v2.fa.gff",
        snp_pos    = "baypass/g10_snp_positions.csv",
        xtx        = "baypass/g10_core_summary_pi_xtx.out",
        betai      = "baypass/g10_trt_summary_betai_reg.out",
        contrast   = "baypass/g10_contrast_summary_contrast.out",
    output:
        png = "final_plots/baypass/baypass_pbs_aligned.png",
        svg = "final_plots/baypass/baypass_pbs_aligned.svg",
    log:
        "logs/baypass/baypass_pbs_aligned.log"
    params:
        script = "scripts/baypass_pbs_aligned.py",
        contrast_names = "B vs T,B vs M,T vs M",
        stat = "pbs",
    conda:
        config["environments"]["polars"]
    resources:
        resources = config.get("default_resources", "")
    shell:
        """
        mkdir -p $(dirname {output.png})
        python {params.script} \
            --divergence-files {input.divergence} \
            --genome-fasta     {input.fasta} \
            --repeat-gff       {input.repeat_gff} \
            --stat             {params.stat} \
            --snp-pos          {input.snp_pos} \
            --xtx              {input.xtx} \
            --betai            {input.betai} \
            --contrast         {input.contrast} \
            --contrast-names  "{params.contrast_names}" \
            --output           {output.png} \
            --output-svg       {output.svg} \
        > {log} 2>&1
        """


rule baypass_wild_manhattan:
    input:
        snp_pos  = "baypass_wild/wild_snp_positions.csv",
        betai    = "baypass_wild/wild_trt_summary_betai_reg.out",
        contrast = "baypass_wild/wild_contrast_summary_contrast.out",
    output:
        png = "final_plots/baypass_wild/baypass_wild_manhattan.png",
        svg = "final_plots/baypass_wild/baypass_wild_manhattan.svg",
    log:
        "logs/baypass/baypass_wild_manhattan.log"
    params:
        script         = "scripts/baypass_manhattan.py",
        contrast_names = "B vs T",
    conda:
        config["environments"]["polars"]
    resources:
        resources = config.get("default_resources", "")
    shell:
        """
        mkdir -p $(dirname {output.png})
        python {params.script} \
            --snp-pos        {input.snp_pos} \
            --betai          {input.betai} \
            --contrast       {input.contrast} \
            --contrast-names "{params.contrast_names}" \
            --output         {output.png} \
            --output-svg     {output.svg} \
        > {log} 2>&1
        """


rule all_baypass_manhattan:
    input:
        "final_plots/baypass/baypass_manhattan.png",
        "final_plots/baypass/baypass_heatmap_manhattan.png",
        "final_plots/baypass/baypass_pbs_aligned.png",
