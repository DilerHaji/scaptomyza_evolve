_WILD_SAMPLES = ["AVB", "AVT", "PSB", "PST", "RMB", "RMT"]
_WILD_VCF     = "fvariants/fvOG_e10fe9w.fixed_no_neff.vcf.gz"
_WILD_PREFIX  = "wild"
_WILD_OUTDIR  = "baypass_wild"
_WILD_BAYPASS = "baypass/sources/g_baypass"


rule baypass_wild_positions_bed:
    input:
        vcf = _WILD_VCF,
    output:
        bed = "baypass_wild/variant_positions.bed",
    log:
        "logs/baypass_wild/positions_bed.log",
    conda:
        config["environments"]["bcftools"]
    resources:
        resources = config["default_resources"]
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools view -H {input.vcf} \
          | awk '{{print $1"\t"($2-1)"\t"$2}}' \
          > {output.bed} 2> {log}
        echo "Sites: $(wc -l < {output.bed})" >> {log}
        """


rule baypass_wild_mpileup:
    input:
        bam = config["input_dir"] + "/{wild_sample}.bam",
        bed = "baypass_wild/variant_positions.bed",
        ref = config["reference_genome"],
    output:
        mp  = "baypass_wild/pileups/{wild_sample}.mpileup",
    wildcard_constraints:
        wild_sample = "|".join(_WILD_SAMPLES),
    log:
        "logs/baypass_wild/mpileup_{wild_sample}.log",
    conda:
        config["environments"]["samtools"]
    resources:
        resources = config["default_resources"]
    shell:
        """
        mkdir -p $(dirname {log})
        samtools mpileup \
            -q 25 -Q 25 \
            -B --ignore-RG \
            -l {input.bed} \
            -f {input.ref} \
            {input.bam} \
            > {output.mp} 2> {log}
        echo "Lines: $(wc -l < {output.mp})" >> {log}
        """


rule baypass_wild_prepare:
    input:
        vcf     = _WILD_VCF,
        pileups = expand("baypass_wild/pileups/{s}.mpileup", s=_WILD_SAMPLES),
    output:
        geno       = "baypass_wild/wild_pooldata.geno",
        omega_geno = "baypass_wild/wild_omega_pooldata.geno",
        poolsize   = "baypass_wild/wild_poolsize.txt",
        cov        = "baypass_wild/wild_treatment.cov",
        contrast   = "baypass_wild/wild_contrasts.con",
        snp_pos    = "baypass_wild/wild_snp_positions.csv",
    log:
        "logs/baypass_wild/prepare.log",
    params:
        script     = "scripts/prepare_baypass_wild.py",
        mpileup_dir = "baypass_wild/pileups",
        output_dir  = "baypass_wild",
        prefix      = _WILD_PREFIX,
        min_cov     = 5,
        min_pop_cov = 1.0,
        thin_step   = 16,
    conda:
        config["environments"]["polars"]
    resources:
        resources = config["default_resources"]
    shell:
        """
        python {params.script} \
            --vcf         {input.vcf} \
            --mpileup-dir {params.mpileup_dir} \
            --output-dir  {params.output_dir} \
            --prefix      {params.prefix} \
            --min-cov     {params.min_cov} \
            --min-pop-cov {params.min_pop_cov} \
            --thin-step   {params.thin_step} \
        > {log} 2>&1
        """


rule baypass_wild_omega:
    input:
        geno     = "baypass_wild/wild_omega_pooldata.geno",
        poolsize = "baypass_wild/wild_poolsize.txt",
    output:
        omega = "baypass_wild/wild_omega_mat_omega.out",
        xtx   = "baypass_wild/wild_omega_summary_pi_xtx.out",
    log:
        "logs/baypass_wild/omega.log",
    params:
        baypass    = _WILD_BAYPASS,
        outprefix  = "baypass_wild/wild_omega",
        nthreads   = 8,
    resources:
        resources = config["default_resources_10cpus"]
    shell:
        """
        {params.baypass} \
            -pooldatafile {input.geno} \
            -poolsizefile {input.poolsize} \
            -outprefix    {params.outprefix} \
            -nthreads     {params.nthreads} \
        > {log} 2>&1
        """


rule baypass_wild_covariate:
    input:
        geno     = "baypass_wild/wild_pooldata.geno",
        poolsize = "baypass_wild/wild_poolsize.txt",
        omega    = "baypass_wild/wild_omega_mat_omega.out",
        cov      = "baypass_wild/wild_treatment.cov",
    output:
        betai = "baypass_wild/wild_trt_summary_betai_reg.out",
    log:
        "logs/baypass_wild/covariate.log",
    params:
        baypass   = _WILD_BAYPASS,
        outprefix = "baypass_wild/wild_trt",
        nthreads  = 8,
    resources:
        resources = config["default_resources_10cpus"]
    shell:
        """
        {params.baypass} \
            -pooldatafile {input.geno} \
            -poolsizefile {input.poolsize} \
            -omegafile    {input.omega} \
            -efile        {input.cov} \
            -outprefix    {params.outprefix} \
            -nthreads     {params.nthreads} \
        > {log} 2>&1
        """


rule baypass_wild_contrast:
    input:
        geno     = "baypass_wild/wild_pooldata.geno",
        poolsize = "baypass_wild/wild_poolsize.txt",
        omega    = "baypass_wild/wild_omega_mat_omega.out",
        contrast = "baypass_wild/wild_contrasts.con",
    output:
        c2 = "baypass_wild/wild_contrast_summary_contrast.out",
    log:
        "logs/baypass_wild/contrast.log",
    params:
        baypass   = _WILD_BAYPASS,
        outprefix = "baypass_wild/wild_contrast",
        nthreads  = 8,
    resources:
        resources = config["default_resources_10cpus"]
    shell:
        """
        {params.baypass} \
            -pooldatafile {input.geno} \
            -poolsizefile {input.poolsize} \
            -omegafile    {input.omega} \
            -contrastfile {input.contrast} \
            -outprefix    {params.outprefix} \
            -nthreads     {params.nthreads} \
        > {log} 2>&1
        """


rule all_baypass_wild:
    input:
        "baypass_wild/wild_contrast_summary_contrast.out",
        "baypass_wild/wild_trt_summary_betai_reg.out",
        "final_plots/baypass_wild/baypass_wild_manhattan.png",
