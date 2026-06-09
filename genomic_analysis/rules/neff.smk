# poolsnp
# poolsnp_correct

rule poolsnp_biallelic:
    input:
        "variants/poolsnp/{poolsnp}.variants"
    output:
        temp("variants/poolsnp_biallelic/{poolsnp}.variants"),
    benchmark: 
        "benchmarks/poolsnp_biallelic/{poolsnp}.log"
    log: 
        "logs/poolsnp_biallelic/{poolsnp}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        
        awk -F '\t' '!($5 ~ /,/) {{print}}' {input} > {output}
                 
        """


rule poolsnp_correct_pool:
    input:
        variants = "variants/poolsnp_biallelic/{poolsnp}.variants",
        header = "variants/poolsnp/{poolsnp}.header"
    output:
        temp("variants/poolsnp_correct_pool/{poolsnp}/{poolsnp_correct}.ncorrect"),
    params: 
        calculate_neff_poolsnp = config["scripts"]["calculate_neff_poolsnp"],
        poolsizes = config["poolsizes"],
        #col = lambda wildcards: POOLSNP_DICT[wildcards.poolsnp].index(wildcards.poolsnp_correct) + 9
    benchmark: 
        "benchmarks/poolsnp_correct_pool/{poolsnp}/{poolsnp_correct}.log"
    log: 
        "logs/poolsnp_correct_pool/{poolsnp}/{poolsnp_correct}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        python {params.calculate_neff_poolsnp} \
            --input {input.variants} \
            --header {input.header} \
            --poolsizes {params.poolsizes} \
            --pool {wildcards.poolsnp_correct} \
            --output {output} \
            2> {log}
                         
        """



def get_correct_pool(wildcards):
    return [f"variants/poolsnp_correct_pool/{wildcards.poolsnp}/" + item + ".ncorrect" for item in PILEUP_DICT[POOLSNP_DICT[wildcards.poolsnp][0]]]


# PASSED ONTO PFILTER.SMK
rule poolsnp_correct:
    input:
        csv = get_correct_pool,
        sample_metadata=config["sample_metadata"]
    output:
        combined = "variants/poolsnp_ncorrect/{poolsnp}.csv",
        chrom = "variants/poolsnp_ncorrect/chrom_{poolsnp}.csv"
    params: 
        prefix = lambda wildcards: "variants/poolsnp_ncorrect/" + wildcards.poolsnp,
        combine_freq = config["scripts"]["combine_freq_poolsnp"]
    benchmark: 
        "benchmarks/poolsnp_correct/{poolsnp}.log"
    log: 
        "logs/poolsnp_correct/{poolsnp}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        
        python {params.combine_freq} {input.csv} --output {output.combined} --metadata {input.sample_metadata}
        
        cat {output.combined} | cut -d ',' -f1 | sort | uniq | grep -v "^CHROM" > {output.chrom}
                                 
        """



