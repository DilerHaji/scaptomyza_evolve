def get_poolsnp_vcf(wildcards):
    entries = get_entries()
    return [f"poolsnp_sorted/{wildcards.poolsnp}/{entry}.sorted.vcf.gz" for entry in entries]

def get_poolsnp_vcf_index(wildcards):
    entries = get_entries()
    return [f"poolsnp_index2/index_{wildcards.poolsnp}/{entry}.done" for entry in entries]


rule aggregate_poolsnp:
    input:
        vcf=get_poolsnp_vcf,
        index_done=get_poolsnp_vcf_index
    output:
        combined="variants/poolsnp/{poolsnp}.combined.vcf.gz",
        #header="variants/poolsnp/{poolsnp}.header",
        #variants="variants/poolsnp/{poolsnp}.variants"
    benchmark: 
        "benchmarks/aggregate_poolsnp/{poolsnp}.log"
    log: 
        "logs/aggregate_poolsnp/{poolsnp}.log"
    resources:
        resources=config["default_resources"]
    conda:
        config["environments"]["bcftools"]
    shell:
        """
        
        bcftools concat -a -Oz -o {output.combined} {input.vcf} 2> {log}

        """


# PASSED ONTO NEFF.SMK 
rule aggregate_poolsnp_preprocessing:
    input:
        vcf="variants/poolsnp/{poolsnp}.combined.vcf.gz"
    output:
        vcf="variants/poolsnp/{poolsnp}.combined.vcf",
        header="variants/poolsnp/{poolsnp}.header",
        variants="variants/poolsnp/{poolsnp}.variants",
    benchmark: 
        "benchmarks/aggregate_poolsnp_preprocessing/{poolsnp}.log"
    log: 
        "logs/aggregate_poolsnp_preprocessing/{poolsnp}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        
        gunzip -c {input.vcf} > {output.vcf}
        grep "#" {output.vcf} > {output.header}
        grep -v "#" {output.vcf} > {output.variants}

        """


rule aggregate_poolsnp_index:
    input:
        "variants/poolsnp/{poolsnp}.combined.vcf.gz"
    output:
        index="variants_index/poolsnp/{poolsnp}.done"
    benchmark: 
        "benchmarks/aggregate_poolsnp_index/{poolsnp}.log"
    log: 
        "logs/aggregate_poolsnp_index/{poolsnp}.log"
    resources:
        resources=config["default_resources"]
    conda:
        config["environments"]["htslib"]
    shell:
        """
        tabix -f -p vcf {input}
        touch {output.index}
        """