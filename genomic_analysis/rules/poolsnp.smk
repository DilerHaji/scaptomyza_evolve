rule poolSNP:
    input: 
        poolsnp="poolsnp/PoolSNP.sh",
        pileup = lambda wildcards: "poolsnp_pileup/" + POOLSNP_DICT[wildcards.poolsnp][0] + "/{entry}.mpileup",
        pileup_done = lambda wildcards: "poolsnp_pileup/" + POOLSNP_DICT[wildcards.poolsnp][0] + "/{entry}.done",
        #pileup = "pileup/{poolsnp}/{entry}.mpileup",
        ref=config["reference_genome"]
    output:
        directory("poolsnp/{poolsnp}/{entry}"),
        "poolsnp/{poolsnp}/{entry}.vcf",
    params: 
        names = lambda wildcards: ",".join(PILEUP_DICT[POOLSNP_DICT[wildcards.poolsnp][0]]),
        pileupdir = lambda wildcards: POOLSNP_DICT[wildcards.poolsnp][0],
        maxcov = lambda wildcards: POOLSNP_DICT[wildcards.poolsnp][1][0],
        mincov = lambda wildcards: POOLSNP_DICT[wildcards.poolsnp][1][1],
        mincount = lambda wildcards: POOLSNP_DICT[wildcards.poolsnp][1][2],
        minfreq = lambda wildcards: POOLSNP_DICT[wildcards.poolsnp][1][3],
        missfrac = lambda wildcards: POOLSNP_DICT[wildcards.poolsnp][1][4],
        jobs = config["poolsnp"]["jobs"],
        badsites = config["poolsnp"]["badsites"],
        allsites = config["poolsnp"]["allsites"]
    conda: 
        config["environments"]["poolsnp"]
    benchmark: 
        "benchmarks/poolSNP/{poolsnp}/{entry}.log"
    log: 
        "logs/poolSNP/{poolsnp}/{entry}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        
        cd poolsnp/
        
        bash PoolSNP.sh \
        mpileup=../poolsnp_pileup/{params.pileupdir}/{wildcards.entry}.mpileup \
        reference={input.ref} \
        names={params.names} \
        max-cov={params.maxcov} \
        min-cov={params.mincov} \
        min-count={params.mincount} \
        min-freq={params.minfreq} \
        miss-frac={params.missfrac} \
        jobs={params.jobs} \
        badsites={params.badsites} \
        allsites={params.allsites} \
        output={wildcards.poolsnp}/{wildcards.entry}
        
        gunzip {wildcards.poolsnp}/{wildcards.entry}.vcf.gz
                          
        """
        


rule index_vcf:
    input:
        "poolsnp/{poolsnp}/{entry}.vcf"
    output:
        compressed="poolsnp_index/{poolsnp}/{entry}.vcf.gz",
        index="poolsnp_index/index_{poolsnp}/{entry}.done"
    benchmark: 
        "benchmarks/index_vcf/{poolsnp}_{entry}.log"
    log: 
        "logs/index_vcf/{poolsnp}_{entry}.log"
    resources:
        resources=config["default_resources"]
    conda:
        config["environments"]["htslib"]
    shell:
        """
        bgzip -c {input} > {output.compressed}
        tabix -p vcf {output.compressed}
        touch {output.index}
        """


rule sort_vcf:
    input:
        vcf="poolsnp_index/{poolsnp}/{entry}.vcf.gz",
        tbi="poolsnp_index/index_{poolsnp}/{entry}.done"
    output:
        "poolsnp_sorted/{poolsnp}/{entry}.sorted.vcf.gz"
    benchmark: 
        "benchmarks/sort_vcf/{poolsnp}_{entry}.log"
    log: 
        "logs/sort_vcf/{poolsnp}_{entry}.log"
    resources:
        resources=config["default_resources"]
    conda:
        config["environments"]["bcftools"]
    shell:
        "bcftools sort -Oz -o {output} {input.vcf}"


rule index2_vcf:
    input:
        "poolsnp_sorted/{poolsnp}/{entry}.sorted.vcf.gz"
    output:
        index="poolsnp_index2/index_{poolsnp}/{entry}.done"
    benchmark: 
        "benchmarks/index2_vcf/{poolsnp}_{entry}.log"
    log: 
        "logs/index2_vcf/{poolsnp}_{entry}.log"
    resources:
        resources=config["default_resources"]
    conda:
        config["environments"]["htslib"]
    shell:
        """
        tabix -f -p vcf {input}
        touch {output.index}
        """


