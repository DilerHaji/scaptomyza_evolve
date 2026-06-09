rule get_poolsnp:
    output:
        "poolsnp/PoolSNP.sh"
    params:
        poolsnp = config["programs"]["poolsnp"],
    benchmark: 
        "benchmarks/get_poolsnp.log"
    log: 
        "logs/get_poolsnp.log"
    shell:
        """
        
        [[ -d poolsnp ]] || mkdir poolsnp
        
        cp {params.poolsnp}/* poolsnp/ -r
		                          
        """

localrules: get_poolsnp

print(PILEUP_DICT)


def get_poolsnp_pileup(wildcards):
    items = PILEUP_DICT[wildcards.pileup]
    if "merged" in wildcards.pileup:
        return expand("msorted/{sample}.bam", sample=items)
    else:
        return expand("sorted/{sample}.bam", sample=items)

def get_poolsnp_pileup_bai(wildcards):
    items = PILEUP_DICT[wildcards.pileup]
    if "merged" in wildcards.pileup:
        return expand("msorted/{sample}.bam.bai", sample=items)
    else:
        return expand("sorted/{sample}.bam.bai", sample=items)




# def get_poolsnp_pileup(wildcards):
#     return ["sorted/" + item + ".bam" for item in PILEUP_DICT[wildcards.pileup]]
# 
# def get_poolsnp_pileup_bai(wildcards):
#     return ["sorted/" + item + ".bam.bai" for item in PILEUP_DICT[wildcards.pileup]]

# def get_poolsnp_pileup(wildcards):
#     return [ "freebayes/bams_sorted/" + item + ".sorted.bam" for item in PILEUP_DICT[wildcards.pileup]]
# 
# def get_poolsnp_pileup_bai(wildcards):
#     return [ "freebayes/bams_indexed/" + item + ".done" for item in PILEUP_DICT[wildcards.pileup]]
# 

rule get_mpileup_poolsnp:
    input:
        bam = get_poolsnp_pileup,
        bai = get_poolsnp_pileup_bai,
        ref=config["reference_genome"],
    output:
        pileup = temp("poolsnp_pileup/{pileup}/{entry}.mpileup"),
        pileup_done = "poolsnp_pileup/{pileup}/{entry}.done"
    params:
        region = lambda wildcards: wildcards.entry,
        q=config["mpileup"]["min_map_quality"],
        Q=config["mpileup"]["min_quality"],
        R=config["mpileup"]["rg_tag"],
    conda: 
        config["environments"]["samtools"]        
    benchmark: 
        "benchmarks/get_poolsnp_mpileup/{pileup}/{entry}.log"
    log: 
        "logs/get_poolsnp_mpileup/{pileup}/{entry}.log"
    resources:
        resources=config["default_resources_10cpus"]
    shell:
        """
        
        [[ -d poolsnp_pileup ]] || mkdir poolsnp_pileup
		[[ -d poolsnp_pileup/{wildcards.pileup} ]] || mkdir poolsnp_pileup/{wildcards.pileup}
        
        samtools mpileup -f {input.ref} -r {params.region} {input.bam} -Q {params.Q} -q {params.q} {params.R} \
        > poolsnp_pileup/{wildcards.pileup}/{wildcards.entry}.mpileup
        
        head -n 5000 poolsnp_pileup/{wildcards.pileup}/{wildcards.entry}.mpileup > poolsnp_pileup/{wildcards.pileup}/sub_{wildcards.entry}.mpileup
        
        touch {output.pileup_done}
        
        """

def get_pileup(wildcards):
    entries = get_entries()
    return [f"poolsnp_pileup/{wildcards.pileup}/{entry}.mpileup" for entry in entries]

def get_pileup_done(wildcards):
    entries = get_entries()
    return [f"poolsnp_pileup/{wildcards.pileup}/{entry}.done" for entry in entries]


rule gather_pileup:
    input:
        pileups = get_pileup,
        pileups_done = get_pileup_done
    output:
        pileup = "gather_pileup/{pileup}.mpileup",
    benchmark: 
        "benchmarks/gather_pileup/{pileup}.log"
    log: 
        "logs/gather_pileup/{pileup}.log"
    resources:
        resources=config["default_resources_10cpus"]
    shell:
        """        
        cat {input.pileups} > {output.pileup}
        """

