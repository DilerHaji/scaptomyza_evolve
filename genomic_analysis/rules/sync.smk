rule get_sites_filtering:
    input:
        variants = expand("variants/poolsnp_biallelic/{poolsnp}.variants", poolsnp=POOLSNP)
    output:
        "sync/sites_for_filtering.txt",
    benchmark: 
        "benchmarks/get_sites_filtering.log"
    log: 
        "logs/get_sites_filtering.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        cat {input.variants} > all.variants

        awk '{{print $1"\t"$2}}' all.variants > all.sites
        sort -k1,1 -k2,2n all.sites | uniq > {output}

        """


#########################################################
# These pileups are used for SNAPE 
#########################################################

rule get_mpileup:
    input:
        bam="sorted/{bam}.bam",
        ref=config["reference_genome"],
    output:
        pileup = temp("bam_pileup/{bam}/{entry}.mpileup"),
        pileup_done = "bam_pileup/{bam}/{entry}.done"
    params:
        region = lambda wildcards: wildcards.entry,
        q=config["mpileup"]["min_map_quality"],
        Q=config["mpileup"]["min_quality"],
        R=config["mpileup"]["rg_tag"],
    conda: 
        config["environments"]["samtools"]        
    benchmark: 
        "benchmarks/get_mpileup/{bam}/{entry}.log"
    log: 
        "logs/get_mpileup/{bam}/{entry}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        
        [[ -d pileup ]] || mkdir pileup
		[[ -d pileup/{wildcards.bam} ]] || mkdir pileup/{wildcards.bam}
        
        samtools mpileup -f {input.ref} -r {params.region} {input.bam} -Q {params.Q} -q {params.q} {params.R} > bam_pileup/{wildcards.bam}/{wildcards.entry}.mpileup
        
        touch {output.pileup}
        
        """






