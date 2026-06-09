rule subset_poolsnp:
    input:
    	combined = "variants/poolsnp_ncorrect/{poolsnp}.csv",
    output:
        "pfilter/tmp/{poolsnp}/{entry}.csv"
    params: 
        subset_ncorrect = config["scripts"]["subset_ncorrect"]
    conda:
        config["environments"]["polars"]
    benchmark: 
        "benchmarks/subset_poolsnp/{poolsnp}/{entry}.log"
    log: 
        "logs/subset_poolsnp/{poolsnp}/{entry}.log"
    resources:
        resources=config["default_resources_5cpus"]
    shell:
        """
        python {params.subset_ncorrect} {input.combined} {output} {wildcards.entry} Bergland2014
            
        """


# rule subset_poolsnp:
#     input:
#     	combined = "variants/poolsnp_ncorrect/{poolsnp}.csv",
#       chrom_batch1 = "chrom_batch1.txt"
#     output:
#         temp("pfilter/tmp/{poolsnp}/{chrom}.csv")
#     benchmark: 
#         "benchmarks/subset_poolsnp/{poolsnp}/{chrom}.log"
#     log: 
#         "logs/subset_poolsnp/{poolsnp}/{chrom}.log"
#     resources:
#         resources=config["default_resources"]
#     shell:
#         """
#         
#         awk -v chr={wildcards.chrom} -F, 'NR == 1 || $1 == chr' {input.combined} > {output} 2> {log}
#             
#         """


print(POOLSNP_DICT)

# This rule needs a lot of memory. Using high CPU resources to allow for greater memory on savio
rule pfilter:
    input:
    	"pfilter/tmp/{poolsnp}/{entry}.csv"
    output:
    	occurance_filtered = "pfilter/{poolsnp}/{entry}_filtered_variants.csv",
    params:
        script = config["scripts"]["pfilter"],
        output_prefix = lambda wildcards: os.path.join("pfilter", wildcards.poolsnp, wildcards.entry),
        freq_column = config["pfilter"]["freq_column"], 
#         freq_cutoff_lower = config["pfilter"]["freq_cutoff_lower"],
#         freq_cutoff_upper = config["pfilter"]["freq_cutoff_upper"],
#         min_occurance_gen = config["pfilter"]["min_occurance_gen"], 
#         min_occurance_pop = config["pfilter"]["min_occurance_pop"], 
        freq_cutoff_lower = lambda wildcards: POOLSNP_DICT[wildcards.poolsnp][2][0],
        freq_cutoff_upper = lambda wildcards: POOLSNP_DICT[wildcards.poolsnp][2][1],
        min_occurance_gen = lambda wildcards: POOLSNP_DICT[wildcards.poolsnp][2][2], 
        min_occurance_pop = lambda wildcards: POOLSNP_DICT[wildcards.poolsnp][2][3],
    benchmark: 
        "benchmarks/pfilter/{poolsnp}/{entry}.log"
    conda: 
        config["environments"]["polars"]
    log: 
        "logs/pfilter/{poolsnp}/{entry}.log"
    resources:
        resources=config["default_resources_24cpus"]
    shell:
        """
        
        python {params.script} {input} {params.freq_column} {params.freq_cutoff_lower} {params.freq_cutoff_upper} {params.min_occurance_gen} {params.min_occurance_pop} {params.output_prefix} > {log} 2>&1
        
        """

def gather_pfilter(wildcards):
    ENTRY = get_entries()
    return [f"pfilter/{variants}/{entry}_filtered_variants.csv" 
            for variants in VARIANTS_DICT[wildcards.variants] 
            for entry in ENTRY] 

rule pfilter_gather:
    input:
        gather_pfilter
    output:
        "{variants}_final_variants.csv"
    params: 
        final_variants_agg=config["scripts"]["final_variants_agg"]
    benchmark: 
        "benchmarks/pfilter_gather/{variants}.log"
    log: 
        "logs/pfilter_gather/{variants}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        
        bash {params.final_variants_agg} {output} {input}
     
        """

# rule pfilter_aggregate:
#     input:
#         "{variants}_final_variants.csv"
#     output:
#         "{variants}_final_variants_agg.csv"
#     params: 
#         aggregate_final_variants=config["scripts"]["aggregate_final_variants"]
#     benchmark: 
#         "benchmarks/pfilter_aggregate/{variants}.log"
#     conda: 
#         config["environments"]["polars"]
#     log: 
#         "logs/pfilter_aggregate/{variants}.log"
#     resources:
#         resources=config["default_resources"]
#     shell:
#         """
#         
#         python {params.aggregate_final_variants} {input} {output}
#      
#         """

rule final_variants:
    input:
        gather_pfilter
    output:
        "{variants}_final_variants.txt"
    params: 
        final_variants=config["scripts"]["final_variants"]
    benchmark: 
        "benchmarks/variants/{variants}_final_variants.log"
    log: 
        "logs/variants/{variants}_final_variants.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        
        bash {params.final_variants} {output} {input}
         
        """