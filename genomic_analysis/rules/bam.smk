rule sort_bam:
    input: 
    	"realigned/{bam}.bam"
    output:
        "sorted/{bam}.bam"
    conda:  
        config["environments"]["samtools"]
    benchmark: 
        "benchmarks/sort_bam/{bam}.log"
    log: 
        "logs/sort_bam/{bam}.log"
    resources:
        resources=config["default_resources"]
    shell:
        "samtools sort -o {output} {input}"



rule index_bam:
    input: 
    	"sorted/{bam}.bam"
    output:
        "sorted/{bam}.bam.bai"
    conda:  
        config["environments"]["samtools"]
    benchmark: 
        "benchmarks/index_bam/{bam}.log"
    log: 
        "logs/index_bam/{bam}.log"
    resources:
        resources=config["default_resources"]
    shell:
        "samtools index {input} {output}"


rule merge_bams:
    input:
    	lambda wildcards: ["sorted/" + item + ".bam" for item in MERGE_DICT[wildcards.merge]],
    	#expand("sorted/{bam}.bam.bai", bam=BAM)
    output:
        "merged/{merge}.bam"
    params: 
        #samples_to_merge = lambda wildcards: [f'sorted/{item}.bam' for item in get_key(MERGE_DICT, wildcards.merge) ],
        extra = config["merge_bams"]["extra"]
    conda: 
        config["environments"]["samtools"]
    benchmark: 
        "benchmarks/merge_bams/{merge}.log"
    log: 
        "logs/merge_bams/{merge}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        
        samtools merge {params.extra} {output} {input}
        
        """



rule sort_merged_bam:
    input: 
    	"merged/{merge}.bam"
    output:
        "msorted/{merge}.bam"
    conda:  
        config["environments"]["samtools"]
    benchmark: 
        "benchmarks/sort_merged_bam/{merge}.log"
    log: 
        "logs/sort_merged_bam/{merge}.log"
    resources:
        resources=config["default_resources"]
    shell:
        "samtools sort -o {output} {input}"



rule index_merged_bam:
    input: 
    	"msorted/{merge}.bam"
    output:
        "msorted/{merge}.bam.bai"
    conda:  
        config["environments"]["samtools"]
    benchmark: 
        "benchmarks/index_merged_bam/{merge}.log"
    log: 
        "logs/index_merged_bam/{merge}.log"
    resources:
        resources=config["default_resources"]
    shell:
        "samtools index {input} {output}"
