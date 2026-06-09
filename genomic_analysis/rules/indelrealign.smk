rule index_bam_RealignerTargetCreator:
    input: 
    	os.path.join(config["input_dir"], "{bam}.bam")
    output:
        os.path.join(config["input_dir"], "{bam}.bam.bai")
    conda:  
        config["environments"]["samtools"]
    benchmark: 
        "benchmarks/index_bam_RealignerTargetCreator/{bam}.log"
    log: 
        "logs/index_bam_RealignerTargetCreator/{bam}.log"
    resources:
        resources=config["default_resources"]
    shell:
        "samtools index {input} {output}"
        

rule RealignerTargetCreator: 
    input: 
        bam = os.path.join(config["input_dir"], "{bam}.bam"),
        idx = os.path.join(config["input_dir"], "{bam}.bam.bai"),
        reference_genome = config["reference_genome"]
    output:
        "realigned/{bam}.indels.intervals"
    params:
        extra = config["indel"]["extra"],
        jar_path = config["indel"]["jar_path"]
    conda: 
        config["environments"]["gatk3"]
    benchmark: 
        "benchmarks/RealignerTargetCreator/{bam}.log"
    log: 
        "logs/RealignerTargetCreator/{bam}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        [[ -f $CONDA_PREFIX/opt/gatk-3.5/GenomeAnalysisTK.jar ]] || gatk-register {params.jar_path}
        
        [[ -d realigned ]] || mkdir realigned
        
        GenomeAnalysisTK -T RealignerTargetCreator \
        -R {input.reference_genome} \
        -I {input.bam} \
        -o {output}
        
        """


rule IndelRealigner: 
    input: 
        target_intervals = "realigned/{bam}.indels.intervals", 
        bam = os.path.join(config["input_dir"], "{bam}.bam"),
        reference_genome = config["reference_genome"]
    output:
        "realigned/{bam}.bam"
    params:
        config["indel"]["extra"],
        jar_path = config["indel"]["jar_path"]
    conda: 
        config["environments"]["gatk3"]
    benchmark: 
        "benchmarks/IndelRealigner/{bam}.log"
    log: 
        "logs/IndelRealigner/{bam}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        
        #gatk-register {params.jar_path}
    	    	  
        GenomeAnalysisTK \
        -T IndelRealigner \
        -R {input.reference_genome} \
        -I {input.bam} \
        -targetIntervals {input.target_intervals} \
        -o {output}
        
        """
