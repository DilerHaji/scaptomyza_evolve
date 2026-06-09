rule fix_fvariants:
    input:
        combined="fvariants/{fvar}.combined.vcf.gz",
    output:
        fixed="fvariants/{fvar}.fixed.vcf.gz",
    params: 
        fix_vcf = config["scripts"]["fix_vcf"]
    benchmark: 
        "benchmarks/fix_fvariants/{fvar}.log"
    log: 
        "logs/fix_fvariants/{fvar}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        bash {params.fix_vcf} {input} {output}
        """


rule fix_fvariants_no_neff:
    input:
        combined="fvariants/{fvar}.combined_no_neff.vcf.gz",
    output:
        fixed="fvariants/{fvar}.fixed_no_neff.vcf.gz",
    params: 
        fix_vcf = config["scripts"]["fix_vcf"]
    benchmark: 
        "benchmarks/fix_fvariants/{fvar}.log"
    log: 
        "logs/fix_fvariants/{fvar}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        bash {params.fix_vcf} {input} {output}
        """


rule fvariants_to_sync:
    input:
        fixed="fvariants/{fvar}.fixed.vcf.gz",
    output:
        done = "fvariants/{fvar}/all.done",
        sync = "fvariants/{fvar}.sync",
    params: 
        grenedalf = config["programs"]["grenedalf"],
        out1 = lambda wildcards: "fvariants/" + wildcards.fvar,
        input_call = config["grenedalf"]["input_call_sync"],
    conda:
        config["environments"]["grenedalf"]
    benchmark: 
        "benchmarks/fvariant_to_sync/{fvar}.log"
    log: 
        "logs/fvariant_to_sync/{fvar}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        
		[[ -d {params.out1} ]] || mkdir {params.out1}
		
		cp {input} {params.out1}

		cd {params.out1}

        {params.grenedalf} sync \
        {params.input_call} . \
        --file-prefix {wildcards.fvar} \
        --out-dir . > all.done
        
        mv {wildcards.fvar}sync.sync ../{wildcards.fvar}.sync
        
        rm *vcf.gz
            		                
        """