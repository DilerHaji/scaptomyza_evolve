rule calculate_ly:
    input:
        af = "delta_tmp/{delta}/afmatsites.csv",
        neff = "delta_tmp/{delta}/neffsites.csv"
    output:
       "ly_input/{delta}.csv"
    params:
        out_prefix = "ly/{delta}",
        script = config["scripts"]["lynch"],
        pairs = lambda wildcards: DELTA_DICT[wildcards.delta][0]
    conda: 
        config["environments"]["polars"]
    benchmark: 
        "benchmarks/calculate_lynch/{delta}.log"
    log: 
        "logs/calculate_lynch/{delta}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        
        if [ -s {input.af} ]; then
            python {params.script} \
            --af_file {input.af} \
            --neff_file {input.neff} \
            --samples {params.pairs} \
            --output {output}
        else
            touch {output}
        fi

        """

rule collect_ly:
    input:
        lambda wildcards: [f"ly_input/{sample}.csv" for sample in LY_DICT[wildcards.ly][1].split("|")]
    output:
       "ly/{ly}.csv"
    output:
        "pca_delta/input/{pca}.csv"
    conda:
        config["environments"]["polars"]
    resources:
        resources = config["default_resources_10cpus"]
    script:
        config["scripts"]["pca_input"]




rule calculate_ly_ne:
    input:
        af = "delta_tmp/{delta}/afmatsites.csv",
        neff = "delta_tmp/{delta}/neffsites.csv"
    output:
       "ne/{delta}.csv"
    params:
        out_prefix = "ne/{delta}",
        script = config["scripts"]["lynch_ne"],
        pairs = lambda wildcards: DELTA_DICT[wildcards.delta]
    conda: 
        config["environments"]["polars"]
    benchmark: 
        "benchmarks/calculate_lynch/{delta}.log"
    log: 
        "logs/calculate_lynch/{delta}.log"
    resources:
        resources=config["default_resources_10cpus"]
    shell:
        """
        
        if [ -s {input.af} ]; then
            python {params.script} \
            --af_file {input.af} \
            --neff_file {input.neff} \
            --sample {params.pairs} \
            --output {output} \
            --af_lower 0.01 \
            --af_upper 0.99 \
            --min_samples_percent 100 \
            --generation_interval 1
        else
            touch {output}
        fi

        """



