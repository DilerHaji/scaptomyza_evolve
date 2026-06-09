import numpy as np

# rule pglm_convert:
#     input:
#         vcf = lambda wildcards: os.path.join("fvariants/sort_filtered_region", GLM_DICT[wildcards.glm][0], wildcards.entry + ".vcf.gz"),
#         metadata = config["sample_metadata"],
#     output:
#         csv = "glm/convert/{glm}/con_{entry}.csv",
#         done = temp("glm/convert/{glm}/con_{entry}.done")
#     params:
#         script = config["scripts"]["convert_glm"],
#         trts = lambda wildcards: GLM_DICT[wildcards.glm][1],
#     conda: 
#         config["environments"]["polars"]
#     benchmark: 
#         "benchmarks/glm_convert/{glm}/{entry}.log"
#     log: 
#         "logs/glm_convert/{glm}/{entry}.log"
#     resources:
#         resources=config["default_resources_5cpus"]
#     shell:
#         """
#         
#         variants=$(zcat {input.vcf} | grep -v '^#' | wc -l) || true
#         
#         if [ "$variants" -gt 0 ]; then
#             python {params.script} \
#             -v {input.vcf} \
#             -m {input.metadata} \
#             -o {output.csv} \
#             -t {params.trts} 2> {log}
#             
#             touch {output.csv}
#        
#         else
#            
#             touch {output.csv}
#        
#         fi
#         
#         touch {output.done}
#         
#         """


rule pglm_convert:
    input:
        vcf = lambda wildcards: os.path.join("fvariants/sort_filtered_region", GLM_DICT[wildcards.glm][0], wildcards.entry + ".vcf.gz"),
        vcf_index = lambda wildcards: os.path.join("fvariants/index_filtered_region", GLM_DICT[wildcards.glm][0], wildcards.entry + ".done"),
        metadata = config["sample_metadata"],
    output:
        csv = "glm/convert/{glm}/con_{entry}.csv",
        done = temp("glm/convert/{glm}/con_{entry}.done")
    params:
        script = config["scripts"]["convert_glm"],
        trts = lambda wildcards: GLM_DICT[wildcards.glm][1],
    conda: 
        config["environments"]["polars"]
    benchmark: 
        "benchmarks/glm_convert/{glm}/{entry}.log"
    log: 
        "logs/glm_convert/{glm}/{entry}.log"
    resources:
        resources=config["default_resources_5cpus"]
    shell:
        """
        # Count variants to avoid running python on empty files
        variants=$(zcat {input.vcf} | grep -v '^#' | wc -l) || true
        
        if [ "$variants" -gt 0 ]; then
            python {params.script} \
            -v {input.vcf} \
            -m {input.metadata} \
            -o {output.csv} \
            -t {params.trts} 2> {log}
            
            touch {output.csv}
        else
            touch {output.csv}
        fi
        
        touch {output.done}
        """


rule pglm_preprocessing:
    input:
        csv = "glm/convert/{glm}/con_{entry}.csv",
        done = "glm/convert/{glm}/con_{entry}.done"
    output:
        neff_mat = "glm/preprocessing/{glm}/pre_{entry}_neff_mat.csv",
        afmat = "glm/preprocessing/{glm}/pre_{entry}_afmat.csv",
        sites = "glm/preprocessing/{glm}/pre_{entry}_sites.csv",
        samps = "glm/preprocessing/{glm}/pre_{entry}_samps.csv",
    params:
        script = config["scripts"]["glm_preprocessing"],
        prefix = lambda wildcards: os.path.join("glm", "preprocessing", wildcards.glm, "pre_" + wildcards.entry),
        freq_column = config["pfilter"]["freq_column"],
        neff_column = "neff", # option to change this moved to VCF processing for PCA
        variables = config["glm"]["variables"],
    conda: 
        config["environments"]["polars"]
    benchmark: 
        "benchmarks/glm1_preprocessing_glm/{glm}/{entry}.log"
    log: 
        "logs/glm1_preprocessing_glm/{glm}/{entry}.log"
    resources:
        resources=config["default_resources_24cpus"]
    shell:
        """
        
        if [ ! -s {input.csv} ]; then
            touch {output.neff_mat}
            touch {output.afmat}
            touch {output.sites}
            touch {output.samps}
        else
            python {params.script} \
            {input.csv} \
            {params.freq_column} \
            {params.neff_column} \
            {params.prefix} \
            {params.variables} > {log} 2>&1
        fi
            
        """



rule pglm_gethaf:
    input:
        neff_mat = "glm/preprocessing/{glm}/pre_{entry}_neff_mat.csv",
        afmat = "glm/preprocessing/{glm}/pre_{entry}_afmat.csv",
        sites = "glm/preprocessing/{glm}/pre_{entry}_sites.csv",
        samps = "glm/preprocessing/{glm}/pre_{entry}_samps.csv",
    output:
        rdata = "glm/haf/{glm}/gethaf_{entry}.Rdata",
        rds = "glm/haf/{glm}/gethaf_{entry}_neff.RDS",
    params:
        glm_gethaf = config["scripts"]["glm_gethaf"],
        prefix = lambda wildcards: os.path.join("glm", "haf", wildcards.glm, "gethaf_" + wildcards.entry),
    conda: 
        config["environments"]["base_r"]
    benchmark: 
        "benchmarks/glm1_gethaf_glm/{glm}/{entry}.log"
    log: 
        "logs/glm1_gethaf_glm/{glm}/{entry}.log"
    resources:
        resources=config["default_resources_2cpus"]
    shell:
        """
        
        if [ ! -s {input.neff_mat} ] || [ ! -s {input.afmat} ] || [ ! -s {input.sites} ] || [ ! -s {input.samps} ]; then
            touch {output.rdata}
            touch {output.rds}
        else
            Rscript {params.glm_gethaf} \
            {input.neff_mat} \
            {input.afmat} \
            {input.sites} \
            {input.samps} \
            {params.prefix} > {log} 2>&1
        fi
        
        if [ ! -f {output.rdata} ]; then
            touch {output.rdata}
        fi

        if [ ! -f {output.rds} ]; then
            touch {output.rds}
        fi

        """


rule pglm_merge_haf:
    input:
        rdata = expand("glm/haf/{{glm}}/gethaf_{entry}.Rdata", entry=get_entries()),
        rds   = expand("glm/haf/{{glm}}/gethaf_{entry}_neff.RDS", entry=get_entries())
    output:
        rdata = "glm/haf/{glm}/merged_gethaf.Rdata",
        rds   = "glm/haf/{glm}/merged_gethaf_neff.RDS"
    params:
        merge_script = "scripts/glm_merge_haf.r"
    log:
        "logs/glm_merge_haf/{glm}/merge.log"
    conda: 
        config["environments"]["base_r"]
    resources:
        resources=config["default_resources_24cpus"]
    shell:
        """
        Rscript {params.merge_script} \
        {output.rdata} \
        {output.rds} \
        {input.rdata} \
        {input.rds} > {log} 2>&1
        """


rule pglm:
    input:
        rdata = "glm/haf/{glm}/gethaf_{entry}.Rdata",
        rds = "glm/haf/{glm}/gethaf_{entry}_neff.RDS",
    output:
        csv = "glm/glm/{glm}/glm_{entry}.csv",
        done = temp("glm/glm/{glm}/glm_{entry}.done")
    params:
        glm = config["scripts"]["glm"],
        prefix = lambda wildcards: os.path.join("glm", "glm", wildcards.glm, "glm_" + wildcards.entry),
        mainEffect = config["glm"]["mainEffect"],
        effect1 = config["glm"]["effect1"],
        effect1type = config["glm"]["effect1type"],
        effect2 = config["glm"]["effect2"],
        effect2type = config["glm"]["effect2type"],
        interaction = config["glm"]["interaction"],
        repName = config["glm"]["repName"],
        repNamecontrast = config["glm"]["repNamecontrast"],
        repNameignore = config["glm"]["repNameignore"],
        repNameinter = config["glm"]["repNameinter"],
        repNameeff = config["glm"]["repNameeff"],
        trtPopGroups = config["glm"]["trtPopGroups"],
        makePlots = config["glm"]["makePlots"],
        saveModels = config["glm"]["saveModels"],
        mixedEffects = config["glm"]["mixedEffects"],
        compareModels = config["glm"]["compareModels"],
        trtContrast = lambda wildcards: GLM_DICT[wildcards.glm][3],
        # NEW: Get treatments from config (default to NA if not set)
        selectTrts = config["glm"].get("selectTrts", "NA")
    conda: 
        config["environments"]["base_r"]
    benchmark: 
        "benchmarks/glm/{glm}/{entry}.log"
    log: 
        "logs/glm/{glm}/{entry}.log"
    resources:
        resources=config["default_resources_2cpus"]
    shell:
        """

        if [ ! -s {input.rdata} ] || [ ! -s {input.rds} ]; then
            touch {output.csv}
        else
            Rscript {params.glm} \
            {input.rdata} \
            --readDepth {input.rds} \
            --mainEffect {params.mainEffect} \
            --effect1 {params.effect1} \
            --effect1type {params.effect1type} \
            --effect2 {params.effect2} \
            --effect2type {params.effect2type} \
            --interaction {params.interaction} \
            --repName {params.repName} \
            --repNamecontrast {params.repNamecontrast} \
            --repNameignore {params.repNameignore} \
            --repNameinter {params.repNameinter} \
            --repNameeff {params.repNameeff} \
            --trtPopGroups {params.trtPopGroups} \
            --trtContrast {params.trtContrast} \
            --mixedEffects {params.mixedEffects} \
            --compareModels {params.compareModels} \
            --randomEffectsType nested \
            --makePlots {params.makePlots} \
            --saveModels {params.saveModels} \
            --selectTrts {params.selectTrts} \
            --saveAs csv \
            --outDir {params.prefix} > {log} 2>&1
        fi
        
        touch {output.done}
        
        """

def pgather_glm(wildcards):
    entries = get_entries()
    return [f"glm/glm/{wildcards.glm}/glm_{entry}.csv" for entry in entries]

def pgather_glm_done(wildcards):
    entries = get_entries()
    return [f"glm/glm/{wildcards.glm}/glm_{entry}.done" for entry in entries]


rule pglm_final:
    input:
        csv = pgather_glm,
        done = pgather_glm_done
    output:
        glm = "glm_final/{glm}.csv",
    benchmark: 
        "benchmarks/glm_final/{glm}.log"
    log: 
        "logs/glm_final/{glm}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        
        header_written=false
        for file in {input.csv}; do
            if [[ -s "$file" ]]; then
                head -n 1 "$file" > {output.glm}
                header_written=true
                break
            fi
        done

        for file in {input.csv}; do
            tail -n +2 $file >> {output.glm}
        done
         
        """








###### Leave-one-out #####

rule pglm_loo:
    input:
        rdata = "glm/haf/{glm}/gethaf_{entry}.Rdata",
        rds = "glm/haf/{glm}/gethaf_{entry}_neff.RDS",
    output:
        csv = "glm/glm_loo/{glm}/glm_{entry}.csv",
        done = temp("glm/glm_loo/{glm}/glm_{entry}.done")
    params:
        # POINT THIS TO THE NEW R SCRIPT PROVIDED BELOW
        glm_script = "scripts/calc_GLM_LOO.R", 
        prefix = lambda wildcards: os.path.join("glm", "glm_loo", wildcards.glm, "glm_" + wildcards.entry),
        mainEffect = config["glm"]["mainEffect"],
        effect1 = config["glm"]["effect1"],
        effect1type = config["glm"]["effect1type"],
        effect2 = config["glm"]["effect2"],
        effect2type = config["glm"]["effect2type"],
        interaction = config["glm"]["interaction"],
        repName = config["glm"]["repName"],
        repNamecontrast = config["glm"]["repNamecontrast"],
        repNameignore = config["glm"]["repNameignore"],
        repNameinter = config["glm"]["repNameinter"],
        repNameeff = config["glm"]["repNameeff"],
        trtPopGroups = config["glm"]["trtPopGroups"],
        mixedEffects = config["glm"]["mixedEffects"],
        randomEffectsType = "nested", # Defaulting to nested based on your script context, adjust if needed
        trtContrast = lambda wildcards: GLM_DICT[wildcards.glm][3],
        selectTrts = config["glm"].get("selectTrts", "NA")
    conda: 
        config["environments"]["base_r"]
    benchmark: 
        "benchmarks/glm_loo/{glm}/{entry}.log"
    log: 
        "logs/glm_loo/{glm}/{entry}.log"
    resources:
        resources=config["default_resources_2cpus"]
    shell:
        """
        if [ ! -s {input.rdata} ] || [ ! -s {input.rds} ]; then
            touch {output.csv}
        else
            Rscript {params.glm_script} \
            {input.rdata} \
            --readDepth {input.rds} \
            --mainEffect {params.mainEffect} \
            --effect1 {params.effect1} \
            --effect1type {params.effect1type} \
            --effect2 {params.effect2} \
            --effect2type {params.effect2type} \
            --interaction {params.interaction} \
            --repName {params.repName} \
            --repNamecontrast {params.repNamecontrast} \
            --repNameignore {params.repNameignore} \
            --repNameinter {params.repNameinter} \
            --repNameeff {params.repNameeff} \
            --trtPopGroups {params.trtPopGroups} \
            --trtContrast {params.trtContrast} \
            --mixedEffects {params.mixedEffects} \
            --randomEffectsType {params.randomEffectsType} \
            --selectTrts {params.selectTrts} \
            --saveAs csv \
            --outDir {params.prefix} > {log} 2>&1
        fi
        
        touch {output.done}
        """


# Helper functions to gather LOO results
def pgather_glm_loo(wildcards):
    entries = get_entries()
    return [f"glm/glm_loo/{wildcards.glm}/glm_{entry}.csv" for entry in entries]

def pgather_glm_loo_done(wildcards):
    entries = get_entries()
    return [f"glm/glm_loo/{wildcards.glm}/glm_{entry}.done" for entry in entries]

# The Aggregation Rule
rule pglm_final_loo:
    input:
        csv = pgather_glm_loo,
        done = pgather_glm_loo_done
    output:
        glm = "glm_final_loo/{glm}.csv",
    benchmark: 
        "benchmarks/glm_final_loo/{glm}.log"
    log: 
        "logs/glm_final_loo/{glm}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        
        header_written=false
        # Loop through files to find the first non-empty one to extract the header
        for file in {input.csv}; do
            if [[ -s "$file" ]]; then
                head -n 1 "$file" > {output.glm}
                header_written=true
                break
            fi
        done

        # Concatenate the content (skipping header) of all files
        for file in {input.csv}; do
            if [[ -s "$file" ]]; then
                tail -n +2 $file >> {output.glm}
            fi
        done
         
        """