def gather_poolsnp(wildcards):
    return ["variants/poolsnp/" + item + ".combined.vcf.gz" for item in FVAR_DICT[wildcards.fvar][0]]

def gather_poolsnp_done(wildcards):
    return ["variants_index/poolsnp/" + item + ".done" for item in FVAR_DICT[wildcards.fvar][0]]


rule poolsnp_gathered: 
    input:
        vcf=gather_poolsnp,
        index_done=gather_poolsnp_done
    output:
        vcf="gathered_variants/{fvar}.vcf.gz",
        index="gathered_variants/{fvar}.vcf.gz.csi"
    benchmark: 
        "benchmarks/poolsnp_gathered/{fvar}.log"
    log: 
        "logs/poolsnp_gathered/{fvar}.log"
    resources:
        resources=config["default_resources_5cpus"]
    conda:
        config["environments"]["bcftools"]
    shell:
        """
        
        input_count=$(echo {input.vcf} | wc -w)
        
        if [ "$input_count" -eq 1 ]; then
            cp {input.vcf} {output.vcf} 2>> {log}
        else
            bcftools merge --threads 24 -Oz -o {output.vcf} {input.vcf} 2>> {log}
        fi

        bcftools index --threads 24 {output.vcf} 2>> {log}
                        
        """

rule fvariants_get_header:
    input:
        "gathered_variants/{fvar}.vcf.gz"
    output:
        "fvariants/{fvar}.header"
    benchmark: 
        "benchmarks/fvariants_get_header/{fvar}.log"
    log: 
        "logs/fvariants_get_header/{fvar}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        (zcat {input} | head -n 10000 | grep ^"#" > {output}) || true
        """

# rule fvariants_subset_region:
#     input:
#         vcf = "gathered_variants/{fvar}.vcf.gz",
#         index = "gathered_variants/{fvar}.vcf.gz.csi",
#         header = "fvariants/{fvar}.header"
#     output:
#         vcf=temp("fvariants/subset_region/{fvar}/{entry}.vcf.gz"),
#         done=temp("fvariants/subset_region/{fvar}/{entry}.done")
#     params:
#         make_single_bed = config["scripts"]["make_single_bed"],
#         region = "{entry}",
#         temp_bed = "fvariants/subset_region/{fvar}/{entry}_region.bed",
#         temp_body = "fvariants/subset_region/{fvar}/{entry}_body.vcf"
#     conda:
#         config["environments"]["htslib"]
#     resources:
#         resources=config["default_resources"]
#     shell:
#         """
#         mkdir -p $(dirname {params.temp_bed})
#         bash {params.make_single_bed} {params.region} {params.temp_bed}
#         tabix -f -R {params.temp_bed} {input.vcf} > {params.temp_body}
#         cat {input.header} {params.temp_body} | bgzip -c > {output.vcf}
#         rm {params.temp_bed} {params.temp_body}
#         touch {output.done}
#         """

# rule fvariants_subset_region:
#     input:
#         vcf = "gathered_variants/{fvar}.vcf.gz",
#         index = "gathered_variants/{fvar}.vcf.gz.csi",
#         header = "fvariants/{fvar}.header",
#         # We need the index to get contig lengths for the BED file
#         fai = config["reference_genome"] + ".fai" 
#     output:
#         vcf=temp("fvariants/subset_region/{fvar}/{entry}.vcf.gz"),
#         done=temp("fvariants/subset_region/{fvar}/{entry}.done")
#     params:
#         # Retrieve the comma-separated string of contigs for this Bin ID
#         contigs = lambda wildcards: CONTIG_BINS[wildcards.entry],
#         temp_bed = "fvariants/subset_region/{fvar}/{entry}_region.bed",
#         temp_body = "fvariants/subset_region/{fvar}/{entry}_body.vcf"
#     conda:
#         config["environments"]["htslib"]
#     resources:
#         resources=config["default_resources"]
#     shell:
#         """
#         # 1. Create a dynamic BED file for all contigs in this bin
#         rm -f {params.temp_bed}
#         
#         # Split comma-separated contigs and create BED format (Name 0 Length) using the .fai file
#         # This approach avoids command-line length limits with tabix
#         IFS=',' read -r -a contig_array <<< "{params.contigs}"
#         for c in "${{contig_array[@]}}"; do
#             grep -P "^$c\t" {input.fai} | awk '{{print $1"\t0\t"$2}}' >> {params.temp_bed}
#         done
# 
#         # 2. Extract regions
#         # tabix -R uses the BED file to extract multiple regions efficiently
#         tabix -h -R {params.temp_bed} {input.vcf} | bgzip -c > {output.vcf}
#         
#         # Cleanup
#         rm {params.temp_bed}
#         touch {output.done}
#         """

rule fvariants_subset_region:
    input:
        vcf = "gathered_variants/{fvar}.vcf.gz",
        index = "gathered_variants/{fvar}.vcf.gz.csi",
        header = "fvariants/{fvar}.header",
        # We restore the .fai dependency to get exact contig lengths
        fai = config["reference_fai"]
    output:
        vcf=temp("fvariants/subset_region/{fvar}/{entry}.vcf.gz"),
        done=temp("fvariants/subset_region/{fvar}/{entry}.done")
    params:
        # Important: Split the comma-separated string into a Python list
        # so Snakemake passes them as individual arguments to the shell loop
        contigs = lambda wildcards: CONTIG_BINS.get(wildcards.entry, "").split(","),
        temp_bed = "fvariants/subset_region/{fvar}/{entry}_region.bed"
    conda:
        config["environments"]["htslib"]
    resources:
        resources=config["default_resources"]
    shell:
        """
        # 1. Create a dynamic BED file for all contigs in this bin
        rm -f {params.temp_bed}
        
        # Iterate over the contig list provided by python
        for c in {params.contigs:q}; do
            # Grep the contig name from the .fai file (column 1) and get the length (column 2)
            # -P "^$c\t" ensures we match the exact contig name at the start of the line
            grep -P "^$c\t" {input.fai} | awk '{{print $1"\t0\t"$2}}' >> {params.temp_bed}
        done

        # 2. Extract regions
        # tabix -R uses the BED file to extract exactly the regions defined
        tabix -h -R {params.temp_bed} {input.vcf} | bgzip -c > {output.vcf}
        
        # Cleanup
        rm {params.temp_bed}
        touch {output.done}
        """


rule fvariants_sort_region:
    input:
        vcf="fvariants/subset_region/{fvar}/{entry}.vcf.gz",
        done="fvariants/subset_region/{fvar}/{entry}.done"
    output:
        vcf = temp("fvariants/sort_region/{fvar}/{entry}.vcf.gz"),
        done = temp("fvariants/sort_region/{fvar}/{entry}.done")
    resources:
        resources=config["default_resources"]
    conda:
        config["environments"]["bcftools"]
    shell:
        """
        bcftools sort -o {output.vcf} -O z {input.vcf}
        touch {output.done}
        """

# KEY VARIANT FILTERING STEP 
rule create_variant_chunks:
    input:
        variants = lambda wildcards: FVAR_DICT[wildcards.fvar][1] + "_final_variants.txt",
    output:
        expand("fvariants/variant_chunks/{{fvar}}/chunk_{i}.txt", i=range(1, config["tabix"]["chunks"] + 1))
    params:
        variant_chunks = config["scripts"]["variant_chunks"],
        chunks = config["tabix"]["chunks"],
        prefix = "fvariants/variant_chunks/{fvar}/chunk_",
        dir0 = "fvariants/variant_chunks/{fvar}"
    shell:
        """
        bash {params.variant_chunks} {input.variants} {params.prefix} {params.chunks} {params.dir0}
        """

localrules: create_variant_chunks


rule pfilterPCA_filter_region:
    input:
        vcf = "fvariants/sort_region/{fvar}/{entry}.vcf.gz",
        done = "fvariants/sort_region/{fvar}/{entry}.done",
        chunks = expand("fvariants/variant_chunks/{{fvar}}/chunk_{i}.txt", i=range(1, config["tabix"]["chunks"] + 1)),
        header = "fvariants/{fvar}.header"
    output:
        vcf = temp("fvariants/filter_region/{fvar}/{entry}.vcf.gz"),
        temp_body = temp("fvariants/filter_region/{fvar}/{entry}_body.vcf"),
        done = temp("fvariants/filter_region/{fvar}/{entry}.done")
    params: 
        dir0 = "fvariants/filter_region/{fvar}/{entry}",
        tabix = config["scripts"]["tabix"],
        sub = "fvariants/filter_region/{fvar}/sub_{entry}.vcf",
    log: 
        "logs/pfilterPCA_filter_region/{fvar}/{entry}.log"
    resources:
        resources = config["default_resources"]
    conda:
        config["environments"]["htslib"]
    shell:
        """      
        tabix -f -p vcf {input.vcf}
    
        for chunk in {input.chunks}; do
            tabix -f -R $chunk {input.vcf} >> {output.temp_body}
        done
        
        touch {output.temp_body}
        
        cat {input.header} {output.temp_body} | bgzip -c > {output.vcf}
        touch {output.vcf}
        zcat {output.vcf} | head -n 5000 > {params.sub} || true
        touch {params.sub}
        touch {output.done}
        """

def get_corrections(wildcards): 
    return(f"variants/poolsnp_ncorrect/{poolsnp}.csv" for poolsnp in FVAR_DICT[wildcards.fvar][0])


rule pfilterPCA_get_corrections:
    input:
        get_corrections
    output:
        corrections = temp("fvariants/get_corrections/{fvar}.csv"),
        header = temp("fvariants/get_corrections/header_{fvar}.csv"),
    params: 
        sub = "fvariants/get_corrections/sub_{fvar}.csv"
    resources:
        resources = config["default_resources"]
    shell:
        """
        awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input} > {output.corrections}
        head -n 5000 {output.corrections} > {params.sub}
        head -n 1 {output.corrections} > {output.header}
        """


# rule fvariants_correct_filter_region:
#     input:
#         vcf = "fvariants/filter_region/{fvar}/{entry}.vcf.gz",
#         done = "fvariants/filter_region/{fvar}/{entry}.done",
#         corrections = "fvariants/get_corrections/{fvar}.csv",
#         header = "fvariants/get_corrections/header_{fvar}.csv"
#     output:
#         vcf = temp("fvariants/corrected_filter_region/{fvar}/{entry}.vcf.gz"),
#         tmp1 = temp("fvariants/get_corrections/{fvar}/tmp_{entry}.csv"),
#         #tmp2 = "fvariants/get_corrections/{fvar}/2{entry}.csv",
#         done = temp("fvariants/corrected_filter_region/{fvar}/{entry}.done")
#     params: 
#         neff_vcf = config["scripts"]["neff_vcf"],
#         subset_ncorrect = config["scripts"]["subset_ncorrect"],
#         sub = "fvariants/get_corrections/{fvar}/sub_{entry}.vcf",
#         correction = lambda wildcards: FVAR_DICT[wildcards.fvar][2]
#     log: 
#         "logs/fvariants_correct_filter_region/{fvar}/{entry}.log"
#     resources:
#         resources = config["default_resources"]
#     conda:
#         config["environments"]["polars"]
#     shell:
#         """
#         
#         python {params.subset_ncorrect} {input.corrections} {output.tmp1} {wildcards.entry} {params.correction}
#         
#         if [ -s {output.tmp1} ]; then
#             python {params.neff_vcf} \
#             --vcf {input.vcf} \
#             --output {output.vcf} \
#             --correction_files {output.tmp1} \
#             --correction_column {params.correction}
#             
#             zcat {output.vcf} | head -n 5000 > {params.sub} || true
#        
#         else
#             touch {output.vcf}
#         fi
#            
#         touch {output.done}
#         """

rule fvariants_correct_filter_region:
    input:
        vcf = "fvariants/filter_region/{fvar}/{entry}.vcf.gz",
        done = "fvariants/filter_region/{fvar}/{entry}.done",
        corrections = "fvariants/get_corrections/{fvar}.csv",
        header = "fvariants/get_corrections/header_{fvar}.csv"
    output:
        vcf = temp("fvariants/corrected_filter_region/{fvar}/{entry}.vcf.gz"),
        tmp1 = temp("fvariants/get_corrections/{fvar}/tmp_{entry}.csv"),
        done = temp("fvariants/corrected_filter_region/{fvar}/{entry}.done")
    params: 
        neff_vcf = config["scripts"]["neff_vcf"],
        subset_ncorrect = config["scripts"]["subset_ncorrect"],
        sub = "fvariants/get_corrections/{fvar}/sub_{entry}.vcf",
        correction = lambda wildcards: FVAR_DICT[wildcards.fvar][2],

        # --- CHANGED: Use the helper function from common.smk ---
        contigs = get_contigs_from_bin
    log: 
        "logs/fvariants_correct_filter_region/{fvar}/{entry}.log"
    resources:
        resources = config["default_resources_24cpus"]
    conda:
        config["environments"]["polars"]
    shell:
        """
        python {params.subset_ncorrect} {input.corrections} {output.tmp1} "{params.contigs}" {params.correction}
        
        if [ -s {output.tmp1} ]; then
            python {params.neff_vcf} \
            --vcf {input.vcf} \
            --output {output.vcf} \
            --correction_files {output.tmp1} \
            --correction_column {params.correction}
            
            zcat {output.vcf} | head -n 5000 > {params.sub} || true
       
        else
            touch {output.vcf}
        fi
           
        touch {output.done}
        """


# rule fvariants_correct_filter_region_bgzip:
#     input:
#         vcf = "fvariants/corrected_filter_region/{fvar}/{entry}.vcf",
#         done = "fvariants/corrected_filter_region/{fvar}/{entry}.done"
#     output:
#         "fvariants/corrected_filter_region/{fvar}/{entry}.vcf.gz",
#     resources:
#         resources = config["default_resources"]
#     conda:
#         config["environments"]["htslib"]
#     shell:
#         """
#         if [ -s {input.vcf} ]; then
#             bgzip -c {input.vcf} > {output}
#         else
#             touch {output}
#         fi
#         """


rule fvariants_sort_filtered_region:
    input:
        vcf = "fvariants/corrected_filter_region/{fvar}/{entry}.vcf.gz",
        done = "fvariants/corrected_filter_region/{fvar}/{entry}.done"
    output:
        vcf = "fvariants/sort_filtered_region/{fvar}/{entry}.vcf.gz",
        done = temp("fvariants/sort_filtered_region/{fvar}/{entry}.done")
    resources:
        resources=config["default_resources"]
    conda:
        config["environments"]["bcftools"]
    shell:
        """
        if [ -s {input.vcf} ]; then
            bcftools sort {input.vcf} | bcftools norm -d all -o {output.vcf} -O z
        else
            touch {output.vcf} || true
        fi
        touch {output.done} || true
        """
        

rule fvariants_index_filtered_region:
    input:
        vcf = "fvariants/sort_filtered_region/{fvar}/{entry}.vcf.gz",
        done = "fvariants/sort_filtered_region/{fvar}/{entry}.done"
    output:
        "fvariants/index_filtered_region/{fvar}/{entry}.done"
    resources:
        resources=config["default_resources"]
    conda:
        config["environments"]["htslib"]
    shell:
        """
        if [ -s {input.vcf} ]; then
            tabix -f -p vcf {input.vcf} || true
        fi
        touch {output} || true
        """


def get_fvariants_vcf(wildcards):
    entries = get_entries()
    return [f"fvariants/sort_filtered_region/{wildcards.fvar}/{entry}.vcf.gz" for entry in entries]

def get_fvariants_vcf_index(wildcards):
    entries = get_entries()
    return [f"fvariants/index_filtered_region/{wildcards.fvar}/{entry}.done" for entry in entries]


rule aggregate_fvariants:
    input:
        vcf=get_fvariants_vcf,
        index_done=get_fvariants_vcf_index
    output:
        combined="fvariants/{fvar}.combined.vcf.gz",
        #header="variants/poolsnp/{poolsnp}.header",
        #variants="variants/poolsnp/{poolsnp}.variants"
    benchmark: 
        "benchmarks/aggregate_fvariants/{fvar}.log"
    log: 
        "logs/aggregate_fvariants/{fvar}.log"
    resources:
        resources=config["default_resources"]
    conda:
        config["environments"]["bcftools"]
    shell:
        """
        non_empty_files=$(for file in {input.vcf}; do [ -s "$file" ] && echo "$file"; done)
        if [ -n "$non_empty_files" ]; then
            bcftools concat -Oz -o {output.combined} $non_empty_files 2> {log}
        else
            touch {output.combined} || true
        fi        
        """


### Getting a list of positions for downstream applications 
### VCF file should be 1-based since poolSNP uses the mpileup file generated by samtools mpileup, which is 1-based 

rule fvariants_positions: 
    input:
        "fvariants/{fvar}.combined.vcf.gz"
    output:
        "fvariants/{fvar}.positions"
    benchmark: 
        "benchmarks/fvariants_positions/{fvar}.log"
    log: 
        "logs/fvariants_positions/{fvar}.log"
    resources:
        resources=config["default_resources"]
    conda:
        config["environments"]["bcftools"]
    shell:
        """
        
        bcftools query -f '%CHROM\\t%POS\\n' {input} | sort -k1,1 -k2,2n > {output}
            
        """





#### This version is used for calculations that internally correct for pooled sequencing bias 
####
####
####

rule fvariants_sort_filtered_region_no_neff:
    input:
        vcf = "fvariants/filter_region/{fvar}/{entry}.vcf.gz",
        done = "fvariants/filter_region/{fvar}/{entry}.done"
    output:
        vcf = "fvariants/sort_filtered_region_no_neff/{fvar}/{entry}.vcf.gz",
        done = temp("fvariants/sort_filtered_region_no_neff/{fvar}/{entry}.done")
    resources:
        resources=config["default_resources"]
    conda:
        config["environments"]["bcftools"]
    shell:
        """
        if [ -s {input.vcf} ]; then
            bcftools sort {input.vcf} | bcftools norm -d all -o {output.vcf} -O z
        else
            touch {output.vcf} || true
        fi
        touch {output.done} || true
        """
        

rule fvariants_index_filtered_region_no_neff:
    input:
        vcf = "fvariants/sort_filtered_region_no_neff/{fvar}/{entry}.vcf.gz",
        done = "fvariants/sort_filtered_region_no_neff/{fvar}/{entry}.done"
    output:
        "fvariants/index_filtered_region_no_neff/{fvar}/{entry}.done"
    resources:
        resources=config["default_resources"]
    conda:
        config["environments"]["htslib"]
    shell:
        """
        if [ -s {input.vcf} ]; then
            tabix -f -p vcf {input.vcf} || true
        fi
        touch {output} || true
        """


def get_fvariants_vcf_no_neff(wildcards):
    entries = get_entries()
    return [f"fvariants/sort_filtered_region_no_neff/{wildcards.fvar}/{entry}.vcf.gz" for entry in entries]

def get_fvariants_vcf_index_no_neff(wildcards):
    entries = get_entries()
    return [f"fvariants/index_filtered_region_no_neff/{wildcards.fvar}/{entry}.done" for entry in entries]


rule aggregate_fvariants_no_neff:
    input:
        vcf=get_fvariants_vcf_no_neff,
        index_done=get_fvariants_vcf_index_no_neff
    output:
        combined="fvariants/{fvar}.combined_no_neff.vcf.gz",
        #header="variants/poolsnp/{poolsnp}.header",
        #variants="variants/poolsnp/{poolsnp}.variants"
    benchmark: 
        "benchmarks/aggregate_fvariants/{fvar}.log"
    log: 
        "logs/aggregate_fvariants/{fvar}.log"
    resources:
        resources=config["default_resources"]
    conda:
        config["environments"]["bcftools"]
    shell:
        """
        non_empty_files=$(for file in {input.vcf}; do [ -s "$file" ] && echo "$file"; done)
        if [ -n "$non_empty_files" ]; then
            bcftools concat -Oz -o {output.combined} $non_empty_files 2> {log}
        else
            touch {output.combined}
        fi        
        """




