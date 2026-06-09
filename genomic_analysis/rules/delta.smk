print(DELTA_DICT)
print(DELTA)


import functools
import shutil


@functools.lru_cache(maxsize=None)
def get_entries_cached():
    try:
        return get_entries()
    except NameError:
        return list(CONTIG_BINS.keys())

def gather_pre_neff(wildcards):
    return expand("glm/preprocessing/{glm}/pre_{entry}_neff_mat.csv", 
                  glm=wildcards.glm, entry=get_entries_cached())

def gather_pre_afmat(wildcards):
    return expand("glm/preprocessing/{glm}/pre_{entry}_afmat.csv", 
                  glm=wildcards.glm, entry=get_entries_cached())

def gather_pre_sites(wildcards):
    return expand("glm/preprocessing/{glm}/pre_{entry}_sites.csv", 
                  glm=wildcards.glm, entry=get_entries_cached())

def gather_pre_samps(wildcards):
    entries = get_entries_cached()
    return f"glm/preprocessing/{wildcards.glm}/pre_{entries[0]}_samps.csv"



rule aggregate_glm_preprocessing:
    wildcard_constraints:
        glm = "[^/]+"
    input:
        neff_mats = gather_pre_neff,
        afmats = gather_pre_afmat,
        sites = gather_pre_sites,
        samps = gather_pre_samps
    output:
        neff_mat = "glm/preprocessing/{glm}_neff_mat.csv",
        afmat = "glm/preprocessing/{glm}_afmat.csv",
        sites = "glm/preprocessing/{glm}_sites.csv",
        samps = "glm/preprocessing/{glm}_samps.csv"
    log:
        "logs/aggregate_glm/{glm}.log"
    resources:
        resources=config.get("default_resources", "")
    run:
        def concatenate_files(input_files, output_file):
            with open(output_file, 'w') as outfile:
                header_found = False
                if input_files:
                    with open(input_files[0], 'r') as infile:
                        header = infile.readline()
                        if header:
                            outfile.write(header)
                            header_found = True
                
                for fname in input_files:
                    with open(fname, 'r') as infile:
                        header_line = infile.readline()
                        if not header_line:
                            continue
                        for line in infile:
                            outfile.write(line)

        concatenate_files(input.neff_mats, output.neff_mat)
        concatenate_files(input.afmats, output.afmat)
        concatenate_files(input.sites, output.sites)
        
        shutil.copyfile(input.samps, output.samps)


rule get_matrix:
    input:
        neff_mat = lambda wildcards: "glm/preprocessing/" + DELTA_DICT[wildcards.delta][1] + "_neff_mat.csv",
        afmat =  lambda wildcards: "glm/preprocessing/" + DELTA_DICT[wildcards.delta][1] + "_afmat.csv",
        sites =  lambda wildcards: "glm/preprocessing/" + DELTA_DICT[wildcards.delta][1] + "_sites.csv",
        samps =  lambda wildcards: "glm/preprocessing/" + DELTA_DICT[wildcards.delta][1] + "_samps.csv",
    output:
       afmatsites = "delta_tmp/{delta}/afmatsites.csv",
       neffsites = "delta_tmp/{delta}/neffsites.csv"
    benchmark: 
        "benchmarks/{delta}/get_matrix.log"
    log: 
        "logs/{delta}/get_matrix.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        paste -d ',' {input.sites} {input.afmat} > {output.afmatsites} || true
        if [ ! -s {output.afmatsites} ]; then touch {output.afmatsites}; fi
        
        paste -d ',' {input.sites} {input.neff_mat} > {output.neffsites} || true
        if [ ! -s {output.neffsites} ]; then touch {output.neffsites}; fi
        """

rule calculate_delta:
    input:
       af = "delta_tmp/{delta}/afmatsites.csv",
    output:
       "delta/{delta}.csv"
    params:
        out_prefix = "delta/{delta}",
        script = config["scripts"]["deltaf"],
        pairs = lambda wildcards: DELTA_DICT[wildcards.delta][0],
        reference = lambda wildcards: DELTA_DICT[wildcards.delta][2]
    conda: 
        config["environments"]["polars"]
    benchmark: 
        "benchmarks/calculate_delta/{delta}.log"
    log: 
        "logs/calculate_delta/{delta}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        if [ -s {input.af} ]; then
            python {params.script} \
            {input.af} \
            {params.pairs} \
            --reference_sample {params.reference} \
            --output_prefix {params.out_prefix}
        else
            touch {output}
        fi
        """
