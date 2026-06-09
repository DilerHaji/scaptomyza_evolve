# ==============================================================================
# SETUP & HELPERS
# ==============================================================================

# 1. Load FAI lengths to determine chromosome sizes
def load_fai_lengths(fai_path):
    lengths = {}
    try:
        with open(fai_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    lengths[parts[0]] = int(parts[1])
    except FileNotFoundError:
        pass # Allow dry-runs to proceed if file missing
    return lengths

# Initialize global variable
FAI_LENGTHS = load_fai_lengths(config["reference_genome"] + ".fai")

# 2. Helper to format region strings from bin dictionary
def get_regions_param(wildcards):
    """
    Returns semicolon-separated regions (chr:start:end) for the python script.
    Handles both comma-separated strings (Whole Chroms) and BED lists (Split Bins).
    """
    if wildcards.entry not in CONTIG_BINS:
        return ""
    
    bin_data = CONTIG_BINS[wildcards.entry]
    regions = []
    
    # Case A: String "chr1,chr2" -> Whole Chromosomes
    if isinstance(bin_data, str):
        chroms = bin_data.split(",")
        for c in chroms:
            c = c.strip()
            # Lookup total length from FAI
            length = FAI_LENGTHS.get(c, 0)
            if length > 0:
                regions.append(f"{c}:0:{length}")

    # Case B: List of BED strings ["chr1\t0\t1000"] -> Specific Intervals
    elif isinstance(bin_data, list):
        for line in bin_data:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                c, s, e = parts[0], parts[1], parts[2]
                regions.append(f"{c}:{s}:{e}")
            
    return ";".join(regions)

# ==============================================================================
# RULES FOR PBS MIXTURE ANALYSIS
# ==============================================================================

# 1. MELT (Raw Calculation)
rule grenfst_div_melt_mixture_bin:
    input:
        csv = lambda wildcards: "grenfst/fst_queue/" + GRENFST_MULTIPLOT_DICT[wildcards.grenfst_multiplot][0] + ".csv"
    output:
        csv = "grenfst/divergence_mixture_lmm/bins/{grenfst_multiplot}/{entry}.csv"
    params:
        script = "scripts/melt_pbs_combinatorial.py",
        target_reps = lambda w: ",".join(GRENFST_MULTIPLOT_DICT[w.grenfst_multiplot][1]),
        ref_reps    = lambda w: ",".join(GRENFST_MULTIPLOT_DICT[w.grenfst_multiplot][2]),
        
        # Outgroup Configuration
        outgroup_reps = lambda w: config["mixture_reps"] if isinstance(config.get("mixture_reps"), str) else ",".join(config.get("mixture_reps", [])),
        outgroup_gen = "dynamic",

        generations = lambda w: ",".join(GRENFST_MULTIPLOT_DICT[w.grenfst_multiplot][3]),
        regions = get_regions_param,

        # Founder generation anchor (optional): F reps used to add PBS≈0 at gen 00
        founder_reps = config.get("founder_reps", ""),
        founder_gen  = config.get("founder_gen",  "00"),
        # Mode: per_rep (rigorous matched-pair with pooled B/M) or combinatorial (legacy)
        mode = config.get("pbs_melt_mode", "per_rep"),
    conda:
        config["environments"]["polars"]
    resources:
        resources=config.get("default_resources", "")
    shell:
        """
        python {params.script} \
        --input {input.csv} \
        --output {output.csv} \
        --target-reps "{params.target_reps}" \
        --ref-reps "{params.ref_reps}" \
        --outgroup-reps "{params.outgroup_reps}" \
        --generations "{params.generations}" \
        --outgroup-gen "{params.outgroup_gen}" \
        --founder-reps "{params.founder_reps}" \
        --founder-gen  "{params.founder_gen}" \
        --regions "{params.regions}" \
        --mode "{params.mode}"
        """

# 2. MEDIAN AGGREGATION (Reduce)
rule grenfst_div_calc_medians_mixture:
    input:
        lambda wildcards: expand("grenfst/divergence_mixture_lmm/bins/{grenfst_multiplot}/{entry}.csv",
                                 grenfst_multiplot=wildcards.grenfst_multiplot,
                                 entry=get_entries())
    output:
        medians = "grenfst/divergence_mixture_lmm/medians/{grenfst_multiplot}.csv"
    params:
        script = "scripts/calc_chrom_medians.py"
    conda:
        config["environments"]["polars"]
    resources:
        resources=config.get("default_resources_40", "")
    shell:
        """
        python {params.script} \
        --inputs {input} \
        --output {output.medians}
        """

# 3. JOIN PBE (Map)
rule grenfst_div_join_pbe:
    input:
        csv = "grenfst/divergence_mixture_lmm/bins/{grenfst_multiplot}/{entry}.csv",
        medians = "grenfst/divergence_mixture_lmm/medians/{grenfst_multiplot}.csv"
    output:
        csv = "grenfst/divergence_mixture_lmm/bins_pbe/{grenfst_multiplot}/{entry}.csv"
    params:
        script = "scripts/join_pbe_medians.py"
    conda:
        config["environments"]["polars"]
    resources:
        resources=config.get("default_resources", "")
    shell:
        """
        python {params.script} \
        --input {input.csv} \
        --medians {input.medians} \
        --output {output.csv}
        """

# 4. FIT LMM (Map)
rule grenfst_div_lmm_fit_mixture_bin:
    input:
        csv = "grenfst/divergence_mixture_lmm/bins_pbe/{grenfst_multiplot}/{entry}.csv"
    output:
        stats = "grenfst/divergence_mixture_lmm/stats_bins/{grenfst_multiplot}/{entry}.csv"
    params:
        script = "scripts/fit_pbs_lmm.py",
        chroms = get_contigs_from_bin
    conda:
        config["environments"]["polars"]
    resources:
        resources=config["default_resources_10cpus"]
    threads: 4
    shell:
        """
        python {params.script} \
        --input {input.csv} \
        --output {output.stats} \
        --chroms "{params.chroms}" \
        --threads {threads} 
        """

# 5. FINAL AGGREGATE
rule grenfst_div_lmm_aggregate_mixture:
    input:
        lambda wildcards: expand("grenfst/divergence_mixture_lmm/stats_bins/{grenfst_multiplot}/{entry}.csv",
                                 grenfst_multiplot=wildcards.grenfst_multiplot,
                                 entry=get_entries())
    output:
        final = "grenfst/divergence_mixture_lmm/{grenfst_multiplot}/divergence_stats.csv"
    params:
        script = "scripts/aggregate_pbs_bins.py"
    conda:
        config["environments"]["polars"]
    resources:
        resources=config.get("default_resources", "")
    shell:
        """
        python {params.script} \
        --inputs {input} \
        --output {output.final} 
        """
        
# # rules/pbs_mixture.smk
# 
# # Load FAI lengths globally to look up chromosome sizes
# def load_fai_lengths(fai_path):
#     lengths = {}
#     try:
#         with open(fai_path, 'r') as f:
#             for line in f:
#                 parts = line.strip().split('\t')
#                 if len(parts) >= 2:
#                     lengths[parts[0]] = int(parts[1])
#     except FileNotFoundError:
#         pass # Allow dry-runs to proceed if file missing
#     return lengths
# 
# # Load lengths immediately
# FAI_LENGTHS = load_fai_lengths(config["reference_genome"] + ".fai")
# 
# def get_regions_param(wildcards):
#     """
#     Converts bin definitions into semicolon-separated regions (chr:start:end).
#     Handles both formats:
#     1. String: "chr1,chr2" (from grenedalf_interval.smk) -> Uses full chrom length from FAI
#     2. List: ["chr1\t0\t100"] (from common.smk)
#     """
#     if wildcards.entry not in CONTIG_BINS:
#         return ""
#     
#     bin_data = CONTIG_BINS[wildcards.entry]
#     regions = []
# 
#     # Case 1: Comma-separated string (Whole Chromosomes)
#     if isinstance(bin_data, str):
#         chroms = bin_data.split(",")
#         for c in chroms:
#             c = c.strip()
#             if not c: continue
#             # Lookup total length from FAI
#             length = FAI_LENGTHS.get(c, 0)
#             if length > 0:
#                 regions.append(f"{c}:0:{length}")
#             else:
#                 # If FAI lookup fails, skip or warn. 
#                 pass
# 
#     # Case 2: List of BED strings (Specific Intervals)
#     elif isinstance(bin_data, list):
#         for line in bin_data:
#             parts = line.strip().split("\t")
#             if len(parts) >= 3:
#                 c, s, e = parts[0], parts[1], parts[2]
#                 regions.append(f"{c}:{s}:{e}")
#             
#     return ";".join(regions)
# 
# def get_chroms_from_bin(wildcards):
#     """
#     Extracts comma-separated list of chromosomes for the bin.
#     """
#     if wildcards.entry not in CONTIG_BINS:
#         return ""
#     
#     bin_data = CONTIG_BINS[wildcards.entry]
#     chroms = set()
#     
#     if isinstance(bin_data, str):
#         # format: "chr1,chr2"
#         for c in bin_data.split(","):
#             if c.strip(): chroms.add(c.strip())
#             
#     elif isinstance(bin_data, list):
#         # format: ["chr1\t0\t100", ...]
#         for line in bin_data:
#             c = line.strip().split("\t")[0]
#             chroms.add(c)
#             
#     return ",".join(sorted(list(chroms)))
# 
# # ==============================================================================
# # TRACK 2 (MIXTURE): Combinatorial Divergence Analysis (Mixture as Outgroup)
# # ==============================================================================
# 
# rule grenfst_div_melt_mixture_bin:
#     input:
#         csv = lambda wildcards: "grenfst/fst_queue/" + GRENFST_MULTIPLOT_DICT[wildcards.grenfst_multiplot][0] + ".csv"
#     output:
#         csv = "grenfst/divergence_mixture_lmm/bins/{grenfst_multiplot}/{entry}.csv"
#     params:
#         script = "scripts/melt_pbs_combinatorial.py",
#         # 1. Target (T)
#         target_reps = lambda w: ",".join(GRENFST_MULTIPLOT_DICT[w.grenfst_multiplot][1]),
#         # 2. Reference (B)
#         ref_reps    = lambda w: ",".join(GRENFST_MULTIPLOT_DICT[w.grenfst_multiplot][2]),
#         # 3. Outgroup (M) -> Passed into the 'founder' argument of the script
#         founder_reps = lambda w: config["mixture_reps"] if isinstance(config.get("mixture_reps"), str) else ",".join(config.get("mixture_reps", [])),
#         # 4. Generations
#         generations = lambda w: ",".join(GRENFST_MULTIPLOT_DICT[w.grenfst_multiplot][3]),
#         # 5. Tell script to use current generation for the outgroup
#         founder_gen = "dynamic",
#         regions = get_regions_param 
#     conda:
#         config["environments"]["polars"]
#     resources:
#         resources=config.get("default_resources", "")
#     log:
#         "logs/divergence_mixture_lmm/melt_bins/{grenfst_multiplot}_{entry}.log"
#     shell:
#         """
#         python {params.script} \
#         --input {input.csv} \
#         --output {output.csv} \
#         --target-reps "{params.target_reps}" \
#         --ref-reps "{params.ref_reps}" \
#         --founder-reps "{params.founder_reps}" \
#         --generations "{params.generations}" \
#         --founder-gen "{params.founder_gen}" \
#         --regions "{params.regions}" \
#         > {log} 2>&1
#         """
# 
# # 2. FIT STEP
# rule grenfst_div_lmm_fit_mixture_bin:
#     input:
#         csv = "grenfst/divergence_mixture_lmm/bins/{grenfst_multiplot}/{entry}.csv"
#     output:
#         stats = "grenfst/divergence_mixture_lmm/stats_bins/{grenfst_multiplot}/{entry}.csv"
#     params:
#         script = "scripts/fit_pbs_lmm.py",
#         chroms = get_chroms_from_bin
#     conda:
#         config["environments"]["polars"]
#     resources:
#         resources=config["default_resources_10cpus"]
#     threads: 4
#     log:
#         "logs/divergence_mixture_lmm/fit_bins/{grenfst_multiplot}_{entry}.log"
#     shell:
#         """
#         python {params.script} \
#         --input {input.csv} \
#         --output {output.stats} \
#         --chroms "{params.chroms}" \
#         --threads {threads} \
#         --plot-diagnostics \
#         > {log} 2>&1
#         """
# 
# # 3. AGGREGATE STEP
# rule grenfst_div_lmm_aggregate_mixture:
#     input:
#         lambda wildcards: expand("grenfst/divergence_mixture_lmm/stats_bins/{{grenfst_multiplot}}/{entry}.csv", entry=get_entries())
#     output:
#         final = "grenfst/divergence_mixture_lmm/{grenfst_multiplot}/divergence_stats.csv"
#     params:
#         script = "scripts/aggregate_pbs_bins.py"
#     conda:
#         config["environments"]["polars"]
#     resources:
#         resources=config.get("default_resources", "")
#     log:
#         "logs/divergence_mixture_lmm/aggregate_{grenfst_multiplot}.log"
#     shell:
#         """
#         python {params.script} \
#         --inputs {input} \
#         --output {output.final} \
#         > {log} 2>&1
#         """
