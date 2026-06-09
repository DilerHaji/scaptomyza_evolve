# 2. Dynamic Binning of Contigs (50Mb chunks)
BIN_SIZE_TARGET = 50 * 1000 * 1000 

def get_contig_bins(fai_path):
    """
    Reads FAI, groups contigs into bins of approx BIN_SIZE_TARGET.
    Returns: Dict { "0": "chr1", "1": "scaffold_1,scaffold_2" }
    """
    bins = {}
    current_bin_idx = 0
    current_bin_size = 0
    current_contigs = []

    try:
        with open(fai_path, "r") as f:
            for line in f:
                parts = line.split("\t")
                chrom = parts[0]
                length = int(parts[1])

                # If adding this contig exceeds target AND we already have content, close bin
                if (current_bin_size + length > BIN_SIZE_TARGET) and current_contigs:
                    bins[str(current_bin_idx)] = ",".join(current_contigs)
                    current_bin_idx += 1
                    current_bin_size = 0
                    current_contigs = []

                current_contigs.append(chrom)
                current_bin_size += length

        # Add the final bin
        if current_contigs:
            bins[str(current_bin_idx)] = ",".join(current_contigs)
            
    except FileNotFoundError:
        # Warning only, so Snakemake can still parse the file during dry-runs even if FAI missing
        print(f"Warning: FAI file not found at {fai_path}", file=sys.stderr)
        return {}
    
    return bins

# Load Bins Globally (separate from common.smk's CONTIG_BINS to avoid overwriting it)
GREN_CONTIG_BINS = get_contig_bins(config["reference_genome"] + ".fai")
BIN_IDS = list(GREN_CONTIG_BINS.keys())

def get_chroms_for_bin(wildcards):
    return GREN_CONTIG_BINS[wildcards.bin_id]



def get_grenedalf_pairs(wildcards):
    vcf_key = wildcards.grenfst
    pairs = []
    
    # 1. Get replicate lists
    founder_reps = config.get("founder_reps", "")
    founder_gen = config.get("founder_gen", "00")
    mixture_reps = config.get("mixture_reps", "") # e.g. M1, M2, M3, M4

    # Helper to split string lists
    def parse_list(s):
        if isinstance(s, str) and s: return [x.strip() for x in s.split(',') if x.strip()]
        return []

    founder_reps = parse_list(founder_reps)
    mixture_reps = parse_list(mixture_reps)

    # 2. Iterate through your Multiplot Dict
    for analysis, params in GRENFST_MULTIPLOT_DICT.items():
        if params[0] == vcf_key:
            targets = params[1] # e.g. T1, T2...
            refs = params[2]    # e.g. B1, B2...
            gens = params[3]    # e.g. 01, 02...

            for gen in gens:
                # A. Target vs Reference (Treatment vs Baseline)
                for t in targets:
                    for r in refs:
                        pairs.append(f"{t}G{gen}\t{r}G{gen}")

                # B. Standard Robustness (vs G00 Founders)
                if founder_reps:
                    for t in targets:
                        for f in founder_reps:
                            pairs.append(f"{t}G{gen}\t{f}G{founder_gen}")
                    for r in refs:
                        for f in founder_reps:
                            pairs.append(f"{r}G{gen}\t{f}G{founder_gen}")

                # C. Mixture Triplet Logic (Need T vs M and B vs M at CURRENT gen)
                if mixture_reps:
                    for t in targets:
                        for m in mixture_reps:
                            pairs.append(f"{t}G{gen}\t{m}G{gen}")

                    for r in refs:
                        for m in mixture_reps:
                            pairs.append(f"{r}G{gen}\t{m}G{gen}")

    # D. Founder vs Founder pairs (needed for gen=0 PBS anchor in melt script)
    if founder_reps:
        for i, f1 in enumerate(founder_reps):
            for f2 in founder_reps[i+1:]:
                pairs.append(f"{f1}G{founder_gen}\t{f2}G{founder_gen}")

    return "\n".join(sorted(list(set(pairs))))

# def get_grenedalf_pairs(wildcards):
#     vcf_key = wildcards.grenfst
#     pairs = []
#     
#     # 1. Get Founder configuration globally
#     founder_reps = config.get("founder_reps", "")
#     founder_gen = config.get("founder_gen", "00") # Default to 00 if not set
#     
#     if isinstance(founder_reps, str) and founder_reps:
#         founder_reps = [f.strip() for f in founder_reps.split(',') if f.strip()]
#     elif not founder_reps:
#         founder_reps = []
# 
#     # 2. Iterate through your Multiplot Dict to find configurations for this VCF
#     # Structure: Key -> (VCF_Key, Targets, Refs, Gens)
#     for analysis, params in GRENFST_MULTIPLOT_DICT.items():
#         if params[0] == vcf_key:
#             targets = params[1]
#             refs = params[2]
#             gens = params[3]
# 
#             # Generate pairs for every generation requested
#             for gen in gens:
#                 # A. Target vs Reference (for all combinations, matching fst_slope logic)
#                 for t in targets:
#                     for r in refs:
#                         # Assumes naming convention RepGGen (e.g., T1G01)
#                         # Ensure this matches your pool_sizes.tsv!
#                         pairs.append(f"{t}G{gen}\t{r}G{gen}")
# 
#                 # B. Robust Analysis (vs Founders)
#                 if founder_reps:
#                     # Target vs Founder
#                     for t in targets:
#                         for f in founder_reps:
#                             pairs.append(f"{t}G{gen}\t{f}G{founder_gen}")
#                     
#                     # Reference vs Founder
#                     for r in refs:
#                         for f in founder_reps:
#                             pairs.append(f"{r}G{gen}\t{f}G{founder_gen}")
# 
#     # Remove duplicates and sort
#     unique_pairs = sorted(list(set(pairs)))
#     
#     # Return as a newline-separated string for writing to file
#     return "\n".join(unique_pairs)
#     
    
rule grenfst_interval_bin:
    input: 
        grenfst = lambda wildcards: "fvariants/" + GRENFST_DICT[wildcards.grenfst][0] + ".fixed_no_neff.vcf.gz",
        # No change needed here
    output:
        csv = "grenfst/fst_queue_bins/{grenfst}/{bin_id}/fst.csv"
    params: 
        run_dir = "grenfst/fst_queue_bins/{grenfst}/{bin_id}",
        method = config["grenedalf"]["method"],
        pool_sizes_file = "pool_sizes.tsv", 
        window_width = lambda wildcards: GRENFST_DICT[wildcards.grenfst][1],
        window_stride = lambda wildcards: GRENFST_DICT[wildcards.grenfst][2],
        filter_total_snp_min_frequency = lambda wildcards: GRENFST_DICT[wildcards.grenfst][3],
        pair_list = get_grenedalf_pairs,
        # Get the comma-separated string, replace comma with newline for the file
        chrom_list_newline = lambda wildcards: GREN_CONTIG_BINS[wildcards.bin_id].replace(",", "\n")
    resources:
        resources=config["default_resources"]
    conda:
        config["environments"]["grenedalf"]
    shell:
        """
        VCF_ABS=$(readlink -f {input.grenfst})
        POOLS_ABS=$(readlink -f {params.pool_sizes_file})

        [[ -d {params.run_dir} ]] || mkdir -p {params.run_dir}
        cd {params.run_dir}

        echo -e "{params.pair_list}" > pairs_raw.txt
        # Filter pairs to only those where both samples exist in the VCF
        zcat "$VCF_ABS" | grep "^#CHROM" | tr '\t' '\n' | tail -n +10 > vcf_samples.txt
        awk 'NR==FNR{{s[$1]=1; next}} ($1 in s) && ($2 in s)' vcf_samples.txt pairs_raw.txt > pairs_to_calculate.txt

        # --- NEW LOGIC: Write regions to file ---
        echo -e "{params.chrom_list_newline}" > region_list.txt

        rm -f fst.csv fst_tmp.csv

        set +e
        grenedalf fst \
        --vcf-path "$VCF_ABS" \
        --filter-total-snp-min-frequency {params.filter_total_snp_min_frequency} \
        --window-type interval \
        --window-interval-width {params.window_width} \
        --window-interval-stride {params.window_stride} \
        --window-average-policy valid-snps \
        --method unbiased-nei \
        --pool-sizes "$POOLS_ABS" \
        --comparand-list pairs_to_calculate.txt \
        --filter-region-list region_list.txt \
        --allow-file-overwriting \
        --file-suffix _tmp
        GREN_EXIT=$?
        set -e

        if [[ -f fst_tmp.csv ]]; then
            if [[ $GREN_EXIT -eq 0 ]]; then
                mv fst_tmp.csv fst.csv  # complete run
            else
                rm -f fst_tmp.csv
                echo "ERROR: grenedalf was interrupted mid-run (exit $GREN_EXIT)" >&2
                exit 1  # partial output — must rerun
            fi
        else
            touch fst.csv  # no output file: empty region or non-fatal error
        fi
        """


# ==============================================================================
# 3. GATHER: Merge Bin FST results
# ==============================================================================
rule grenfst_queue:
    input: 
        # EXPAND over BIN_IDS instead of get_chromosomes()
        lambda wildcards: expand(
            "grenfst/fst_queue_bins/{grenfst}/{bin_id}/fst.csv", 
            grenfst=wildcards.grenfst, 
            bin_id=BIN_IDS
        )
    output:
        "grenfst/fst_queue/{grenfst}.csv"
    params: 
        output_dir="grenfst/fst_queue"
    resources:
        resources=config["default_resources"]
    resources:
        resources=config["default_resources"]
    benchmark: 
        "benchmarks/grenfst/{grenfst}.log"
    log: 
        "logs/grenfst/{grenfst}.log"
    shell:
        """
        [[ -d {params.output_dir} ]] || mkdir -p {params.output_dir}

        # --- ROBUST MERGE LOGIC ---
        
        # 1. Find the first file with content (size > 0) to get the header
        HEADER_FOUND=0
        for f in {input}; do
            if [[ -s "$f" ]]; then
                head -n 1 "$f" > {output}
                HEADER_FOUND=1
                break
            fi
        done
        
        # Guard: If no files had content (entire genome empty?), touch output
        if [[ $HEADER_FOUND -eq 0 ]]; then
            touch {output}
            exit 0
        fi

        # 2. Append content from ALL files (skipping their headers)
        for f in {input}; do
            if [[ -s "$f" ]]; then
                tail -n +2 "$f" >> {output}
            fi
        done
        """






##########################################
# Plot the slope of FST change as a Manhattan plot
# (Unchanged, assuming input CSV format is consistent)
##########################################

        
rule grenfst_queue_slope_plot:
    input: 
        "grenfst/lmm_results/{grenfst_multiplot}/consensus_lmm.csv"
    output: 
        "grenfst/fst_queue_slope/{grenfst_multiplot}/plot.png"
    params:
        script = "scripts/fst_slope_plot_lmm.py"
    conda: 
        config["environments"]["polars"]
    benchmark: 
        "benchmarks/grenfst_queue_slope/{grenfst_multiplot}.log"
    log: 
        "logs/grenfst_queue_slope/{grenfst_multiplot}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        python {params.script} \
        --input {input} \
        --output {output} \
        --top-n 10 \
        --p-cutoff 0.05 \
        --ymax 0.0025
        """
        


##########################################
# Find outliers from the FST slope
# (Unchanged)
##########################################
rule grenfst_queue_slope_outliers:
    input: 
        csv = "grenfst/lmm_results/{grenfst_multiplot}/consensus_lmm.csv",
        repeats = "../popgen/repeat_content/masked_sfla_v2/sfla_v2.fa.out",
        map = "chrom_map.tsv"
    output:
        csv = "grenfst/fst_queue_slope/{grenfst_multiplot}/outliers.csv",
        plot = "grenfst/fst_queue_slope/{grenfst_multiplot}/manhattan.png"
    params: 
        script = config["scripts"]["fst_slope_peaks"],
    conda:
        config["environments"]["polars"]
    benchmark: 
        "benchmarks/grenfst_queue_slope/{grenfst_multiplot}.log"
    log: 
        "logs/grenfst_queue_slope/{grenfst_multiplot}.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        python {params.script} \
        --input {input.csv} \
        --repeats {input.repeats} \
        --chrom-map {input.map} \
        --output {output.plot} \
        --out-csv {output.csv} \
        --slope-col slope_median \
        --error-col z_score_sd \
        --max-error 1 \
        --z-score 3.5 \
        --extend-bp 500 \
        --top-n 22 \
        --max-repeat 0 \
        --fst-window-size 1000 \
        --neighbor-window 100
        """

##########################################
# Plot the FSTs themselves for a series of samples 
# UPDATED: Plots Mean + Variance Ribbon
##########################################
rule grenfst_queue_plot:
    input: 
        csv = lambda wildcards: "grenfst/fst_queue/" + GRENFST_MULTIPLOT_DICT[wildcards.grenfst_multiplot][0] + ".csv",
        #slope_peaks = "grenfst/fst_queue_slope/{grenfst_multiplot}/outliers.csv"
        peaks = "final_aligned_landscape_v4_candidate_sites_detailed.tsv"
    output:
       "grenfst/fst_queue_plot/{grenfst_multiplot}/multiplot.png"
    params: 
        target_reps = lambda wildcards: ",".join(GRENFST_MULTIPLOT_DICT[wildcards.grenfst_multiplot][1]), 
        ref_reps    = lambda wildcards: ",".join(GRENFST_MULTIPLOT_DICT[wildcards.grenfst_multiplot][2]), 
        script = config["scripts"]["multiFST_plot"]
    conda:
        config["environments"]["polars"]
    benchmark: 
        "benchmarks/grenfst_queue_plot/{grenfst_multiplot}.log"
    log: 
        "logs/grenfst_queue_plot/{grenfst_multiplot}.log"
    resources:
        resources=config["default_resources_24cpus"]
    shell:
        """
        python {params.script} \
        --input {input.csv} \
        --output {output} \
        --target-reps {params.target_reps} \
        --ref-reps {params.ref_reps} \
        --peaks {input.peaks} \
        --top-n 10 \
        --ymax 0.02 \
        --show-xlabels
        """




# ##########################################
# # SCATTER: Run FST per chromosome (INTERVAL MODE)
# # Optimization: Uses --comparand-list to calculate ONLY necessary pairs
# # CHANGED: Uses --window-type interval (bp) to fix "diagonal tilt" artifacts
# ##########################################
# rule grenfst_interval_chrom:
#     input: 
#         grenfst=lambda wildcards: "fvariants/" + GRENFST_DICT[wildcards.grenfst][0] + ".fixed_no_neff.vcf.gz",
#         scaffolds='chromosome_names.txt' 
#     output:
#         # We keep the path "fst_queue_chroms" for compatibility with downstream rules,
#         # even though we are technically doing intervals now.
#         csv = temp("grenfst/fst_queue_chroms/{grenfst}/{chrom}/fst.csv")
#     params: 
#         run_dir = "grenfst/fst_queue_chroms/{grenfst}/{chrom}",
#         method=config["grenedalf"]["method"],
#         pool_sizes_file="pool_sizes.tsv", 
#         # Map the dictionary values to Width (bp) and Stride (bp)
#         window_width = lambda wildcards: GRENFST_DICT[wildcards.grenfst][1],
#         window_stride = lambda wildcards: GRENFST_DICT[wildcards.grenfst][2],
#         filter_total_snp_min_frequency = lambda wildcards: GRENFST_DICT[wildcards.grenfst][3],
#         pair_list = get_grenedalf_pairs
#     conda:
#         config["environments"]["grenedalf"]
#     resources:
#         resources=config["default_resources"]
#     shell:
#         """
#         VCF_ABS=$(readlink -f {input.grenfst})
#         POOLS_ABS=$(readlink -f {params.pool_sizes_file})
# 
#         [[ -d {params.run_dir} ]] || mkdir -p {params.run_dir}
#         cd {params.run_dir}
# 
#         echo -e "{params.pair_list}" > pairs_to_calculate.txt
# 
#         # --- GUARD LOGIC ---
#         # 1. Run grenedalf with INTERVAL windows
#         # 2. || true : prevents snakemake from failing if no SNPs found in chrom
#         grenedalf fst \
#         --vcf-path "$VCF_ABS" \
#         --filter-total-snp-min-frequency {params.filter_total_snp_min_frequency} \
#         --window-type interval \
#         --window-interval-width {params.window_width} \
#         --window-interval-stride {params.window_stride} \
#         --window-average-policy valid-snps \
#         --method unbiased-nei \
#         --pool-sizes "$POOLS_ABS" \
#         --comparand-list pairs_to_calculate.txt \
#         --filter-region {wildcards.chrom} || true
# 
#         # 3. Ensure the output file exists so Snakemake doesn't complain "Missing output file"
#         if [[ ! -f fst.csv ]]; then
#             touch fst.csv
#         fi
#         """
# 
# ##########################################
# # GATHER: Merge chromosomal FST results
# # Combines the individual CSVs into the final one at the original location
# # (Unchanged, compatible with interval output)
# ##########################################
# rule grenfst_queue:
#     input: 
#         lambda wildcards: expand(
#             "grenfst/fst_queue_chroms/{grenfst}/{chrom}/fst.csv", 
#             grenfst=wildcards.grenfst, 
#             chrom=get_chromosomes()
#         )
#     output:
#         "grenfst/fst_queue/{grenfst}.csv"
#     params: 
#         output_dir="grenfst/fst_queue"
#     resources:
#         resources=config["default_resources"]
#     benchmark: 
#         "benchmarks/grenfst/{grenfst}.log"
#     log: 
#         "logs/grenfst/{grenfst}.log"
#     shell:
#         """
#         [[ -d {params.output_dir} ]] || mkdir -p {params.output_dir}
# 
#         # --- ROBUST MERGE LOGIC ---
#         
#         # 1. Find the first file with content (size > 0) to get the header
#         HEADER_FOUND=0
#         for f in {input}; do
#             if [[ -s "$f" ]]; then
#                 head -n 1 "$f" > {output}
#                 HEADER_FOUND=1
#                 break
#             fi
#         done
#         
#         # Guard: If no files had content (entire genome empty?), touch output
#         if [[ $HEADER_FOUND -eq 0 ]]; then
#             touch {output}
#             exit 0
#         fi
# 
#         # 2. Append content from ALL files (skipping their headers)
#         for f in {input}; do
#             if [[ -s "$f" ]]; then
#                 tail -n +2 "$f" >> {output}
#             fi
#         done
#         """