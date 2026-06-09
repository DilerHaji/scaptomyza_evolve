import os
import glob
import functools
import math

# ==============================================================================
# CONFIG & SETTINGS
# ==============================================================================

# Target: 10Mb
BIN_SIZE_TARGET = 1 * 1000 * 1000 

# ==============================================================================
# HYBRID BINNING LOGIC (Descriptive IDs)
# ==============================================================================

def get_hybrid_bins_bed(fai_path, target_size, map_output_file="bin_key_map.tsv"):
    """
    1. Reads FAI.
    2. Splits large contigs -> IDs: {Chrom}_part{N}
    3. Merges small contigs -> IDs: multi_{N}
    4. Writes a TSV map file for user verification (only if not already present,
       to avoid Lustre lock contention when many jobs start simultaneously).
    """
    bins = {}

    # buffers for small contigs
    current_bed_lines = []
    current_size = 0
    multi_bin_counter = 1

    try:
        with open(fai_path, "r") as f_in:
            for line in f_in:
                parts = line.split("\t")
                chrom = parts[0]
                length = int(parts[1])

                # --- CASE 1: Large Contig (Split it) ---
                if length > target_size:
                    # First: Flush any pending small buffer
                    if current_bed_lines:
                        bin_id = f"multi_{multi_bin_counter:03d}"
                        bins[bin_id] = current_bed_lines
                        multi_bin_counter += 1
                        current_bed_lines, current_size = [], 0

                    # Now split the large one
                    num_chunks = math.ceil(length / target_size)
                    for i in range(num_chunks):
                        start = i * target_size
                        end = min((i + 1) * target_size, length)

                        # ID format: chr1_part01
                        bin_id = f"{chrom}_part{i+1:02d}"
                        bed_line = f"{chrom}\t{start}\t{end}"

                        bins[bin_id] = [bed_line]

                # --- CASE 2: Small Contig (Group it) ---
                else:
                    # If full, flush buffer
                    if current_size + length > target_size:
                        bin_id = f"multi_{multi_bin_counter:03d}"
                        bins[bin_id] = current_bed_lines
                        multi_bin_counter += 1
                        current_bed_lines, current_size = [], 0

                    # Add to buffer
                    current_bed_lines.append(f"{chrom}\t0\t{length}")
                    current_size += length

            # Flush final buffer
            if current_bed_lines:
                bin_id = f"multi_{multi_bin_counter:03d}"
                bins[bin_id] = current_bed_lines

    except FileNotFoundError:
        # Warning only (for dry runs before indexing)
        return {}

    # Write the map file only if it doesn't already exist, to avoid concurrent
    # write contention on Lustre when many SLURM jobs start simultaneously.
    if bins and not os.path.exists(map_output_file):
        try:
            with open(map_output_file, "w") as f_out:
                f_out.write("Bin_ID\tType\tRegion_Count\tTotal_Length\tRegions (First 3...)\n")
                for bin_id, bed_lines in bins.items():
                    total_len = sum(int(l.split("\t")[2]) - int(l.split("\t")[1]) for l in bed_lines)
                    btype = "Split" if len(bed_lines) == 1 and "_part" in bin_id else "Multi"
                    f_out.write(f"{bin_id}\t{btype}\t{len(bed_lines)}\t{total_len}\t{';'.join(bed_lines[:3])}...\n")
        except OSError:
            pass  # Non-fatal: map file is only for debugging

    return bins

# Load Bins Globally
# This writes 'bin_key_map.tsv' to your folder every time Snakemake parses this file
CONTIG_BINS = get_hybrid_bins_bed(config["reference_genome"] + ".fai", BIN_SIZE_TARGET)

# Update the global entry list function
@functools.lru_cache(maxsize=None)
def get_entries():
    return list(CONTIG_BINS.keys())



################################################################################
##### General Functions ####
################################################################################

def get_key(mapping, key):
    return mapping.get(key, [])

################################################################################
##### Unpacking mpileup ####
################################################################################

rule unpack: 
    output:
        "unpacked/{pileup}.mpileup"
    params: 
        pileup=lambda wildcards: os.path.join(config["input_dir"], "{wildcards.pileup}.mpileup.gz")
    benchmark: 
        "benchmarks/fst/{pileup}_unpack.log"
    log: 
        "logs/fst/{pileup}_unpack.log"
    resources:
        resources=config["default_resources"]
    shell:
        """
        if [[ -f {output} ]]; then 
            echo "pileup already unpacked"
        else
            gunzip -c {params.pileup} > {output}
        fi
        """

################################################################################
##### MAKE BED ####
################################################################################

checkpoint create_bed:
    input:
        config["reference_genome"]
    output:
        bed = "reference.bed",
        scaffolds = 'chromosome_names.txt'
    params:
        chunk_size = config["create_bed"]["chunk_size"],
        top_n_chrom = config["create_bed"]["top_n_chrom"],
        get_chroms = config["scripts"]["get_chroms"]
    benchmark: 
        "benchmarks/create_bed.log"
    log: 
        "logs/create_bed.log"
    shell:
        """
        python {params.get_chroms} {input} {params.top_n_chrom} {output.scaffolds}
        python scripts/generate_bed.py --ref {input} --output reference_all.bed --chunk_size {params.chunk_size}
        grep -F -f {output.scaffolds} reference_all.bed > {output.bed}
        """

localrules: create_bed




# Fallback for scaffold-based parallelism if needed
def get_chromosomes():
    checkpoint_output = checkpoints.create_bed.get().output.scaffolds
    with open(checkpoint_output) as f:
        chrom = [line.strip() for line in f]
    return chrom
    

# Helper to safely retrieve contigs for a bin
def get_contigs_from_bin(wildcards):
    try:
        # Get the list of bed lines for this bin ID
        lines = CONTIG_BINS[wildcards.entry]
        # Extract just the chromosome names (first column)
        names = [line.split("\t")[0] for line in lines]
        # Return comma-separated string
        return ",".join(names)
    except KeyError:
        # This will print to your terminal if the error happens again, showing us what keys DO exist
        print(f"\n[ERROR] Could not find bin '{wildcards.entry}' in CONTIG_BINS.")
        print(f"Total bins loaded: {len(CONTIG_BINS)}")
        print(f"Sample keys available: {list(CONTIG_BINS.keys())[:5]}\n")
        raise