#!/bin/bash
set -euo pipefail

usage() {
    cat << EOF
Usage: $0 [OPTIONS] <input.vcf.gz> <output.vcf>

Fix VCF file AD field format for grenedalf compatibility.

ARGUMENTS:
    input.vcf.gz    Input VCF file (can be gzipped or uncompressed)
    output.vcf      Output fixed VCF file (will be uncompressed)

OPTIONS:
    -h, --help      Show this help message
    -k, --keep-multiallelic  Keep multi-allelic sites (default: skip them)
    -v, --verbose   Verbose output

DESCRIPTION:
    This script fixes VCF files to be compatible with grenedalf by:
    1. Changing AD header from Number=1 to Number=R
    2. Converting RD and AD fields to proper AD format (ref,alt)
    3. Skipping multi-allelic sites by default (use -k to keep them)

EXAMPLES:
    $0 input.vcf.gz output.vcf
    $0 -v input.vcf.gz output.vcf
    $0 --keep-multiallelic input.vcf.gz output.vcf

EOF
}

KEEP_MULTIALLELIC=false
VERBOSE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            exit 0
            ;;
        -k|--keep-multiallelic)
            KEEP_MULTIALLELIC=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -*)
            echo "Error: Unknown option $1" >&2
            usage >&2
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

if [[ $# -ne 2 ]]; then
    echo "Error: Exactly 2 arguments required" >&2
    usage >&2
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' does not exist" >&2
    exit 1
fi


if [[ "$INPUT_FILE" == *.gz ]]; then
    CAT_CMD="zcat"
else
    CAT_CMD="cat"
fi

log() {
    if [[ "$VERBOSE" == true ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
    fi
}

TOTAL_VARIANTS=0
SKIPPED_MULTIALLELIC=0
PROCESSED_VARIANTS=0
HEADER_LINES=0

log "Starting VCF processing..."
log "Input: $INPUT_FILE"
log "Output: $OUTPUT_FILE"
log "Keep multi-allelic: $KEEP_MULTIALLELIC"

$CAT_CMD "$INPUT_FILE" | awk -v keep_multiallelic="$KEEP_MULTIALLELIC" -v verbose="$VERBOSE" '
BEGIN {
    FS = OFS = "\t"
    total_variants = 0
    skipped_multiallelic = 0
    processed_variants = 0
    header_lines = 0
}

# Function to log messages (only if verbose)
function log_msg(msg) {
    if (verbose == "true") {
        print "[" strftime("%Y-%m-%d %H:%M:%S") "] " msg > "/dev/stderr"
    }
}

# Process header lines
/^##/ {
    header_lines++
    
    # Fix the AD FORMAT header
    if (/^##FORMAT=<ID=AD/) {
        # Replace Number=1 with Number=R and update description
        gsub(/Number=1/, "Number=R")
        gsub(/Description="[^"]*"/, "Description=\"Allelic depths for the ref and alt alleles in the order listed\"")
        log_msg("Fixed AD header line")
    }
    print
    next
}

# Process column header line
/^#CHROM/ {
    header_lines++
    print
    next
}

# Process variant lines
!/^#/ {
    total_variants++
    
    # Check for multi-allelic sites (ALT field contains comma)
    if ($5 ~ /,/) {
        skipped_multiallelic++
        if (keep_multiallelic == "false") {
            log_msg("Skipping multi-allelic variant at " $1 ":" $2)
            next
        } else {
            log_msg("WARNING: Keeping multi-allelic variant at " $1 ":" $2 " - may cause issues with grenedalf")
        }
    }
    
    processed_variants++
    
    # Process each sample (columns 10 and beyond)
    for (i = 10; i <= NF; i++) {
        # Split the sample data by colons
        split($i, fields, ":")
        
        # Check if we have the expected number of fields (GT:RD:AD:DP:FREQ = 5)
        if (length(fields) >= 5) {
            gt = fields[1]
            rd = fields[2] 
            ad = fields[3]
            dp = fields[4]
            freq = fields[5]
            
            # Create new AD field by combining RD and AD
            new_ad = rd "," ad
            
            # Reconstruct the sample data with the new AD format
            $i = gt ":" rd ":" new_ad ":" dp ":" freq
        } else {
            # If sample has wrong structure, keep as is but log warning
            log_msg("WARNING: Sample " (i-9) " at " $1 ":" $2 " has unexpected format: " $i)
        }
    }
    
    print
    
    # Progress reporting for large files
    if (total_variants % 10000 == 0) {
        log_msg("Processed " total_variants " variants...")
    }
}

END {
    print "# Processing complete" > "/dev/stderr"
    print "# Header lines: " header_lines > "/dev/stderr"
    print "# Total variants: " total_variants > "/dev/stderr"
    print "# Multi-allelic skipped: " skipped_multiallelic > "/dev/stderr"
    print "# Variants processed: " processed_variants > "/dev/stderr"
    
    if (total_variants == 0) {
        print "ERROR: No variants found in input file" > "/dev/stderr"
        exit 1
    }
    
    if (processed_variants == 0) {
        print "ERROR: No variants were processed - check input file format" > "/dev/stderr"
        exit 1
    }
}
' > "$OUTPUT_FILE"

if [[ $? -ne 0 ]]; then
    rm -f "$OUTPUT_FILE"
    exit 1
fi

if command -v du >/dev/null 2>&1; then
    INPUT_SIZE=$(du -h "$INPUT_FILE" | cut -f1)
    OUTPUT_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
fi