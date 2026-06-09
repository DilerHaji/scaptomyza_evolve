"""
Extract allele frequencies from a samtools mpileup file.

For each site in the pileup, counts reference and alternate allele reads
and computes allele frequency. Only biallelic SNPs are retained.

Input:  snakemake.input.mpileup  — samtools mpileup output (single sample)
Output: snakemake.output.af      — TSV: chrom, pos, ref, depth, ref_count, alt_count, alt_af
"""

import sys
import re


def parse_pileup_bases(ref, bases_str):
    """Parse mpileup bases string, return counts of A, C, G, T."""
    # Remove read start/end markers and mapping quality
    bases = re.sub(r"\^.", "", bases_str)  # remove read start + mapq char
    bases = bases.replace("$", "")         # remove read end
    # Remove indel sequences: +3ACG or -2AT
    bases = re.sub(r"[+-](\d+)", lambda m: f"[{m.group(1)}]", bases)
    # Now remove the bracketed number and that many following chars
    cleaned = []
    i = 0
    while i < len(bases):
        if bases[i] == "[":
            end = bases.index("]", i)
            n = int(bases[i+1:end])
            i = end + 1 + n  # skip the indel bases
        else:
            cleaned.append(bases[i])
            i += 1
    bases = "".join(cleaned)

    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    for b in bases:
        if b in ".,":
            counts[ref.upper()] += 1
        elif b.upper() in counts:
            counts[b.upper()] += 1
        # skip * (deletion) and N

    return counts


def main():
    input_file = snakemake.input.mpileup
    output_file = snakemake.output.af

    with open(input_file) as fin, open(output_file, "w") as fout:
        fout.write("chrom\tpos\tref\tdepth\tref_count\talt_allele\talt_count\talt_af\n")

        for line in fin:
            fields = line.strip().split("\t")
            if len(fields) < 6:
                continue

            chrom = fields[0]
            pos = fields[1]
            ref = fields[2].upper()
            depth = int(fields[3])

            if depth == 0 or ref not in "ACGT":
                continue

            bases_str = fields[4]
            counts = parse_pileup_bases(ref, bases_str)

            # Find the most common non-ref allele
            ref_count = counts[ref]
            alt_alleles = {k: v for k, v in counts.items() if k != ref and v > 0}

            if not alt_alleles:
                # Monomorphic at this site — still record it
                fout.write(f"{chrom}\t{pos}\t{ref}\t{depth}\t{ref_count}\t.\t0\t0.0\n")
                continue

            alt_allele = max(alt_alleles, key=alt_alleles.get)
            alt_count = alt_alleles[alt_allele]
            total = ref_count + alt_count

            if total == 0:
                continue

            alt_af = alt_count / total
            fout.write(f"{chrom}\t{pos}\t{ref}\t{depth}\t{ref_count}\t{alt_allele}\t{alt_count}\t{alt_af:.6f}\n")


main()
