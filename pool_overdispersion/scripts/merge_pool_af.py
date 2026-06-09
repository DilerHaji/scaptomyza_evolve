"""
Merge per-sample AF tables into one wide table.

Each sample's AF file has: chrom, pos, ref, depth, ref_count, alt_allele, alt_count, alt_af

Output: chrom, pos, ref, {sample}_depth, {sample}_ref_count, {sample}_alt_count, {sample}_af
for all samples.

Sites are joined on (chrom, pos). Sites missing in a sample get NaN.
"""

import polars as pl
import gzip
import os


def main():
    sample_files = snakemake.input
    samples = snakemake.params.samples
    output_file = snakemake.output[0]

    dfs = []
    for sample, fpath in zip(samples, sample_files):
        df = pl.read_csv(
            fpath,
            separator="\t",
            has_header=True,
            dtypes={"chrom": pl.Utf8, "pos": pl.Int64},
        )
        # Rename columns with sample prefix
        df = df.rename({
            "depth": f"{sample}_depth",
            "ref_count": f"{sample}_ref_count",
            "alt_count": f"{sample}_alt_count",
            "alt_af": f"{sample}_af",
        })
        # Keep chrom, pos, ref, and the renamed columns
        df = df.select([
            "chrom", "pos", "ref",
            f"{sample}_depth", f"{sample}_ref_count",
            f"{sample}_alt_count", f"{sample}_af",
        ])
        dfs.append(df)

    # Join all on (chrom, pos, ref)
    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.join(df, on=["chrom", "pos", "ref"], how="outer")

    # Sort by chrom, pos
    merged = merged.sort(["chrom", "pos"])

    # Write compressed
    merged.write_csv(output_file, separator="\t")


main()
