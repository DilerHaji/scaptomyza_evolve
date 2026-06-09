#!/usr/bin/env python3

import polars as pl
import sys
import os

def check_empty(file_path):
    pass

def create_empty_file(file_path):
    with open(file_path, 'w') as f:
        pass

def create_empty_outputs(output_prefix):
    create_empty_file(f"{output_prefix}_filtered_freq.csv")
    create_empty_file(f"{output_prefix}_filtered_counts.csv")
    create_empty_file(f"{output_prefix}_filtered_variants.csv")


def main(df, freq_column, freq_cutoff_lower, freq_cutoff_upper, percent_threshold, min_pop_count, output_prefix):
    required_cols = ["pop", "trt", "gen"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Input file must contain '{col}' column.")

    df_filtered = df.filter(
        (pl.col(freq_column) > float(freq_cutoff_lower)) &
        (pl.col(freq_column) < float(freq_cutoff_upper))
    )

    tmp_filtered_freq = f"{output_prefix}_filtered_freq.csv"
    df_filtered.write_csv(tmp_filtered_freq)
    check_empty(tmp_filtered_freq)

    if len(df_filtered) == 0:
        create_empty_file(f"{output_prefix}_filtered_counts.csv")
        create_empty_file(f"{output_prefix}_filtered_variants.csv")
        return

    rep_counts = (
        df_filtered.group_by(["CHROM", "POS", "pop", "trt"])
        .agg(pl.n_unique("gen").alias("gens_polymorphic"))
    )

    rep_totals = (
        df.group_by(["pop", "trt"])
        .agg(pl.n_unique("gen").alias("gens_total"))
    )

    rep_stats = rep_counts.join(rep_totals, on=["pop", "trt"], how="left")
    rep_stats = rep_stats.with_columns(
        ((pl.col("gens_polymorphic") / pl.col("gens_total")) * 100)
        .round(2)
        .alias("percent_generations_polymorphic")
    )

    rep_summary = (
        rep_stats
        .with_columns(
            (pl.col("percent_generations_polymorphic") >= percent_threshold)
            .alias("rep_meets_threshold")
        )
        .group_by(["CHROM", "POS"])
        .agg(pl.sum("rep_meets_threshold").alias("replicates_meeting_threshold"))
    )

    filtered_counts_file = f"{output_prefix}_filtered_counts.csv"
    rep_summary.write_csv(filtered_counts_file)
    check_empty(filtered_counts_file)

    variants_to_keep = rep_summary.filter(
        pl.col("replicates_meeting_threshold") >= int(min_pop_count)
    ).select(["CHROM", "POS"])

    filtered_df = df_filtered.join(variants_to_keep, on=["CHROM", "POS"], how="inner").sort(["CHROM", "POS"])

    filtered_file = f"{output_prefix}_filtered_variants.csv"
    filtered_df.write_csv(filtered_file)
    check_empty(filtered_file)

if __name__ == "__main__":
    input_file, freq_column, freq_cutoff_lower, freq_cutoff_upper, percent_threshold, min_pop_count, output_prefix = sys.argv[1:]

    freq_cutoff_lower = float(freq_cutoff_lower)
    freq_cutoff_upper = float(freq_cutoff_upper)
    percent_threshold = float(percent_threshold)
    min_pop_count = int(min_pop_count)

    try:
        df = pl.read_csv(input_file)
    except pl.exceptions.NoDataError:
        create_empty_outputs(output_prefix)
        sys.exit(0)

    if len(df) == 0:
        create_empty_outputs(output_prefix)
        sys.exit(0)

    main(df, freq_column, freq_cutoff_lower, freq_cutoff_upper, percent_threshold, min_pop_count, output_prefix)
