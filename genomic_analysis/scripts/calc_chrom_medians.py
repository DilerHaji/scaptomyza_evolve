#!/usr/bin/env python3
import polars as pl
import argparse
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True, help="List of all bin CSV files")
    parser.add_argument("--output", required=True, help="Output medians CSV")
    args = parser.parse_args()

    schema_overrides = {
        "chrom": pl.String,
        "gen": pl.Int64,       # Force Integer
        "lineage": pl.String,
        "pbs_target": pl.Float64,
        "pbs_ref": pl.Float64,
        "T_bm": pl.Float64,
        "T_tm": pl.Float64
    }

    cols_needed = list(schema_overrides.keys())

    try:
        q = pl.concat([
            pl.scan_csv(f, schema_overrides=schema_overrides).select(cols_needed) 
            for f in args.inputs
        ])

        medians = q.group_by(["gen", "lineage", "chrom"]).agg([
            pl.col("pbs_target").median().alias("med_pbs_target"),
            pl.col("pbs_ref").median().alias("med_pbs_ref"),
            pl.col("T_bm").median().alias("med_T_bm"),
            pl.col("T_tm").median().alias("med_T_tm")
        ]).collect()

        if medians.height == 0:
            print("Warning: Aggregation resulted in 0 rows.", file=sys.stderr)

        medians.write_csv(args.output)
        
    except Exception as e:
        with open(args.output, "w") as f:
            f.write("gen,lineage,chrom,med_pbs_target,med_pbs_ref,med_T_bm,med_T_tm\n")
        sys.exit(1) # Force exit code 1 so Snakemake marks it as failed and shows the log

if __name__ == "__main__":
    main()