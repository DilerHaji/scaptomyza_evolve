#!/usr/bin/env python3
import polars as pl
import argparse
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input bin CSV")
    parser.add_argument("--medians", required=True, help="Medians CSV")
    parser.add_argument("--output", required=True, help="Output PBE CSV")
    args = parser.parse_args()

    bin_schema = {
        "chrom": pl.String,
        "gen": pl.Int64,
        "lineage": pl.String,
        "total_divergence": pl.Float64,
        "T_tm": pl.Float64,
        "T_bm": pl.Float64,
        "pbs_target": pl.Float64,
        "pbs_ref": pl.Float64,
        "pbsn1_target": pl.Float64,
        "pbsn1_ref": pl.Float64
    }

    median_schema = {
        "chrom": pl.String,
        "gen": pl.Int64,
        "lineage": pl.String,
        "med_pbs_target": pl.Float64,
        "med_pbs_ref": pl.Float64,
        "med_T_bm": pl.Float64,
        "med_T_tm": pl.Float64
    }

    try:
        lf = pl.scan_csv(args.input, schema_overrides=bin_schema)
        med_lf = pl.scan_csv(args.medians, schema_overrides=median_schema)

        lf = lf.join(med_lf, on=["gen", "lineage", "chrom"], how="left")
        
        lf = lf.with_columns([
            (pl.col("pbs_target") - (pl.col("T_bm") * (pl.col("med_pbs_target") / pl.col("med_T_bm")))).alias("pbe_target"),
            (pl.col("pbs_ref") - (pl.col("T_tm") * (pl.col("med_pbs_ref") / pl.col("med_T_tm")))).alias("pbe_ref")
        ])
        

        lf.collect().write_csv(args.output)
        
    except Exception as e:
        header = "chrom,start,end,gen,lineage,total_divergence,T_tm,T_bm,pbs_target,pbs_ref,pbsn1_target,pbsn1_ref,pbe_target,pbe_ref\n"
        with open(args.output, "w") as f:
            f.write(header)
        sys.exit(1) # Return error so Snakemake sees it

if __name__ == "__main__":
    main()