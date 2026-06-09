# scripts/aggregate_permutations.py
import polars as pl
import argparse
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()


    lf = pl.concat([
        pl.scan_csv(f).with_columns(pl.lit(i).alias("perm_id")) 
        for i, f in enumerate(args.inputs)
    ])


    consensus = lf.group_by(["chrom", "start", "end"]).agg([
        pl.col("slope_lmm").median().alias("slope_median"),
        pl.col("z_score").median().alias("z_score_median"),
        pl.col("z_score").std().alias("z_score_sd"),
        pl.col("p_value").median().alias("p_value_median"),
        pl.count().alias("n_perms")
    ])


    df = consensus.collect().sort(["chrom", "start"])
    df.write_csv(args.output)

if __name__ == "__main__":
    main()