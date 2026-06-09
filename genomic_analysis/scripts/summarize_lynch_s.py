import polars as pl
import argparse
import os
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mean_s_files", nargs="+", required=True)
    parser.add_argument("--labels", nargs="+", required=True)
    parser.add_argument("--window_sizes", nargs="+", type=int, default=[200000, 500000])
    parser.add_argument("--candidate_regions", help="TSV with chrom, start, end columns")
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--max_abs_s", type=float, default=1.0,
                        help="Filter SNPs with |mean_s| above this value (default 1.0)")
    return parser.parse_args()


def load_and_merge(files, labels):
    dfs = []
    for f, label in zip(files, labels):
        if os.path.getsize(f) == 0:
            continue
        df = pl.read_csv(f).rename({"mean_s": f"s_{label}"})
        dfs.append(df)

    if not dfs:
        raise ValueError("No valid input files")

    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.join(df, on=["CHROM", "POS"], how="inner")

    return merged


def window_summary(df, s_cols, window_size):
    windowed = df.with_columns(
        (pl.col("POS") // window_size * window_size).alias("win_start")
    )

    agg_exprs = []
    for col in s_cols:
        agg_exprs.extend([
            pl.col(col).abs().median().alias(f"median_abs_{col}"),
            pl.col(col).abs().mean().alias(f"mean_abs_{col}"),
            pl.col(col).median().alias(f"median_{col}"),
            pl.col(col).count().alias(f"n_{col}"),
        ])

    result = windowed.group_by(["CHROM", "win_start"]).agg(agg_exprs).sort(["CHROM", "win_start"])

    median_abs_cols = [f"median_abs_{c}" for c in s_cols]
    result = result.with_columns(
        pl.mean_horizontal(median_abs_cols).alias("grand_median_abs_s")
    )

    return result


def flag_candidate_regions(windowed_df, candidate_file, window_size):
    if candidate_file is None or not os.path.exists(candidate_file):
        return windowed_df.with_columns(pl.lit(False).alias("is_candidate"))

    cands = pl.read_csv(candidate_file, separator="\t")

    chrom_col = [c for c in cands.columns if "chrom" in c.lower()][0]
    start_col = [c for c in cands.columns if "start" in c.lower()][0]
    end_col = [c for c in cands.columns if "end" in c.lower()][0]

    candidate_windows = set()
    for row in cands.iter_rows(named=True):
        chrom = row[chrom_col]
        region_start = row[start_col]
        region_end = row[end_col]
        win_s = (region_start // window_size) * window_size
        win_e = (region_end // window_size) * window_size
        for w in range(win_s, win_e + window_size, window_size):
            candidate_windows.add((chrom, w))

    flags = [
        (row["CHROM"], row["win_start"]) in candidate_windows
        for row in windowed_df.iter_rows(named=True)
    ]
    return windowed_df.with_columns(pl.Series("is_candidate", flags))


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.output_prefix) or ".", exist_ok=True)

    merged = load_and_merge(args.mean_s_files, args.labels)
    s_cols = [c for c in merged.columns if c.startswith("s_")]

    for col in s_cols:
        merged = merged.with_columns(
            pl.when(pl.col(col).is_finite() & (pl.col(col).abs() <= args.max_abs_s))
            .then(pl.col(col))
            .otherwise(None)
            .alias(col)
        )

    for col in s_cols:
        n_valid = merged[col].drop_nulls().len()

    abs_vals = []
    for col in s_cols:
        vals = merged[col].drop_nulls().to_numpy()
        abs_vals.extend(np.abs(vals).tolist())

    abs_arr = np.array(abs_vals)

    clipped = abs_arr[abs_arr <= 1.0]

    summary_lines = ["metric,value"]
    summary_lines.append(f"n_snps_finite,{len(abs_arr)}")
    summary_lines.append(f"n_snps_below_1,{len(clipped)}")
    summary_lines.append(f"pct_below_1,{len(clipped)/len(abs_arr)*100:.1f}")
    for q in [0.25, 0.50, 0.75, 0.90, 0.95, 0.99]:
        summary_lines.append(f"p{int(q*100)}_abs_s,{np.quantile(abs_arr, q):.6f}")

    with open(f"{args.output_prefix}_genomewide_s.csv", "w") as f:
        f.write("\n".join(summary_lines) + "\n")

    for wsize in args.window_sizes:
        windowed = window_summary(merged, s_cols, wsize)
        windowed = flag_candidate_regions(windowed, args.candidate_regions, wsize)

        outfile = f"{args.output_prefix}_windows_{wsize}.csv"
        windowed.write_csv(outfile)

        grand = windowed["grand_median_abs_s"].drop_nulls().to_numpy()
        grand = grand[np.isfinite(grand)]

        if "is_candidate" in windowed.columns:
            cand = windowed.filter(pl.col("is_candidate"))
            bg = windowed.filter(~pl.col("is_candidate"))
            if len(cand) > 0:
                cand_s = cand["grand_median_abs_s"].drop_nulls().to_numpy()
                bg_s = bg["grand_median_abs_s"].drop_nulls().to_numpy()
                cand_s = cand_s[np.isfinite(cand_s)]
                bg_s = bg_s[np.isfinite(bg_s)]

                # Levene dir_thr implied
                s_cand = np.median(cand_s)
                s_bg = np.median(bg_s)


if __name__ == "__main__":
    main()
