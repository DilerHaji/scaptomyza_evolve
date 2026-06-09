#!/usr/bin/env python3
import argparse
from typing import List, Optional, Tuple, Dict
import sys
import polars as pl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

def find_col(cols: List[str], keys: List[str], default: Optional[str] = None) -> str:
    low = [c.lower() for c in cols]
    for k in keys:
        k = k.lower()
        for c, lc in zip(cols, low):
            if k in lc: return c
    if default is not None and default in cols: return default
    return None

def parse_chrom_map(map_file: str) -> Dict[str, str]:
    mapping = {}
    try:
        with open(map_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 2: continue
                rep_name = parts[0].strip()
                csv_name = parts[1].strip()
                mapping[csv_name] = rep_name
    except Exception as e:
        sys.exit(f"Error reading chrom map: {e}")
    return mapping

def parse_repeat_masker(repeat_file: str) -> Dict[str, List[Tuple[int, int]]]:
    repeats = {}
    try:
        with open(repeat_file, 'r') as f:
            for line in f:
                if not line.strip() or "position" in line or "score" in line: continue
                parts = line.split()
                if len(parts) < 7: continue
                chrom = parts[4]
                try:
                    start = int(parts[5])
                    end = int(parts[6])
                    if chrom not in repeats: repeats[chrom] = []
                    repeats[chrom].append((start, end))
                except ValueError: continue
    except Exception as e:
        sys.exit(f"Error reading repeat file: {e}")
    return repeats

def filter_by_window_content(df: pl.DataFrame, repeat_file: str, 
                             chrom_col: str, start_col: str, end_col: str,
                             max_frac: float, neighbor_window: int,
                             chrom_map_file: Optional[str] = None) -> pl.DataFrame:
    
    rep_dict = parse_repeat_masker(repeat_file)
    if not rep_dict: return df

    chrom_map = {}
    if chrom_map_file:
        chrom_map = parse_chrom_map(chrom_map_file)

    unique_chroms = df[chrom_col].unique().to_list()
    all_fractions = np.zeros(len(df), dtype=np.float32)
    
    mapped_count = 0
    missing_count = 0

    for chrom in unique_chroms:
        target_rep_chrom = chrom
        if chrom_map and chrom in chrom_map:
            target_rep_chrom = chrom_map[chrom]
        
        if target_rep_chrom not in rep_dict:
            missing_count += 1
            continue
        
        mapped_count += 1
        
        chrom_df = df.with_row_index("orig_idx").filter(pl.col(chrom_col) == chrom)
        if chrom_df.is_empty(): continue
        
        starts = chrom_df[start_col].to_numpy()
        ends = chrom_df[end_col].to_numpy()
        orig_indices = chrom_df["orig_idx"].to_numpy()
        
        chrom_repeats = rep_dict[target_rep_chrom]
        max_pos = max(ends.max(), max(end for _, end in chrom_repeats))
        
        genome_mask = np.zeros(max_pos + 1, dtype=np.int8)
        for s, e in chrom_repeats:
            s = max(0, s)
            e = min(e, max_pos)
            genome_mask[s:e] = 1 
        
        integral = np.zeros(len(genome_mask) + 1, dtype=np.int32)
        integral[1:] = np.cumsum(genome_mask)
        
        starts = np.clip(starts, 0, max_pos).astype(int)
        ends = np.clip(ends, 0, max_pos).astype(int)
        
        repeat_counts = integral[ends] - integral[starts]
        window_lens = ends - starts
        
        with np.errstate(divide='ignore', invalid='ignore'):
            fractions = repeat_counts / window_lens
            fractions = np.nan_to_num(fractions)
            
        all_fractions[orig_indices] = fractions

    df = df.with_columns(pl.Series(name="raw_rep_frac", values=all_fractions))

    target_col = "raw_rep_frac"
    if neighbor_window > 0:
        roll_size = (neighbor_window * 2) + 1
        df = df.sort([chrom_col, start_col])
        df = df.with_columns(
            pl.col("raw_rep_frac")
            .rolling_mean(window_size=roll_size, center=True, min_periods=1)
            .over(chrom_col)
            .alias("neighborhood_rep_frac")
        )
        target_col = "neighborhood_rep_frac"

    filtered_df = df.filter(pl.col(target_col) <= max_frac)
    
    n_removed = len(df) - len(filtered_df)
    return filtered_df

def per_contig_layout(df: pl.DataFrame, chrom: str, start: str, top_n: int, spacing: int):
    lens = (
        df.select([pl.col(chrom), pl.col(start).alias("pos")])
        .group_by(chrom)
        .agg([pl.col("pos").min().alias("minp"), pl.col("pos").max().alias("maxp")])
        .with_columns((pl.col("maxp") - pl.col("minp")).alias("len"))
        .sort("len", descending=True)
    )
    order = lens.select(chrom).head(top_n).to_series().to_list()
    lengths = {r[chrom]: int(r["len"]) for r in lens.iter_rows(named=True)}
    offsets, acc = {}, 0
    for i, c in enumerate(order):
        offsets[c] = acc
        acc += lengths[c] + spacing
    return order, lengths, offsets


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--repeats", required=False)
    ap.add_argument("--chrom-map", required=False, help="TSV mapping: RepeatName [tab] CSVName")
    
    ap.add_argument("--slope-col", default="slope_median")
    ap.add_argument("--pval-col", default="p_value_median")
    
    ap.add_argument("--error-col", required=False, help="Column name for standard error (e.g. slope_se)")
    ap.add_argument("--max-error", type=float, required=False, help="Discard windows with error > this value")

    ap.add_argument("--max-repeat", type=float, default=0.20)
    ap.add_argument("--neighbor-window", type=int, default=0)
    ap.add_argument("--fst-window-size", type=int, default=None)

    ap.add_argument("--output", required=True)
    ap.add_argument("--out-csv", required=False)

    ap.add_argument("--z-score", type=float, default=3.0)
    ap.add_argument("--extend-bp", type=int, default=0)
    ap.add_argument("--top-n", type=int, default=22)
    ap.add_argument("--ymax", type=float, default=None)
    
    ap.add_argument("--p-cutoff", type=float, default=0.01)
    ap.add_argument("--col-sig", default="#E69F00")
    ap.add_argument("--col-bg", default="#808080")
    ap.add_argument("--spacing", type=int, default=5_000_000)
    
    args = ap.parse_args()

    df = pl.read_csv(args.input)
    cols = df.columns
    
    chrom_col = find_col(cols, ["chrom", "chr"])
    start_col = find_col(cols, ["start", "pos", "begin"])

    if chrom_col is None or start_col is None:
        sys.exit("Could not identify Chromosome or Start/Pos columns.")

    if args.fst_window_size:
        end_col = "calculated_end"
        df = df.with_columns((pl.col(start_col) + args.fst_window_size).alias(end_col))
    else:
        end_col = find_col(cols, ["end", "stop"])
        if end_col is None:
            sys.exit("No 'end' column and no --fst-window-size provided.")

    if args.slope_col not in cols: sys.exit(f"Missing column {args.slope_col}")

    if args.error_col and args.max_error is not None:
        if args.error_col not in cols:
            sys.exit(f"Error column '{args.error_col}' not found in CSV.")

        before_count = len(df)
        df = df.filter(pl.col(args.error_col) <= args.max_error)
        removed = before_count - len(df)

        if df.is_empty(): sys.exit("All data removed after error filtering.")

    if args.repeats:
        df = filter_by_window_content(
            df, args.repeats, chrom_col, start_col, end_col,
            max_frac=args.max_repeat,
            neighbor_window=args.neighbor_window,
            chrom_map_file=args.chrom_map
        )
        if df.is_empty(): sys.exit("All data removed after filtering.")

    slopes = df[args.slope_col].to_numpy()
    mean_val = np.mean(slopes)
    std_val = np.std(slopes)
    threshold = mean_val + (args.z_score * std_val)

    if args.out_csv:
        outliers = df.filter(pl.col(args.slope_col) >= threshold)
        outliers = outliers.with_columns([
            (pl.col(start_col) - args.extend_bp).clip(0).alias("region_start"),
            (pl.col(end_col) + args.extend_bp).alias("region_end")
        ])
        outliers.write_csv(args.out_csv)

    order, lengths, offsets = per_contig_layout(df, chrom_col, start_col, args.top_n, args.spacing)
    
    mins = (
        df.select([pl.col(chrom_col), pl.col(start_col)])
        .filter(pl.col(chrom_col).is_in(order))
        .group_by(chrom_col).agg(pl.col(start_col).min().alias("__min__"))
    )
    off_map = pl.DataFrame({chrom_col: list(offsets.keys()), "__off__": list(offsets.values())})

    plot_df = (
        df.lazy()
        .filter(pl.col(chrom_col).is_in(order))
        .join(mins.lazy(), on=chrom_col, how="left")
        .join(off_map.lazy(), on=chrom_col, how="left")
        .with_columns([
            ((pl.col(start_col) - pl.col("__min__") + pl.col("__off__")) / 1e6).alias("x_mb"),
            pl.col(args.slope_col).alias("y_pos"),
            pl.col(args.pval_col).alias("pval")
        ])
        .collect()
    )

    if plot_df.is_empty(): sys.exit("No data for plotting.")

    x_vals = plot_df["x_mb"].to_numpy()
    y_vals = plot_df["y_pos"].to_numpy()
    p_vals = plot_df["pval"].to_numpy()
    mask_sig = p_vals < args.p_cutoff
    mask_bg = ~mask_sig
    
    plt.rcParams.update({"figure.dpi": 300})
    fig, ax = plt.subplots(figsize=(18, 6))

    for i in range(len(order) - 1):
        cname = order[i]
        mid = offsets[cname] + lengths[cname] + (args.spacing / 2)
        ax.axvline(mid / 1e6, color='#E0E0E0', linewidth=1.5, zorder=0)

    if np.any(mask_bg):
        ax.scatter(x_vals[mask_bg], y_vals[mask_bg], c=args.col_bg, s=10, alpha=0.3, lw=0, rasterized=True)
    if np.any(mask_sig):
        ax.scatter(x_vals[mask_sig], y_vals[mask_sig], c=args.col_sig, s=10, alpha=1.0, lw=0, rasterized=True)

    ax.axhline(0, color="black", lw=1)
    ax.axhline(mean_val, color="gray", ls="--", lw=1)
    ax.axhline(threshold, color="black", ls="--", lw=1.2)

    if args.ymax: ax.set_ylim(top=args.ymax)
    ax.set_ylabel("Slope Median", fontweight='bold')
    ax.set_xlabel("Genomic Position (Mb)")
    ax.set_xticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=args.col_sig, label=f'Sig (P < {args.p_cutoff})'),
        Line2D([0], [0], color='black', lw=1, linestyle='--', label=f'{args.z_score}σ Limit')
    ]
    ax.legend(handles=legend_elements, loc="upper left")

    fig.tight_layout()
    fig.savefig(args.output)

if __name__ == "__main__":
    main()