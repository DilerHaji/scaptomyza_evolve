#!/usr/bin/env python3
import argparse
from typing import List, Optional
import sys

import polars as pl
import matplotlib.pyplot as plt
import numpy as np

def find_col(cols: List[str], keys: List[str], default: Optional[str] = None) -> str:
    low = [c.lower() for c in cols]
    for k in keys:
        k = k.lower()
        for c, lc in zip(cols, low):
            if k in lc: return c
    if default is not None and default in cols: return default
    raise SystemExit(f"Could not find any of: {keys}")

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
    ap.add_argument("--output", required=True)
    ap.add_argument("--top-n", type=int, default=22, help="Number of chromosomes to plot")
    
    ap.add_argument("--p-cutoff", type=float, default=0.01, 
                    help="P-value threshold for coloring points Gold (default: 0.01)")
    
    ap.add_argument("--ymax", type=float, default=None, help="Force Y-axis max")
    ap.add_argument("--spacing", type=int, default=5_000_000)
    
    ap.add_argument("--col-sig", default="#E69F00", help="Significant (Gold)")
    ap.add_argument("--col-bg", default="#808080", help="Non-significant (Grey)")
    
    args = ap.parse_args()

    try:
        df = pl.read_csv(args.input)
    except Exception as e:
        sys.exit(f"Error reading CSV: {e}")

    cols = df.columns
    
    chrom_col = find_col(cols, ["chrom", "chr"])
    start_col = find_col(cols, ["start", "pos"])
    
    slope_col = "slope_median"
    pval_col = "p_value_median"

    if slope_col not in cols or pval_col not in cols:
        sys.exit(f"Input must contain '{slope_col}' and '{pval_col}'")

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
            pl.col(slope_col).alias("y_pos"),
            pl.col(pval_col).alias("pval")
        ])
        .collect()
    )

    if plot_df.is_empty(): sys.exit("No data after filtering.")

    x_vals = plot_df["x_mb"].to_numpy()
    y_vals = plot_df["y_pos"].to_numpy()
    p_vals = plot_df["pval"].to_numpy()

    mask_sig = p_vals < args.p_cutoff
    mask_bg = ~mask_sig
    
    plt.rcParams.update({"figure.dpi": 300})
    fig, ax = plt.subplots(figsize=(18, 6))

    point_size = 10 

    for i in range(len(order) - 1):
        chrom_name = order[i]
        chrom_end = offsets[chrom_name] + lengths[chrom_name]
        mid_gap = chrom_end + (args.spacing / 2)
        ax.axvline(mid_gap / 1e6, color='#E0E0E0', linewidth=1.5, zorder=0)

    if np.any(mask_bg):
        ax.scatter(x_vals[mask_bg], y_vals[mask_bg], c=args.col_bg, s=point_size, 
                   alpha=0.3, lw=0, rasterized=True, zorder=1, label="ns")

    if np.any(mask_sig):
        ax.scatter(x_vals[mask_sig], y_vals[mask_sig], c=args.col_sig, s=point_size, 
                   alpha=1.0, lw=0, rasterized=True, zorder=2, label=f"P < {args.p_cutoff}")

    ax.axhline(0, color="black", lw=1)
    
    mean_val = np.mean(y_vals)
    std_val = np.std(y_vals)
    z3 = mean_val + (3 * std_val)
    
    ax.axhline(mean_val, color="gray", ls="--", lw=1, label="Mean Slope")
    ax.axhline(z3, color="black", ls="--", lw=1.2, label="3σ Outlier")

    if args.ymax:
        ax.set_ylim(bottom=-0.0001, top=args.ymax)

    ax.set_ylabel("Slope Median", fontsize=11, fontweight='bold')
    ax.set_xlabel("Genomic Position (Mb)", fontsize=10)
    
    ax.set_xticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=args.col_sig, label=f'Sig (P < {args.p_cutoff})', markersize=10),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=args.col_bg, label='Non-sig', markersize=10),
        Line2D([0], [0], color='black', lw=1, linestyle='--', label='3σ Slope')
    ]
    ax.legend(handles=legend_elements, loc="upper left", frameon=True)

    fig.tight_layout()
    fig.savefig(args.output)

if __name__ == "__main__":
    main()