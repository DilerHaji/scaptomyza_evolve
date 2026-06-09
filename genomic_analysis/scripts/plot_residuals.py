#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import List, Optional, Tuple, Dict
import sys

import polars as pl
import matplotlib
# Force headless backend for clusters
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


def find_col(cols: List[str], keys: List[str], default: Optional[str] = None) -> str:
    low = [c.lower() for c in cols]
    for k in keys:
        k = k.lower()
        for c, lc in zip(cols, low):
            if k in lc:
                return c
    if default is not None and default in cols:
        return default
    return None 

def per_contig_layout(
    df: pl.DataFrame, chrom: str, start: str, top_n: int, spacing: int
) -> Tuple[List[str], Dict[str, int], Dict[str, int]]:
    lens = (
        df.select([pl.col(chrom), pl.col(start).alias("pos")])
        .group_by(chrom)
        .agg([pl.col("pos").min().alias("minp"), pl.col("pos").max().alias("maxp")])
        .with_columns((pl.col("maxp") - pl.col("minp")).alias("len"))
        .sort("len", descending=True)
    )

    order = lens.select(chrom).head(top_n).to_series().to_list()
    if not order:
        raise SystemExit("No chromosomes found to plot.")

    lengths = {r[chrom]: int(r["len"]) for r in lens.iter_rows(named=True)}
    offsets, acc = {}, 0
    for c in order:
        offsets[c] = acc
        acc += lengths[c] + spacing
    return order, lengths, offsets

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--input", required=True, help="CSV file (fst_slopes_corrected.csv)")
    ap.add_argument("--output", required=True, help="Output image path (e.g. manhattan.png)")
    
    ap.add_argument("--slope-col", default="slope_residual", help="Column to plot on Y-axis")
    ap.add_argument("--r2-col", default="r2", help="Column for shading")
    
    ap.add_argument("--width", type=float, default=20, help="Figure width (inches)")
    ap.add_argument("--height", type=float, default=6, help="Figure height (inches)")
    ap.add_argument("--point-size", type=float, default=4, help="Scatter point size")
    ap.add_argument("--alpha", type=float, default=0.8, help="Point opacity (0.0 - 1.0)")
    ap.add_argument("--title", default=None, help="Custom plot title")

    ap.add_argument("--min-r2", type=float, default=0.1, help="Hide points with very low R2")
    ap.add_argument("--min-snps", type=float, default=0, help="Filter out windows with mean density < N")
    ap.add_argument("--top-n", type=int, default=20, help="Number of chromosomes to plot")
    
    ap.add_argument("--ymin", type=float, default=None, help="Force Y-axis min")
    ap.add_argument("--ymax", type=float, default=None, help="Force Y-axis max")
    
    args = ap.parse_args()

    path = Path(args.input)
    if not path.exists():
        sys.exit(f"Input file not found: {path}")

    try:
        head = pl.read_csv(path, n_rows=1)
        cols = head.columns
        chrom_col = find_col(cols, ["chrom", "chr", "chromosome"]) or sys.exit("No chrom col found")
        start_col = find_col(cols, ["start", "pos", "position"]) or sys.exit("No start col found")

        slope_col = find_col(cols, [args.slope_col, "slope_raw", "slope_wls"]) or sys.exit(f"Slope col '{args.slope_col}' not found")
        r2_col = find_col(cols, [args.r2_col, "r2_wls", "r_squared"]) or sys.exit(f"R2 col '{args.r2_col}' not found")

        density_col = find_col(cols, ["mean_density", "snps", "density"])

        if args.min_snps > 0 and not density_col:
            sys.exit(f"Error: --min-snps {args.min_snps} requested, but no density column found in CSV.")

    except Exception as e:
        sys.exit(f"Error reading header: {e}")

    coord_df = pl.read_csv(path, columns=[chrom_col, start_col])
    order, lengths, offsets = per_contig_layout(coord_df, chrom_col, start_col, args.top_n, 5_000_000)
    
    cols_to_select = [chrom_col, start_col, slope_col, r2_col]
    if density_col: cols_to_select.append(density_col)

    q = pl.scan_csv(path).select(cols_to_select)
    q = q.filter(pl.col(chrom_col).is_in(order))
    
    if args.min_r2 > 0:
        q = q.filter(pl.col(r2_col) >= args.min_r2)

    if args.min_snps > 0 and density_col:
        q = q.filter(pl.col(density_col) >= args.min_snps)
        
    q = q.filter(
        pl.col(slope_col).is_not_null() & 
        pl.col(slope_col).is_finite()
    )

    mins = (
        coord_df.select([pl.col(chrom_col), pl.col(start_col)])
        .filter(pl.col(chrom_col).is_in(order))
        .group_by(chrom_col)
        .agg(pl.col(start_col).min().alias("__min__"))
    )
    
    off_map = pl.DataFrame({chrom_col: list(offsets.keys()), "__off__": list(offsets.values())})

    plot_df = (
        q.join(mins.lazy(), on=chrom_col, how="left")
        .join(off_map.lazy(), on=chrom_col, how="left")
        .with_columns([
            ((pl.col(start_col) - pl.col("__min__") + pl.col("__off__")) / 1e6).alias("x_mb")
        ])
        .collect()
    )

    if plot_df.height == 0:
        sys.exit("Error: No data points remained after filtering.")

    y_vals = plot_df[slope_col].to_numpy()
    mean_slope = np.mean(y_vals)
    std_slope = np.std(y_vals)
    z3 = mean_slope + (3 * std_slope)
    z_min3 = mean_slope - (3 * std_slope)

    bounds = [0.0, 0.50, 0.75, 0.90, 0.95, 1.0]
    hex_colors = ["#D3D3D3", "#A9A9A9", "#808080", "#404040", "#000000"]
    cmap = mcolors.ListedColormap(hex_colors)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
        
    plt.rcParams.update({"figure.dpi": 300})
    
    fig, ax = plt.subplots(figsize=(args.width, args.height))

    x = plot_df["x_mb"].to_numpy()
    c_vals = plot_df[r2_col].to_numpy()
    
    sc = ax.scatter(x, y_vals, c=c_vals, s=args.point_size, cmap=cmap, norm=norm, 
                    alpha=args.alpha, linewidths=0, rasterized=True)

    ax.axhline(0, color="gray", linestyle="-", linewidth=0.5, alpha=0.5)
    ax.axhline(mean_slope, color="cornflowerblue", linestyle="--", linewidth=1, label="Mean")
    ax.axhline(z3, color="red", linestyle="--", linewidth=1, label=f"+3σ ({z3:.1e})")
    
    if z_min3 > np.min(y_vals):
        ax.axhline(z_min3, color="red", linestyle="--", linewidth=1, label=f"-3σ ({z_min3:.1e})")

    sorted_offsets = sorted(offsets.items(), key=lambda i: i[1])
    for chrom_name, off_val in sorted_offsets:
        ax.axvline(off_val / 1e6, color='black', linewidth=0.5, alpha=0.1)
        mid_point = (off_val + (lengths[chrom_name]/2)) / 1e6
        clean_name = chrom_name.replace("chr_", "").replace("_HRSCAF", "")
        ax.text(mid_point, ax.get_ylim()[0], clean_name, 
                ha='center', va='bottom', fontsize=6, rotation=90, alpha=0.5)

    cbar = plt.colorbar(sc, ax=ax, pad=0.01, ticks=bounds)
    cbar.set_label(f"Reliability ($R^2$)", rotation=270, labelpad=15)

    if args.ymin is not None: ax.set_ylim(bottom=args.ymin)
    if args.ymax is not None: ax.set_ylim(top=args.ymax)
    
    ax.set_xlabel("Genomic Position (Mb)")
    ax.set_ylabel(f"Trend ({slope_col})")
    
    title_txt = args.title if args.title else f"Genome-Wide Divergence Trends (Top {args.top_n} Contigs)"
    if args.min_snps > 0:
        title_txt += f"\n(Filtered: SNP Density >= {args.min_snps})"

    ax.set_title(title_txt, loc="left", fontweight="bold")
    ax.legend(loc="upper left", frameon=True)

    fig.tight_layout()
    fig.savefig(args.output)

if __name__ == "__main__":
    main()