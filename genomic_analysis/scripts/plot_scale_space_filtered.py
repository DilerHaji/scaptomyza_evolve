import polars as pl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import sys
import re
import numpy as np
import copy

def parse_window_size(filename):
    match = re.search(r'_(\d+)_(\d+)\.csv$', filename)
    if match: return int(match.group(1))
    return 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs='+', required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--bin-size", type=int, default=5000)
    parser.add_argument("--column", type=str, default="z_score_median")
    parser.add_argument("--top-n", type=int, default=50)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--normalize", action="store_true", help="Z-score normalization")
    group.add_argument("--percentile", action="store_true", help="Convert to Percentile Rank (0-1)")
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument("--vmin", type=float, default=None)
    parser.add_argument("--gap", type=int, default=15)
    parser.add_argument("--width", type=float, default=20.0)
    parser.add_argument("--height", type=float, default=4.0)
    args = parser.parse_args()

    lazy_dfs = []
    window_sizes_found = set()
    for f in args.inputs:
        w_size = parse_window_size(f)
        window_sizes_found.add(w_size)
        ldf = pl.scan_csv(f).select([
            pl.col("chrom"), pl.col("end"), pl.col(args.column).alias("val")
        ]).with_columns([
            pl.lit(w_size).alias("window_size"),
            (pl.col("end") // args.bin_size).alias("pos_bin")
        ])
        lazy_dfs.append(ldf)

    if not lazy_dfs: sys.exit("No inputs.")
    
    df = pl.concat(lazy_dfs).group_by(["chrom", "window_size", "pos_bin"]).agg(pl.mean("val")).collect()

    if args.normalize:
        stats = df.group_by("window_size").agg([
            pl.col("val").median().alias("med"),
            (pl.col("val") - pl.col("val").median()).abs().median().alias("mad")
        ])
        df = df.join(stats, on="window_size").with_columns(
            ((pl.col("val") - pl.col("med")) / (pl.col("mad") * 1.4826)).alias("val_norm")
        )
        df = df.drop("val").rename({"val_norm": "val"})
        
    elif args.percentile:
        df = df.with_columns(
            (pl.col("val").rank() / pl.col("val").count()).over("window_size").alias("val_rank")
        )
        df = df.drop("val").rename({"val_rank": "val"})

    top_scaffolds = df.group_by("chrom").agg(pl.max("pos_bin").alias("max_bin")).sort("max_bin", descending=True).head(args.top_n)
    sorted_chroms = top_scaffolds["chrom"].to_list()
    sorted_windows = sorted(list(window_sizes_found))
    w_map = {size: i for i, size in enumerate(sorted_windows)}
    
    matrix_parts = [] 
    chrom_ticks = []
    chrom_labels = []
    current_x_offset = 0
    spacer_col = np.full((len(sorted_windows), args.gap), np.nan)

    for i, chrom in enumerate(sorted_chroms):
        c_data = df.filter(pl.col("chrom") == chrom)
        max_bin = c_data["pos_bin"].max()
        matrix = np.full((len(sorted_windows), max_bin + 1), np.nan)
        
        rows = [w_map[w] for w in c_data["window_size"].to_list()]
        cols = c_data["pos_bin"].to_list()
        vals = c_data["val"].to_list()
        
        valid_mask = [c <= max_bin for c in cols]
        if all(valid_mask): matrix[rows, cols] = vals
        else:
            r, c_idx, v = np.array(rows)[valid_mask], np.array(cols)[valid_mask], np.array(vals)[valid_mask]
            matrix[r, c_idx] = v

        matrix_parts.append(matrix)
        width = matrix.shape[1]
        
        if width > 5: 
            chrom_ticks.append(current_x_offset + (width / 2))
            parts = chrom.split("_")
            label = "Sc" + parts[2] if (len(parts) > 2 and parts[1].startswith("Sc")) else chrom[:10]
            chrom_labels.append(label)
        
        current_x_offset += width
        if i < len(sorted_chroms) - 1:
            matrix_parts.append(spacer_col)
            current_x_offset += args.gap

    full_matrix = np.hstack(matrix_parts)

    fig, ax = plt.subplots(figsize=(args.width, args.height), dpi=300)
    cmap = copy.copy(plt.cm.viridis)
    cmap.set_bad(color='#333333')
    
    im = ax.imshow(full_matrix, aspect='auto', origin='lower', cmap=cmap, 
                   vmin=args.vmin, vmax=args.vmax, interpolation='nearest')
    
    ax.set_xticks(chrom_ticks)
    ax.set_xticklabels(chrom_labels, rotation=90, fontsize=8)
    ax.set_yticks(range(len(sorted_windows)))
    ax.set_yticklabels(sorted_windows, fontsize=7)
    ax.set_ylabel("Window Size (bp)", fontsize=9)
    for s in ax.spines.values(): s.set_visible(False)
    
    cbar = plt.colorbar(im, ax=ax, fraction=0.015, pad=0.01)
    
    if args.percentile:
        lbl = "Percentile Rank"
    elif args.normalize:
        lbl = "Relative Z (Sigma)"
    else:
        lbl = args.column
    cbar.set_label(lbl, fontsize=8)

    plt.tight_layout()
    plt.savefig(args.output)

if __name__ == "__main__":
    main()