import polars as pl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import matplotlib.colors as mcolors
import argparse
import sys
import os
import numpy as np
import copy
import re

PRETTY_NAMES = {
    "btwTB": "Triculture vs Biculture",
    "btwTM": "Triculture vs Monoculture",
    "btwBM": "Biculture vs Monoculture",
    "btwMV": "Monoculture vs Variant"
}

def get_pretty_name(group_raw):
    for k, v in PRETTY_NAMES.items():
        if k in group_raw:
            return v
    return group_raw

def parse_metadata(filepath):
    folder_name = os.path.basename(os.path.dirname(filepath))
    parts = folder_name.split('_')
    if len(parts) >= 3 and parts[-1].isdigit() and parts[-2].isdigit():
        window_size = int(parts[-2])
        group_name = "_".join(parts[:-2]) 
        return group_name, window_size
    return "Unknown", 0

def get_mako_cmap():
    colors = ["#0B0405", "#28192F", "#3E3253", "#4A5076", "#467091", "#3E91A8", "#35B4AB", "#59D7A4", "#9BEBA7", "#DEF5E5"]
    return mcolors.LinearSegmentedColormap.from_list("mako_custom", colors)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs='+', required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--common-only", action="store_true")
    parser.add_argument("--bin-size", type=int, default=5000)
    parser.add_argument("--column", type=str, default="z_score_median")
    parser.add_argument("--top-n", type=int, default=50)
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--normalize", action="store_true")
    group.add_argument("--percentile", action="store_true")

    parser.add_argument("--cmap", type=str, default="mako") 
    parser.add_argument("--vmin", type=float, default=None)
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument("--gap", type=int, default=25)
    parser.add_argument("--width", type=float, default=20.0)
    parser.add_argument("--panel-height", type=float, default=2.5)
    parser.add_argument("--font-scale", type=float, default=1.0, help="Scaling factor for all text elements")
    
    args = parser.parse_args()

    lazy_dfs = []
    windows_per_group = {} 
    
    for f in args.inputs:
        group, w_size = parse_metadata(f)
        if group not in windows_per_group: windows_per_group[group] = set()
        windows_per_group[group].add(w_size)
        
        ldf = pl.scan_csv(f).select([
            pl.col("chrom"), pl.col("end"), pl.col(args.column).alias("val")
        ])
        ldf = ldf.with_columns([
            pl.lit(group).alias("group"),
            pl.lit(w_size).alias("window_size"),
            (pl.col("end") // args.bin_size).alias("pos_bin")
        ])
        lazy_dfs.append(ldf)

    if not lazy_dfs: sys.exit("No inputs.")
    all_data = pl.concat(lazy_dfs).group_by(["group", "window_size", "chrom", "pos_bin"]).agg(pl.mean("val")).collect()

    if args.percentile:
        all_data = all_data.with_columns(
            (pl.col("val").rank() / pl.col("val").count()).over(["group", "window_size"]).alias("val_rank")
        )
        all_data = all_data.drop("val").rename({"val_rank": "val"})
    elif args.normalize:
        stats = all_data.group_by(["group", "window_size"]).agg([
            pl.col("val").median().alias("med"),
            (pl.col("val") - pl.col("val").median()).abs().median().alias("mad")
        ])
        all_data = all_data.join(stats, on=["group", "window_size"]).with_columns(
            ((pl.col("val") - pl.col("med")) / (pl.col("mad") * 1.4826)).alias("val_norm")
        )
        all_data = all_data.drop("val").rename({"val_norm": "val"})

    all_groups = sorted(list(windows_per_group.keys()))
    if args.common_only:
        common = windows_per_group[all_groups[0]]
        for g in all_groups[1:]: common = common.intersection(windows_per_group[g])
        sorted_windows = sorted(list(common))
        all_data = all_data.filter(pl.col("window_size").is_in(sorted_windows))
    else:
        all_wins = set()
        for s in windows_per_group.values(): all_wins.update(s)
        sorted_windows = sorted(list(all_wins))

    w_map = {size: i for i, size in enumerate(sorted_windows)}

    scaf_stats = all_data.group_by("chrom").agg(pl.max("pos_bin").alias("max_bin")).sort("max_bin", descending=True).head(args.top_n)
    sorted_chroms = scaf_stats["chrom"].to_list()
    chrom_sizes = {row["chrom"]: row["max_bin"] for row in scaf_stats.iter_rows(named=True)}

    FS_AXIS = 10 * args.font_scale
    FS_LABEL = 12 * args.font_scale
    FS_TICK = 9 * args.font_scale
    FS_FOOTER_BIG = 10 * args.font_scale
    FS_FOOTER_SMALL = 8 * args.font_scale

    num_panels = len(all_groups)
    height_ratios = [10] * num_panels + [1.5] 
    total_h = (args.panel_height * num_panels) + 1.2
    
    fig = plt.figure(figsize=(args.width, total_h), dpi=300)
    gs = GridSpec(num_panels + 1, 2, width_ratios=[50, 1], height_ratios=height_ratios, wspace=0.02, hspace=0.08)

    if args.cmap == "mako":
        base_cmap = get_mako_cmap()
    else:
        try: base_cmap = copy.copy(plt.get_cmap(args.cmap))
        except: base_cmap = get_mako_cmap()
    
    cmap = copy.copy(base_cmap)
    
    cmap.set_bad(color='white')

    axes = []
    x_layout = [] 
    
    if args.vmin is None: args.vmin = all_data["val"].min()
    if args.vmax is None: args.vmax = all_data["val"].max()

    for i, group in enumerate(all_groups):
        ax = fig.add_subplot(gs[i, 0])
        axes.append(ax)
        
        group_data = all_data.filter(pl.col("group") == group)
        matrix_parts = []
        spacer = np.full((len(sorted_windows), args.gap), np.nan) # NaN makes it White via set_bad
        curr_x = 0

        for j, chrom in enumerate(sorted_chroms):
            c_data = group_data.filter(pl.col("chrom") == chrom)
            max_bin = chrom_sizes[chrom]
            matrix = np.full((len(sorted_windows), max_bin + 1), args.vmin)
            
            if not c_data.is_empty():
                rows = [w_map[w] for w in c_data["window_size"].to_list()]
                cols = c_data["pos_bin"].to_list()
                vals = c_data["val"].to_list()
                valid = [c <= max_bin for c in cols]
                if all(valid): matrix[rows, cols] = vals
                else: matrix[np.array(rows)[valid], np.array(cols)[valid]] = np.array(vals)[valid]

            matrix_parts.append(matrix)
            width = matrix.shape[1]

            if i == 0:
                size_mb = (width * args.bin_size) / 1e6
                parts = chrom.split("_")
                # Format: "Scaffold 7"
                clean = f"Scaffold {parts[2]}" if (len(parts)>2 and parts[1].startswith("Sc")) else chrom
                
                x_layout.append({
                    "label": clean,
                    "mb": size_mb,
                    "start": curr_x,
                    "width": width
                })
            
            curr_x += width + args.gap
            if j < len(sorted_chroms) - 1:
                matrix_parts.append(spacer)

        full_matrix = np.hstack(matrix_parts)
        total_width_bins = full_matrix.shape[1]

        im = ax.imshow(full_matrix, aspect='auto', origin='lower', cmap=cmap, 
                       vmin=args.vmin, vmax=args.vmax, interpolation='nearest')
        
        ax.set_yticks(range(0, len(sorted_windows), 2))
        yticklabels = [f"{w//1000}k" for w in sorted_windows[::2]]
        ax.set_yticklabels(yticklabels, fontsize=FS_TICK)
        ax.set_ylabel(get_pretty_name(group), fontsize=FS_LABEL, fontweight='bold', labelpad=15)
        
        ax.set_xticks([])
        for spine in ax.spines.values(): spine.set_visible(False)

    ax_foot = fig.add_subplot(gs[num_panels, 0])
    ax_foot.set_xlim(0, total_width_bins)
    ax_foot.set_ylim(0, 1)
    ax_foot.axis('off')

    for item in x_layout:
        rect = mpatches.Rectangle((item["start"], 0.7), item["width"], 0.3, 
                                  facecolor='#555555', edgecolor='white', linewidth=0.5)
        ax_foot.add_patch(rect)
        
        center_x = item["start"] + (item["width"] / 2)
        if item["width"] > (total_width_bins * 0.015): 
            ax_foot.text(center_x, 0.45, item["label"], 
                         ha='center', va='center', color='black', fontsize=FS_FOOTER_BIG, fontweight='bold')
            ax_foot.text(center_x, 0.1, f"{item['mb']:.1f} Mb", 
                         ha='center', va='center', color='#444444', fontsize=FS_FOOTER_SMALL)

    ax_cbar = fig.add_subplot(gs[0:num_panels, 1])
    cbar = plt.colorbar(im, cax=ax_cbar, orientation='vertical')
    
    if args.percentile:
        lbl = "Percentile Rank"
    elif args.normalize:
        lbl = "Relative Z (Sigma)"
    else:
        lbl = args.column
    cbar.set_label(lbl, fontsize=FS_LABEL, fontweight='bold')
    cbar.ax.tick_params(labelsize=FS_TICK)

    plt.savefig(args.output, bbox_inches='tight')

if __name__ == "__main__":
    main()