#!/usr/bin/env python3
import argparse
import polars as pl
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
import sys
import os
import copy


matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams['axes.titlesize'] = 8
matplotlib.rcParams['axes.labelsize'] = 7
matplotlib.rcParams['xtick.labelsize'] = 6
matplotlib.rcParams['ytick.labelsize'] = 6
matplotlib.rcParams['legend.fontsize'] = 6
matplotlib.rcParams['svg.fonttype'] = 'none' 
matplotlib.rcParams['pdf.fonttype'] = 42


def parse_window_size(filepath):
    folder_name = os.path.basename(os.path.dirname(filepath))
    parts = folder_name.split('_')
    ints = [int(p) for p in parts if p.isdigit()]
    if len(ints) >= 2: return ints[-2]
    elif len(ints) == 1: return ints[0]
    return 0

def load_bed_annotations(bed_file, target_chrom, start, end):
    if not bed_file: return [], []
    try:
        q = pl.scan_csv(bed_file, separator='\t', has_header=False, truncate_ragged_lines=True)
        df = q.filter(
            (pl.col("column_1") == target_chrom) & 
            (pl.col("column_3") >= start) & 
            (pl.col("column_2") <= end)
        ).collect()
        
        data = []
        labels = set()
        if df.height > 0:
            for row in df.iter_rows():
                c, s, e, lbl = row[0], row[1], row[2], row[3]
                data.append((s, e, lbl))
                labels.add(lbl)
        return data, sorted(list(labels))
    except Exception as e:
        return [], []

def process_heatmap_region(files, col, chrom, start, end, bin_size):
    file_map = [(parse_window_size(f), f) for f in files]
    file_map.sort(key=lambda x: x[0])
    sorted_sizes = [x[0] for x in file_map]
    sorted_files = [x[1] for x in file_map]
    
    length = end - start
    n_bins = int(length / bin_size) + 1
    matrix = np.full((len(sorted_sizes), n_bins), np.nan)
    
    for r_idx, fpath in enumerate(sorted_files):
        q = pl.scan_csv(fpath).filter(
            (pl.col("chrom") == chrom) &
            (pl.col("end") >= start) &
            (pl.col("start") <= end)
        )
        df = q.collect()
        if df.height == 0: continue
        
        df = df.with_columns([
            ((pl.col("start") + pl.col("end")) / 2.0).alias("midpoint")
        ])
        df = df.with_columns([
            (((pl.col("midpoint") - start) / bin_size).floor().cast(pl.Int64)).alias("bin_idx")
        ])
        
        agg = df.group_by("bin_idx").agg(pl.col(col).mean())
        b_indices = agg["bin_idx"].to_numpy()
        b_vals = agg[col].to_numpy()
        
        valid = (b_indices >= 0) & (b_indices < n_bins)
        if np.any(valid):
            matrix[r_idx, b_indices[valid]] = b_vals[valid]
            
    return matrix, sorted_sizes

def smooth_and_aggregate(pos, vals, ods, window_size, use_median=True):
    sort_idx = np.argsort(pos)
    s_pos = pos[sort_idx]
    s_vals = vals[sort_idx]
    s_ods = ods[sort_idx]
    
    if window_size <= 0:
        return s_pos, s_vals, s_ods

    y_smooth = np.zeros_like(s_vals)
    od_driver = np.zeros_like(s_vals) 
    
    radius = window_size // 2
    left_idxs = np.searchsorted(s_pos, s_pos - radius, side='left')
    right_idxs = np.searchsorted(s_pos, s_pos + radius, side='right')
    
    for i in range(len(s_pos)):
        l, r = left_idxs[i], right_idxs[i]
        chunk_v = s_vals[l:r]
        chunk_o = s_ods[l:r]
        
        if len(chunk_v) > 0:
            y_smooth[i] = np.median(chunk_v) if use_median else np.mean(chunk_v)
            od_driver[i] = np.max(chunk_o)
        else:
            y_smooth[i] = s_vals[i]
            od_driver[i] = s_ods[i]
            
    return s_pos, y_smooth, od_driver


def save_colorbar_legend(cmap, norm, label, filename, ticks=None, ticklabels=None):
    fig_leg = plt.figure(figsize=(1.5, 3), dpi=300)
    ax_leg = fig_leg.add_axes([0.3, 0.1, 0.2, 0.8])
    cb = matplotlib.colorbar.ColorbarBase(ax_leg, cmap=cmap, norm=norm, orientation='vertical')
    cb.set_label(label, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabels is not None:
        cb.set_ticklabels(ticklabels)
    plt.savefig(filename, bbox_inches='tight', transparent=True)
    plt.savefig(filename.replace('.png', '.svg'), bbox_inches='tight', transparent=True)
    plt.close(fig_leg)

def save_elements_legend(handles, title, filename):
    if not handles: return
    n = len(handles)
    fig_leg = plt.figure(figsize=(2, 0.3 * n + 0.5), dpi=300)
    fig_leg.legend(handles=handles, title=title, loc='center', frameon=False, fontsize=7, title_fontsize=8)
    plt.savefig(filename, bbox_inches='tight', transparent=True)
    plt.savefig(filename.replace('.png', '.svg'), bbox_inches='tight', transparent=True)
    plt.close(fig_leg)


def render_zoom_plot(region_id, chrom, start, end, 
                     lrt_df, 
                     heatmap_matrix, heatmap_sizes, 
                     annot_data, annot_labels, 
                     args, output_base):

    filename = f"{output_base}.png"
    
    if "singular_fit" in lrt_df.columns:
        lrt_df = lrt_df.filter(pl.col("singular_fit") == False)
    if "converged" in lrt_df.columns:
        lrt_df = lrt_df.filter(pl.col("converged") == True)
        
    if lrt_df.height == 0:
        return

    FIG_W = args.width
    FIG_H = args.height
    fs_tick, fs_label = 6, 8
    CAT_COLORS = ["#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00"] 

    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=300)
    
    ratios = args.ratios if args.ratios else [1.0, 2.0, 0.4, 0.2]
    gs = GridSpec(4, 1, height_ratios=ratios, 
                  hspace=0.1, left=0.12, right=0.95, top=0.95, bottom=0.05)
    
    ax_heat = fig.add_subplot(gs[0])
    ax_glm = fig.add_subplot(gs[1], sharex=ax_heat)
    ax_annot = fig.add_subplot(gs[2], sharex=ax_heat)
    ax_foot = fig.add_subplot(gs[3])

    cmap = copy.copy(plt.get_cmap("viridis"))
    cmap.set_bad(color='white', alpha=1.0)
    
    hm_valid = heatmap_matrix[~np.isnan(heatmap_matrix)]
    if len(hm_valid) > 0:
        vmin, vmax = np.percentile(hm_valid, 5), np.percentile(hm_valid, 99)
    else:
        vmin, vmax = 0, 1

    ax_heat.imshow(heatmap_matrix, aspect='auto', origin='lower', cmap=cmap,
                   extent=[start, end, 0, len(heatmap_sizes)], 
                   vmin=vmin, vmax=vmax, interpolation='nearest')

    n_rows = len(heatmap_sizes)
    if n_rows > 1:
        tick_indices = np.linspace(0, n_rows - 1, min(5, n_rows), dtype=int)
    else:
        tick_indices = [0]

    ax_heat.set_yticks([i + 0.5 for i in tick_indices])
    ax_heat.set_yticklabels([f"{int(heatmap_sizes[i]/1000)}" for i in tick_indices], fontsize=fs_tick)
    ax_heat.set_ylabel("Window\n(kbp)", fontsize=fs_tick, rotation=0, ha='right', va='center', labelpad=5)
    ax_heat.tick_params(axis='x', labelbottom=False)

    save_colorbar_legend(cmap, mcolors.Normalize(vmin=vmin, vmax=vmax), 
                         args.heatmap_col, f"{output_base}_legend_heatmap.png")

    raw_pos = lrt_df["pos"].to_numpy()
    
    if args.pval_col not in lrt_df.columns:
        print(f"Error: Column '{args.pval_col}' not found in CSV. Available: {lrt_df.columns}")
        return
        
    raw_pval = lrt_df[args.pval_col].to_numpy()
    if "Overdispersion" in lrt_df.columns:
        raw_od = lrt_df["Overdispersion"].fill_null(1.0).to_numpy()
    else:
        raw_od = np.ones(len(raw_pos))

    raw_logp = -np.log10(np.where(raw_pval <= 1e-300, 1e-300, raw_pval))
    
    use_median = (args.smooth_type == 'median')
    windows = sorted(args.smooth_window, reverse=True)
    
    global_max_y = 0

    for w_idx, win_size in enumerate(windows):
        s_pos, s_y, s_od_driver = smooth_and_aggregate(
            raw_pos, raw_logp, raw_od, win_size, use_median=use_median
        )
        
        od_metric = np.abs(s_od_driver - 1.0)
        sort_order = np.argsort(od_metric)[::-1] 
        
        s_pos_sorted = s_pos[sort_order]
        s_y_sorted = s_y[sort_order]
        
        c_solid = "#222222" 
        ax_glm.scatter(s_pos_sorted, s_y_sorted, c=c_solid, s=args.point_size, 
                       alpha=0.9, linewidths=0, zorder=10 + w_idx)
        
        if len(s_y) > 0:
            global_max_y = max(global_max_y, np.max(s_y))

    y_top = global_max_y * 1.1 if global_max_y > 1 else 2.0
    ax_glm.set_ylim(0, y_top)
    
    smooth_lbl = "Median" if use_median else "Mean"
    label_y = r"%s $-\log_{10}(P_{LRT})$" % smooth_lbl
    ax_glm.set_ylabel(label_y, fontsize=fs_label)
    
    ax_glm.tick_params(axis='y', labelsize=fs_tick)
    ax_glm.tick_params(axis='x', labelbottom=False)

    win_handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#222222', 
               label=f'{w}bp' if w>0 else 'Raw', markersize=6) 
        for w in sorted(windows)
    ]
    global_dir = os.path.dirname(output_base)
    save_elements_legend(win_handles, "Window", os.path.join(global_dir, "GLOBAL_legend_windows.png"))
    
    ax_annot.set_facecolor('none')
    gene_handles = []
    if annot_data and annot_labels:
        cat_map = {cat: CAT_COLORS[i % len(CAT_COLORS)] for i, cat in enumerate(annot_labels)}
        
        bar_h = args.gene_bar_height
        bar_y = (1.0 - bar_h) / 2.0 
        
        min_w = (end - start) * args.min_gene_width
        
        for (g_start, g_end, label) in annot_data:
            w = max(g_end - g_start, min_w)
            rect = Rectangle((g_start, bar_y), w, bar_h, 
                             facecolor=cat_map.get(label, 'k'), edgecolor=None, linewidth=0)
            ax_annot.add_patch(rect)
        gene_handles = [Patch(facecolor=cat_map[c], label=c) for c in annot_labels]

    ax_annot.set_ylabel("Genes", fontsize=fs_tick, rotation=0, ha='right', va='center', labelpad=5)
    ax_annot.set_yticks([])
    ax_annot.set_ylim(0, 1)

    if gene_handles:
        save_elements_legend(gene_handles, "Genes", f"{output_base}_legend_genes.png")

    ax_foot.axis('off')
    
    ax_annot.tick_params(axis='x', labelbottom=True, labelsize=fs_tick, rotation=45)
    
    def x_fmt(x, pos):
        return f'{x/1e6:.2f}M' if x > 1e6 else f'{x/1000:.0f}k'
    ax_annot.xaxis.set_major_formatter(mticker.FuncFormatter(x_fmt))
    
    for label in ax_annot.get_xticklabels():
        label.set_horizontalalignment('right')

    for ax in [ax_heat, ax_glm, ax_annot]:
        ax.set_xlim(start, end)

    plt.savefig(filename, bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_base}.svg", bbox_inches='tight')
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description="Generate Zoom plots from LRT Results CSV")
    parser.add_argument("--lrt-csv", required=True, help="Path to lrt_results.csv")
    parser.add_argument("--heatmap-files", nargs='+', required=True, help="List of FST/Heatmap csv files")
    parser.add_argument("--output-dir", required=True, help="Folder to save PNGs")
    parser.add_argument("--bed-file", type=str, help="Gene annotation BED")
    
    parser.add_argument("--pval-col", default="LRT_p", help="Column name for P-value (default: LRT_p)")
    
    parser.add_argument("--heatmap-col", default="z_score_median", help="Column in heatmap files")
    parser.add_argument("--buffer", type=int, default=10000, help="bp padding around the region")
    parser.add_argument("--bin-size", type=int, default=2000, help="Heatmap bin size")
    
    parser.add_argument("--smooth-type", choices=['median', 'mean'], default='median', 
                        help="Smoothing method for Y-axis points (default: median)")
    
    parser.add_argument("--smooth-window", type=int, nargs='+', default=[1000, 5000, 10000], 
                        help="List of window sizes (bp) to plot overlayed. 0 = Raw Points.")

    parser.add_argument("--point-size", type=float, default=15.0, help="Size of scatter points (default: 15)")
    parser.add_argument("--gene-bar-height", type=float, default=0.8, help="Height of gene bars (0.0-1.0, default: 0.8)")
    parser.add_argument("--min-gene-width", type=float, default=0.005, help="Min width of gene bar as fraction of plot width (default: 0.005)")
    
    parser.add_argument("--width", type=float, default=5.0, help="Plot width in inches")
    parser.add_argument("--height", type=float, default=6.0, help="Plot height in inches")
    parser.add_argument("--ratios", type=float, nargs='+', help="Height ratios: Heatmap GLM Genes Footer (4 values)")
    
    args = parser.parse_args()

    df = pl.read_csv(args.lrt_csv)
    if not os.path.exists(args.output_dir): os.makedirs(args.output_dir)

    if "Region_ID" not in df.columns:
        print("Warning: 'Region_ID' column not found. Creating dummy ID based on chrom/pos chunks if possible or exiting.")

    region_ids = df["Region_ID"].unique().to_list()
    print(f"--- Found {len(region_ids)} regions to plot ---")

    for rid in region_ids:
        focal_df = df.filter(pl.col("Region_ID") == rid)
        if focal_df.height == 0: continue
        
        chrom = focal_df["chrom"][0]
        min_pos = focal_df["pos"].min()
        max_pos = focal_df["pos"].max()
        
        start = max(0, min_pos - args.buffer)
        end = max_pos + args.buffer
        
        plot_df = df.filter(
            (pl.col("chrom") == chrom) &
            (pl.col("pos") >= start) &
            (pl.col("pos") <= end)
        )

        hm_matrix, hm_sizes = process_heatmap_region(
            args.heatmap_files, args.heatmap_col, chrom, start, end, args.bin_size
        )
        annot_data, annot_labels = load_bed_annotations(args.bed_file, chrom, start, end)
        
        out_base = os.path.join(args.output_dir, f"{rid}_{chrom}")
        
        render_zoom_plot(
            rid, chrom, start, end,
            plot_df,
            hm_matrix, hm_sizes,
            annot_data, annot_labels,
            args, out_base
        )

if __name__ == "__main__":
    main()