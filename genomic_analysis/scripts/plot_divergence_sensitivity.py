#!/usr/bin/env python3
import argparse
import polars as pl
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import os
import sys
import numpy as np
from scipy.ndimage import gaussian_filter1d

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
mpl.rcParams['font.size'] = 7
mpl.rcParams['axes.titlesize'] = 8
mpl.rcParams['axes.labelsize'] = 7
mpl.rcParams['xtick.labelsize'] = 6
mpl.rcParams['ytick.labelsize'] = 6
mpl.rcParams['axes.linewidth'] = 0.8
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['legend.fontsize'] = 6

mpl.rcParams['svg.fonttype'] = 'none' 

def parse_args():
    parser = argparse.ArgumentParser(description="Focusing Divergence with Mixture Treatment")
    parser.add_argument("--candidates", required=True, help="Path to *_candidate_sites_detailed.tsv")
    parser.add_argument("--af-data", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output-dir", default="plots/divergence_focus_mix")
    parser.add_argument("--top-n", type=int, default=50)
    parser.add_argument("--max-dist", type=int, default=1000)
    parser.add_argument("--divergence-ylim", type=float, help="Symmetric Y-limit for Divergence plot.")
    parser.add_argument("--num-steps", type=int, default=20)
    parser.add_argument("--min-step", type=int, default=10)
    parser.add_argument("--ribbon-scale", type=float, default=0.5)
    parser.add_argument("--smoothness", type=float, default=1.0)
    parser.add_argument("--white-fade-start", type=float, default=0.0)
    parser.add_argument("--white-fade-end", type=float, default=0.85)
    parser.add_argument("--outline-width", type=float, default=0.5)
    parser.add_argument("--trt-1", default="B", help="Bottom Plot")
    parser.add_argument("--trt-2", default="T", help="Top Plot")
    parser.add_argument("--trt-mix", default="M", help="Middle Plot (Mixture)")    
    parser.add_argument("--col-1", default="#2E86F0", help="Color for Trt 1 (Blue)")
    parser.add_argument("--col-2", default="#E8A60C", help="Color for Trt 2 (Gold)")
    parser.add_argument("--col-mix", default="#67002e", help="Color for Mixture (Maroon/Purple)")
    parser.add_argument("--width", type=float, default=6.0, help="Figure width in inches")
    parser.add_argument("--height", type=float, default=4.5, help="Figure height in inches")
    parser.add_argument("--ratios", nargs=3, type=float, default=[1, 1, 1],
                        help="Height ratios for [Trt2, Mix, Trt1].")
    return parser.parse_args()

def normalize_cols(df_or_schema):
    ren = {}
    cols = df_or_schema if isinstance(df_or_schema, list) else df_or_schema.columns
    for c in cols:
        if c.lower() == "chrom": ren[c] = "CHROM"
        if c.lower() == "pos": ren[c] = "POS"
    return ren

def get_broad_data(af_lazy, meta, chrom, center, max_dist, t1, t2, t_mix):
    start, end = center - max_dist, center + max_dist
    ren = normalize_cols(af_lazy.collect_schema().names())
    
    af_df = af_lazy.rename(ren).filter(
        (pl.col("CHROM").cast(pl.String) == str(chrom)) & 
        (pl.col("POS") >= start) & (pl.col("POS") <= end)
    ).collect()
    
    if af_df.height == 0: return None

    id_vars = ["CHROM", "POS"]
    val_vars = [c for c in af_df.columns if c not in id_vars]
    long_df = af_df.unpivot(index=id_vars, on=val_vars, variable_name="sample", value_name="af")
    
    joined = long_df.join(meta, on="sample", how="inner").filter(
        pl.col("trt").is_in([t1, t2, t_mix])
    ).sort("gen")
    
    pdf = joined.to_pandas()
    if pdf.empty: return None

    pdf = pdf.sort_values(["POS", "replicate", "gen"])
    start_vals = pdf.groupby(["POS", "replicate"])["af"].transform("first")
    pdf["delta_af"] = pdf["af"] - start_vals
    
    return pdf

def apply_smoothing(x, y, sigma):
    if sigma <= 0: return x, y
    x_high = np.linspace(x.min(), x.max(), 200)
    y_linear = np.interp(x_high, x, y)
    scaled_sigma = sigma * 2.0 
    y_smooth = gaussian_filter1d(y_linear, scaled_sigma, mode='nearest')
    return x_high, y_smooth

def calculate_ribbon_empirical(df, scale_factor, sigma):
    stats_df = df.groupby('gen')['delta_af'].agg(['mean', 'sem']).reset_index().sort_values('gen')
    x = stats_df['gen'].values
    y_mean = stats_df['mean'].values
    width = stats_df['sem'].values * scale_factor
    
    y_l_raw = y_mean - width
    y_h_raw = y_mean + width
    
    x_s, y_l_smooth = apply_smoothing(x, y_l_raw, sigma)
    _, y_h_smooth = apply_smoothing(x, y_h_raw, sigma)
    
    return x_s, y_l_smooth, y_h_smooth

def calculate_diff_empirical_signed(sub, top_trt, bottom_trt, scale_factor, sigma):
    stats = sub.groupby(['gen', 'trt'])['delta_af'].agg(['mean', 'sem']).unstack()
    if top_trt not in stats['mean'] or bottom_trt not in stats['mean']: return None
    
    gens = stats.index.values
    m_top = stats['mean'][top_trt].values
    m_bot = stats['mean'][bottom_trt].values
    s_top = stats['sem'][top_trt].values
    s_bot = stats['sem'][bottom_trt].values
    
    diff_mean = m_top - m_bot
    sem_comb = np.sqrt(s_top**2 + s_bot**2)
    width = sem_comb * scale_factor
    
    y_l_raw = diff_mean - width
    y_h_raw = diff_mean + width
    
    x_s, y_l_smooth = apply_smoothing(gens, y_l_raw, sigma)
    _, y_h_smooth = apply_smoothing(gens, y_h_raw, sigma)
    
    return x_s, y_l_smooth, y_h_smooth

def mix_white(hex_color, amount_white):
    rgb = np.array(mcolors.hex2color(hex_color))
    white = np.array([1.0, 1.0, 1.0])
    mixed = rgb * (1 - amount_white) + white * amount_white
    return mixed

def generate_log_steps(min_val, max_val, num_steps):
    steps = np.geomspace(min_val, max_val, num_steps).astype(int)
    steps = np.unique(np.insert(steps, 0, 0))
    return steps

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    dists = generate_log_steps(args.min_step, args.max_dist, args.num_steps)
    dists_ordered = np.sort(dists) 
    n_bins = len(dists_ordered)
    fade_factors = np.linspace(args.white_fade_start, args.white_fade_end, n_bins)

    cand_df = pl.read_csv(args.candidates, separator='\t')
    peaks_df = (
        cand_df.sort("LogP", descending=True)
               .unique(subset=["Region_ID"], keep="first")
               .sort("Region_ID")
               .head(args.top_n)
    )
    peaks = peaks_df.select([
        pl.col("Chrom").alias("CHROM"), 
        pl.col("Pos").alias("POS"),
        pl.col("Region_ID")
    ]).to_dicts()
    
    af_lazy = pl.scan_csv(args.af_data)
    meta = pl.read_csv(args.metadata)
    if "replicate" not in meta.columns:
        meta = meta.with_columns(pl.lit("1").alias("replicate"))

    for i, p in enumerate(peaks):
        chrom = str(p['CHROM'])
        pos = int(p['POS'])
        rid = p['Region_ID']
        
        master_df = get_broad_data(af_lazy, meta, chrom, pos, args.max_dist, 
                                   args.trt_1, args.trt_2, args.trt_mix)
        if master_df is None: continue
        
        fig = plt.figure(figsize=(args.width, args.height))
        
        gs = fig.add_gridspec(3, 3, 
                              height_ratios=args.ratios, 
                              width_ratios=[1, 1, 0.05], 
                              wspace=0.55, hspace=0.15)
        
        ax_top = fig.add_subplot(gs[0, 0])
        ax_mid = fig.add_subplot(gs[1, 0], sharex=ax_top, sharey=ax_top)
        ax_bot = fig.add_subplot(gs[2, 0], sharex=ax_top, sharey=ax_top)
        
        ax_diff = fig.add_subplot(gs[:, 1]) 
        cbar_ax = fig.add_subplot(gs[:, 2])
        
        for step_idx, dist in enumerate(dists_ordered):
            z_order = 1000 + (step_idx * 10)
            if dist == 0: sub = master_df[master_df['POS'] == pos]
            else: mask = np.abs(master_df['POS'].values - pos) <= dist; sub = master_df[mask]
            if sub.empty: continue

            t_sub = sub[sub['trt'] == args.trt_2]
            if not t_sub.empty:
                fill_color = mix_white(args.col_2, fade_factors[step_idx])
                x, y_l, y_h = calculate_ribbon_empirical(t_sub, args.ribbon_scale, args.smoothness)
                ax_top.fill_between(x, y_l, y_h, color=fill_color, alpha=1.0, 
                                   edgecolor='white', linewidth=args.outline_width, zorder=z_order)

            m_sub = sub[sub['trt'] == args.trt_mix]
            if not m_sub.empty:
                fill_color = mix_white(args.col_mix, fade_factors[step_idx])
                x, y_l, y_h = calculate_ribbon_empirical(m_sub, args.ribbon_scale, args.smoothness)
                ax_mid.fill_between(x, y_l, y_h, color=fill_color, alpha=1.0, 
                                   edgecolor='white', linewidth=args.outline_width, zorder=z_order)

            b_sub = sub[sub['trt'] == args.trt_1]
            if not b_sub.empty:
                fill_color = mix_white(args.col_1, fade_factors[step_idx])
                x, y_l, y_h = calculate_ribbon_empirical(b_sub, args.ribbon_scale, args.smoothness)
                ax_bot.fill_between(x, y_l, y_h, color=fill_color, alpha=1.0, 
                                   edgecolor='white', linewidth=args.outline_width, zorder=z_order)

        z_start = 5000
        for step_idx, dist in enumerate(dists_ordered):
            col_diff = mix_white("#111111", fade_factors[step_idx])
            z_order = z_start + (step_idx * 10)
            if dist == 0: sub = master_df[master_df['POS'] == pos]
            else: mask = np.abs(master_df['POS'].values - pos) <= dist; sub = master_df[mask]
            if sub.empty: continue
            
            res = calculate_diff_empirical_signed(sub, args.trt_2, args.trt_1, args.ribbon_scale, args.smoothness)
            if res:
                x, y_l, y_h = res
                ax_diff.fill_between(x, y_l, y_h, color=col_diff, alpha=1.0, 
                                     edgecolor='white', linewidth=args.outline_width, zorder=z_order)

        def add_refs(ax):
            ax.axhline(0, color='gray', lw=0.8, ls=':', zorder=5000)
            ax.axhline(0.1, color='gray', lw=0.8, ls='--', alpha=0.5, zorder=5000)
            ax.axhline(-0.1, color='gray', lw=0.8, ls='--', alpha=0.5, zorder=5000)

        add_refs(ax_top)
        ax_top.set_title(r"A. $\Delta$ Allele Frequency", fontweight='bold', loc='left', fontsize=8)
        ax_top.text(0.05, 0.85, args.trt_2, transform=ax_top.transAxes, color=args.col_2, fontweight='bold', fontsize=8, zorder=6000)
        plt.setp(ax_top.get_xticklabels(), visible=False)
        ax_top.set_ylabel(r"$\Delta$ AF", fontsize=7)

        add_refs(ax_mid)
        ax_mid.text(0.05, 0.85, args.trt_mix, transform=ax_mid.transAxes, color=args.col_mix, fontweight='bold', fontsize=8, zorder=6000)
        plt.setp(ax_mid.get_xticklabels(), visible=False)
        ax_mid.set_ylabel(r"$\Delta$ AF", fontsize=7)

        add_refs(ax_bot)
        ax_bot.text(0.05, 0.85, args.trt_1, transform=ax_bot.transAxes, color=args.col_1, fontweight='bold', fontsize=8, zorder=6000)
        ax_bot.set_xlabel("Generation", fontsize=7)
        ax_bot.set_ylabel(r"$\Delta$ AF", fontsize=7)

        ax_diff.axhline(0, color='gray', lw=1, ls='-', zorder=5000)
        ax_diff.set_title(f"B. Divergence ({args.trt_2} vs {args.trt_1})", fontweight='bold', loc='left', fontsize=8)
        ax_diff.set_xlabel("Generation", fontsize=7)
        ax_diff.set_ylabel(f"Diff ({args.trt_2} - {args.trt_1})", fontsize=7, labelpad=1)
        if args.divergence_ylim: ax_diff.set_ylim(-args.divergence_ylim, args.divergence_ylim)
        
        ax_diff.text(0.02, 0.98, f"Favors {args.trt_2}", transform=ax_diff.transAxes, 
                     color=args.col_2, fontweight='bold', ha='left', va='top', zorder=10000, fontsize=7)
        ax_diff.text(0.02, 0.02, f"Favors {args.trt_1}", transform=ax_diff.transAxes, 
                     color=args.col_1, fontweight='bold', ha='left', va='bottom', zorder=10000, fontsize=7)

        cbar_colors = [mix_white("#111111", f) for f in fade_factors]
        cmap_discrete = mcolors.ListedColormap(cbar_colors)
        norm = mpl.colors.LogNorm(vmin=args.min_step, vmax=args.max_dist)
        cb = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap_discrete, norm=norm, orientation='vertical')
        cb.set_label("Pooling Radius (bp)", fontsize=7)
        cb.ax.tick_params(labelsize=6)
        
        clean_name = chrom.split("_")[-1] if len(chrom) > 15 else chrom
        plt.suptitle(f"{rid} | {clean_name}:{pos}", y=0.98, fontsize=9)
        
        outfile_png = os.path.join(args.output_dir, f"divergence_{rid}_{chrom}_{pos}.png")
        plt.savefig(outfile_png, bbox_inches='tight', dpi=300)

        outfile_svg = os.path.join(args.output_dir, f"divergence_{rid}_{chrom}_{pos}.svg")
        plt.savefig(outfile_svg, bbox_inches='tight')
        
        plt.close()

if __name__ == "__main__":
    main()