#!/usr/bin/env python3
import argparse
import re
import os
from pathlib import Path
from typing import List, Tuple, Dict, Any, Optional
import polars as pl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy import stats

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--target-reps", required=True)
    ap.add_argument("--ref-reps", required=True)
    ap.add_argument("--generations", required=False)
    ap.add_argument("--peaks", help="Candidate regions file (Region_ID, Chrom, Pos)")
    
    ap.add_argument("--spacing", type=int, default=3_200_000)
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--ymax", type=float, default=None)
    ap.add_argument("--font-scale", type=float, default=1.2)
    ap.add_argument("--show-xlabels", action="store_true")
    
    ap.add_argument("--format", default="png")
    ap.add_argument("--threshold", default=0.99)
    ap.add_argument("--window", default=1)
    
    return ap.parse_args()

def clean_contig_label(name: str) -> str:
    parts = re.split(r'[_]+', name)
    if parts:
        if parts[-1].isdigit(): return parts[-1]
        return parts[-1]
    return name

def per_contig_layout(df: pl.DataFrame, chrom_col: str, start_col: str, top_n: int, spacing: int):
    lens = (
        df.select([pl.col(chrom_col), pl.col(start_col).alias("pos")])
        .group_by(chrom_col)
        .agg([pl.col("pos").min().alias("minp"), pl.col("pos").max().alias("maxp")])
        .with_columns((pl.col("maxp") - pl.col("minp")).alias("len"))
        .sort("len", descending=True)
    )
    order = lens.select(chrom_col).head(top_n).to_series().to_list()
    lengths = {r[chrom_col]: int(r["len"]) for r in lens.iter_rows(named=True)}
    mins = {r[chrom_col]: int(r["minp"]) for r in lens.iter_rows(named=True)}
    
    offsets, acc = {}, 0
    for c in order:
        offsets[c] = acc
        acc += lengths[c] + spacing
    return order, lengths, offsets, mins

def detect_generations(cols: List[str], targets: List[str]) -> List[str]:
    gens = set()
    targets_clean = [t.strip() for t in targets]
    for t in targets_clean:
        pattern = re.compile(rf"{re.escape(t)}[_.]?G(\d+)", re.IGNORECASE)
        for c in cols:
            match = pattern.search(c)
            if match: gens.add(match.group(1))
    return sorted(list(gens), key=lambda x: int(x) if x.isdigit() else x)

def find_col_for_rep_gen(cols: List[str], target: str, ref: str, gen: str) -> Optional[str]:
    t = target.lower(); r = ref.lower(); g = gen.lower()
    for c in cols:
        cl = c.lower()
        if "fst" not in cl: continue
        if t in cl and r in cl:
            if f"g{g}" in cl or f".{g}" in cl or f"_{g}" in cl: return c
    return None

def plot_stacked_figure(outfile: str, row_configs: List[Dict], df: pl.DataFrame, 
                        peak_mask: np.ndarray, peak_rects: List, chrom_colors: np.ndarray, 
                        layout_info: Tuple, ymax_arg: float, font_scale: float):
    n_rows = len(row_configs)
    if n_rows == 0: return

    fig_h = 2.5 * n_rows
    fig = plt.figure(figsize=(14, fig_h))
    outer_grid = gridspec.GridSpec(n_rows, 1, figure=fig, hspace=0.3)
    xticks, xticklabels = layout_info
    x_vals_global = df["x_mb"].to_numpy()
    
    COL_BG = "#999999"
    COL_OUT = "#D55E00"

    for i, config in enumerate(row_configs):
        col_name = config["col"]
        title = config["title"]
        if col_name not in df.columns: continue
        y_vals = df[col_name].to_numpy()
        
        inner_grid = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer_grid[i], width_ratios=[4, 1], wspace=0.05)
        
        ax_man = fig.add_subplot(inner_grid[0])
        if ymax_arg: plot_ymax = ymax_arg
        else:
            if len(y_vals) > 0 and not np.all(np.isnan(y_vals)):
                local_max = np.nanquantile(y_vals, 0.999)
                plot_ymax = max(0.1, local_max * 1.1)
            else: plot_ymax = 1.0
        
        ax_man.scatter(x_vals_global, y_vals, c=chrom_colors, s=5, alpha=0.8, lw=0, rasterized=True)
        ax_man.set_ylim(0, plot_ymax)
        
        if peak_rects:
            bar_y = plot_ymax * 0.95
            for (p_start, p_end) in peak_rects:
                w = max(p_end - p_start, 0.5)
                m = (p_start + p_end) / 2
                ax_man.hlines(y=bar_y, xmin=m-(w/2), xmax=m+(w/2), colors=COL_OUT, linewidth=3, zorder=10)

        ax_man.set_ylabel("FST", fontsize=9)
        ax_man.set_title(title, loc="left", fontweight="bold", fontsize=10)
        ax_man.grid(axis="y", linestyle="-", linewidth=0.5, alpha=0.2)
        
        if i == n_rows - 1:
            final_ticks = xticks[:5]
            final_labels = xticklabels[:5]
            ax_man.set_xticks(final_ticks)
            ax_man.set_xticklabels(final_labels, rotation=0, fontsize=9)
            ax_man.set_xlabel("Chromosome (Top 5)", fontsize=9)
        else: ax_man.set_xticks([])

        ax_hist = fig.add_subplot(inner_grid[1], sharey=ax_man)
        bg_data = y_vals[~peak_mask]; out_data = y_vals[peak_mask]
        bg_data = bg_data[np.isfinite(bg_data)]; out_data = out_data[np.isfinite(out_data)]
        mu_bg = np.mean(bg_data) if len(bg_data) > 0 else 0
        mu_out = np.mean(out_data) if len(out_data) > 0 else 0
        
        if len(bg_data) > 5 and np.std(bg_data) > 1e-9:
            try:
                kde_bg = stats.gaussian_kde(bg_data)
                y_grid = np.linspace(0, plot_ymax, 100)
                x_grid = kde_bg(y_grid)
                ax_hist.fill_betweenx(y_grid, 0, x_grid, color=COL_BG, alpha=0.5, lw=0)
                ax_hist.plot(x_grid, y_grid, color=COL_BG, lw=1)
                ax_hist.hlines(mu_bg, 0, np.max(x_grid), color="black", linestyle="-", linewidth=1, alpha=0.7)
            except: pass

        if len(out_data) > 5 and np.std(out_data) > 1e-9:
            try:
                kde_out = stats.gaussian_kde(out_data)
                y_grid = np.linspace(0, plot_ymax, 100)
                x_grid = kde_out(y_grid)
                ax_hist.fill_betweenx(y_grid, 0, x_grid, color=COL_OUT, alpha=0.6, lw=0)
                ax_hist.plot(x_grid, y_grid, color=COL_OUT, lw=1)
                ax_hist.hlines(mu_out, 0, np.max(x_grid), color="black", linestyle="-", linewidth=1)
            except: pass
            
        ax_hist.axis("off")
        ax_hist.spines['left'].set_visible(True)
        ax_hist.spines['left'].set_linewidth(0.5)

    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)

def plot_region_history(
    region_id: str,
    outfile: str,
    median_col_names: List[Tuple[str, str]], # [(gen, colname)]
    df: pl.DataFrame,
    region_mask: np.ndarray,
    ymax_arg: float
):
    """Generates a separate plot for ONE region across generations."""
    
    n_gens = len(median_col_names)
    if n_gens == 0: return

    # Vertical stack of histograms
    fig, axes = plt.subplots(nrows=n_gens, ncols=1, figsize=(6, 2 * n_gens), sharex=True)
    if n_gens == 1: axes = [axes]
    
    COL_BG = "#999999" # Genome
    COL_REG = "#D55E00" # Region
    
    max_fst_global = 0.0
    for _, col in median_col_names:
        vals = df[col].to_numpy()
        vals = vals[np.isfinite(vals)]
        if len(vals) > 0:
            max_fst_global = max(max_fst_global, np.max(vals))
            
    if ymax_arg: max_fst_global = ymax_arg
    else: max_fst_global = max(0.1, max_fst_global * 1.1)

    for i, (gen, col_name) in enumerate(median_col_names):
        ax = axes[i]
        
        y_vals = df[col_name].to_numpy()
        
        bg_data = y_vals # Genome (Plotting all vs region is standard)
        reg_data = y_vals[region_mask] # Specific Region
        
        bg_data = bg_data[np.isfinite(bg_data)]
        reg_data = reg_data[np.isfinite(reg_data)]
        
        mu_bg = np.mean(bg_data) if len(bg_data) > 0 else 0
        mu_reg = np.mean(reg_data) if len(reg_data) > 0 else 0
        diff = mu_reg - mu_bg
        
        if len(bg_data) > 5:
            try:
                kde_bg = stats.gaussian_kde(bg_data)
                x_grid = np.linspace(0, max_fst_global, 200)
                y_grid = kde_bg(x_grid)
                ax.fill_between(x_grid, 0, y_grid, color=COL_BG, alpha=0.3, label="Genome")
                ax.axvline(mu_bg, color=COL_BG, linestyle="--", lw=1)
            except: pass
            
        if len(reg_data) > 1:
            try:
                if len(reg_data) < 5:
                    ax.hist(reg_data, bins=10, density=True, color=COL_REG, alpha=0.6, label=region_id)
                else:
                    kde_reg = stats.gaussian_kde(reg_data)
                    x_grid = np.linspace(0, max_fst_global, 200)
                    y_grid = kde_reg(x_grid)
                    ax.fill_between(x_grid, 0, y_grid, color=COL_REG, alpha=0.6, label=region_id)
                ax.axvline(mu_reg, color=COL_REG, linestyle="-", lw=1.5)
            except: pass
        
        ax.set_title(f"Gen {gen} | ΔFST: {diff:.4f}", loc="left", fontsize=10, fontweight="bold")
        ax.set_ylabel("Density")
        ax.grid(axis="x", linestyle=":", alpha=0.3)
        ax.set_xlim(0, max_fst_global)
        
        if i == 0:
            ax.legend(loc="upper right", fontsize=8)
            
    axes[-1].set_xlabel("FST")
    fig.suptitle(f"History of {region_id}", y=0.99)
    plt.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)

def main():
    args = parse_args()
    base = 10 * args.font_scale
    plt.rcParams.update({"font.size": base, "figure.dpi": 300})

    df_raw = pl.read_csv(args.input)
    cols = df_raw.columns
    chrom_col = next((c for c in cols if c.lower() in ["chrom", "chr", "chromosome"]), cols[0])
    start_col = next((c for c in cols if c.lower() in ["start", "pos", "window_start"]), cols[1])

    targets = args.target_reps.split(",")
    refs = args.ref_reps.split(",")
    
    if args.generations: gens = args.generations.split(",")
    else:
        gens = detect_generations(cols, targets)
        if not gens: raise SystemExit("No generations detected.")

    gen_to_cols = {}
    pair_to_gen_col = {} 
    for t in targets:
        for r in refs:
            pair_to_gen_col[f"{t}_vs_{r}"] = {}

    for g in gens:
        found_in_gen = []
        for t in targets:
            for r in refs:
                col = find_col_for_rep_gen(cols, t, r, g)
                if col:
                    found_in_gen.append(col)
                    pair_to_gen_col[f"{t}_vs_{r}"][g] = col
        gen_to_cols[g] = found_in_gen

    df_calc = df_raw
    median_col_names = []
    for g in gens:
        g_cols = gen_to_cols[g]
        if not g_cols: continue
        med_name = f"Median_G{g}"
        df_calc = df_calc.with_columns(pl.concat_list(g_cols).list.median().alias(med_name))
        median_col_names.append((g, med_name))

    order, lengths, offsets, mins = per_contig_layout(df_calc, chrom_col, start_col, args.top_n, args.spacing)
    df = df_calc.filter(pl.col(chrom_col).is_in(order))
    
    off_df = pl.DataFrame({chrom_col: list(offsets.keys()), "__off__": list(offsets.values())})
    min_df = pl.DataFrame({chrom_col: list(mins.keys()), "__min__": list(mins.values())})
    df = (
        df.join(off_df, on=chrom_col).join(min_df, on=chrom_col)
          .with_columns([
              ((pl.col(start_col) - pl.col("__min__") + pl.col("__off__")) / 1e6).alias("x_mb"),
              pl.col(start_col).alias("pos_bp")
          ])
    )
    
    peak_rects = []
    peak_mask = np.zeros(df.height, dtype=bool)
    
    region_masks = {} # { "Region_1": boolean_array, ... }
    
    if args.peaks and Path(args.peaks).exists():
        try:
            is_tsv = args.peaks.endswith(".tsv")
            separator = "\t" if is_tsv else ","
            p_df = pl.read_csv(args.peaks, separator=separator)
            
            if "Region_ID" in p_df.columns and "Pos" in p_df.columns:
                c_col_cand = next((c for c in p_df.columns if c.lower() in ["chrom", "chr"]), None)
                if c_col_cand:
                    p_df = p_df.filter(pl.col(c_col_cand).is_in(order))
                    
                    p_agg = (p_df.group_by([c_col_cand, "Region_ID"])
                            .agg([
                                pl.col("Pos").min().alias("peak_start"),
                                pl.col("Pos").max().alias("peak_end")
                            ])
                            .rename({c_col_cand: "chrom"}))
                    
                    d_chrom = df[chrom_col].to_numpy()
                    d_pos = df["pos_bp"].to_numpy()
                    
                    for r in p_agg.iter_rows(named=True):
                        c_mask = (d_chrom == r["chrom"])
                        s_mask = (d_pos >= r["peak_start"]) & (d_pos <= r["peak_end"])
                        peak_mask |= (c_mask & s_mask)
                    
                    for r in p_agg.iter_rows(named=True):
                        rid = r["Region_ID"]
                        c_mask = (d_chrom == r["chrom"])
                        s_mask = (d_pos >= r["peak_start"]) & (d_pos <= r["peak_end"])
                        region_masks[rid] = (c_mask & s_mask)

                    p_mapped = (
                        p_agg.join(off_df, on="chrom").join(min_df, on="chrom")
                            .with_columns([
                                ((pl.col("peak_start") - pl.col("__min__") + pl.col("__off__")) / 1e6).alias("x_start"),
                                ((pl.col("peak_end") - pl.col("__min__") + pl.col("__off__")) / 1e6).alias("x_end")
                            ])
                    )
                    peak_rects = p_mapped.select(["x_start", "x_end"]).rows()
        except Exception as e:
            pass

    bounds = []
    for c in order:
        a = offsets[c] / 1e6; b = (offsets[c] + lengths[c]) / 1e6
        bounds.append((c, a, b))
    xticks = [(a + b) / 2 for _, a, b in bounds]
    xticklabels = [clean_contig_label(c) for c, _, _ in bounds]
    
    try:
        mapping = {c: k for k, c in enumerate(order)}
        c_series = df[chrom_col].map_dict(mapping, default=-1)
    except: c_series = pl.Series([-1]*df.height)

    chrom_indices = c_series.to_numpy()
    colors_gray = ["#4D4D4D", "#B3B3B3"]
    chrom_colors = np.array([colors_gray[ci % 2] for ci in chrom_indices])

    median_rows = []
    for g, col in median_col_names:
        median_rows.append({"col": col, "title": f"Median FST - Gen {g}"})
    
    if median_rows:
        plot_stacked_figure(args.output, median_rows, df, peak_mask, peak_rects, 
                            chrom_colors, (xticks, xticklabels), args.ymax, args.font_scale)

    out_path = Path(args.output)
    parent_dir = out_path.parent
    stem = out_path.stem
    ext = out_path.suffix

    region_dir = parent_dir / "region_plots"
    region_dir.mkdir(exist_ok=True)
    
    for rid, mask in region_masks.items():
        if np.sum(mask) == 0: continue
        
        reg_outfile = str(region_dir / f"{stem}_{rid}_history{ext}")
        plot_region_history(rid, reg_outfile, median_col_names, df, mask, args.ymax)

if __name__ == "__main__":
    main()