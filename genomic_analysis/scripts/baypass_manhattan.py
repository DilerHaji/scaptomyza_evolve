#!/usr/bin/env python3
import argparse
import os
import re
import sys
import numpy as np
import polars as pl
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.font_manager as font_manager

font_dirs = ['/usr/share/fonts/truetype/msttcorefonts', '/usr/share/fonts/dejavu']
for d in font_dirs:
    for ff in font_manager.findSystemFonts(fontpaths=[d]):
        try:
            font_manager.fontManager.addfont(ff)
        except Exception:
            pass
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.rcParams['axes.spines.top'] = False
matplotlib.rcParams['axes.spines.right'] = False
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['pdf.fonttype'] = 42


def clean_label(name):
    parts = re.split(r'[_]+', name)
    if parts and parts[-1].isdigit():
        return parts[-1]
    return parts[-1] if parts else name


def build_layout(df, top_n, spacing):
    lens = (df.select(['chrom', 'pos'])
            .group_by('chrom')
            .agg([pl.col('pos').min().alias('minp'),
                  pl.col('pos').max().alias('maxp')])
            .with_columns((pl.col('maxp') - pl.col('minp')).alias('len'))
            .sort('len', descending=True))

    order = lens.head(top_n)['chrom'].to_list()
    lengths = {r['chrom']: int(r['len']) for r in lens.iter_rows(named=True)}
    mins = {r['chrom']: int(r['minp']) for r in lens.iter_rows(named=True)}

    offsets, acc = {}, 0
    for c in order:
        offsets[c] = acc
        acc += lengths[c] + spacing
    return order, lengths, offsets, mins


def manhattan_panel(ax, x_mb, y_vals, chrom_colors, ylabel, title,
                    ymax=None, threshold=None, threshold_label=None,
                    point_size=2, alpha=0.7, highlight_mask=None,
                    highlight_color='#D55E00'):
    sizes = np.full(len(x_mb), point_size)
    colors = chrom_colors.copy()

    if highlight_mask is not None:
        colors[highlight_mask] = highlight_color
        sizes[highlight_mask] = point_size * 2.5

    ax.scatter(x_mb, y_vals, c=colors, s=sizes,
               alpha=alpha, linewidths=0, rasterized=True)

    if ymax is not None:
        ax.set_ylim(0, ymax)

    ax.set_ylabel(ylabel, fontsize=7)
    ax.set_title(title, fontsize=8, fontweight='bold', loc='left')
    ax.grid(axis='y', linestyle='-', linewidth=0.3, alpha=0.2)

    if threshold is not None:
        ax.axhline(threshold, color='#D55E00', linestyle='--',
                    linewidth=0.7, alpha=0.8, zorder=1)
        if threshold_label:
            ax.text(0.99, threshold / ax.get_ylim()[1] + 0.02,
                    threshold_label, transform=ax.get_yaxis_transform(),
                    ha='right', va='bottom', fontsize=5, color='#D55E00')

    # Median line
    finite = np.isfinite(y_vals)
    if finite.any():
        med = np.nanmedian(y_vals[finite])
        ax.axhline(med, color='#999999', linestyle=':', linewidth=0.5, zorder=1)
        ax.text(0.99, 0.95, f'median={med:.2f}', transform=ax.transAxes,
                ha='right', va='top', fontsize=5, color='#666666')

    if highlight_mask is not None:
        n_high = np.sum(highlight_mask)
        ax.text(0.01, 0.95, f'{n_high} outliers',
                transform=ax.transAxes, ha='left', va='top', fontsize=5,
                color=highlight_color)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--snp-pos', required=True,
                    help='SNP position CSV (mrk,chrom,pos)')
    ap.add_argument('--xtx', default=None,
                    help='Core model pi_xtx summary file')
    ap.add_argument('--betai', default=None,
                    help='Covariate model betai_reg summary file')
    ap.add_argument('--contrast', default=None,
                    help='Contrast model summary file')
    ap.add_argument('--contrast-names', default='B vs T,B vs M,T vs M',
                    help='Comma-separated contrast labels')
    ap.add_argument('--output', required=True,
                    help='Output Manhattan PNG')
    ap.add_argument('--output-svg', default=None)
    ap.add_argument('--spacing', type=int, default=3_200_000)
    ap.add_argument('--top-n', type=int, default=10)
    ap.add_argument('--fig-width', type=float, default=14)
    ap.add_argument('--panel-height', type=float, default=2.0)
    ap.add_argument('--xtx-threshold', type=float, default=None,
                    help='XtX significance threshold (e.g. from POD simulations)')
    ap.add_argument('--bf-threshold', type=float, default=20,
                    help='BF(dB) threshold for "decisive" evidence (default 20 dB)')
    args = ap.parse_args()

    if not any([args.xtx, args.betai, args.contrast]):
        sys.exit("ERROR: Provide at least one of --xtx, --betai, --contrast")

    contrast_names = [x.strip() for x in args.contrast_names.split(',')]


    pos_df = pl.read_csv(args.snp_pos)
    n_snps = len(pos_df)

    panels = []  # list -> title, ylabel, values, threshold, threshold_label

    if args.xtx:
        xtx_df = pl.from_pandas(pd.read_csv(args.xtx, delim_whitespace=True))
        if 'XtXst' in xtx_df.columns:
            xtx_vals = xtx_df['XtXst'].to_numpy().astype(float)
            panels.append(('XtX — Core Differentiation (standardized)',
                           'XtXst', xtx_vals,
                           args.xtx_threshold, 'threshold' if args.xtx_threshold else None))
        elif 'M_XtX' in xtx_df.columns:
            xtx_vals = xtx_df['M_XtX'].to_numpy().astype(float)
            panels.append(('XtX — Core Differentiation',
                           'XtX', xtx_vals, None, None))
 
    if args.betai:
        bf_df = pl.from_pandas(pd.read_csv(args.betai, delim_whitespace=True))
        bf_col = None
        for c in bf_df.columns:
            if 'BF' in c:
                bf_col = c
                break
        if bf_col is None:
            bf_col = bf_df.columns[-1]
            print(f"  Using last column as BF: {bf_col}")
        bf_vals = bf_df[bf_col].to_numpy().astype(float)
        if len(bf_vals) == n_snps:
            panels.append((f'BF(dB) — Treatment Association',
                           'BF (dB)', bf_vals,
                           args.bf_threshold, f'{args.bf_threshold} dB'))
        else:
            n_cov = len(bf_vals) // n_snps
            if 'COVARIABLE' in bf_df.columns:
                bf_sub = bf_df.filter(pl.col('COVARIABLE') == 1)
                bf_vals = bf_sub[bf_col].to_numpy().astype(float)
            else:
                bf_vals = bf_vals[:n_snps]
            panels.append((f'BF(dB) — Treatment Association',
                           'BF (dB)', bf_vals,
                           args.bf_threshold, f'{args.bf_threshold} dB'))

    if args.contrast:
        c2_df = pl.from_pandas(pd.read_csv(args.contrast, delim_whitespace=True))
        c2_col = None
        for c in c2_df.columns:
            if 'C2' in c:
                c2_col = c
                break
        if c2_col is None:
            c2_col = c2_df.columns[-1]

        if 'CONTRAST' in c2_df.columns:
            contrast_ids = c2_df['CONTRAST'].unique().sort().to_list()
            for ci, cid in enumerate(contrast_ids):
                label = contrast_names[ci] if ci < len(contrast_names) else f"Contrast {cid}"
                c2_sub = c2_df.filter(pl.col('CONTRAST') == cid)
                c2_vals = c2_sub[c2_col].to_numpy().astype(float)
                if len(c2_vals) == n_snps:
                    panels.append((f'C2 — {label}',
                                   'C2', c2_vals, None, None))
        else:
            c2_vals = c2_df[c2_col].to_numpy().astype(float)
            if len(c2_vals) == n_snps:
                panels.append(('C2 — Contrast', 'C2', c2_vals, None, None))

    if not panels:
        sys.exit("ERROR: No valid panels to plot.")

    order, lengths, offsets, mins = build_layout(pos_df, args.top_n, args.spacing)
    df = pos_df.filter(pl.col('chrom').is_in(order))

    off_df = pl.DataFrame({'chrom': list(offsets.keys()),
                           '__off__': list(offsets.values())})
    min_df = pl.DataFrame({'chrom': list(mins.keys()),
                           '__min__': list(mins.values())})
    df = (df.join(off_df, on='chrom').join(min_df, on='chrom')
            .with_columns([
                ((pl.col('pos') - pl.col('__min__') + pl.col('__off__')) / 1e6)
                .alias('x_mb')
            ]))

    keep_mask = pos_df['chrom'].is_in(order).to_numpy()

    x_mb = df['x_mb'].to_numpy()

    chrom_map = {c: i for i, c in enumerate(order)}
    c_series = df['chrom'].to_list()
    chrom_idx = np.array([chrom_map.get(c, -1) for c in c_series])
    colors_pair = ['#4D4D4D', '#B3B3B3']
    chrom_colors = np.array([colors_pair[ci % 2] for ci in chrom_idx])

    bounds = []
    for c in order:
        a = offsets[c] / 1e6
        b = (offsets[c] + lengths[c]) / 1e6
        bounds.append((c, a, b))
    xticks = [(a + b) / 2 for _, a, b in bounds]
    xlabels = [clean_label(c) for c, _, _ in bounds]

    n_panels = len(panels)
    panel_heights = [args.panel_height] * n_panels
    fig_h = sum(panel_heights) + 0.8
    fig = plt.figure(figsize=(args.fig_width, fig_h), dpi=300)
    gs = gridspec.GridSpec(n_panels, 1, figure=fig, hspace=0.4,
                           height_ratios=panel_heights)

    for i, (title, ylabel, values, threshold, threshold_label) in enumerate(panels):
        ax = fig.add_subplot(gs[i])

        vals = values[keep_mask] if len(values) == n_snps else values
        if len(vals) != len(x_mb):
            print(f"  WARNING: Panel '{title}' has {len(vals)} values, "
                  f"expected {len(x_mb)}. Skipping.")
            ax.text(0.5, 0.5, 'Data length mismatch', transform=ax.transAxes,
                    ha='center', va='center', fontsize=10, color='red')
            continue

        finite_vals = vals[np.isfinite(vals)]
        if len(finite_vals) > 0:
            ymax = np.nanquantile(finite_vals, 0.999) * 1.15
            ymax = max(ymax, 1.0)
        else:
            ymax = 1.0

        if len(finite_vals) > 0:
            q995 = np.nanquantile(finite_vals, 0.995)
            highlight = np.isfinite(vals) & (vals >= q995)
        else:
            highlight = None

        manhattan_panel(ax, x_mb, vals, chrom_colors, ylabel, title,
                        ymax=ymax, threshold=threshold,
                        threshold_label=threshold_label,
                        highlight_mask=highlight)

        if i == n_panels - 1:
            final_ticks = xticks[:min(len(xticks), 8)]
            final_labels = xlabels[:min(len(xlabels), 8)]
            ax.set_xticks(final_ticks)
            ax.set_xticklabels(final_labels, rotation=0, fontsize=7)
            ax.set_xlabel('Chromosome', fontsize=8)
        else:
            ax.set_xticks([])

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    fig.savefig(args.output, bbox_inches='tight', dpi=300)

    if args.output_svg:
        fig.savefig(args.output_svg, bbox_inches='tight', format='svg')

    plt.close(fig)


if __name__ == '__main__':
    main()
