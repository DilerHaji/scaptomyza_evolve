#!/usr/bin/env python3
import argparse
import os
import re
import sys
import numpy as np
import polars as pl

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

def find_fst_cols_for_gen(cols, targets, refs, gen):
    found = []
    gen_str = f"g{gen}"
    for t in targets:
        for r in refs:
            t_gen = f"{t}g{gen}".lower()
            r_gen = f"{r}g{gen}".lower()
            for c in cols:
                cl = c.lower()
                if 'fst' not in cl:
                    continue
                if ((t_gen in cl and r_gen in cl) or
                    (r_gen in cl and t_gen in cl)):
                    found.append(c)
                    break
    return found


def detect_generations(cols, targets):
    gens = set()
    for t in targets:
        pattern = re.compile(rf"{re.escape(t)}G(\d+)", re.IGNORECASE)
        for c in cols:
            m = pattern.search(c)
            if m:
                gens.add(m.group(1))
    return sorted(gens, key=lambda x: int(x))


def clean_label(name):
    parts = re.split(r'[_]+', name)
    if parts and parts[-1].isdigit():
        return parts[-1]
    return parts[-1] if parts else name


def build_layout(df, top_n, spacing):
    lens = (df.select(['chrom', 'start'])
            .group_by('chrom')
            .agg([pl.col('start').min().alias('minp'),
                  pl.col('start').max().alias('maxp')])
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

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', required=True,
                    help='grenedalf FST queue CSV')
    ap.add_argument('--output', required=True,
                    help='Output Manhattan PNG')
    ap.add_argument('--output-svg', default=None)
    ap.add_argument('--output-csv', default=None,
                    help='Output per-generation median FST CSV')
    ap.add_argument('--target-reps', required=True,
                    help='Comma-separated target replicate prefixes (e.g. T1,T2,T3,T4)')
    ap.add_argument('--ref-reps', required=True,
                    help='Comma-separated reference replicate prefixes (e.g. B1,B2,B3,B4)')
    ap.add_argument('--generations', default=None,
                    help='Comma-separated generations (e.g. 01,02,06,07,08,09,10). '
                         'Auto-detected if omitted.')
    ap.add_argument('--spacing', type=int, default=3_200_000)
    ap.add_argument('--top-n', type=int, default=10)
    ap.add_argument('--ymax', type=float, default=None)
    ap.add_argument('--fig-width', type=float, default=14)
    ap.add_argument('--panel-height', type=float, default=1.8,
                    help='Height of each generation panel in inches')
    ap.add_argument('--cv-panel', action='store_true', default=True,
                    help='Add temporal CV panel at top (default: True)')
    ap.add_argument('--no-cv-panel', action='store_true',
                    help='Disable temporal CV panel')
    ap.add_argument('--cv-highlight-pct', type=float, default=90,
                    help='Percentile of mean FST above which to highlight windows')
    args = ap.parse_args()

    targets = [x.strip() for x in args.target_reps.split(',')]
    refs = [x.strip() for x in args.ref_reps.split(',')]

    with open(args.input, 'r') as fh:
        header_line = fh.readline().strip()
    cols = header_line.split(',')

    if args.generations:
        gens = [g.strip() for g in args.generations.split(',')]
    else:
        gens = detect_generations(cols, targets)
    if not gens:
        sys.exit("ERROR: No generations detected in column names.")

    chrom_col = next((c for c in cols if c.lower() in ['chrom', 'chr', 'chromosome']), cols[0])
    start_col = next((c for c in cols if c.lower() in ['start', 'pos', 'window_start']), cols[1])
    end_col = next((c for c in cols if c.lower() in ['end', 'window_end']), None)

    keep_cols = [chrom_col, start_col]
    if end_col:
        keep_cols.append(end_col)
    for g in gens:
        keep_cols.extend(find_fst_cols_for_gen(cols, targets, refs, g))
    keep_cols = list(dict.fromkeys(keep_cols))  # deduplicate, preserve order

    df_raw = pl.read_csv(args.input, columns=keep_cols)
    cols = df_raw.columns

    gen_median_cols = []
    for g in gens:
        fst_cols = find_fst_cols_for_gen(cols, targets, refs, g)
        if not fst_cols:
            continue
        med_name = f"median_fst_G{g}"
        arr = df_raw.select([pl.col(c).cast(pl.Float64) for c in fst_cols]).to_numpy()
        medians = np.nanmedian(arr, axis=1)
        df_raw = df_raw.with_columns(pl.Series(name=med_name, values=medians))
        gen_median_cols.append((g, med_name))

    if not gen_median_cols:
        sys.exit("ERROR: No FST data found for any generation.")

    order, lengths, offsets, mins = build_layout(df_raw, args.top_n, args.spacing)
    df = df_raw.filter(pl.col(chrom_col).is_in(order))

    off_df = pl.DataFrame({chrom_col: list(offsets.keys()),
                           '__off__': list(offsets.values())})
    min_df = pl.DataFrame({chrom_col: list(mins.keys()),
                           '__min__': list(mins.values())})
    df = (df.join(off_df, on=chrom_col).join(min_df, on=chrom_col)
            .with_columns([
                ((pl.col(start_col) - pl.col('__min__') + pl.col('__off__')) / 1e6)
                .alias('x_mb')
            ]))

    chrom_map = {c: i for i, c in enumerate(order)}
    c_series = df[chrom_col].to_list()
    chrom_idx = np.array([chrom_map.get(c, -1) for c in c_series])
    colors_pair = ['#4D4D4D', '#B3B3B3']
    chrom_colors = np.array([colors_pair[ci % 2] for ci in chrom_idx])

    x_mb = df['x_mb'].to_numpy()

    do_cv = args.cv_panel and not args.no_cv_panel and len(gen_median_cols) >= 3
    cv_vals = None
    mean_fst_vals = None
    if do_cv:
        all_gen_arrays = []
        for g, col_name in gen_median_cols:
            all_gen_arrays.append(df[col_name].to_numpy().astype(float))
        stacked = np.stack(all_gen_arrays, axis=0)  # (n_gens, n_windows)
        mean_fst_vals = np.nanmean(stacked, axis=0)
        std_fst_vals = np.nanstd(stacked, axis=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            cv_vals = np.where(mean_fst_vals > 0,
                               std_fst_vals / mean_fst_vals, np.nan)

    if args.output_csv:
        export_cols = [chrom_col, start_col]
        if end_col:
            export_cols.append(end_col)
        export_cols += [mc for _, mc in gen_median_cols]
        os.makedirs(os.path.dirname(os.path.abspath(args.output_csv)), exist_ok=True)
        df.select(export_cols).write_csv(args.output_csv)

    n_gen_panels = len(gen_median_cols)
    n_panels = n_gen_panels + (2 if do_cv else 0)  # CV track + scatter
    panel_heights = []
    panel_names = []

    if do_cv:
        panel_heights.append(args.panel_height * 0.8)  # CV Manhattan
        panel_names.append('cv')
        panel_heights.append(args.panel_height * 0.8)  # mean vs CV scatter
        panel_names.append('scatter')

    for g, _ in gen_median_cols:
        panel_heights.append(args.panel_height)
        panel_names.append(f'G{g}')

    fig_h = sum(panel_heights) + 0.5
    fig = plt.figure(figsize=(args.fig_width, fig_h), dpi=300)
    gs = gridspec.GridSpec(n_panels, 1, figure=fig, hspace=0.35,
                          height_ratios=panel_heights)

    if args.ymax:
        ymax = args.ymax
    else:
        all_vals = []
        for _, col_name in gen_median_cols:
            v = df[col_name].to_numpy()
            v = v[np.isfinite(v)]
            if len(v) > 0:
                all_vals.append(np.nanquantile(v, 0.999))
        ymax = max(0.1, max(all_vals) * 1.1) if all_vals else 0.5

    bounds = []
    for c in order:
        a = offsets[c] / 1e6
        b = (offsets[c] + lengths[c]) / 1e6
        bounds.append((c, a, b))
    xticks = [(a + b) / 2 for _, a, b in bounds]
    xlabels = [clean_label(c) for c, _, _ in bounds]

    panel_idx = 0

    if do_cv:
        ax_cv = fig.add_subplot(gs[panel_idx])
        panel_idx += 1

        high_mean_mask = mean_fst_vals > np.nanpercentile(
            mean_fst_vals[np.isfinite(mean_fst_vals)], args.cv_highlight_pct)
        low_cv_mask = cv_vals < np.nanpercentile(
            cv_vals[np.isfinite(cv_vals)], 25)  # bottom 25% CV
        candidate_mask = high_mean_mask & low_cv_mask & np.isfinite(cv_vals)

        cv_colors = np.where(candidate_mask, '#D55E00', chrom_colors)
        cv_sizes = np.where(candidate_mask, 8, 4)

        ax_cv.scatter(x_mb, cv_vals, c=cv_colors, s=cv_sizes,
                      alpha=0.7, linewidths=0, rasterized=True)
        cv_median = np.nanmedian(cv_vals[np.isfinite(cv_vals)])
        ax_cv.axhline(cv_median, color='#999999', linestyle=':', linewidth=0.5)
        ax_cv.set_ylabel('CV of FST(B,T)\nacross gens', fontsize=7)
        ax_cv.set_title('Temporal Stability: low CV + high mean FST = balancing selection candidate',
                        fontsize=7, fontweight='bold', loc='left')
        ax_cv.set_ylim(0, min(np.nanpercentile(cv_vals[np.isfinite(cv_vals)], 99.5), 5))
        ax_cv.set_xticks([])
        n_cand = np.sum(candidate_mask)
        ax_cv.text(0.99, 0.95, f'{n_cand} candidates (high FST, low CV)',
                   transform=ax_cv.transAxes, ha='right', va='top', fontsize=6,
                   color='#D55E00')

        ax_sc = fig.add_subplot(gs[panel_idx])
        panel_idx += 1

        finite = np.isfinite(mean_fst_vals) & np.isfinite(cv_vals)
        sc_colors = np.where(candidate_mask[finite], '#D55E00', '#999999')
        sc_sizes = np.where(candidate_mask[finite], 6, 2)
        ax_sc.scatter(mean_fst_vals[finite], cv_vals[finite],
                      c=sc_colors, s=sc_sizes, alpha=0.5,
                      linewidths=0, rasterized=True)
        ax_sc.set_xlabel('Mean FST(B,T)', fontsize=7)
        ax_sc.set_ylabel('CV', fontsize=7)
        ax_sc.set_title('Mean FST vs temporal variability', fontsize=7,
                        fontweight='bold', loc='left')
        ax_sc.set_ylim(0, min(np.nanpercentile(cv_vals[np.isfinite(cv_vals)], 99.5), 5))

        ax_sc.axhline(cv_median, color='#cccccc', linestyle=':', linewidth=0.5)
        fst_med = np.nanpercentile(mean_fst_vals[np.isfinite(mean_fst_vals)],
                                   args.cv_highlight_pct)
        ax_sc.axvline(fst_med, color='#cccccc', linestyle=':', linewidth=0.5)
        ax_sc.text(0.98, 0.02, 'Balanced\ncandidates', transform=ax_sc.transAxes,
                   ha='right', va='bottom', fontsize=6, color='#D55E00',
                   fontstyle='italic')

    for i, (g, col_name) in enumerate(gen_median_cols):
        ax = fig.add_subplot(gs[panel_idx])
        panel_idx += 1

        y_vals = df[col_name].to_numpy()

        ax.scatter(x_mb, y_vals, c=chrom_colors, s=4,
                   alpha=0.8, linewidths=0, rasterized=True)
        ax.set_ylim(0, ymax)
        ax.set_ylabel(f'FST\n(B–T)', fontsize=6)
        ax.set_title(f'Generation {int(g):02d}', fontsize=8,
                     fontweight='bold', loc='left')
        ax.grid(axis='y', linestyle='-', linewidth=0.3, alpha=0.2)

        med = np.nanmedian(y_vals[np.isfinite(y_vals)])
        ax.axhline(med, color='#999999', linestyle=':', linewidth=0.5, zorder=1)
        ax.text(0.99, 0.95, f'median={med:.4f}', transform=ax.transAxes,
                ha='right', va='top', fontsize=5, color='#666666')

        if i == n_gen_panels - 1:
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
