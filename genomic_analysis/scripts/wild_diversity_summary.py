#!/usr/bin/env python3
import argparse
import os
import re
from pathlib import Path

import numpy as np
import polars as pl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats as sp_stats


matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.rcParams['axes.spines.top'] = False
matplotlib.rcParams['axes.spines.right'] = False
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['pdf.fonttype'] = 42


WILD_POPS = ['AV', 'PS', 'RM']
HOSTS = ['B', 'T']
WILD_SAMPLES = [f'{p}{h}' for p in WILD_POPS for h in HOSTS]    # AVB, AVT, PSB, ...
FOUNDERS = ['F1G00', 'F2G00', 'F3G00', 'F4G00']
ALL_SAMPLES = WILD_SAMPLES + FOUNDERS

METRICS = ['theta_pi', 'theta_watterson', 'tajimas_d']
METRIC_LABEL = {
    'theta_pi':        r'$\theta_\pi$',
    'theta_watterson': r'$\theta_W$',
    'tajimas_d':       r"Tajima's D",
}

# colourblind-safe; B = orange, T = blue, founder = grey
HOST_COLOR = {'B': '#D55E00', 'T': '#0072B2', 'F': '#7f7f7f'}


def col(sample, metric):
    return f'{sample}.1.{metric}'


def col_passed(sample):
    return f'{sample}.1.passed'


def load_matrix(df, samples, metric):
    arrs = []
    for s in samples:
        vals = df[col(s, metric)].to_numpy().astype(float)
        passed = df[col_passed(s)].to_numpy().astype(float)
        vals = np.where(passed > 0, vals, np.nan)
        arrs.append(vals)
    return np.column_stack(arrs)


def summary_row(vals):
    v = vals[np.isfinite(vals)]
    if v.size == 0:
        return dict(n=0, mean=np.nan, se=np.nan, median=np.nan, q1=np.nan, q3=np.nan)
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    return dict(
        n=int(v.size),
        mean=float(np.mean(v)),
        se=float(np.std(v, ddof=1) / np.sqrt(v.size)),
        median=float(med),
        q1=float(q1),
        q3=float(q3),
    )


def scaffold_label(name):
    m = re.search(r'(\d+)_HRSCAF', name)
    return m.group(1) if m else name


def build_layout(df, top_n, gap):
    chroms = df['chrom'].to_numpy()
    ends = df['end'].to_numpy().astype(np.int64)
    unique_chroms = np.unique(chroms)
    chrom_maxend = {}
    for ch in unique_chroms:
        chrom_maxend[ch] = int(ends[chroms == ch].max())
    sorted_chroms = sorted(chrom_maxend.keys(), key=lambda c: chrom_maxend[c], reverse=True)
    order = sorted_chroms[:top_n]
    offsets, items, cur = {}, [], 0
    for ch in order:
        L = chrom_maxend[ch]
        offsets[ch] = cur
        items.append({'start': cur, 'width': L, 'mid': cur + L / 2,
                      'label': scaffold_label(ch), 'chrom': ch})
        cur += L + gap
    return order, offsets, items, cur



def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--diversity', required=True,
                    help='combined grenedalf diversity CSV')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--top-n', type=int, default=10,
                    help='top N scaffolds for Manhattan (by length)')
    ap.add_argument('--gap', type=int, default=1_000_000)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pl.read_csv(args.diversity)
    rows = []
    for s in ALL_SAMPLES:
        pop = s[:2] if s in WILD_SAMPLES else 'F'
        host = s[2] if s in WILD_SAMPLES else 'F'
        group = 'wild' if s in WILD_SAMPLES else 'founder'
        passed = df[col_passed(s)].to_numpy().astype(float)
        npassed_windows = int(np.sum(passed > 0))
        total_snps = int(np.nansum(passed))
        for m in METRICS:
            stats = summary_row(load_matrix(df, [s], m)[:, 0])
            rows.append({
                'sample': s, 'group': group, 'population': pop, 'host': host,
                'metric': m,
                'windows_nonzero': npassed_windows,
                'total_snps': total_snps,
                **stats,
            })
    tbl = pl.DataFrame(rows)
    tbl.write_csv(outdir / 'wild_diversity_table.tsv', separator='\t')

    paired_rows = []
    for pop in WILD_POPS:
        for m in METRICS:
            b = df[col(f'{pop}B', m)].to_numpy().astype(float)
            t = df[col(f'{pop}T', m)].to_numpy().astype(float)
            bp = df[col_passed(f'{pop}B')].to_numpy().astype(float)
            tp = df[col_passed(f'{pop}T')].to_numpy().astype(float)
            b = np.where(bp > 0, b, np.nan)
            t = np.where(tp > 0, t, np.nan)
            both = np.isfinite(b) & np.isfinite(t)
            if both.sum() < 10:
                stat = pval = np.nan
            else:
                stat, pval = sp_stats.wilcoxon(b[both], t[both])
            paired_rows.append({
                'population': pop, 'metric': m,
                'n_windows': int(both.sum()),
                'median_B': float(np.nanmedian(b)),
                'median_T': float(np.nanmedian(t)),
                'median_diff_B_minus_T': float(np.nanmedian(b[both] - t[both])),
                'wilcoxon_W': float(stat) if np.isfinite(stat) else np.nan,
                'p_value': float(pval) if np.isfinite(pval) else np.nan,
            })
    pl.DataFrame(paired_rows).write_csv(outdir / 'wild_paired_tests.tsv', separator='\t')

    wvf_rows = []
    for m in METRICS:
        wild_mat = load_matrix(df, WILD_SAMPLES, m)
        fnd_mat  = load_matrix(df, FOUNDERS, m)
        w = np.nanmean(wild_mat, axis=1)
        f = np.nanmean(fnd_mat,  axis=1)
        both = np.isfinite(w) & np.isfinite(f)
        if both.sum() >= 10:
            stat, pval = sp_stats.wilcoxon(w[both], f[both])
            r, _ = sp_stats.spearmanr(w[both], f[both])
        else:
            stat = pval = r = np.nan
        wvf_rows.append({
            'metric': m,
            'n_windows': int(both.sum()),
            'median_wild': float(np.nanmedian(w)),
            'median_founder': float(np.nanmedian(f)),
            'ratio_founder_over_wild': (float(np.nanmedian(f)) / float(np.nanmedian(w))
                                         if np.nanmedian(w) != 0 else np.nan),
            'median_diff_founder_minus_wild': float(np.nanmedian(f[both] - w[both])),
            'wilcoxon_W': float(stat) if np.isfinite(stat) else np.nan,
            'p_value': float(pval) if np.isfinite(pval) else np.nan,
            'spearman_r': float(r) if np.isfinite(r) else np.nan,
        })
    pl.DataFrame(wvf_rows).write_csv(outdir / 'wild_vs_founder_tests.tsv', separator='\t')

    fig, axes = plt.subplots(1, 3, figsize=(10, 3.2), dpi=200)
    for ax, m in zip(axes, METRICS):
        data, colors, labels = [], [], []
        for s in ALL_SAMPLES:
            v = load_matrix(df, [s], m)[:, 0]
            data.append(v[np.isfinite(v)])
            h = s[2] if s in WILD_SAMPLES else 'F'
            colors.append(HOST_COLOR[h])
            labels.append(s)
        bp = ax.boxplot(data, positions=np.arange(len(ALL_SAMPLES)),
                        widths=0.6, patch_artist=True, showfliers=False,
                        medianprops=dict(color='black', linewidth=0.8))
        for patch, c in zip(bp['boxes'], colors):
            patch.set_facecolor(c); patch.set_alpha(0.7); patch.set_linewidth(0.4)
        for element in ['whiskers', 'caps']:
            for ln in bp[element]:
                ln.set_linewidth(0.4)
        ax.set_xticks(np.arange(len(ALL_SAMPLES)))
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=6)
        ax.set_ylabel(METRIC_LABEL[m], fontsize=8)
        ax.axvline(len(WILD_SAMPLES) - 0.5, color='k', linestyle=':', linewidth=0.5, alpha=0.4)
        if m == 'tajimas_d':
            ax.axhline(0, color='k', linewidth=0.4, alpha=0.4)
    fig.suptitle('Wild vs founder diversity (per 390 kb window)', fontsize=9, y=1.02)
    fig.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(outdir / f'wild_diversity_boxplots.{ext}', bbox_inches='tight', dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(9, 3.0), dpi=200)
    for ax, m in zip(axes, METRICS):
        wild_mat = load_matrix(df, WILD_SAMPLES, m)
        fnd_mat  = load_matrix(df, FOUNDERS, m)
        w = np.nanmean(wild_mat, axis=1)
        f = np.nanmean(fnd_mat,  axis=1)
        ok = np.isfinite(w) & np.isfinite(f)
        ax.scatter(w[ok], f[ok], s=8, alpha=0.4, c='#333333', linewidths=0, rasterized=True)
        lo = np.nanmin(np.concatenate([w[ok], f[ok]]))
        hi = np.nanmax(np.concatenate([w[ok], f[ok]]))
        pad = 0.05 * (hi - lo) if hi > lo else 0.1
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                'k--', linewidth=0.5, alpha=0.5)
        r, _ = sp_stats.spearmanr(w[ok], f[ok])
        ax.set_xlabel(f'Wild mean {METRIC_LABEL[m]}', fontsize=8)
        ax.set_ylabel(f'Founder mean {METRIC_LABEL[m]}', fontsize=8)
        ax.set_title(f'{METRIC_LABEL[m]}  (ρ={r:.3f})', fontsize=9)
        ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
    fig.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(outdir / f'wild_vs_founder_scatter.{ext}', bbox_inches='tight', dpi=200)
    plt.close(fig)

    order, offsets, items, total = build_layout(df, args.top_n, args.gap)
    sub = df.filter(pl.col('chrom').is_in(order))
    chrom_arr = sub['chrom'].to_numpy()
    mid_arr = ((sub['start'].to_numpy() + sub['end'].to_numpy()) / 2).astype(float)
    x = np.array([offsets[c] + m for c, m in zip(chrom_arr, mid_arr)])

    MIN_SNPS = 500
    D_LO, D_HI = -1.0, 1.0


    fig, axes = plt.subplots(3, 1, figsize=(14, 7), dpi=200, sharex=True)
    for ax, pop in zip(axes, WILD_POPS):
        s_b, s_t = f'{pop}B', f'{pop}T'
        for s, color, label in [(s_b, HOST_COLOR['B'], f'{s_b} (B)'),
                                 (s_t, HOST_COLOR['T'], f'{s_t} (T)')]:
            passed = sub[col_passed(s)].to_numpy().astype(float)
            y_raw = sub[col(s, 'tajimas_d')].to_numpy().astype(float)
            y_raw = np.where(passed > 0, y_raw, np.nan)
            for ch in order:
                mask_ch = chrom_arr == ch
                xi = x[mask_ch]
                yi = y_raw[mask_ch]
                pi = passed[mask_ch]
                sort_idx = np.argsort(xi)
                xi, yi, pi = xi[sort_idx], yi[sort_idx], pi[sort_idx]
                ok = np.isfinite(yi)
                low = ok & (pi < MIN_SNPS)
                if low.any():
                    yi_low = np.clip(yi[low], D_LO, D_HI)
                    ax.scatter(xi[low], yi_low, s=2, c='#cccccc', alpha=0.3,
                               linewidths=0, zorder=1, rasterized=True)
                hi = ok & (pi >= MIN_SNPS)
                if hi.any():
                    in_range = hi & (yi >= D_LO) & (yi <= D_HI)
                    if in_range.any():
                        ax.plot(xi[in_range], yi[in_range], color=color,
                                linewidth=0.6, alpha=0.6, rasterized=True, zorder=2)
                        ax.scatter(xi[in_range], yi[in_range], s=3, c=color,
                                   alpha=0.5, linewidths=0, zorder=3, rasterized=True)
        ax.axhline(0, color='k', linewidth=0.5, alpha=0.5)
        ax.set_ylabel(f"{pop}\nTajima's D", fontsize=9, rotation=0,
                      ha='right', va='center', labelpad=45)
        ax.set_ylim(D_LO - 0.05, D_HI + 0.05)
        ax.tick_params(axis='y', labelsize=6)
        from matplotlib.lines import Line2D
        handles = [Line2D([0],[0], color=HOST_COLOR['B'], lw=1.5, label=f'{pop}B (B-host)'),
                   Line2D([0],[0], color=HOST_COLOR['T'], lw=1.5, label=f'{pop}T (T-host)')]
        ax.legend(handles=handles, fontsize=6, loc='upper right', framealpha=0.7)
    last = axes[-1]
    last.set_xticks([it['mid'] for it in items])
    last.set_xticklabels([it['label'] for it in items], rotation=0, fontsize=7)
    last.set_xlabel('Top-%d scaffolds (by length)' % args.top_n, fontsize=9)
    for ax in axes:
        for i, it in enumerate(items):
            if i % 2 == 1:
                ax.axvspan(it['start'], it['start'] + it['width'],
                           alpha=0.04, color='k', zorder=0)
        for it in items[:-1]:
            ax.axvline(it['start'] + it['width'] + args.gap / 2,
                       color='k', linewidth=0.3, alpha=0.15)
    fig.suptitle("Tajima's D across wild pools (390 kb windows)", fontsize=11, y=0.995)
    fig.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(outdir / f'wild_tajimasD_manhattan.{ext}', bbox_inches='tight', dpi=200)
    plt.close(fig)


if __name__ == '__main__':
    main()
