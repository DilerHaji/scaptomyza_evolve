#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
matplotlib.rcParams['font.size'] = 7
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.rcParams['axes.spines.top'] = False
matplotlib.rcParams['axes.spines.right'] = False
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['pdf.fonttype'] = 42

GROUPS = {
    'Wild B':    ['AVB', 'PSB', 'RMB'],
    'Wild T':    ['AVT', 'PST', 'RMT'],
    'Founders':  ['F1G00', 'F2G00', 'F3G00', 'F4G00'],
    'G10 B':     ['B1G10', 'B2G10', 'B3G10', 'B4G10'],
    'G10 T':     ['T1G10', 'T2G10', 'T3G10', 'T4G10'],
    'G10 M':     ['M1G10', 'M2G10', 'M3G10', 'M4G10'],
}

GROUP_ORDER = ['Wild B', 'Wild T', 'Founders', 'G10 B', 'G10 T', 'G10 M']
GROUP_COLOR = {
    'Wild B': '#D55E00', 'Wild T': '#0072B2', 'Founders': '#7f7f7f',
    'G10 B': '#D55E00', 'G10 T': '#0072B2', 'G10 M': '#009E73',
}

METRICS = ['theta_pi', 'theta_watterson', 'tajimas_d']
METRIC_LABEL = {
    'theta_pi': r'$\theta_\pi$',
    'theta_watterson': r'$\theta_W$',
    'tajimas_d': "Tajima's D",
}

MIN_SNPS = 500


def get_vals(df, samples, metric):
    vals = []
    for s in samples:
        p = df[f'{s}.1.passed'].values.astype(float)
        v = df[f'{s}.1.{metric}'].values.astype(float)
        vals.append(v[(p >= MIN_SNPS) & np.isfinite(v)])
    return np.concatenate(vals)


def get_per_pool_median(df, samples, metric):
    meds = []
    for s in samples:
        p = df[f'{s}.1.passed'].values.astype(float)
        v = df[f'{s}.1.{metric}'].values.astype(float)
        v = v[(p >= MIN_SNPS) & np.isfinite(v)]
        meds.append(float(np.median(v)) if len(v) > 0 else np.nan)
    return meds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--diversity', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.diversity)



    rows = []
    for g in GROUP_ORDER:
        for m in METRICS:
            pool_meds = get_per_pool_median(df, GROUPS[g], m)
            pooled = get_vals(df, GROUPS[g], m)
            rows.append({
                'group': g,
                'metric': m,
                'n_pools': len(GROUPS[g]),
                'n_windows_pooled': len(pooled),
                'pooled_median': float(np.median(pooled)),
                'pooled_q1': float(np.percentile(pooled, 25)),
                'pooled_q3': float(np.percentile(pooled, 75)),
                'per_pool_median_mean': float(np.nanmean(pool_meds)),
                'per_pool_median_sd': float(np.nanstd(pool_meds, ddof=1)) if len(pool_meds) > 1 else np.nan,
            })
    tbl = pd.DataFrame(rows)
    tbl.to_csv(outdir / 'attrition_summary_table.tsv', sep='\t', index=False,
               float_format='%.4f')


    for m in METRICS:
        for g in GROUP_ORDER:
            sub = tbl[(tbl['group'] == g) & (tbl['metric'] == m)]
            med = sub['pooled_median'].values[0]



    ratio_rows = []
    for m in ['theta_pi', 'theta_watterson']:
        wild_med = np.median(get_vals(df, GROUPS['Wild B'] + GROUPS['Wild T'], m))
        fnd_med = np.median(get_vals(df, GROUPS['Founders'], m))
        for g10_name in ['G10 B', 'G10 T', 'G10 M']:
            g10_med = np.median(get_vals(df, GROUPS[g10_name], m))
            ratio_rows.append({
                'metric': m,
                'comparison': f'founder / wild',
                'ratio': fnd_med / wild_med if wild_med != 0 else np.nan,
            })
            ratio_rows.append({
                'metric': m,
                'comparison': f'{g10_name} / founder',
                'ratio': g10_med / fnd_med if fnd_med != 0 else np.nan,
            })
            ratio_rows.append({
                'metric': m,
                'comparison': f'{g10_name} / wild',
                'ratio': g10_med / wild_med if wild_med != 0 else np.nan,
            })
    ratios = pd.DataFrame(ratio_rows).drop_duplicates()
    ratios.to_csv(outdir / 'attrition_ratios.tsv', sep='\t', index=False,
                  float_format='%.4f')

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=200)

    for ax, m in zip(axes, METRICS):
        data_list = []
        for g in GROUP_ORDER:
            data_list.append(get_vals(df, GROUPS[g], m))

        positions = np.arange(len(GROUP_ORDER))
        parts = ax.violinplot(data_list, positions=positions,
                              showmedians=False, showextrema=False)
        for pc, g in zip(parts['bodies'], GROUP_ORDER):
            pc.set_facecolor(GROUP_COLOR[g])
            pc.set_alpha(0.25)
            pc.set_edgecolor(GROUP_COLOR[g])
            pc.set_linewidth(0.8)


        rng = np.random.default_rng(42)
        for i, (data, g) in enumerate(zip(data_list, GROUP_ORDER)):
            jitter = rng.uniform(-0.3, 0.3, size=len(data))
            ax.scatter(i + jitter, data, s=1.5, c=GROUP_COLOR[g], alpha=0.15,
                       linewidths=0, rasterized=True, zorder=1)


        for i, (data, g) in enumerate(zip(data_list, GROUP_ORDER)):
            q1, med, q3 = np.percentile(data, [25, 50, 75])
            ax.vlines(i, q1, q3, color=GROUP_COLOR[g], linewidth=2.5, zorder=4)
            ax.scatter(i, med, color='white', s=15, zorder=5,
                       edgecolors=GROUP_COLOR[g], linewidths=0.8)

        if m == 'tajimas_d':
            ax.axhline(0, color='k', linewidth=0.5, linestyle='--', alpha=0.5)


        ax.axvline(1.5, color='k', linewidth=0.5, linestyle=':', alpha=0.3)
        ax.axvline(2.5, color='k', linewidth=0.5, linestyle=':', alpha=0.3)

        ax.set_xticks(positions)
        ax.set_xticklabels([g.replace(' ', '\n') for g in GROUP_ORDER],
                           fontsize=6, rotation=0)
        ax.set_ylabel(METRIC_LABEL[m], fontsize=9)
        ax.tick_params(axis='y', labelsize=6)

    fig.suptitle('Diversity attrition: wild → founders → G10  (≥500 SNPs/window)',
                 fontsize=10, y=1.01)
    fig.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(outdir / f'attrition_violins.{ext}', bbox_inches='tight', dpi=200)
    plt.close(fig)


if __name__ == '__main__':
    main()
