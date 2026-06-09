#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
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


def scaffold_label(name):
    m = re.search(r'(\d+)_HRSCAF', str(name))
    return m.group(1) if m else str(name)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--betai',    required=True)
    ap.add_argument('--contrast', required=True)
    ap.add_argument('--positions', required=True)
    ap.add_argument('--outdir',   required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    betai = pd.read_csv(args.betai, sep=r'\s+')
    contrast = pd.read_csv(args.contrast, sep=r'\s+')
    pos = pd.read_csv(args.positions)

    df = (pos.merge(
              betai[['MRK', 'BF(dB)', 'Beta_is', 'SD_Beta_is', 'eBPis']]
                   .rename(columns={'MRK': 'mrk'}),
              on='mrk', how='inner')
             .merge(
              contrast[['MRK', 'C2', 'log10(1/pval)']]
                   .rename(columns={'MRK': 'mrk'}),
              on='mrk', how='inner'))

    bf = df['BF(dB)'].values
    beta = df['Beta_is'].values
    c2 = df['C2'].values
    lp = df['log10(1/pval)'].values

    row = {
        'n_snps':          int(len(df)),
        'bf_min':          float(np.nanmin(bf)),
        'bf_max':          float(np.nanmax(bf)),
        'bf_median':       float(np.nanmedian(bf)),
    }
    for p in [50, 75, 90, 95, 99, 99.9]:
        row[f'bf_pct_{p}'] = float(np.nanpercentile(bf, p))
    row['n_bf_gt_10']     = int(np.sum(bf > 10))
    row['n_bf_gt_15']     = int(np.sum(bf > 15))
    row['n_bf_gt_20']     = int(np.sum(bf > 20))
    row['frac_bf_gt_10']  = float(np.mean(bf > 10))

    row['frac_beta_positive'] = float(np.mean(beta > 0))
    row['frac_beta_negative'] = float(np.mean(beta < 0))
    top_bf_mask = bf > 10
    row['mean_abs_beta_top_bf'] = (float(np.mean(np.abs(beta[top_bf_mask])))
                                    if top_bf_mask.sum() else float('nan'))

    row['c2_median']     = float(np.nanmedian(c2))
    row['c2_pct_95']     = float(np.nanpercentile(c2, 95))
    row['c2_pct_99']     = float(np.nanpercentile(c2, 99))
    row['n_log10p_gt_3'] = int(np.sum(lp > 3))
    row['n_log10p_gt_5'] = int(np.sum(lp > 5))

    rho_bf_c2, _ = sp_stats.spearmanr(bf, lp)
    pearson_bf_c2 = float(np.corrcoef(bf, lp)[0, 1])
    row['spearman_bf_vs_c2logpval'] = float(rho_bf_c2)
    row['pearson_bf_vs_c2logpval']  = pearson_bf_c2


    cutoff = max(1, int(0.05 * len(df)))
    top_bf = set(np.argpartition(-bf, cutoff)[:cutoff])
    top_lp = set(np.argpartition(-lp, cutoff)[:cutoff])
    row['top5pct_bf_c2_overlap_frac'] = len(top_bf & top_lp) / len(top_bf)

    pd.DataFrame([row]).to_csv(outdir / 'wild_baypass_summary.tsv',
                                sep='\t', index=False)


    top100 = (df.nlargest(100, 'BF(dB)')
                .assign(scaffold=lambda d: d['chrom'].map(scaffold_label))
                [['mrk', 'scaffold', 'chrom', 'pos',
                  'BF(dB)', 'Beta_is', 'eBPis', 'C2', 'log10(1/pval)']])
    top100.to_csv(outdir / 'wild_baypass_top_hits.tsv', sep='\t', index=False)

    bf_thr = float(np.nanpercentile(bf, 95))
    df['top5pct'] = df['BF(dB)'] > bf_thr
    per_scaf = (df.groupby('chrom')
                  .agg(n_snps=('top5pct', 'size'),
                       n_top5pct=('top5pct', 'sum'))
                  .assign(frac_top5pct=lambda d: d['n_top5pct'] / d['n_snps'])
                  .sort_values('n_top5pct', ascending=False)
                  .head(15)
                  .reset_index())
    per_scaf['scaffold'] = per_scaf['chrom'].map(scaffold_label)
    per_scaf.to_csv(outdir / 'wild_baypass_by_chrom.tsv',
                    sep='\t', index=False)

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), dpi=200)
    ax = axes[0]
    ax.hist(bf, bins=200, color='#555555', edgecolor='none', alpha=0.8)
    ax.axvline(10, color='#D55E00', linestyle='--', linewidth=0.8,
               label=f'BF(dB)=10  n={row["n_bf_gt_10"]:,}')
    ax.axvline(20, color='#0072B2', linestyle=':', linewidth=0.8,
               label=f'BF(dB)=20  n={row["n_bf_gt_20"]:,}')
    ax.set_xlabel('BF(dB)  (β_i covariate model)')
    ax.set_ylabel('SNP count')
    ax.set_title(f'Wild BayPass BF(dB) distribution (n={len(df):,})', fontsize=9)
    ax.legend(fontsize=6)
    ax.set_xlim(np.nanpercentile(bf, 0.1), np.nanpercentile(bf, 99.99))

    ax = axes[1]
    n = len(df)
    idx = np.random.default_rng(0).choice(n, size=min(n, 50000), replace=False)
    ax.scatter(bf[idx], lp[idx], s=1, alpha=0.2, c='#333333',
               linewidths=0, rasterized=True)
    ax.axvline(10, color='#D55E00', linestyle='--', linewidth=0.5, alpha=0.6)
    ax.set_xlabel('BF(dB)   (β_i model)')
    ax.set_ylabel(r'$-\log_{10}(p)$   (C2 contrast)')
    ax.set_title(fr'β_i vs C2 (Spearman ρ={rho_bf_c2:.3f}, Pearson r={pearson_bf_c2:.3f})',
                 fontsize=9)

    fig.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(outdir / f'wild_baypass_diagnostics.{ext}',
                    bbox_inches='tight', dpi=200)
    plt.close(fig)


if __name__ == '__main__':
    main()
