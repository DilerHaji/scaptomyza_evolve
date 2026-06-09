#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys, os, re
from collections import OrderedDict
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'cvtkpy'))
from cvtk.cov import temporal_replicate_cov, stack_temporal_covariances, total_variance
from cvtk.utils import sort_samples, process_samples, reshape_matrix, validate_diploids

OUTDIR = 'variance_analysis/cvtkpy_final'
POWDIR = 'variance_analysis/cvtkpy_power'

# Project colors
COL_B = '#4EA2FF'
COL_T = '#EDB72F'
COL_M = '#9BAB96'
COL_ANTAG = '#D76161'
COL_CONCORDANT = '#4EA2FF'  # concordant uses B color since it's shared
COL_NS = '#bdbdbd'
COL_MEAN = '#333333'

samples_list = None
freqs_all = None
depths_all = None
chroms_arr = None
pos_arr = None


def load_all():
    global samples_list, freqs_all, depths_all, chroms_arr, pos_arr
    if freqs_all is not None:
        return
    samples_list = [l.strip() for l in open('variance_analysis/sample_list.txt') if l.strip()]
    freqs_all, depths_all, chroms_arr, pos_arr = [], [], [], []
    with open('variance_analysis/merged_ad.tsv') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            chroms_arr.append(parts[0])
            pos_arr.append(int(parts[1]))
            tots, alts = [], []
            for ad in parts[4:]:
                if ad in ('.', '.,.'):
                    tots.append(0); alts.append(0)
                else:
                    r, a = ad.split(',')[:2]
                    r, a = int(r), int(a)
                    alts.append(a); tots.append(r + a)
            tot = np.array(tots)
            freqs_all.append(np.where(tot > 0, np.array(alts) / tot, np.nan))
            depths_all.append(tot)
    freqs_all = np.array(freqs_all).T
    depths_all = np.array(depths_all).T
    chroms_arr = np.array(chroms_arr)
    pos_arr = np.array(pos_arr)


def savefig_svg(fig, name):
    path = os.path.join(OUTDIR, name)
    fig.savefig(path, format='svg', bbox_inches='tight')
    with open(path, 'r') as f:
        svg = f.read()
    svg = re.sub(r'<clipPath\b[^>]*>.*?</clipPath>', '', svg, flags=re.DOTALL)
    svg = re.sub(r'<mask\b[^>]*>.*?</mask>', '', svg, flags=re.DOTALL)
    svg = re.sub(r'\s*clip-path="[^"]*"', '', svg)
    svg = re.sub(r'\s*mask="[^"]*"', '', svg)
    with open(path, 'w') as f:
        f.write(svg)


def get_treatment_data(trt):
    full_gens = [1, 2, 6, 7, 8, 9]
    indices, cvtk_samples = [], []
    for rep in range(1, 5):
        indices.append(samples_list.index(f'F{rep}G00'))
        cvtk_samples.append((f'{trt}{rep}', 0))
        for gen in full_gens:
            indices.append(samples_list.index(f'{trt}{rep}G{gen:02d}'))
            cvtk_samples.append((f'{trt}{rep}', gen))
    fr = freqs_all[indices, :]
    dp = depths_all[indices, :]
    ok = np.all(dp >= 20, axis=0) & np.all(np.isfinite(fr), axis=0)
    mf = np.nanmean(fr[:, ok], axis=0)
    ok2 = (mf >= 0.05) & (mf <= 0.95)
    return fr[:, ok][:, ok2], dp[:, ok][:, ok2], cvtk_samples


def get_scaffold_order_and_offsets(df):
    scaf_sizes = df.groupby('chrom')['end'].max().sort_values(ascending=False)
    offsets = {}
    cum = 0
    ordered_scaffolds = []
    for scaf, size in scaf_sizes.items():
        offsets[scaf] = cum
        ordered_scaffolds.append(scaf)
        cum += size + 2e6  # 2Mb gap between scaffolds
    return ordered_scaffolds, offsets, cum


def fig_G_k2():
    load_all()
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
    timepoint_labels = list(range(1, 7))
    gen_labels = ['G01', 'G02', 'G06', 'G07', 'G08', 'G09']
    trt_colors = {'B': COL_B, 'T': COL_T, 'M': COL_M}

    for ax, trt in zip(axes, ['B', 'T', 'M']):
        fr_f, dp_f, cvtk_samples = get_treatment_data(trt)
        sorted_samples, sorted_i = sort_samples(cvtk_samples)
        _, _, R, ntp = process_samples(fr_f, sorted_samples)
        T = ntp - 1
        fr_3d = reshape_matrix(fr_f[sorted_i, :], R)
        dp_3d = reshape_matrix(dp_f[sorted_i, :], R)
        dip = validate_diploids(np.full(28, 80)[sorted_i], R, ntp)
        cov = temporal_replicate_cov(fr_3d, dp_3d, dip, bias_correction=True,
                                      standardize=False, share_first=False)
        temp_covs = stack_temporal_covariances(cov, R, T)
        tv = np.stack([total_variance(fr_3d, dp_3d, dip, t=t, standardize=False,
                                       bias_correction=True) for t in range(1, T+1)])
        G_k2 = np.zeros((T, R))
        for t_end in range(1, T + 1):
            for r in range(R):
                tc = temp_covs[:t_end, :t_end, r]
                offdiag_k2 = np.tril(tc, -2) + np.triu(tc, 2)
                G_k2[t_end-1, r] = np.nansum(offdiag_k2) / tv[t_end-1, r]

        col = trt_colors[trt]
        for r in range(R):
            ax.plot(timepoint_labels, G_k2[:, r], 'o-', color=col,
                    alpha=0.35, markersize=5, linewidth=1.2,
                    label=f'Rep {r+1}' if trt == 'B' else None)
        G_k2_mean = np.nanmean(G_k2, axis=1)
        ax.plot(timepoint_labels, G_k2_mean, 's-', color=COL_MEAN,
                linewidth=2.5, markersize=8, label='Mean' if trt == 'B' else None, zorder=5)
        ax.axhline(0, color='gray', ls='--', alpha=0.5, linewidth=1)
        ax.set_xlabel('Cumulative timepoint', fontsize=11)
        ax.set_title(f'{trt} treatment', fontsize=13, color=col)
        ax.set_xticks(timepoint_labels)
        ax.set_xticklabels(gen_labels, fontsize=9)
        if ax == axes[0]:
            ax.set_ylabel('G(t)  [k ≥ 2 off-diagonal only]', fontsize=11)
            ax.legend(fontsize=9)

    fig.suptitle('Cumulative G statistic excluding adjacent-interval covariance (k ≥ 2)',
                 fontsize=14, y=1.02)
    fig.tight_layout()
    savefig_svg(fig, 'fig_G_k2_by_replicate.svg')
    plt.close(fig)


def fig_convergence():
    summary = pd.read_csv(f'{OUTDIR}/summary.tsv', sep='\t')
    TREATMENTS = ['B', 'T', 'M']
    trt_colors = [COL_B, COL_T, COL_M]

    fig, ax = plt.subplots(figsize=(7, 5))
    x_positions = np.arange(3)
    x_offsets = [-0.12, 0.12]
    markers = ['o', 's']
    marker_sizes = [12, 10]
    labels_corr = ['N = 80 (nominal)', 'N_eff = 29 (founder-validated)']

    for j, (corr, marker, ms, label_c) in enumerate(
            zip(['N=80', 'N_eff=29'], markers, marker_sizes, labels_corr)):
        sub = summary[summary['correction'] == corr]
        for i, (trt, col) in enumerate(zip(TREATMENTS, trt_colors)):
            row = sub[sub['treatment'] == trt].iloc[0]
            cc = float(row['cc_estimate'])
            lo = float(row['cc_ci_lower'])
            hi = float(row['cc_ci_upper'])
            x = x_positions[i] + x_offsets[j]

            ax.plot([x, x], [lo, hi], color=col, linewidth=2, zorder=4)
            cap_w = 0.04
            ax.plot([x - cap_w, x + cap_w], [lo, lo], color=col, linewidth=1.5, zorder=4)
            ax.plot([x - cap_w, x + cap_w], [hi, hi], color=col, linewidth=1.5, zorder=4)

            if j == 0:  # open markers for N=80
                ax.plot(x, cc, marker=marker, markersize=ms,
                        markerfacecolor='white', markeredgecolor=col,
                        markeredgewidth=2, zorder=5, linestyle='none',
                        label=label_c if i == 0 else None)
            else:  # filled for N_eff
                ax.plot(x, cc, marker=marker, markersize=ms,
                        markerfacecolor=col, markeredgecolor='black',
                        markeredgewidth=0.5, zorder=5, linestyle='none',
                        label=label_c if i == 0 else None)

    ax.axhline(0, color='gray', ls='--', alpha=0.5)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(TREATMENTS, fontsize=13)
    ax.set_xlabel('Treatment', fontsize=12)
    ax.set_ylabel('Convergence correlation', fontsize=12)
    ax.set_title('Convergence correlation with 95% block bootstrap CIs', fontsize=13)
    ax.legend(fontsize=9, loc='upper left')
    fig.tight_layout()
    savefig_svg(fig, 'fig_convergence_comparison.svg')
    plt.close(fig)


def fig_bt_distribution():
    df = pd.read_csv(f'{OUTDIR}/window_BT_correlations.tsv', sep='\t')
    fig, ax = plt.subplots(figsize=(8, 5))
    sig_neg = df[df['p_BT_neg'] < 0.05]['r_BT']
    sig_pos = df[(1 - df['p_BT_neg']) < 0.05]['r_BT']
    nonsig = df[(df['p_BT_neg'] >= 0.05) & ((1 - df['p_BT_neg']) >= 0.05)]['r_BT']
    bins = np.linspace(-0.6, 0.8, 60)
    ax.hist(nonsig, bins=bins, color=COL_NS, alpha=0.8, label='Not significant')
    ax.hist(sig_pos, bins=bins, color=COL_CONCORDANT, alpha=0.7,
            label=f'Concordant ({len(sig_pos)}, {100*len(sig_pos)/len(df):.0f}%)')
    ax.hist(sig_neg, bins=bins, color=COL_ANTAG, alpha=0.8,
            label=f'Antagonistic ({len(sig_neg)}, {100*len(sig_neg)/len(df):.0f}%)')
    ax.axvline(0, color='black', ls='--', linewidth=1)
    ax.axvline(df['r_BT'].mean(), color=COL_MEAN, ls='-', linewidth=1.5,
               label=f'Mean = {df["r_BT"].mean():.3f}')
    ax.set_xlabel('B–T correlation (r) per 100 kb window', fontsize=12)
    ax.set_ylabel('Number of windows', fontsize=12)
    ax.set_title('Genome-wide distribution of B–T allele frequency concordance', fontsize=13)
    ax.legend(fontsize=10)
    fig.tight_layout()
    savefig_svg(fig, 'fig_BT_distribution.svg')
    plt.close(fig)



def fig_manhattan():
    df = pd.read_csv(f'{OUTDIR}/window_BT_correlations.tsv', sep='\t')

    scaf_sizes = df.groupby('chrom')['end'].max().sort_values(ascending=False)
    big_scaffolds = scaf_sizes[scaf_sizes >= 500000].index.tolist()
    df = df[df['chrom'].isin(big_scaffolds)].copy()

    ordered = scaf_sizes[scaf_sizes >= 500000].index.tolist()
    offsets = {}
    cum = 0
    for scaf in ordered:
        offsets[scaf] = cum
        cum += float(scaf_sizes[scaf]) + 2e6  # 2Mb gap

    df['genome_pos'] = df.apply(
        lambda row: offsets[row['chrom']] + (row['start'] + row['end']) / 2, axis=1)
    df = df.sort_values('genome_pos')

    fig, ax = plt.subplots(figsize=(18, 5))

    for i, scaf in enumerate(ordered):
        x_min = offsets[scaf] / 1e6
        x_max = (offsets[scaf] + float(scaf_sizes[scaf])) / 1e6
        if i % 2 == 0:
            ax.axvspan(x_min, x_max, alpha=0.06, color='#cccccc', zorder=0)

    sig_neg = df[df['p_BT_neg'] < 0.05]
    sig_pos = df[(1 - df['p_BT_neg']) < 0.05]
    nonsig = df[(df['p_BT_neg'] >= 0.05) & ((1 - df['p_BT_neg']) >= 0.05)]

    ax.scatter(nonsig['genome_pos'] / 1e6, nonsig['r_BT'],
               s=6, color=COL_NS, alpha=0.4, zorder=1,
               label=f'Not significant ({len(nonsig)})')
    ax.scatter(sig_pos['genome_pos'] / 1e6, sig_pos['r_BT'],
               s=8, color=COL_CONCORDANT, alpha=0.5, zorder=2,
               label=f'Concordant ({len(sig_pos)}, {100*len(sig_pos)/len(df):.0f}%)')
    ax.scatter(sig_neg['genome_pos'] / 1e6, sig_neg['r_BT'],
               s=18, color=COL_ANTAG, alpha=0.9, zorder=3,
               edgecolors='darkred', linewidths=0.5,
               label=f'Antagonistic ({len(sig_neg)}, {100*len(sig_neg)/len(df):.0f}%)')

    ax.axhline(0, color='black', ls='-', linewidth=0.5, alpha=0.3)
    ax.set_xlabel('Genomic position (Mb, scaffolds ordered by size)', fontsize=12)
    ax.set_ylabel('B–T correlation (r)', fontsize=12)
    ax.set_title('Genome-wide B–T allele frequency correlation (100 kb windows)', fontsize=13)
    ax.legend(fontsize=10, loc='upper right', framealpha=0.9)

    for i, scaf in enumerate(ordered):
        scaf_mid = (offsets[scaf] + float(scaf_sizes[scaf]) / 2) / 1e6
        short = scaf.split('HRSCAF_')[-1] if 'HRSCAF' in scaf else scaf.split('_')[-1]
        if float(scaf_sizes[scaf]) > 5e6:  # only label scaffolds > 5Mb
            ax.text(scaf_mid, ax.get_ylim()[0] + 0.02, short,
                    ha='center', fontsize=8, alpha=0.5, va='bottom')

    fig.tight_layout()
    savefig_svg(fig, 'fig_manhattan_BT.svg')
    plt.close(fig)

def fig_power():
    power_data = pd.read_csv(f'{POWDIR}/power_grid.tsv', sep='\t')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    # Left: varying N_pool
    npool_colors = {24: COL_ANTAG, 40: COL_T, 60: COL_M, 80: COL_B}
    for np_val in [24, 40, 60, 80]:
        sub = power_data[(power_data['N_pool'] == np_val) & (power_data['Ne'] == 250)]
        if len(sub) > 0:
            label = f'N_pool = {np_val}'
            if np_val == 24: label += ' (ours)'
            elif np_val == 80: label += ' (nominal)'
            ax1.plot(sub['R'], sub['power'], 'o-', color=npool_colors[np_val],
                     markersize=8, linewidth=2, label=label)
    ax1.axhline(0.80, color='gray', ls='--', alpha=0.4, label='80% power')
    ax1.set_xlabel('Number of replicates', fontsize=12)
    ax1.set_ylabel('Power to detect cc = 0.05', fontsize=12)
    ax1.set_title('Varying effective pool size (Ne = 250)', fontsize=13)
    ax1.legend(fontsize=9)
    ax1.set_ylim(-0.05, 1.05)

    ne_colors = {100: COL_B, 250: COL_M, 500: COL_T}
    for ne_val in [100, 250, 500]:
        sub = power_data[(power_data['Ne'] == ne_val) & (power_data['N_pool'] == 24)]
        if len(sub) > 0:
            label = f'Ne = {ne_val}'
            if ne_val == 250: label += ' (ours)'
            elif ne_val == 500: label += ' (census)'
            ax2.plot(sub['R'], sub['power'], 'o-', color=ne_colors[ne_val],
                     markersize=8, linewidth=2, label=label)
    ax2.axhline(0.80, color='gray', ls='--', alpha=0.4)
    ax2.set_xlabel('Number of replicates', fontsize=12)
    ax2.set_title('Varying population Ne (N_pool = 24)', fontsize=13)
    ax2.legend(fontsize=9)

    fig.suptitle('Power to detect convergence correlation = 0.05', fontsize=14, y=1.02)
    fig.tight_layout()
    savefig_svg(fig, 'fig_power_combined.svg')
    plt.close(fig)


if __name__ == '__main__':


    fig_G_k2()
    fig_convergence()
    fig_bt_distribution()
    fig_manhattan()
    fig_power()