#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

OUTDIR = 'variance_analysis/cvtkpy_final'


def fig_bt_distribution():
    df = pd.read_csv(f'{OUTDIR}/window_BT_correlations.tsv', sep='\t')

    fig, ax = plt.subplots(figsize=(8, 5))
    sig_neg = df[df['p_BT_neg'] < 0.05]['r_BT']
    sig_pos = df[(1 - df['p_BT_neg']) < 0.05]['r_BT']
    nonsig = df[(df['p_BT_neg'] >= 0.05) & ((1 - df['p_BT_neg']) >= 0.05)]['r_BT']

    bins = np.linspace(-0.6, 0.8, 60)
    ax.hist(nonsig, bins=bins, color='#bdbdbd', alpha=0.8, label='Not significant')
    ax.hist(sig_pos, bins=bins, color='#4393c3', alpha=0.8,
            label=f'Sig. positive ({len(sig_pos)}, {100*len(sig_pos)/len(df):.0f}%)')
    ax.hist(sig_neg, bins=bins, color='#d6604d', alpha=0.8,
            label=f'Sig. negative ({len(sig_neg)}, {100*len(sig_neg)/len(df):.0f}%)')

    ax.axvline(0, color='black', ls='--', linewidth=1)
    ax.axvline(df['r_BT'].mean(), color='black', ls='-', linewidth=1.5,
               label=f'Mean = {df["r_BT"].mean():.3f}')
    ax.set_xlabel('B-T correlation (r) per 100kb window', fontsize=12)
    ax.set_ylabel('Number of windows', fontsize=12)
    ax.set_title('Genome-wide distribution of B-T allele frequency concordance',
                 fontsize=13)
    ax.legend(fontsize=10)
    fig.tight_layout()
    fig.savefig(f'{OUTDIR}/fig_supp_BT_distribution.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

def fig_power_combined():
    for f in ['power_curve.png', 'power_combined.png']:
        path = f'variance_analysis/cvtkpy_power/{f}'

def fig_summary_composite():
    fig = plt.figure(figsize=(16, 12))

    ax_g = fig.add_subplot(2, 2, 1)
    img = plt.imread(f'{OUTDIR}/fig_G_k2_by_replicate.png')
    ax_g.imshow(img)
    ax_g.axis('off')
    ax_g.set_title('A', fontsize=16, fontweight='bold', loc='left', x=-0.02)

    ax_cc = fig.add_subplot(2, 2, 2)
    img2 = plt.imread(f'{OUTDIR}/fig_convergence_comparison.png')
    ax_cc.imshow(img2)
    ax_cc.axis('off')
    ax_cc.set_title('B', fontsize=16, fontweight='bold', loc='left', x=-0.02)

    ax_bt = fig.add_subplot(2, 2, 3)
    df = pd.read_csv(f'{OUTDIR}/window_BT_correlations.tsv', sep='\t')
    sig_neg = df[df['p_BT_neg'] < 0.05]['r_BT']
    sig_pos = df[(1 - df['p_BT_neg']) < 0.05]['r_BT']
    nonsig = df[(df['p_BT_neg'] >= 0.05) & ((1 - df['p_BT_neg']) >= 0.05)]['r_BT']
    bins = np.linspace(-0.6, 0.8, 60)
    ax_bt.hist(nonsig, bins=bins, color='#bdbdbd', alpha=0.8, label='ns')
    ax_bt.hist(sig_pos, bins=bins, color='#4393c3', alpha=0.8,
               label=f'Concordant ({len(sig_pos)})')
    ax_bt.hist(sig_neg, bins=bins, color='#d6604d', alpha=0.8,
               label=f'Antagonistic ({len(sig_neg)})')
    ax_bt.axvline(0, color='black', ls='--', linewidth=1)
    ax_bt.axvline(df['r_BT'].mean(), color='black', ls='-', linewidth=1.5)
    ax_bt.set_xlabel('B-T correlation per 100kb window')
    ax_bt.set_ylabel('Windows')
    ax_bt.legend(fontsize=8)
    ax_bt.set_title('C', fontsize=16, fontweight='bold', loc='left', x=-0.02)

    ax_an = fig.add_subplot(2, 2, 4)
    img3 = plt.imread(f'{OUTDIR}/fig_antagonistic_scaffolds.png')
    ax_an.imshow(img3)
    ax_an.axis('off')
    ax_an.set_title('D', fontsize=16, fontweight='bold', loc='left', x=-0.02)

    fig.tight_layout()
    fig.savefig(f'{OUTDIR}/fig_supp_composite.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':

    fig_bt_distribution()
    fig_power_combined()
    fig_summary_composite()