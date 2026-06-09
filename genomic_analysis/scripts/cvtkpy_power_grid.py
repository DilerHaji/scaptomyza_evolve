#!/usr/bin/env python3

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'cvtkpy'))
from cvtk.cvtk import TiledTemporalFreqs
from cvtk.gintervals import GenomicIntervals
from collections import OrderedDict


def simulate_and_test(R, timepoints, L, Ne, N_pool, depth, cc_true, seed):
    rng = np.random.default_rng(seed)
    ntp = len(timepoints)
    n_chroms = 20
    loci_per_chrom = L // n_chroms

    p0 = rng.uniform(0.15, 0.85, size=L)
    gen_intervals = [timepoints[i+1] - timepoints[i] for i in range(ntp - 1)]

    freqs_true = np.zeros((R, ntp, L))
    for r in range(R):
        freqs_true[r, 0, :] = p0

    for t_idx, dt in enumerate(gen_intervals):
        p_current = freqs_true[:, t_idx, :]
        total_drift_var = np.mean(p_current * (1 - p_current), axis=0) * dt / (2 * Ne)
        if cc_true > 0:
            shared_var = cc_true * total_drift_var
            shared_delta = rng.normal(0, np.sqrt(np.maximum(shared_var, 0)), size=L)
        else:
            shared_delta = np.zeros(L)

        for r in range(R):
            p_r = freqs_true[r, t_idx, :]
            drift_var_r = p_r * (1 - p_r) * dt / (2 * Ne)
            indep_var = drift_var_r * (1 - cc_true)
            drift_delta = rng.normal(0, np.sqrt(np.maximum(indep_var, 0)), size=L)
            freqs_true[r, t_idx + 1, :] = np.clip(p_r + drift_delta + shared_delta, 0.001, 0.999)

    freqs_obs = np.zeros((R, ntp, L))
    for r in range(R):
        for t in range(ntp):
            p = freqs_true[r, t, :]
            pool_counts = rng.binomial(2 * N_pool, p)
            p_pool = pool_counts / (2 * N_pool)
            d_site = np.maximum(rng.poisson(depth, size=L), 1)
            seq_counts = rng.binomial(d_site, p_pool)
            freqs_obs[r, t, :] = seq_counts / d_site

    freqs_flat = freqs_obs.reshape(R * ntp, L)
    depths_flat = np.full_like(freqs_flat, depth, dtype=int)
    samples = []
    for r in range(R):
        for t in timepoints:
            samples.append((f"R{r+1}", t))

    all_fin = np.all(np.isfinite(freqs_flat), axis=0)
    mf = np.nanmean(freqs_flat[:, all_fin], axis=0)
    maf_ok = (mf >= 0.10) & (mf <= 0.90)
    fr = freqs_flat[:, all_fin][:, maf_ok]
    dp = depths_flat[:, all_fin][:, maf_ok]
    L_filt = fr.shape[1]
    if L_filt < 1000:
        return None

    gi = GenomicIntervals()
    for i in range(L_filt):
        c_idx = min(i // (L_filt // n_chroms), n_chroms - 1)
        pos = (i % (L_filt // n_chroms)) * 200
        gi.append(f'chr{c_idx+1}', pos, pos + 1)

    seqlens = OrderedDict()
    for c in range(1, n_chroms + 1):
        seqlens[f'chr{c}'] = (L_filt // n_chroms) * 200 + 100000
    tiles = GenomicIntervals(seqlens=seqlens)
    for c in range(1, n_chroms + 1):
        for start in range(0, (L_filt // n_chroms) * 200, 100000):
            tiles.append(f'chr{c}', start, start + 100000)

    ttf = TiledTemporalFreqs(tiles, fr, samples, depths=dp,
                              diploids=N_pool, gintervals=gi,
                              swap=True, share_first=False)
    try:
        ci_lo, est, ci_hi = ttf.bootstrap_convergence_corr(B=200, alpha=0.05)
        ci_lo_val = float(np.nanmean(ci_lo))
        ci_hi_val = float(np.nanmean(ci_hi))
        lo, hi = min(ci_lo_val, ci_hi_val), max(ci_lo_val, ci_hi_val)
        return lo > 0
    except Exception:
        return None


def compute_power(R, timepoints, L, Ne, N_pool, depth, cc_true, n_sims):
    n_sig = 0
    n_valid = 0
    for sim in range(n_sims):
        result = simulate_and_test(R, timepoints, L, Ne, N_pool, depth,
                                    cc_true, seed=sim * 1000 + R * 100 + N_pool + Ne)
        if result is not None:
            n_valid += 1
            if result:
                n_sig += 1
    return n_sig / max(n_valid, 1)


def main():
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--n-sims", type=int, default=30)
    p.add_argument("--n-loci", type=int, default=200000)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    timepoints = [0, 1, 2, 6, 7, 8, 9]
    depth = 83
    cc_true = 0.05
    L = args.n_loci
    n_sims = args.n_sims

    rep_values = [4, 6, 8, 10, 15, 20]

    Ne_fixed = 250
    npool_values = [24, 40, 60, 80]

    npool_results = {}
    for np_val in npool_values:
        npool_results[np_val] = []
        for R in rep_values:
            pwr = compute_power(R, timepoints, L, Ne_fixed, np_val, depth,
                                cc_true, n_sims)
            npool_results[np_val].append(pwr)

    fig1, ax1 = plt.subplots(figsize=(8, 6))
    colors = ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4']
    for np_val, color in zip(npool_values, colors):
        label = f'N_pool={np_val}'
        if np_val == 24:
            label += ' (our data)'
        elif np_val == 80:
            label += ' (nominal)'
        ax1.plot(rep_values, npool_results[np_val], 'o-', color=color,
                 markersize=8, linewidth=2, label=label)
    ax1.axhline(0.80, color='gray', ls='--', alpha=0.5, label='80% power')
    ax1.set_xlabel('Number of replicates', fontsize=12)
    ax1.set_ylabel('Power to detect cc=0.05', fontsize=12)
    ax1.set_title(f'Effect of effective pool size (Ne={Ne_fixed}, depth={depth})',
                  fontsize=13)
    ax1.legend(fontsize=10)
    ax1.set_ylim(-0.05, 1.05)
    ax1.set_xticks(rep_values)
    fig1.tight_layout()
    fig1.savefig(os.path.join(args.output_dir, 'power_vs_npool.png'),
                 dpi=150, bbox_inches='tight')
    plt.close(fig1)

    npool_fixed = 24
    ne_values = [100, 250, 500]

    ne_results = {}
    for ne_val in ne_values:
        ne_results[ne_val] = []
        for R in rep_values:
            pwr = compute_power(R, timepoints, L, ne_val, npool_fixed, depth,
                                cc_true, n_sims)
            ne_results[ne_val].append(pwr)

    fig2, ax2 = plt.subplots(figsize=(8, 6))
    colors2 = ['#9467bd', '#8c564b', '#e377c2']
    for ne_val, color in zip(ne_values, colors2):
        label = f'Ne={ne_val}'
        if ne_val == 250:
            label += ' (our estimate)'
        elif ne_val == 500:
            label += ' (census)'
        ax2.plot(rep_values, ne_results[ne_val], 'o-', color=color,
                 markersize=8, linewidth=2, label=label)
    ax2.axhline(0.80, color='gray', ls='--', alpha=0.5, label='80% power')
    ax2.set_xlabel('Number of replicates', fontsize=12)
    ax2.set_ylabel('Power to detect cc=0.05', fontsize=12)
    ax2.set_title(f'Effect of population Ne (N_pool={npool_fixed}, depth={depth})',
                  fontsize=13)
    ax2.legend(fontsize=10)
    ax2.set_ylim(-0.05, 1.05)
    ax2.set_xticks(rep_values)
    fig2.tight_layout()
    fig2.savefig(os.path.join(args.output_dir, 'power_vs_ne.png'),
                 dpi=150, bbox_inches='tight')
    plt.close(fig2)

    fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    for np_val, color in zip(npool_values, colors):
        label = f'N_pool={np_val}'
        if np_val == 24:
            label += ' (ours)'
        elif np_val == 80:
            label += ' (nominal)'
        ax3a.plot(rep_values, npool_results[np_val], 'o-', color=color,
                  markersize=7, linewidth=2, label=label)
    ax3a.axhline(0.80, color='gray', ls='--', alpha=0.4)
    ax3a.set_xlabel('Number of replicates', fontsize=12)
    ax3a.set_ylabel('Power to detect cc=0.05', fontsize=12)
    ax3a.set_title(f'Varying pool size (Ne={Ne_fixed})', fontsize=13)
    ax3a.legend(fontsize=9)
    ax3a.set_ylim(-0.05, 1.05)
    ax3a.set_xticks(rep_values)

    for ne_val, color in zip(ne_values, colors2):
        label = f'Ne={ne_val}'
        if ne_val == 250:
            label += ' (ours)'
        elif ne_val == 500:
            label += ' (census)'
        ax3b.plot(rep_values, ne_results[ne_val], 'o-', color=color,
                  markersize=7, linewidth=2, label=label)
    ax3b.axhline(0.80, color='gray', ls='--', alpha=0.4)
    ax3b.set_xlabel('Number of replicates', fontsize=12)
    ax3b.set_title(f'Varying population Ne (N_pool={npool_fixed})', fontsize=13)
    ax3b.legend(fontsize=9)
    ax3b.set_xticks(rep_values)

    fig3.suptitle(f'Power to detect convergence correlation = {cc_true}\n'
                  f'(timepoints={timepoints}, depth={depth}, {L//1000}k loci)',
                  fontsize=14, y=1.03)
    fig3.tight_layout()
    fig3.savefig(os.path.join(args.output_dir, 'power_combined.png'),
                 dpi=150, bbox_inches='tight')
    plt.close(fig3)

    import pandas as pd
    rows = []
    for np_val in npool_values:
        for i, R in enumerate(rep_values):
            rows.append({'Ne': Ne_fixed, 'N_pool': np_val, 'R': R,
                         'power': npool_results[np_val][i], 'cc_true': cc_true})
    for ne_val in ne_values:
        for i, R in enumerate(rep_values):
            rows.append({'Ne': ne_val, 'N_pool': npool_fixed, 'R': R,
                         'power': ne_results[ne_val][i], 'cc_true': cc_true})
    pd.DataFrame(rows).to_csv(os.path.join(args.output_dir, 'power_grid.tsv'),
                               sep='\t', index=False)

if __name__ == "__main__":
    main()
