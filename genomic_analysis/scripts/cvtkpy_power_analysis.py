#!/usr/bin/env python3

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'cvtkpy'))


def simulate_temporal_freqs(R, timepoints, L, Ne, N_pool, depth,
                            convergence_corr_true=0.0, seed=None):
    rng = np.random.default_rng(seed)
    ntp = len(timepoints)
    T = ntp - 1

    p0 = rng.uniform(0.15, 0.85, size=L)

    freqs_true = np.zeros((R, ntp, L))
    freqs_observed = np.zeros((R, ntp, L))

    for r in range(R):
        freqs_true[r, 0, :] = p0

    gen_intervals = [timepoints[i+1] - timepoints[i] for i in range(T)]

    for t_idx in range(T):
        dt = gen_intervals[t_idx]

        p_current = freqs_true[:, t_idx, :]  # (R, L)

        total_drift_var = np.mean(p_current * (1 - p_current), axis=0) * dt / (2 * Ne)
        if convergence_corr_true > 0:
            shared_var = convergence_corr_true * total_drift_var
            shared_delta = rng.normal(0, np.sqrt(np.maximum(shared_var, 0)), size=L)
        else:
            shared_delta = np.zeros(L)

        for r in range(R):
            p_r = freqs_true[r, t_idx, :]
            drift_var_r = p_r * (1 - p_r) * dt / (2 * Ne)

            indep_var = drift_var_r * (1 - convergence_corr_true)
            drift_delta = rng.normal(0, np.sqrt(np.maximum(indep_var, 0)), size=L)

            p_new = p_r + drift_delta + shared_delta
            p_new = np.clip(p_new, 0.001, 0.999)
            freqs_true[r, t_idx + 1, :] = p_new

    for r in range(R):
        for t in range(ntp):
            p_true = freqs_true[r, t, :]
            n_chrom = 2 * N_pool
            pool_counts = rng.binomial(n_chrom, p_true)
            p_pool = pool_counts / n_chrom
            d_per_site = rng.poisson(depth, size=L)
            d_per_site = np.maximum(d_per_site, 1)
            seq_counts = rng.binomial(d_per_site, p_pool)
            p_obs = seq_counts / d_per_site
            freqs_observed[r, t, :] = p_obs

    freqs_flat = freqs_observed.reshape(R * ntp, L)
    depths_flat = np.full_like(freqs_flat, depth, dtype=int)
    samples = []
    for r in range(R):
        for t in timepoints:
            samples.append((f"R{r+1}", t))

    return freqs_flat, depths_flat, samples


def compute_convergence_corr(freqs, depths, samples, N_pool, n_bootstrap=200):
    from cvtk.cvtk import TiledTemporalFreqs
    from cvtk.gintervals import GenomicIntervals

    L = freqs.shape[1]

    all_finite = np.all(np.isfinite(freqs), axis=0)
    mf = np.nanmean(freqs[:, all_finite], axis=0)
    maf_ok = (mf >= 0.10) & (mf <= 0.90)
    fr = freqs[:, all_finite][:, maf_ok]
    dp = depths[:, all_finite][:, maf_ok]
    L_filt = fr.shape[1]

    if L_filt < 1000:
        return np.nan, np.nan, np.nan

    from collections import OrderedDict
    n_chroms = 20
    loci_per_chrom = L_filt // n_chroms
    gi = GenomicIntervals()
    for i in range(L_filt):
        chrom_idx = i // loci_per_chrom
        if chrom_idx >= n_chroms:
            chrom_idx = n_chroms - 1
        pos = (i % loci_per_chrom) * 200
        gi.append(f'chr{chrom_idx+1}', pos, pos + 1)

    seqlens = OrderedDict()
    for c in range(1, n_chroms + 1):
        seqlens[f'chr{c}'] = loci_per_chrom * 200 + 100000
    tiles = GenomicIntervals(seqlens=seqlens)
    for c in range(1, n_chroms + 1):
        for start in range(0, loci_per_chrom * 200, 100000):
            tiles.append(f'chr{c}', start, start + 100000)

    ttf = TiledTemporalFreqs(tiles, fr, samples, depths=dp,
                              diploids=N_pool, gintervals=gi,
                              swap=True, share_first=False)

    cc = ttf.convergence_corr(bias_correction=True)
    cc_val = float(np.nanmean(cc))

    try:
        ci_lo, est, ci_hi = ttf.bootstrap_convergence_corr(
            B=n_bootstrap, alpha=0.05)
        ci_lo_val = float(np.nanmean(ci_lo))
        ci_hi_val = float(np.nanmean(ci_hi))
        return cc_val, min(ci_lo_val, ci_hi_val), max(ci_lo_val, ci_hi_val)
    except Exception:
        return cc_val, np.nan, np.nan


def main():
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--n-sims", type=int, default=50,
                   help="Simulations per parameter combination")
    p.add_argument("--n-loci", type=int, default=50000,
                   help="Number of loci to simulate (default 50k for speed)")
    p.add_argument("--n-bootstrap", type=int, default=200,
                   help="Bootstrap resamples per simulation")
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    R = 4
    timepoints = [0, 1, 2, 6, 7, 8, 9]
    Ne = 250       # approximate Ne from Session 1
    N_pool = 24    # effective pool size (de-duped depths, overdispersion-corrected)
    depth = 83     # de-duplicated depth

    cc_values = [0.0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20]
    n_sims = args.n_sims

    results = []
    for cc_true in cc_values:
        n_significant = 0
        cc_estimates = []
        ci_widths = []

        for sim in range(n_sims):

            freqs, depths_arr, samples = simulate_temporal_freqs(
                R, timepoints, args.n_loci, Ne, N_pool, depth,
                convergence_corr_true=cc_true, seed=sim * 1000 + int(cc_true * 10000))

            cc_est, ci_lo, ci_hi = compute_convergence_corr(
                freqs, depths_arr, samples, N_pool, args.n_bootstrap)

            cc_estimates.append(cc_est)
            if not np.isnan(ci_lo) and not np.isnan(ci_hi):
                ci_widths.append(ci_hi - ci_lo)
                if ci_lo > 0:
                    n_significant += 1

        power = n_significant / n_sims
        mean_est = np.nanmean(cc_estimates)
        mean_ci_width = np.nanmean(ci_widths) if ci_widths else np.nan

        results.append({
            'cc_true': cc_true,
            'power': power,
            'mean_estimate': mean_est,
            'mean_ci_width': mean_ci_width,
            'n_sims': n_sims,
        })

    import pandas as pd
    res_df = pd.DataFrame(results)
    res_df.to_csv(os.path.join(args.output_dir, 'power_analysis.tsv'),
                  sep='\t', index=False)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot([r['cc_true'] for r in results],
             [r['power'] for r in results],
             'ko-', markersize=8, linewidth=2)
    ax1.axhline(0.05, color='red', ls='--', alpha=0.5, label='α = 0.05')
    ax1.axhline(0.80, color='blue', ls='--', alpha=0.5, label='80% power')
    ax1.set_xlabel('True convergence correlation')
    ax1.set_ylabel('Power (fraction of sims with CI > 0)')
    ax1.set_title(f'Power to detect convergence > 0\n'
                  f'(R={R}, T={len(timepoints)}, N_pool={N_pool}, Ne={Ne}, '
                  f'L={args.n_loci})')
    ax1.legend()
    ax1.set_ylim(-0.05, 1.05)

    ax2.errorbar([r['cc_true'] for r in results],
                 [r['mean_estimate'] for r in results],
                 yerr=[r['mean_ci_width']/2 for r in results],
                 fmt='ko-', markersize=8, capsize=5)
    ax2.plot([0, max(cc_values)], [0, max(cc_values)], 'r--', alpha=0.5,
             label='1:1 line')
    ax2.set_xlabel('True convergence correlation')
    ax2.set_ylabel('Mean estimated convergence correlation')
    ax2.set_title('Estimation accuracy')
    ax2.legend()

    fig.tight_layout()
    fig.savefig(os.path.join(args.output_dir, 'power_curve.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)


    rep_results = []
    for R_test in [4, 6, 8, 10, 15, 20]:
        n_sig = 0
        for sim in range(n_sims):
            freqs, depths_arr, samples = simulate_temporal_freqs(
                R_test, timepoints, args.n_loci, Ne, N_pool, depth,
                convergence_corr_true=0.05, seed=sim * 100 + R_test)
            _, ci_lo, ci_hi = compute_convergence_corr(
                freqs, depths_arr, samples, N_pool, args.n_bootstrap)
            if not np.isnan(ci_lo) and ci_lo > 0:
                n_sig += 1
        power = n_sig / n_sims
        rep_results.append({'R': R_test, 'power': power})

    fig2, ax = plt.subplots(figsize=(6, 5))
    ax.plot([r['R'] for r in rep_results],
            [r['power'] for r in rep_results],
            'ko-', markersize=8, linewidth=2)
    ax.axhline(0.80, color='blue', ls='--', alpha=0.5, label='80% power')
    ax.set_xlabel('Number of replicates')
    ax.set_ylabel('Power')
    ax.set_title(f'Power vs replicates (cc_true=0.05, N_pool={N_pool}, Ne={Ne})')
    ax.legend()
    ax.set_ylim(-0.05, 1.05)
    fig2.tight_layout()
    fig2.savefig(os.path.join(args.output_dir, 'power_vs_replicates.png'),
                 dpi=150, bbox_inches='tight')
    plt.close(fig2)

    pd.DataFrame(rep_results).to_csv(
        os.path.join(args.output_dir, 'power_vs_replicates.tsv'),
        sep='\t', index=False)

if __name__ == "__main__":
    main()
