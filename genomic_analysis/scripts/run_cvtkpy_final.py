#!/usr/bin/env python3
import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import OrderedDict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'cvtkpy'))
from cvtk.cvtk import TiledTemporalFreqs
from cvtk.gintervals import GenomicIntervals
from cvtk.cov import stack_temporal_covariances

MIN_DEPTH = 20
MIN_MAF = 0.05  
TILE_SIZE = 100_000 
N_BOOTSTRAP = 1000
TREATMENTS = ["B", "T", "M"]
INTERVAL_LABELS = ["G00→G01", "G01→G02", "G02→G06", "G06→G07", "G07→G08", "G08→G09"]


def load_data(ad_path, sample_list_path):
    samples = [l.strip() for l in open(sample_list_path) if l.strip()]
    chroms, positions, freqs_list, depths_list = [], [], [], []
    with open(ad_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chroms.append(parts[0])
            positions.append(int(parts[1]))
            tots, alts = [], []
            for ad in parts[4:]:
                if ad in (".", ".,."):
                    tots.append(0); alts.append(0)
                else:
                    r, a = ad.split(",")[:2]
                    r, a = int(r), int(a)
                    alts.append(a); tots.append(r + a)
            tot = np.array(tots)
            freqs_list.append(np.where(tot > 0, np.array(alts) / tot, np.nan))
            depths_list.append(tot)
    freqs = np.array(freqs_list).T
    depths = np.array(depths_list).T
    return samples, np.array(chroms), np.array(positions), freqs, depths


def get_treatment_samples(all_samples, treatment):
    full_gens = [1, 2, 6, 7, 8, 9]
    indices, cvtk_samples = [], []
    for rep in range(1, 5):
        fname = f"F{rep}G00"
        if fname in all_samples:
            indices.append(all_samples.index(fname))
            cvtk_samples.append((f"{treatment}{rep}", 0))
        for gen in full_gens:
            sname = f"{treatment}{rep}G{gen:02d}"
            if sname in all_samples:
                indices.append(all_samples.index(sname))
                cvtk_samples.append((f"{treatment}{rep}", gen))
    return indices, cvtk_samples


def make_tiles(chroms, positions, tile_size):
    chrom_bounds = OrderedDict()
    for c, p in zip(chroms, positions):
        if c not in chrom_bounds:
            chrom_bounds[c] = [p, p]
        else:
            chrom_bounds[c][0] = min(chrom_bounds[c][0], p)
            chrom_bounds[c][1] = max(chrom_bounds[c][1], p)
    seqlens = OrderedDict()
    for c, (mn, mx) in chrom_bounds.items():
        seqlens[c] = mx + tile_size
    tiles = GenomicIntervals(seqlens=seqlens)
    for chrom, (min_p, max_p) in chrom_bounds.items():
        for ws in range(int(min_p) - int(min_p) % tile_size,
                        int(max_p) + 1, tile_size):
            tiles.append(chrom, ws, ws + tile_size)
    return tiles


def run_one(treatment, all_samples, chroms, positions, freqs, depths,
            output_dir, diploids_val, label):
    indices, cvtk_samples = get_treatment_samples(all_samples, treatment)
    freq_sub = freqs[indices, :]
    depth_sub = depths[indices, :]

    ok = (np.all(depth_sub >= MIN_DEPTH, axis=0) &
          np.all(np.isfinite(freq_sub), axis=0))
    mf = np.nanmean(freq_sub[:, ok], axis=0)
    ok2 = (mf >= MIN_MAF) & (mf <= 1 - MIN_MAF)
    fr = freq_sub[:, ok][:, ok2]
    dp = depth_sub[:, ok][:, ok2]
    ch = chroms[ok][ok2]
    po = positions[ok][ok2]

    gi = GenomicIntervals()
    for c, p in zip(ch, po):
        gi.append(c, p, p + 1)
    tiles = make_tiles(ch, po, TILE_SIZE)

    ttf = TiledTemporalFreqs(tiles, fr, cvtk_samples, depths=dp,
                              diploids=diploids_val, gintervals=gi,
                              swap=True, share_first=False)


    G = ttf.calc_G()  # shape (T, R)

    cc = ttf.convergence_corr(bias_correction=True)  # shape (1, T, T) or similar
    cc_scalar = float(np.nanmean(cc))

    gw_cov = ttf.calc_cov(bias_correction=True, standardize=True)

    gw_cov_unstd = ttf.calc_cov(bias_correction=True, standardize=False)
    temp_covs = stack_temporal_covariances(gw_cov_unstd, ttf.R, ttf.T)
    avg_temp = np.mean(temp_covs, axis=2)
    k1_residual = np.mean(np.diag(avg_temp, 1))
    k2_mean = np.mean([np.mean(np.diag(avg_temp, k)) for k in range(2, ttf.T)])

    cc_ci_lo, cc_est, cc_ci_hi = ttf.bootstrap_convergence_corr(
        B=N_BOOTSTRAP, alpha=0.05, progress_bar=True)
    cc_lo = float(np.nanmean(cc_ci_lo))
    cc_hi = float(np.nanmean(cc_ci_hi))

    G_ci_lo, G_est, G_ci_hi = ttf.bootstrap_G(
        B=N_BOOTSTRAP, alpha=0.05, average_replicates=True, progress_bar=True)

    return {
        'treatment': treatment,
        'diploids': diploids_val,
        'label': label,
        'n_loci': fr.shape[1],
        'n_tiles': ttf.ntiles,
        'R': ttf.R,
        'T': ttf.T,
        'G_matrix': G,
        'G_final_per_rep': G[-1, :] if G.ndim == 2 else G,
        'G_final_mean': float(np.nanmean(G[-1, :])) if G.ndim == 2 else float(G),
        'G_bs_trajectory': G_est,
        'G_bs_ci_lo': G_ci_lo,
        'G_bs_ci_hi': G_ci_hi,
        'cc_matrix': np.array(cc).squeeze(),
        'cc_scalar': cc_scalar,
        'cc_ci_lo': min(cc_lo, cc_hi),
        'cc_ci_hi': max(cc_lo, cc_hi),
        'gw_cov': gw_cov,
        'k1_residual': k1_residual,
        'k2_mean': k2_mean,
        'avg_temp_cov': avg_temp,
    }


def plot_G_trajectories(results_n80, results_neff, output_dir):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    timepoint_labels = list(range(1, 7))

    for ax, trt in zip(axes, TREATMENTS):
        r80 = results_n80[trt]
        rnf = results_neff[trt]

        G = r80['G_matrix']
        for rep in range(G.shape[1]):
            ax.plot(timepoint_labels, G[:, rep], 'o-', alpha=0.3, color='C0',
                    markersize=4, linewidth=1)

        G_bs = r80['G_bs_trajectory']
        G_lo = r80['G_bs_ci_lo']
        G_hi = r80['G_bs_ci_hi']
        ax.plot(timepoint_labels, G_bs, 'o-', color='C0', linewidth=2,
                markersize=6, label=f'N=80', zorder=5)
        ax.fill_between(timepoint_labels,
                         np.minimum(G_lo, G_hi),
                         np.maximum(G_lo, G_hi),
                         alpha=0.2, color='C0')

        ax.axhline(0, color='gray', ls='--', alpha=0.5)
        ax.set_xlabel('Cumulative timepoint')
        ax.set_title(f'{trt} treatment')
        ax.set_xticks(timepoint_labels)
        ax.set_xticklabels([f't={t}' for t in timepoint_labels])
        if ax == axes[0]:
            ax.set_ylabel('G(t)')
            ax.legend(fontsize=9)

    fig.suptitle('Cumulative G statistic (fraction of AF change from selection)',
                 fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, 'fig_G_trajectories.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_convergence_summary(results_n80, results_neff, output_dir):
    fig, ax = plt.subplots(figsize=(8, 5))

    x_positions = np.arange(len(TREATMENTS))
    width = 0.35

    # N=80 bars
    cc_80 = [results_n80[t]['cc_scalar'] for t in TREATMENTS]
    ci_lo_80 = [results_n80[t]['cc_ci_lo'] for t in TREATMENTS]
    ci_hi_80 = [results_n80[t]['cc_ci_hi'] for t in TREATMENTS]
    err_lo_80 = [cc - lo for cc, lo in zip(cc_80, ci_lo_80)]
    err_hi_80 = [hi - cc for cc, hi in zip(cc_80, ci_hi_80)]
    ax.bar(x_positions - width/2, cc_80, width, label='N=80 (nominal)',
           color='C0', alpha=0.7)
    ax.errorbar(x_positions - width/2, cc_80,
                yerr=[err_lo_80, err_hi_80],
                fmt='none', color='black', capsize=5)

    # N_eff bars
    cc_nf = [results_neff[t]['cc_scalar'] for t in TREATMENTS]
    ci_lo_nf = [results_neff[t]['cc_ci_lo'] for t in TREATMENTS]
    ci_hi_nf = [results_neff[t]['cc_ci_hi'] for t in TREATMENTS]
    err_lo_nf = [cc - lo for cc, lo in zip(cc_nf, ci_lo_nf)]
    err_hi_nf = [hi - cc for cc, hi in zip(cc_nf, ci_hi_nf)]
    ax.bar(x_positions + width/2, cc_nf, width,
           label='N_eff=29 (founder-validated)',
           color='C1', alpha=0.7)
    ax.errorbar(x_positions + width/2, cc_nf,
                yerr=[err_lo_nf, err_hi_nf],
                fmt='none', color='black', capsize=5)

    ax.axhline(0, color='gray', ls='--', alpha=0.5)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(TREATMENTS)
    ax.set_xlabel('Treatment')
    ax.set_ylabel('Convergence correlation')
    ax.set_title('Convergence correlation with 95% block bootstrap CIs')
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, 'fig_convergence_comparison.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_temporal_cov_heatmap(results, output_dir, label_suffix):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, trt in zip(axes, TREATMENTS):
        r = results[trt]
        cov = r['gw_cov']
        im = ax.imshow(cov, cmap='RdBu_r', aspect='auto',
                       vmin=-np.percentile(np.abs(cov), 95),
                       vmax=np.percentile(np.abs(cov), 95))
        ax.set_title(f'{trt} ({r["label"]})')
        ax.set_xlabel('sample index')
        ax.set_ylabel('sample index')
        plt.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle(f'Standardized temporal-replicate covariance matrix ({label_suffix})',
                 fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f'fig_cov_heatmap_{label_suffix}.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--merged-ad", required=True)
    p.add_argument("--sample-list", required=True)
    p.add_argument("--output-dir", required=True)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    samples, chroms, positions, freqs, depths = load_data(
        args.merged_ad, args.sample_list)

    results_n80 = {}
    for trt in TREATMENTS:
        results_n80[trt] = run_one(trt, samples, chroms, positions, freqs, depths,
                                    args.output_dir, diploids_val=80, label='N=80')

    results_neff = {}
    for trt in TREATMENTS:
        results_neff[trt] = run_one(trt, samples, chroms, positions, freqs, depths,
                                     args.output_dir, diploids_val=29, label='N_eff=29')

    rows = []
    for label, results in [('N=80', results_n80), ('N_eff=29', results_neff)]:
        for trt in TREATMENTS:
            r = results[trt]
            rows.append({
                'treatment': trt,
                'correction': label,
                'n_loci': r['n_loci'],
                'n_tiles_1Mb': r['n_tiles'],
                'G_final_mean': r['G_final_mean'],
                'cc_estimate': r['cc_scalar'],
                'cc_ci_lower': r['cc_ci_lo'],
                'cc_ci_upper': r['cc_ci_hi'],
                'k1_residual': r['k1_residual'],
                'k2_mean': r['k2_mean'],
            })
    pd.DataFrame(rows).to_csv(os.path.join(args.output_dir, 'summary.tsv'),
                               sep='\t', index=False)

    plot_G_trajectories(results_n80, results_neff, args.output_dir)
    plot_convergence_summary(results_n80, results_neff, args.output_dir)
    plot_temporal_cov_heatmap(results_n80, args.output_dir, 'N80')
    plot_temporal_cov_heatmap(results_neff, args.output_dir, 'Neff29')

if __name__ == "__main__":
    main()
