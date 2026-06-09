#!/usr/bin/env python3

import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'cvtkpy'))
from cvtk.cvtk import TiledTemporalFreqs
from cvtk.gintervals import GenomicIntervals
from cvtk.cov import (temporal_replicate_cov, stack_temporal_covariances,
                       calc_hets)
from cvtk.utils import (sort_samples, process_samples, reshape_matrix,
                         validate_diploids, integerize)

MIN_DEPTH = 20
MIN_MAF = 0.10
TILE_SIZE = 100_000  # 100kb, matching Buffalo & Coop for Drosophila


def load_data(ad_path, sample_list_path):
    samples = [l.strip() for l in open(sample_list_path) if l.strip()]
    chroms, positions, freqs_list, depths_list = [], [], [], []
    with open(ad_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chroms.append(parts[0])
            positions.append(int(parts[1]))
            refs, alts, tots = [], [], []
            for ad in parts[4:]:
                if ad in (".", ".,."):
                    refs.append(0); alts.append(0); tots.append(0)
                else:
                    r, a = ad.split(",")[:2]
                    r, a = int(r), int(a)
                    refs.append(r); alts.append(a); tots.append(r + a)
            tot = np.array(tots)
            alt = np.array(alts)
            freqs_list.append(np.where(tot > 0, alt / tot, np.nan))
            depths_list.append(tot)
    freqs = np.array(freqs_list).T
    depths = np.array(depths_list).T
    chroms = np.array(chroms)
    positions = np.array(positions)
    return samples, chroms, positions, freqs, depths


def get_treatment_indices(all_samples, treatment):
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


def estimate_neff(freq_sub, depth_sub, cvtk_samples, nominal_n=80):
    sorted_samples, sorted_i = sort_samples(cvtk_samples)
    _, _, nrep, ntp = process_samples(freq_sub, sorted_samples)
    fr_3d = reshape_matrix(freq_sub[sorted_i, :], nrep)
    dp_3d = reshape_matrix(depth_sub[sorted_i, :], nrep)
    dip = validate_diploids(np.full(len(cvtk_samples), nominal_n)[sorted_i],
                            nrep, ntp)
    raw_cov = temporal_replicate_cov(fr_3d, dp_3d, dip,
                                     bias_correction=False, standardize=False,
                                     share_first=False)
    R, T = nrep, ntp - 1
    avg_temp = np.mean(stack_temporal_covariances(raw_cov, R, T), axis=2)
    raw_k1 = np.mean(np.diag(avg_temp, 1))
    hets = calc_hets(fr_3d, dp_3d, dip)
    mean_het = float(np.nanmean(hets))
    mean_depth = float(np.nanmean(depth_sub))
    noise_total = -raw_k1 / (0.5 * mean_het)
    x = (noise_total - 1 / mean_depth) / (1 + 1 / mean_depth)
    n_eff = 1 / (2 * x) if x > 0 else nominal_n
    return n_eff


def make_tiles(chroms, positions, tile_size):
    from collections import OrderedDict
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
        for win_start in range(int(min_p) - int(min_p) % tile_size,
                               int(max_p) + 1, tile_size):
            tiles.append(chrom, win_start, win_start + tile_size)
    return tiles


def make_diagnostic_data(ttf, offdiag_k=1):
    tile_depths = ttf.depth_by_tile()
    mean_hets = ttf.calc_het_by_tile()
    seqids = ttf.tile_df['seqid'].values

    results = {}
    for use_correction in [False, True]:
        covs = ttf.calc_cov_by_tile(bias_correction=use_correction,
                                     standardize=False)
        offdiag = [np.nanmean(np.diag(c, k=offdiag_k)) for c in covs]
        diag = [np.nanmean(np.diag(c, k=0)) for c in covs]
        df = pd.DataFrame({
            'depth': tile_depths,
            'diag': diag,
            'offdiag': offdiag,
            'mean_het': mean_hets,
            'seqid': seqids,
        })
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        results[use_correction] = df.dropna()

    return results


def plot_correction_diagnostic(results_dict, treatment, n_label, ax_array):
    before = results_dict[False]
    after = results_dict[True]

    for col, (df, label) in enumerate([(before, 'before correction'),
                                        (after, 'after correction')]):
        ax = ax_array[0, col]
        ax.scatter(df['depth'], df['diag'], s=3, alpha=0.4,
                   c=integerize(df['seqid']), cmap='tab20')
        z = np.polyfit(df['depth'], df['diag'], 1)
        xr = np.linspace(df['depth'].min(), df['depth'].max(), 100)
        ax.plot(xr, np.polyval(z, xr), 'r-', lw=1.5)
        ax.set_ylabel('variance (diagonal)')
        ax.set_xlabel('mean tile depth')
        ax.set_title(f'{treatment} {n_label} — {label}')
        slope_str = f"slope={z[0]:.2e}"
        ax.annotate(slope_str, xy=(0.05, 0.90), xycoords='axes fraction',
                    fontsize=8)
        ax = ax_array[1, col]
        ax.scatter(df['depth'], df['offdiag'], s=3, alpha=0.4,
                   c=integerize(df['seqid']), cmap='tab20')
        z2 = np.polyfit(df['depth'], df['offdiag'], 1)
        ax.plot(xr, np.polyval(z2, xr), 'r-', lw=1.5)
        ax.axhline(0, color='0.6', ls='--', zorder=0)
        ax.set_ylabel('covariance (k=1 off-diag)')
        ax.set_xlabel('mean tile depth')
        slope_str2 = f"slope={z2[0]:.2e}"
        ax.annotate(slope_str2, xy=(0.05, 0.90), xycoords='axes fraction',
                    fontsize=8)


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

    for trt in ["B", "T", "M"]:

        indices, cvtk_samples = get_treatment_indices(samples, trt)
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

        n_eff = estimate_neff(fr, dp, cvtk_samples)
        n_eff_int = int(round(n_eff))

        gi = GenomicIntervals()
        for c, pp in zip(ch, po):
            gi.append(c, pp, pp + 1)
        tiles = make_tiles(ch, po, TILE_SIZE)

        ttf_80 = TiledTemporalFreqs(
            tiles, fr, cvtk_samples, depths=dp, diploids=80,
            gintervals=gi, swap=True, share_first=False)
        diag_80 = make_diagnostic_data(ttf_80)

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        plot_correction_diagnostic(diag_80, trt, 'N=80', axes)
        fig.suptitle(f'{trt} treatment — Correction diagnostic (N=80, nominal)',
                     fontsize=13, y=1.02)
        fig.tight_layout()
        fig.savefig(os.path.join(args.output_dir, f'{trt}_diagnostic_N80.png'),
                    dpi=150, bbox_inches='tight')
        plt.close(fig)

        ttf_neff = TiledTemporalFreqs(
            tiles, fr, cvtk_samples, depths=dp, diploids=n_eff_int,
            gintervals=gi, swap=True, share_first=False)
        diag_neff = make_diagnostic_data(ttf_neff)

        fig2, axes2 = plt.subplots(2, 2, figsize=(12, 8))
        plot_correction_diagnostic(diag_neff, trt, f'N_eff={n_eff_int}', axes2)
        fig2.suptitle(f'{trt} treatment — Correction diagnostic (N_eff={n_eff_int})',
                      fontsize=13, y=1.02)
        fig2.tight_layout()
        fig2.savefig(os.path.join(args.output_dir, f'{trt}_diagnostic_Neff.png'),
                     dpi=150, bbox_inches='tight')
        plt.close(fig2)

        fig3, axes3 = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
        for ax, (label, diag_data, n_val) in zip(axes3,
                [('N=80', diag_80, 80), (f'N_eff={n_eff_int}', diag_neff, n_eff_int)]):
            after = diag_data[True]
            ax.scatter(after['depth'], after['offdiag'], s=3, alpha=0.4,
                       c=integerize(after['seqid']), cmap='tab20')
            z = np.polyfit(after['depth'], after['offdiag'], 1)
            xr = np.linspace(after['depth'].min(), after['depth'].max(), 100)
            ax.plot(xr, np.polyval(z, xr), 'r-', lw=1.5)
            ax.axhline(0, color='0.6', ls='--', zorder=0)
            ax.set_xlabel('mean tile depth')
            ax.set_ylabel('corrected k=1 covariance')
            median_resid = np.nanmedian(after['offdiag'])
            ax.set_title(f'{trt} {label} (median={median_resid:.5f})')
        fig3.suptitle(f'{trt} — k=1 residual: N=80 vs N_eff',
                      fontsize=13, y=1.02)
        fig3.tight_layout()
        fig3.savefig(os.path.join(args.output_dir, f'{trt}_k1_comparison.png'),
                     dpi=150, bbox_inches='tight')
        plt.close(fig3)

        diag_80_k2 = make_diagnostic_data(ttf_80, offdiag_k=2)
        fig4, ax4 = plt.subplots(1, 1, figsize=(6, 4))
        after_k2 = diag_80_k2[True]
        ax4.scatter(after_k2['depth'], after_k2['offdiag'], s=3, alpha=0.4,
                    c=integerize(after_k2['seqid']), cmap='tab20')
        z4 = np.polyfit(after_k2['depth'], after_k2['offdiag'], 1)
        xr4 = np.linspace(after_k2['depth'].min(), after_k2['depth'].max(), 100)
        ax4.plot(xr4, np.polyval(z4, xr4), 'r-', lw=1.5)
        ax4.axhline(0, color='0.6', ls='--', zorder=0)
        ax4.set_xlabel('mean tile depth')
        ax4.set_ylabel('corrected k=2 covariance')
        ax4.set_title(f'{trt} — k=2 off-diagonal (N=80, should be flat at 0)')
        ax4.annotate(f"slope={z4[0]:.2e}", xy=(0.05, 0.90),
                     xycoords='axes fraction', fontsize=9)
        fig4.tight_layout()
        fig4.savefig(os.path.join(args.output_dir, f'{trt}_k2_diagnostic.png'),
                     dpi=150, bbox_inches='tight')
        plt.close(fig4)

if __name__ == "__main__":
    main()
