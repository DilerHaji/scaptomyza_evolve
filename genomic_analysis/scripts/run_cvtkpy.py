#!/usr/bin/env python3
import sys
import os
import numpy as np
import pandas as pd
from collections import OrderedDict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'cvtkpy'))
from cvtk.cvtk import TemporalFreqs, TiledTemporalFreqs
from cvtk.gintervals import GenomicIntervals
from cvtk.cov import (temporal_replicate_cov, stack_temporal_covariances,
                       calc_hets)
from cvtk.utils import (sort_samples, process_samples, reshape_matrix,
                         validate_diploids)

MIN_DEPTH = 20
MIN_MAF = 0.10


def load_data(ad_path, sample_list_path):
    samples = [line.strip() for line in open(sample_list_path) if line.strip()]

    chroms, positions = [], []
    freqs_list = []
    depths_list = []

    with open(ad_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chroms.append(parts[0])
            positions.append(int(parts[1]))
            refs, alts, tots = [], [], []
            for ad_field in parts[4:]:
                if ad_field == "." or ad_field == ".,.":
                    refs.append(0); alts.append(0); tots.append(0)
                else:
                    r, a = ad_field.split(",")[:2]
                    r, a = int(r), int(a)
                    refs.append(r); alts.append(a); tots.append(r + a)
            tot = np.array(tots)
            alt = np.array(alts)
            freq = np.where(tot > 0, alt / tot, np.nan)
            freqs_list.append(freq)
            depths_list.append(tot)

    freqs = np.array(freqs_list).T    # samples × loci
    depths = np.array(depths_list).T  # samples × loci
    chroms = np.array(chroms)
    positions = np.array(positions)

    return samples, chroms, positions, freqs, depths


def get_treatment_samples(all_samples, treatment):
    full_gens = [1, 2, 6, 7, 8, 9]
    indices = []
    cvtk_samples = []

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

    k2plus = []
    for k in range(2, T):
        k2plus.extend(np.diag(avg_temp, k).tolist())
    k2_mean = np.mean(k2plus) if k2plus else 0.0

    return n_eff, raw_k1, k2_mean, mean_het, mean_depth


def make_gintervals(chroms, positions):
    gi = GenomicIntervals()
    for c, p in zip(chroms, positions):
        gi.append(c, p, p + 1)
    return gi


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
        for win_start in range(int(min_p) - int(min_p) % tile_size,
                               int(max_p) + 1, tile_size):
            win_end = win_start + tile_size
            tiles.append(chrom, win_start, win_end)

    return tiles


def run_analysis(treatment, all_samples, chroms, positions, freqs, depths,
                 output_dir, tile_size=200000, n_bootstrap=1000):

    indices, cvtk_samples = get_treatment_samples(all_samples, treatment)

    timepoints = sorted(set(s[1] for s in cvtk_samples))

    freq_sub = freqs[indices, :]
    depth_sub = depths[indices, :]

    min_depth_ok = np.all(depth_sub >= MIN_DEPTH, axis=0)
    mean_freq = np.nanmean(freq_sub, axis=0)
    maf_ok = (mean_freq >= MIN_MAF) & (mean_freq <= 1 - MIN_MAF)
    all_finite = np.all(np.isfinite(freq_sub), axis=0)
    keep = min_depth_ok & maf_ok & all_finite

    freq_filtered = freq_sub[:, keep]
    depth_filtered = depth_sub[:, keep]
    chroms_filtered = chroms[keep]
    positions_filtered = positions[keep]

    n_eff, raw_k1, k2_mean, mean_het, mean_depth = estimate_neff(
        freq_filtered, depth_filtered, cvtk_samples, nominal_n=80)
    n_eff_int = int(round(n_eff))

    gi = make_gintervals(chroms_filtered, positions_filtered)

    tiles = make_tiles(chroms_filtered, positions_filtered, tile_size)

    ttf = TiledTemporalFreqs(
        tiles,
        freq_filtered,
        cvtk_samples,
        depths=depth_filtered,
        diploids=n_eff_int,
        gintervals=gi,
        swap=True,
        share_first=False,
    )

    gw_cov = ttf.calc_cov(bias_correction=True, standardize=True)
    np.savetxt(os.path.join(output_dir, f"{treatment}_gw_cov.tsv"),
               gw_cov, delimiter="\t", fmt="%.8f")

    G = ttf.calc_G()
    G_final = G[-1, :] if G.ndim == 2 else G
    G_mean = float(np.nanmean(G_final))
    gw_cov_unstd = ttf.calc_cov(bias_correction=True, standardize=False)
    temp_covs = stack_temporal_covariances(gw_cov_unstd, ttf.R, ttf.T)
    T = ttf.T
    G_k2_per_rep = []
    for r in range(ttf.R):
        tc = temp_covs[:, :, r]
        offdiag_k2 = np.tril(tc, -2) + np.triu(tc, 2)
        total_cov_k2 = np.nansum(offdiag_k2)
        tv = ttf.calc_var(t=T, standardize=False, bias_correction=True)
        G_k2_per_rep.append(total_cov_k2 / tv[r])
    G_k2_per_rep = np.array(G_k2_per_rep)
    G_k2_mean = float(np.nanmean(G_k2_per_rep))

    k1_corrected = np.mean(np.diag(np.mean(temp_covs, axis=2), 1))

    g_results = pd.DataFrame({
        "replicate": [f"rep{r+1}" for r in range(ttf.R)] + ["mean"],
        "G_neff_corrected": list(G_final) + [G_mean],
        "G_k2_validation": list(G_k2_per_rep) + [G_k2_mean],
    })
    g_results.to_csv(os.path.join(output_dir, f"{treatment}_G.tsv"),
                     sep="\t", index=False)
    if G.ndim == 2:
        g_df = pd.DataFrame(G, columns=[f"rep{r+1}" for r in range(G.shape[1])],
                            index=[f"t{t+1}" for t in range(G.shape[0])])
        g_df.to_csv(os.path.join(output_dir, f"{treatment}_G_matrix.tsv"),
                    sep="\t")


    conv_corr = ttf.convergence_corr(bias_correction=True)
    conv_corr_val = float(np.nanmean(conv_corr)) if hasattr(conv_corr, '__len__') else float(conv_corr)

    tile_conv_matrix = ttf.convergence_corr_by_tile(bias_correction=True)
    triu_idx = np.triu_indices(T, k=1)
    tile_conv_scalar = []
    for i in range(tile_conv_matrix.shape[0]):
        tile_mat = tile_conv_matrix[i]
        tile_conv_scalar.append(float(np.nanmean(tile_mat[triu_idx])))

    tile_df = ttf.tile_df.copy()
    tile_df["convergence_corr"] = tile_conv_scalar
    tile_df["n_snps"] = [len(idx) for idx in ttf.tile_indices]
    tile_df.to_csv(os.path.join(output_dir, f"{treatment}_tiles.tsv"),
                   sep="\t", index=False)


    G_ci_lower, G_bs_est, G_ci_upper = ttf.bootstrap_G(
        B=n_bootstrap, alpha=0.05, average_replicates=True, progress_bar=True)
    G_bs_final = float(G_bs_est[-1]) if hasattr(G_bs_est, '__len__') else float(G_bs_est)
    G_ci_lo = float(G_ci_lower[-1]) if hasattr(G_ci_lower, '__len__') else float(G_ci_lower)
    G_ci_hi = float(G_ci_upper[-1]) if hasattr(G_ci_upper, '__len__') else float(G_ci_upper)

    cc_ci_lower, cc_bs_est, cc_ci_upper = ttf.bootstrap_convergence_corr(
        B=n_bootstrap, alpha=0.05, progress_bar=True)
    cc_bs_val = float(np.nanmean(cc_bs_est)) if hasattr(cc_bs_est, '__len__') else float(cc_bs_est)
    cc_ci_lo = float(np.nanmean(cc_ci_lower)) if hasattr(cc_ci_lower, '__len__') else float(cc_ci_lower)
    cc_ci_hi = float(np.nanmean(cc_ci_upper)) if hasattr(cc_ci_upper, '__len__') else float(cc_ci_upper)

    bs_df = pd.DataFrame({
        "statistic": ["G", "convergence_corr"],
        "estimate": [G_bs_final, cc_bs_val],
        "ci_lower": [min(G_ci_lo, G_ci_hi), min(cc_ci_lo, cc_ci_hi)],
        "ci_upper": [max(G_ci_lo, G_ci_hi), max(cc_ci_lo, cc_ci_hi)],
    })
    bs_df.to_csv(os.path.join(output_dir, f"{treatment}_bootstrap.tsv"),
                 sep="\t", index=False)

    summary = {
        "treatment": treatment,
        "n_replicates": ttf.R,
        "n_timepoints": ttf.ntimepoints,
        "timepoints": str(timepoints),
        "n_loci": ttf.L,
        "n_tiles": ttf.ntiles,
        "N_eff": round(n_eff, 1),
        "overdispersion": round(80 / n_eff, 2),
        "G_neff_corrected": G_mean,
        "G_k2_validation": G_k2_mean,
        "residual_k1": k1_corrected,
        "G_bs_estimate": G_bs_final,
        "G_bs_ci_lower": min(G_ci_lo, G_ci_hi),
        "G_bs_ci_upper": max(G_ci_lo, G_ci_hi),
        "convergence_corr": conv_corr_val,
        "cc_bs_estimate": cc_bs_val,
        "cc_bs_ci_lower": min(cc_ci_lo, cc_ci_hi),
        "cc_bs_ci_upper": max(cc_ci_lo, cc_ci_hi),
    }

    return summary


def main():
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--merged-ad", required=True)
    p.add_argument("--sample-list", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--tile-size", type=int, default=200000,
                   help="Window size for tiled analysis (default: 200kb)")
    p.add_argument("--n-bootstrap", type=int, default=1000,
                   help="Number of block bootstrap resamples (default: 1000)")
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    samples, chroms, positions, freqs, depths = load_data(
        args.merged_ad, args.sample_list)

    summaries = []
    for trt in ["B", "T", "M"]:
        summary = run_analysis(trt, samples, chroms, positions, freqs, depths,
                               args.output_dir, args.tile_size, args.n_bootstrap)
        summaries.append(summary)

    pd.DataFrame(summaries).to_csv(
        os.path.join(args.output_dir, "cvtkpy_summary.tsv"),
        sep="\t", index=False)


if __name__ == "__main__":
    main()
