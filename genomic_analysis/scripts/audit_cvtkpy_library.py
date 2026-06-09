#!/usr/bin/env python3

import sys, os
import numpy as np
import pandas as pd
from collections import OrderedDict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'cvtkpy'))
from cvtk.cvtk import TiledTemporalFreqs
from cvtk.gintervals import GenomicIntervals
from cvtk.cov import stack_temporal_covariances

MIN_DEPTH = 20
MIN_MAF = 0.05
TILE_SIZE = 100_000
TREATMENTS = ["B", "T", "M"]


def load_data():
    samples = [l.strip() for l in open("variance_analysis/sample_list.txt") if l.strip()]
    chroms, positions, freqs_list, depths_list = [], [], [], []
    with open("variance_analysis/merged_ad.tsv") as f:
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
    freqs = np.array(freqs_list).T  # samples × sites
    depths = np.array(depths_list).T
    return samples, np.array(chroms), np.array(positions), freqs, depths


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


def get_samples_full(all_samples, treatment):
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


def get_samples_2tp(all_samples, treatment):
    indices, cvtk_samples = [], []
    for rep in range(1, 5):
        for gen in [1, 9]:
            sname = f"{treatment}{rep}G{gen:02d}"
            if sname in all_samples:
                indices.append(all_samples.index(sname))
                cvtk_samples.append((f"{treatment}{rep}", gen))
    return indices, cvtk_samples


def run_cc_analysis(all_samples, chroms, positions, freqs, depths,
                    treatment, sample_getter, diploids_val, label):
    indices, cvtk_samples = sample_getter(all_samples, treatment)
    freq_sub = freqs[indices, :]
    depth_sub = depths[indices, :]

    # Filter
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

    n_sites = fr.shape[1]

    cc = ttf.convergence_corr(bias_correction=True)
    cc_arr = np.array(cc).squeeze()

    result = {
        'treatment': treatment,
        'label': label,
        'diploids': diploids_val,
        'n_sites': n_sites,
        'R': ttf.R,
        'T': ttf.T,
    }

    if cc_arr.ndim == 2:
        result['cc_full_mean'] = float(np.nanmean(cc_arr))
        result['cc_diag_mean'] = float(np.nanmean(np.diag(cc_arr)))
        n = cc_arr.shape[0]
        mask_offdiag = ~np.eye(n, dtype=bool)
        result['cc_offdiag_mean'] = float(np.nanmean(cc_arr[mask_offdiag]))
        result['cc_matrix'] = cc_arr
    elif cc_arr.ndim == 0:
        result['cc_full_mean'] = float(cc_arr)
        result['cc_diag_mean'] = float(cc_arr)
        result['cc_offdiag_mean'] = np.nan
        result['cc_matrix'] = None
    else:
        result['cc_full_mean'] = float(np.nanmean(cc_arr))
        result['cc_diag_mean'] = float(cc_arr)
        result['cc_offdiag_mean'] = np.nan
        result['cc_matrix'] = cc_arr


    if ttf.T > 1:
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
        result['G_k2_per_rep'] = np.array(G_k2_per_rep)
        result['G_k2_mean'] = float(np.mean(G_k2_per_rep))

    return result


def main():
    all_samples, chroms, positions, freqs, depths = load_data()

    full_results = {}
    for trt in TREATMENTS:
        full_results[trt] = run_cc_analysis(
            all_samples, chroms, positions, freqs, depths,
            trt, get_samples_full, diploids_val=29, label="full_7tp")


    full_n80 = {}
    for trt in TREATMENTS:
        full_n80[trt] = run_cc_analysis(
            all_samples, chroms, positions, freqs, depths,
            trt, get_samples_full, diploids_val=80, label="full_7tp_N80")



    two_tp = {}
    for trt in TREATMENTS:
        two_tp[trt] = run_cc_analysis(
            all_samples, chroms, positions, freqs, depths,
            trt, get_samples_2tp, diploids_val=29, label="2tp_totalDp")


    for trt in TREATMENTS:
        r = full_results[trt]
        print(f"{trt:>4} | {'7tp N_eff=29':>20} | {r['cc_full_mean']:8.4f} | "
              f"{r['cc_diag_mean']:8.4f} | {r['cc_offdiag_mean']:10.4f}")
        r80 = full_n80[trt]
        print(f"{trt:>4} | {'7tp N=80':>20} | {r80['cc_full_mean']:8.4f} | "
              f"{r80['cc_diag_mean']:8.4f} | {r80['cc_offdiag_mean']:10.4f}")
        r2 = two_tp[trt]
        print(f"{trt:>4} | {'2tp total Δp N=29':>20} | {r2['cc_full_mean']:8.4f} | "
              f"{'n/a':>8} | {'n/a':>10}")



    for trt in TREATMENTS:
        m = full_results[trt].get('cc_matrix')
        if m is not None and m.ndim == 2:
            print(f"\n  {trt} (T={full_results[trt]['T']}, shape={m.shape}):")
            for i in range(m.shape[0]):
                row = "  " + "".join(f"{m[i,j]:8.4f}" for j in range(m.shape[1]))
                print(row)

    # G statistic
    print("\n\n--- G statistic (k≥2 off-diagonal only, N_eff=29) ---")
    print(f"{'Trt':>4} | {'G_k2_mean':>10} | per-rep values")
    for trt in TREATMENTS:
        r = full_results[trt]
        if 'G_k2_mean' in r:
            per_rep = ", ".join(f"{v:.4f}" for v in r['G_k2_per_rep'])
            print(f"{trt:>4} | {r['G_k2_mean']:10.4f} | [{per_rep}]")

    # Save summary
    rows = []
    for trt in TREATMENTS:
        for label, res in [("7tp_N29", full_results[trt]),
                           ("7tp_N80", full_n80[trt]),
                           ("2tp_N29", two_tp[trt])]:
            rows.append({
                'treatment': trt, 'config': label,
                'cc_full': res['cc_full_mean'],
                'cc_diag': res['cc_diag_mean'],
                'cc_offdiag': res.get('cc_offdiag_mean', np.nan),
                'G_k2': res.get('G_k2_mean', np.nan),
            })
    pd.DataFrame(rows).to_csv(
        "variance_analysis/section1_rigorous/cvtkpy_audit_library.tsv",
        sep="\t", index=False)
    print("\nSaved: variance_analysis/section1_rigorous/cvtkpy_audit_library.tsv")


if __name__ == "__main__":
    main()
