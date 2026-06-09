#!/usr/bin/env python3
import sys
import os
import numpy as np
import pandas as pd
from scipy import stats
from collections import Counter
from itertools import combinations

OUTDIR = "variance_analysis/section1_rigorous"


def load_data():
    samples = [l.strip() for l in open("variance_analysis/sample_list.txt") if l.strip()]
    chroms, positions, ad_ref, ad_alt = [], [], [], []
    with open("variance_analysis/merged_ad.tsv") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chroms.append(parts[0])
            positions.append(int(parts[1]))
            refs, alts = [], []
            for ad in parts[4:]:
                if ad == "." or ad == ".,.":
                    refs.append(0); alts.append(0)
                else:
                    r, a = ad.split(",")[:2]
                    refs.append(int(r)); alts.append(int(a))
            ad_ref.append(refs); ad_alt.append(alts)
    ad_ref = np.array(ad_ref, dtype=np.int32)
    ad_alt = np.array(ad_alt, dtype=np.int32)
    total = ad_ref + ad_alt
    freq = np.where(total > 0, ad_alt / total, np.nan)
    return samples, np.array(chroms), np.array(positions), freq, total


def compute_F_ratio_of_averages(freq, total, rep_indices, min_depth=10, min_maf=0.05,
                                 site_subset=None):
    k = len(rep_indices)
    if k < 2:
        return np.nan, 0

    f_r = freq[:, rep_indices]
    d_r = total[:, rep_indices]
    mask = (np.all(d_r >= min_depth, axis=1) & np.all(np.isfinite(f_r), axis=1))
    pbar = np.nanmean(f_r, axis=1)
    maf_ok = np.minimum(pbar, 1 - pbar) >= min_maf
    mask = mask & maf_ok

    if site_subset is not None:
        mask = mask & site_subset

    if mask.sum() < 100:
        return np.nan, 0

    pbar_v = pbar[mask]
    f_v = f_r[mask]

    ss = np.nansum((f_v - pbar_v[:, None])**2, axis=1) / (k - 1)
    het = pbar_v * (1 - pbar_v)

    valid = het > 1e-10
    if valid.sum() == 0:
        return np.nan, 0

    F = np.mean(ss[valid]) / np.mean(het[valid])
    return F, int(valid.sum())



def block_bootstrap_ne(freq, total, samples, chroms, positions,
                       treatment, n_bootstraps=1000, block_size=100000, seed=42):
    rng = np.random.RandomState(seed)

    blocks = []  # list of arrays of site indices
    for chrom in np.unique(chroms):
        chrom_idx = np.where(chroms == chrom)[0]
        if len(chrom_idx) == 0:
            continue
        chrom_pos = positions[chrom_idx]
        min_p, max_p = chrom_pos.min(), chrom_pos.max()
        block_starts = np.arange(min_p, max_p + block_size, block_size)
        for bs in block_starts:
            in_block = chrom_idx[(chrom_pos >= bs) & (chrom_pos < bs + block_size)]
            if len(in_block) >= 10:
                blocks.append(in_block)

    n_blocks = len(blocks)

    def compute_ne_from_subset(site_mask):
        gens, fs = [], []
        for gen in range(1, 10):
            rep_names = [f"{treatment}{r}G{gen:02d}" for r in range(1, 5)]
            rep_idx = [samples.index(s) for s in rep_names if s in samples]
            if len(rep_idx) != 4:
                continue
            F_val, n = compute_F_ratio_of_averages(freq, total, rep_idx,
                                                    site_subset=site_mask)
            if not np.isnan(F_val):
                gens.append(gen); fs.append(F_val)
        if len(gens) < 3:
            return np.nan
        sl, _, _, _, _ = stats.linregress(gens, fs)
        return 1 / (2 * sl) if sl > 0 else np.inf

    all_sites = np.ones(freq.shape[0], dtype=bool)
    ne_point = compute_ne_from_subset(all_sites)

    ne_boots = []
    for b in range(n_bootstraps):
        sampled_blocks = rng.choice(n_blocks, n_blocks, replace=True)
        site_mask = np.zeros(freq.shape[0], dtype=bool)
        for bi in sampled_blocks:
            site_mask[blocks[bi]] = True
        ne_b = compute_ne_from_subset(site_mask)
        if not np.isnan(ne_b) and ne_b < 1e6:
            ne_boots.append(ne_b)

    ne_boots = np.array(ne_boots)
    return {
        "Ne_point": ne_point,
        "Ne_median": np.median(ne_boots),
        "Ne_2.5": np.percentile(ne_boots, 2.5),
        "Ne_97.5": np.percentile(ne_boots, 97.5),
        "n_boot": len(ne_boots),
    }



def lynch_temporal_ne(freq, total, samples, treatment,
                      min_depth=20, min_maf=0.05):
    n_pool = 80
    n_hap = 2 * n_pool
    timepoints = [1, 2, 6, 7, 8, 9]

    ne_per_pair = {}

    for rep in range(1, 5):
        rep_sample_idx = []
        rep_tps = []
        for tp in timepoints:
            sname = f"{treatment}{rep}G{tp:02d}"
            if sname in samples:
                rep_sample_idx.append(samples.index(sname))
                rep_tps.append(tp)

        if len(rep_sample_idx) < 2:
            continue

        for i in range(len(rep_tps) - 1):
            t1, t2 = rep_tps[i], rep_tps[i+1]
            idx1, idx2 = rep_sample_idx[i], rep_sample_idx[i+1]

            f1 = freq[:, idx1]
            f2 = freq[:, idx2]
            d1 = total[:, idx1]
            d2 = total[:, idx2]

            mask = ((d1 >= min_depth) & (d2 >= min_depth) &
                    np.isfinite(f1) & np.isfinite(f2))
            pbar = (f1 + f2) / 2
            mask = mask & (np.minimum(pbar, 1-pbar) >= min_maf)

            if mask.sum() < 100:
                continue

            f1m = f1[mask]; f2m = f2[mask]
            pbar_m = pbar[mask]
            d1m = d1[mask]; d2m = d2[mask]

            het = pbar_m * (1 - pbar_m)
            valid = het > 1e-10

            Fhat = (f1m[valid] - f2m[valid])**2 / het[valid]

            inv_neff1 = 1.0/n_hap + 1.0/d1m[valid] - 1.0/(n_hap * d1m[valid])
            inv_neff2 = 1.0/n_hap + 1.0/d2m[valid] - 1.0/(n_hap * d2m[valid])

            phi_per_site = Fhat - inv_neff1 - inv_neff2
            phi_mean = np.mean(phi_per_site)

            t_diff = t2 - t1
            if phi_mean > 0:
                ne_pair = t_diff / (2 * phi_mean)
            else:
                ne_pair = np.inf

            ne_per_pair[(rep, t1, t2)] = ne_pair

    ne_vals = [v for v in ne_per_pair.values() if np.isfinite(v) and v > 0]
    if not ne_vals:
        return np.nan, ne_per_pair
    harmonic_mean = len(ne_vals) / np.sum([1.0/v for v in ne_vals])
    return harmonic_mean, ne_per_pair




def simulate_drift_with_selection(Ne, n_reps, n_sites, timepoints, n_pool, depth,
                                  s_strength=0.0, frac_selected=0.0, phi=2.78,
                                  seed=None):
    rng = np.random.RandomState(seed)
    n_tp = len(timepoints)
    max_gen = max(timepoints)

    starting_freqs = rng.uniform(0.1, 0.9, n_sites)

    s_per_site = np.zeros(n_sites)
    if s_strength > 0 and frac_selected > 0:
        n_selected = int(n_sites * frac_selected)
        selected_idx = rng.choice(n_sites, n_selected, replace=False)
        signs = rng.choice([-1, 1], n_selected)
        s_per_site[selected_idx] = s_strength * signs

    freq_obs = np.zeros((n_reps, n_tp, n_sites))

    for rep in range(n_reps):
        founder_draws = rng.binomial(1000, starting_freqs)
        true_p = founder_draws / 1000

        current_p = true_p.copy()
        tp_idx = 0

        for g in range(max_gen + 1):
            if tp_idx < n_tp and g == timepoints[tp_idx]:
                n_eff = max(1, int(n_pool / phi))
                n_haploid_pool = 2 * n_eff
                for site in range(n_sites):
                    p = current_p[site]
                    alt_pool = rng.binomial(n_haploid_pool, p)
                    p_pool = alt_pool / n_haploid_pool
                    alt_reads = rng.binomial(depth, p_pool)
                    freq_obs[rep, tp_idx, site] = alt_reads / depth
                tp_idx += 1

            if g < max_gen:
                if s_strength > 0:
                    delta = s_per_site * current_p * (1 - current_p)
                    current_p = current_p + delta
                n_haploid = 2 * Ne
                new_counts = rng.binomial(n_haploid, np.clip(current_p, 0.001, 0.999))
                current_p = new_counts / n_haploid
                current_p = np.clip(current_p, 0.001, 0.999)

    return freq_obs


def compute_F_regression_from_sim(freq_obs, timepoints):
    n_reps, n_tp, n_sites = freq_obs.shape

    F_vals = []
    valid_tps = []
    for ti, tp in enumerate(timepoints):
        if tp == 0:
            continue
        p_reps = freq_obs[:, ti, :]
        pbar = np.mean(p_reps, axis=0)
        maf = np.minimum(pbar, 1 - pbar)
        ok = maf >= 0.05

        if ok.sum() < 100:
            continue

        ss = np.sum((p_reps[:, ok] - pbar[ok][None, :])**2, axis=0) / (n_reps - 1)
        het = pbar[ok] * (1 - pbar[ok])
        valid = het > 1e-10
        F = np.mean(ss[valid]) / np.mean(het[valid])
        F_vals.append(F)
        valid_tps.append(tp)

    if len(F_vals) < 3:
        return np.nan, np.nan

    sl, _, _, _, _ = stats.linregress(valid_tps, F_vals)
    Ne_est = 1 / (2 * sl) if sl > 0 else np.inf
    return sl, Ne_est


def compute_temporal_cov_from_sim(freq_obs, timepoints):
    n_reps, n_tp, n_sites = freq_obs.shape
    T = n_tp - 1

    delta_p = np.diff(freq_obs, axis=1)
    het = np.zeros_like(delta_p)
    for t in range(T):
        p_mid = (freq_obs[:, t, :] + freq_obs[:, t+1, :]) / 2
        het[:, t, :] = 2 * p_mid * (1 - p_mid)

    grand_mean = np.mean(freq_obs, axis=(0, 1))
    ok = (grand_mean >= 0.10) & (grand_mean <= 0.90)

    delta_f = delta_p[:, :, ok]
    het_f = het[:, :, ok]

    tcov = np.zeros((T, T))
    for r in range(n_reps):
        for s in range(T):
            for t in range(T):
                h_prod = np.sqrt(het_f[r, s, :] * het_f[r, t, :])
                h_prod[h_prod < 1e-10] = np.nan
                tcov[s, t] += np.nanmean(delta_f[r, s, :] * delta_f[r, t, :] / h_prod)
    tcov /= n_reps

    diag_mean = np.mean(np.diag(tcov))
    adj_mean = np.mean([tcov[i, i+1] for i in range(T-1)])

    V_noise = abs(adj_mean)
    V_drift = diag_mean - 2 * V_noise
    noise_fraction = (2 * V_noise) / diag_mean if diag_mean > 0 else np.nan

    return {
        "diag": diag_mean,
        "adj": adj_mean,
        "V_noise": V_noise,
        "V_drift": V_drift,
        "noise_fraction": noise_fraction,
        "Ne_temporal": 1/(2*V_drift) if V_drift > 0 else np.inf,
    }




def main():
    os.makedirs(OUTDIR, exist_ok=True)
    samples, chroms, positions, freq, total = load_data()

    f_results = []
    for trt in ["B", "T", "M"]:
        for gen in range(1, 10):
            rep_names = [f"{trt}{r}G{gen:02d}" for r in range(1, 5)]
            rep_idx = [samples.index(s) for s in rep_names if s in samples]
            k = len(rep_idx)
            F_val, n_sites = compute_F_ratio_of_averages(freq, total, rep_idx)
            f_results.append({
                "treatment": trt, "generation": gen, "k": k,
                "F_ratio_avg": F_val, "n_sites": n_sites
            })
    f_df = pd.DataFrame(f_results)
    f_df.to_csv(os.path.join(OUTDIR, "F_ratio_of_averages.tsv"), sep="\t", index=False)

    ne_point = {}
    for trt in ["B", "T", "M"]:
        sub = f_df[(f_df["treatment"] == trt) & (f_df["k"] == 4)]
        sub = sub.dropna()
        if len(sub) >= 3:
            sl, ic, r, p, se = stats.linregress(sub["generation"].values, sub["F_ratio_avg"].values)
            Ne = 1/(2*sl) if sl > 0 else np.inf
            ne_point[trt] = Ne

    boot_results = []
    for trt in ["B", "T", "M"]:
        boot = block_bootstrap_ne(freq, total, samples, chroms, positions,
                                   trt, n_bootstraps=300, block_size=100000)
        boot["treatment"] = trt
        boot_results.append(boot)
    boot_df = pd.DataFrame(boot_results)
    boot_df.to_csv(os.path.join(OUTDIR, "ne_bootstrap_cis.tsv"), sep="\t", index=False)


    lynch_results = []
    for trt in ["B", "T", "M"]:
        ne_lynch, pairs = lynch_temporal_ne(freq, total, samples, trt)
        lynch_results.append({"treatment": trt, "Ne_lynch": ne_lynch, "n_pairs": len(pairs)})
    lynch_df = pd.DataFrame(lynch_results)
    lynch_df.to_csv(os.path.join(OUTDIR, "ne_lynch.tsv"), sep="\t", index=False)


    timepoints_sim = [0, 1, 2, 6, 7, 8, 9]

    sim_results = []
    n_sims = 100  # per parameter combo

    Ne_sim = 500  # use census Ne

    drift_only = []
    for s in range(n_sims):
        freq_obs = simulate_drift_with_selection(
            Ne=Ne_sim, n_reps=4, n_sites=5000,
            timepoints=timepoints_sim, n_pool=80, depth=83,
            s_strength=0.0, frac_selected=0.0, phi=2.78, seed=s)
        sl, ne_b = compute_F_regression_from_sim(freq_obs, timepoints_sim)
        tcov = compute_temporal_cov_from_sim(freq_obs, timepoints_sim)
        drift_only.append({
            "scenario": "drift_only",
            "s": 0, "frac": 0,
            "ne_between": ne_b,
            "ne_temporal": tcov["Ne_temporal"],
            "noise_frac": tcov["noise_fraction"],
            "diag": tcov["diag"],
            "adj": tcov["adj"],
        })
    sim_results.extend(drift_only)

    no_sel = [d for d in drift_only if not np.isnan(d["noise_frac"])]
    median_noise_frac = np.median([d["noise_frac"] for d in no_sel])

    for s_strength in [0.05, 0.1, 0.2]:
        for frac in [0.01, 0.05]:
            for sim in range(n_sims):
                freq_obs = simulate_drift_with_selection(
                    Ne=Ne_sim, n_reps=4, n_sites=5000,
                    timepoints=timepoints_sim, n_pool=80, depth=83,
                    s_strength=s_strength, frac_selected=frac, phi=2.78, seed=sim+10000)
                sl, ne_b = compute_F_regression_from_sim(freq_obs, timepoints_sim)
                tcov = compute_temporal_cov_from_sim(freq_obs, timepoints_sim)
                sim_results.append({
                    "scenario": f"s={s_strength}_frac={frac}",
                    "s": s_strength, "frac": frac,
                    "ne_between": ne_b,
                    "ne_temporal": tcov["Ne_temporal"],
                    "noise_frac": tcov["noise_fraction"],
                    "diag": tcov["diag"],
                    "adj": tcov["adj"],
                })

    sim_df = pd.DataFrame(sim_results)
    sim_df.to_csv(os.path.join(OUTDIR, "simulation_results.tsv"), sep="\t", index=False)

    # Summary
    for scenario in sim_df["scenario"].unique():
        sub = sim_df[sim_df["scenario"] == scenario]
        ne_med = np.nanmedian(sub["ne_between"])
        ne_lo = np.nanpercentile(sub["ne_between"], 25)
        ne_hi = np.nanpercentile(sub["ne_between"], 75)

    for trt in ["B", "T", "M"]:
        boot = next(b for b in boot_results if b["treatment"] == trt)
        lynch = next(l for l in lynch_results if l["treatment"] == trt)

if __name__ == "__main__":
    main()
