#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd
from scipy import stats

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
    return samples, freq, total


def main():
    samples, freq, total = load_data()
    for trt in ["B", "T", "M"]:
        rep_indices_g01 = [samples.index(f"{trt}{r}G01") for r in range(1, 5) if f"{trt}{r}G01" in samples]
        rep_indices_g09 = [samples.index(f"{trt}{r}G09") for r in range(1, 5) if f"{trt}{r}G09" in samples]

        if len(rep_indices_g01) != 4 or len(rep_indices_g09) != 4:
            continue

        f01 = freq[:, rep_indices_g01]  # n_sites × 4
        f09 = freq[:, rep_indices_g09]
        d01 = total[:, rep_indices_g01]
        d09 = total[:, rep_indices_g09]

        mask = (np.all(d01 >= 10, axis=1) & np.all(d09 >= 10, axis=1) &
                np.all(np.isfinite(f01), axis=1) & np.all(np.isfinite(f09), axis=1))
        pbar01 = np.nanmean(f01, axis=1)
        mask = mask & (np.minimum(pbar01, 1 - pbar01) >= 0.05)

        f01_v = f01[mask]
        f09_v = f09[mask]
        delta_p = f09_v - f01_v  # n_sites × 4

        within_var = np.mean(np.var(delta_p, axis=0, ddof=0))  # variance of all Δp values within each rep, averaged across sites... wait
        within_var_per_rep = np.mean(delta_p**2, axis=0)  # one value per rep
        within_var = np.mean(within_var_per_rep)

        pbar09 = np.nanmean(f09_v, axis=1)
        between_var_per_site = np.var(f09_v, axis=1, ddof=1)  # variance across reps at each site
        between_var = np.mean(between_var_per_site)

        mean_delta = np.mean(delta_p, axis=1)  # n_sites
        between_delta_var = np.var(delta_p, axis=1, ddof=1)  # n_sites
        between_delta_mean = np.mean(between_delta_var)

        ratio = between_delta_mean / within_var if within_var > 0 else np.nan
        parallel = 1 - ratio

    for trt in ["B", "T", "M"]:
        rep_indices_g01 = [samples.index(f"{trt}{r}G01") for r in range(1, 5) if f"{trt}{r}G01" in samples]
        rep_indices_g09 = [samples.index(f"{trt}{r}G09") for r in range(1, 5) if f"{trt}{r}G09" in samples]
        if len(rep_indices_g01) != 4 or len(rep_indices_g09) != 4:
            continue

        f01 = freq[:, rep_indices_g01]
        f09 = freq[:, rep_indices_g09]
        d01 = total[:, rep_indices_g01]
        d09 = total[:, rep_indices_g09]

        mask = (np.all(d01 >= 10, axis=1) & np.all(d09 >= 10, axis=1) &
                np.all(np.isfinite(f01), axis=1) & np.all(np.isfinite(f09), axis=1))
        pbar01 = np.nanmean(f01, axis=1)
        mask = mask & (np.minimum(pbar01, 1 - pbar01) >= 0.05)

        f01_v = f01[mask]
        f09_v = f09[mask]
        delta_p = f09_v - f01_v

        total_var = np.var(delta_p.flatten(), ddof=1)

        mean_delta_per_site = np.mean(delta_p, axis=1)  # one per site
        shared_var = np.var(mean_delta_per_site, ddof=1)

        residual = delta_p - mean_delta_per_site[:, None]
        rep_var = np.var(residual.flatten(), ddof=1)

        pct_shared = 100 * shared_var / (shared_var + rep_var)


if __name__ == "__main__":
    main()
