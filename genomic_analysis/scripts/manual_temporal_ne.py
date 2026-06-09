#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd

N_EFF = 29
MIN_PHI = 1e-12
MAF_LOW = 0.01
MAF_HIGH = 0.99
MIN_DEPTH = 10
TREATMENTS = ["B", "T", "M"]

GENS_ALL = {
    1: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    2: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    3: [1, 2, 6, 7, 8, 9, 10],
    4: [1, 2, 6, 7, 8, 9, 10],
}


def load_data():
    samples = [l.strip() for l in open("variance_analysis/sample_list.txt") if l.strip()]
    ad_ref, ad_alt = [], []
    with open("variance_analysis/merged_ad.tsv") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
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


def ne_one_rep(samples, freq, total, trt, rep, n_eff=N_EFF,
               include_g10=False):
    gens = GENS_ALL[rep].copy()
    if not include_g10 and 10 in gens:
        gens.remove(10)

    sample_idx = {}
    for g in gens:
        sname = f"{trt}{rep}G{g:02d}"
        if sname in samples:
            sample_idx[g] = samples.index(sname)
    gens = sorted(sample_idx.keys())

    if len(gens) < 2:
        return None

    phi_sum = None
    phi_count = None
    T_total = 0
    pair_info = []

    for g1, g2 in zip(gens[:-1], gens[1:]):
        T = g2 - g1
        i1, i2 = sample_idx[g1], sample_idx[g2]
        p1 = freq[:, i1]; p2 = freq[:, i2]
        d1 = total[:, i1]; d2 = total[:, i2]

        ok = ((d1 >= MIN_DEPTH) & (d2 >= MIN_DEPTH)
              & np.isfinite(p1) & np.isfinite(p2)
              & (p1 >= MAF_LOW) & (p1 <= MAF_HIGH)
              & (p2 >= MAF_LOW) & (p2 <= MAF_HIGH))

        if ok.sum() < 1000:
            continue

        p1v = p1[ok]; p2v = p2[ok]
        d1v = d1[ok]; d2v = d2[ok]

        pbar = (p1v + p2v) / 2
        het = pbar * (1 - pbar)
        valid = het > 1e-10
        F_ijk = np.where(valid, (p1v - p2v)**2 / np.maximum(het, 1e-10), 0)

        inv_n1 = 1/(2*n_eff) + 1/d1v - 1/(2*n_eff*d1v)
        inv_n2 = 1/(2*n_eff) + 1/d2v - 1/(2*n_eff*d2v)

        phi_ijk = F_ijk - inv_n1 - inv_n2

        keep = valid 
        n_keep = int(keep.sum())

        if n_keep < 1000:
            continue

        phi_mean_pair = float(np.mean(phi_ijk[keep]))
        Ne_pair = T / (2 * phi_mean_pair) if phi_mean_pair > 0 else np.inf

        pair_info.append(dict(
            g1=g1, g2=g2, T=T,
            phi_mean=phi_mean_pair,
            n_sites=n_keep,
            Ne_pair=Ne_pair,
        ))
        T_total += T

    if not pair_info:
        return None

    Ne_values = [p["Ne_pair"] for p in pair_info
                 if np.isfinite(p["Ne_pair"]) and p["Ne_pair"] > 0]
    if not Ne_values:
        return None

    Ne_harmonic = len(Ne_values) / sum(1/v for v in Ne_values)
    Ne_arith = np.mean(Ne_values)

    return dict(
        treatment=trt, replicate=rep,
        n_pairs=len(pair_info),
        T_total=T_total,
        Ne_harmonic=Ne_harmonic,
        Ne_arith=Ne_arith,
        pairs=pair_info,
    )


def main():
    samples, freq, total = load_data()
    all_rows = []
    per_rep_summary = []
    for trt in TREATMENTS:
        rep_results = []
        for rep in [1, 2, 3, 4]:
            r = ne_one_rep(samples, freq, total, trt, rep)
            if r is None:
                print(f"  {trt}{rep}: no valid data", file=sys.stderr)
                continue
            rep_results.append(r)
            per_rep_summary.append(dict(
                treatment=trt, replicate=rep, n_pairs=r["n_pairs"],
                Ne_harmonic=r["Ne_harmonic"], Ne_arith=r["Ne_arith"],
            ))
            for p in r["pairs"]:
                all_rows.append(dict(
                    treatment=trt, replicate=rep,
                    g1=p["g1"], g2=p["g2"], T=p["T"],
                    phi_mean=p["phi_mean"], n_sites=p["n_sites"],
                    Ne_pair=p["Ne_pair"],
                ))

        ne_harm_values = [r["Ne_harmonic"] for r in rep_results
                          if np.isfinite(r["Ne_harmonic"])]
        ne_arith_values = [r["Ne_arith"] for r in rep_results
                           if np.isfinite(r["Ne_arith"])]

    pd.DataFrame(all_rows).to_csv(
        "variance_analysis/section1_rigorous/manual_temporal_ne_pairs.tsv",
        sep="\t", index=False)
    pd.DataFrame(per_rep_summary).to_csv(
        "variance_analysis/section1_rigorous/manual_temporal_ne_per_rep.tsv",
        sep="\t", index=False)

    df = pd.DataFrame(per_rep_summary)
    trt_summary = df.groupby("treatment").agg(
        n_reps=("replicate", "count"),
        Ne_mean_of_harmonic=("Ne_harmonic", "mean"),
        Ne_median_of_harmonic=("Ne_harmonic", "median"),
        Ne_min=("Ne_harmonic", "min"),
        Ne_max=("Ne_harmonic", "max"),
    ).reset_index()
    trt_summary.to_csv(
        "variance_analysis/section1_rigorous/manual_temporal_ne_summary.tsv",
        sep="\t", index=False)

    f_reg = {"B": 242, "T": 347, "M": 219}
    pool = {"B": 208, "T": 211, "M": 226}
    for trt in TREATMENTS:
        row = trt_summary[trt_summary["treatment"] == trt]
        if len(row):
            manual = row["Ne_mean_of_harmonic"].values[0]
            print(f"{trt:<10} {f_reg[trt]:>11.0f} {pool[trt]:>12.0f} {manual:>10.0f}")

if __name__ == "__main__":
    main()
