#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd

TREATMENTS = ["B", "T", "M"]
GENS_4REP = [1, 2, 6, 7, 8, 9]
N_EFF = 29
MIN_DEPTH = 10
MIN_MAF = 0.05


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


def build_common_filter(samples, freq, total):
    founder_idx = [samples.index(f"F{i}G00") for i in range(1, 5)]
    founder_freq = np.nanmean(freq[:, founder_idx], axis=1)
    maf_ok = (np.minimum(founder_freq, 1 - founder_freq) >= MIN_MAF) & np.isfinite(founder_freq)
    depth_ok = np.ones(freq.shape[0], dtype=bool)
    for trt in TREATMENTS:
        for gen in GENS_4REP:
            for rep in range(1, 5):
                sname = f"{trt}{rep}G{gen:02d}"
                if sname in samples:
                    depth_ok &= (total[:, samples.index(sname)] >= MIN_DEPTH)
    return maf_ok & depth_ok


def pair_mean(mat, func):
    k = mat.shape[1]
    vals = []
    for a in range(k):
        for b in range(a+1, k):
            vals.append(func(mat[:, a], mat[:, b]))
    return float(np.mean(vals))


def cc_total_delta(samples, freq, total, mask, bias_correct, n_eff=N_EFF):
    out = {}
    keep = np.where(mask)[0]
    for trt in TREATMENTS:
        delta_list = []
        noise_var_list = []
        for rep in range(1, 5):
            s0 = f"{trt}{rep}G01"; s1 = f"{trt}{rep}G09"
            if s0 not in samples or s1 not in samples:
                continue
            i0 = samples.index(s0); i1 = samples.index(s1)
            d = freq[keep, i1] - freq[keep, i0]
            delta_list.append(d)
            if bias_correct:
                p0 = freq[keep, i0]; p1 = freq[keep, i1]
                d0 = total[keep, i0]; d1 = total[keep, i1]
                het0 = p0 * (1 - p0); het1 = p1 * (1 - p1)
                nv = (het0 * (1/(2*n_eff) + 1/d0) +
                      het1 * (1/(2*n_eff) + 1/d1))
                noise_var_list.append(nv)
        D = np.stack(delta_list, axis=1)  # n_sites × 4
        D0 = D - D.mean(axis=0, keepdims=True)
        n_sites = D.shape[0]


        def cov_pair(x, y):
            return np.sum(x * y) / (n_sites - 1)
        num = pair_mean(D0, cov_pair)


        if bias_correct:
            var_per_rep = []
            for r in range(D.shape[1]):
                raw_var = np.sum(D0[:, r]**2) / (n_sites - 1)
                noise_contrib = noise_var_list[r].mean()
                var_per_rep.append(max(raw_var - noise_contrib, 1e-12))
            vr = np.array(var_per_rep)
        else:
            vr = np.array([np.sum(D0[:, r]**2) / (n_sites - 1)
                           for r in range(D.shape[1])])


        denom_vals = []
        for a in range(len(vr)):
            for b in range(a+1, len(vr)):
                denom_vals.append(np.sqrt(vr[a] * vr[b]))
        denom = float(np.mean(denom_vals))

        cc = num / denom if denom > 0 else np.nan
        out[trt] = dict(cc=cc, num=num, denom=denom, var_per_rep=vr.tolist())
    return out


def cc_per_interval(samples, freq, total, mask, bias_correct, n_eff=N_EFF,
                    mode="diagonal"):
    keep = np.where(mask)[0]
    intervals = list(zip(GENS_4REP[:-1], GENS_4REP[1:]))
    n_int = len(intervals)

    out = {}
    for trt in TREATMENTS:
        R = 4
        deltas = np.zeros((R, n_int, len(keep)))
        noise_var = np.zeros((R, n_int, len(keep)))
        for r in range(R):
            for ti, (g0, g1) in enumerate(intervals):
                s0 = f"{trt}{r+1}G{g0:02d}"; s1 = f"{trt}{r+1}G{g1:02d}"
                if s0 not in samples or s1 not in samples:
                    continue
                i0 = samples.index(s0); i1 = samples.index(s1)
                d = freq[keep, i1] - freq[keep, i0]
                deltas[r, ti, :] = d
                if bias_correct:
                    p0 = freq[keep, i0]; p1 = freq[keep, i1]
                    d0 = total[keep, i0]; d1 = total[keep, i1]
                    het0 = p0*(1-p0); het1 = p1*(1-p1)
                    noise_var[r, ti, :] = (het0*(1/(2*n_eff) + 1/d0)
                                            + het1*(1/(2*n_eff) + 1/d1))


        cc_mat = np.full((n_int, n_int), np.nan)
        for t in range(n_int):
            for s in range(n_int):
                num_vals = []
                denom_vals = []
                for a in range(R):
                    for b in range(a+1, R):
                        x = deltas[a, t, :]; y = deltas[b, s, :]
                        x0 = x - x.mean(); y0 = y - y.mean()
                        num_vals.append(np.sum(x0*y0) / (len(x0)-1))
                        if bias_correct:
                            var_a = max(np.sum(x0**2)/(len(x0)-1)
                                        - noise_var[a, t, :].mean(), 1e-12)
                            var_b = max(np.sum(y0**2)/(len(y0)-1)
                                        - noise_var[b, s, :].mean(), 1e-12)
                        else:
                            var_a = np.sum(x0**2)/(len(x0)-1)
                            var_b = np.sum(y0**2)/(len(y0)-1)
                        denom_vals.append(np.sqrt(var_a * var_b))
                num = float(np.mean(num_vals))
                denom = float(np.mean(denom_vals))
                cc_mat[t, s] = num / denom if denom > 0 else np.nan

        if mode == "diagonal":
            scalar = float(np.nanmean(np.diag(cc_mat)))
        else:
            scalar = float(np.nanmean(cc_mat))
        out[trt] = dict(cc_scalar=scalar, cc_matrix=cc_mat)
    return out


def main():
    samples, freq, total = load_data()
    mask = build_common_filter(samples, freq, total)

    # 1 and 2: total-Δp cc
    raw_tot = cc_total_delta(samples, freq, total, mask, bias_correct=False)
    cor_tot = cc_total_delta(samples, freq, total, mask, bias_correct=True)
    # 3 and 4: per-interval cc
    intv_diag = cc_per_interval(samples, freq, total, mask,
                                bias_correct=True, mode="diagonal")
    intv_full = cc_per_interval(samples, freq, total, mask,
                                bias_correct=True, mode="full")

    for trt in TREATMENTS:
        print(f"{trt:<10} {raw_tot[trt]['cc']:15.4f} "
              f"{cor_tot[trt]['cc']:16.4f} "
              f"{intv_diag[trt]['cc_scalar']:18.4f} "
              f"{intv_full[trt]['cc_scalar']:18.4f}")

    for trt in TREATMENTS:
        print(f"\n  {trt}:")
        m = intv_diag[trt]['cc_matrix']
        print("       " + "".join(f"{i:>8d}" for i in range(m.shape[0])))
        for i in range(m.shape[0]):
            row = "  " + f"{i:>4d}  "
            for j in range(m.shape[1]):
                v = m[i, j]
                row += f" {v:7.4f}" if np.isfinite(v) else "    nan "
            print(row)

    df = pd.DataFrame([
        dict(treatment=t,
             cc_total_raw=raw_tot[t]['cc'],
             cc_total_corrected=cor_tot[t]['cc'],
             cc_interval_diag=intv_diag[t]['cc_scalar'],
             cc_interval_full=intv_full[t]['cc_scalar'])
        for t in TREATMENTS
    ])
    df.to_csv("variance_analysis/section1_rigorous/cc_audit.tsv",
              sep="\t", index=False)


if __name__ == "__main__":
    main()
