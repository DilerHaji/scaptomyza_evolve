#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from collections import Counter
import sys
import os

C_B = "#4C86A8"
C_T = "#E07B54"
C_M = "#6AAB6E"
C_F = "#8B6BAE"
C_DRIFT = "#999999"
C_SIM = "#333333"
OUTDIR = "variance_analysis/figures"


def load_freq_data():
    samples = [l.strip() for l in open("variance_analysis/sample_list.txt") if l.strip()]
    ad_ref, ad_alt, chroms, positions = [], [], [], []
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


def compute_F(freq, total, rep_indices):
    k = len(rep_indices)
    if k < 2:
        return np.nan, np.nan
    f_r = freq[:, rep_indices]
    d_r = total[:, rep_indices]
    mask = (np.all(d_r >= 10, axis=1) & np.all(np.isfinite(f_r), axis=1))
    pbar = np.nanmean(f_r, axis=1)
    maf_ok = np.minimum(pbar, 1 - pbar) >= 0.05
    mask = mask & maf_ok
    if mask.sum() < 1000:
        return np.nan, np.nan
    pbar_v = pbar[mask]
    f_v = f_r[mask]
    denom = pbar_v * (1 - pbar_v)
    denom[denom < 1e-10] = np.nan
    ss = np.nansum((f_v - pbar_v[:, None])**2, axis=1)
    F_site = (ss / (k - 1)) / denom
    return np.nanmean(F_site), np.nanstd(F_site) / np.sqrt(mask.sum())


def compute_autocorrelation(delta, pos, dist_bins, n_sample=1000, seed=42):
    ok = np.isfinite(delta)
    d_ok = delta[ok]
    p_ok = pos[ok]
    var_d = np.nanvar(d_ok)
    if var_d < 1e-15 or len(d_ok) < 100:
        return [np.nan] * len(dist_bins)
    n_s = len(d_ok)
    rng = np.random.RandomState(seed)
    samp = rng.choice(n_s, min(n_sample, n_s), replace=False)
    r_vals = []
    for lo, hi in dist_bins:
        prods = []
        for si in samp:
            dists = np.abs(p_ok - p_ok[si])
            in_range = (dists >= lo) & (dists < hi) & (np.arange(n_s) != si)
            if in_range.sum() > 0:
                prods.append(np.mean(d_ok[si] * d_ok[in_range]))
        r_vals.append(np.mean(prods) / var_d if prods else np.nan)
    return r_vals


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    samples, chroms, positions, freq, total = load_freq_data()
    trt_colors = {"B": C_B, "T": C_T, "M": C_M}

    sim_df = pd.read_csv("variance_analysis/drift_null/drift_null_results.tsv", sep="\t")
    sim_phi278 = sim_df[sim_df["phi"] == 2.78]

    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.35,
                           left=0.06, right=0.97, top=0.93, bottom=0.08)

    ax_a = fig.add_subplot(gs[0, 0])

    f_data = {}
    for trt in ["B", "T", "M"]:
        gens, Fs, SEs = [], [], []
        for gen in range(1, 10):
            rep_names = [f"{trt}{r}G{gen:02d}" for r in range(1, 5)]
            rep_idx = [samples.index(s) for s in rep_names if s in samples]
            F_val, F_se = compute_F(freq, total, rep_idx)
            if not np.isnan(F_val):
                gens.append(gen); Fs.append(F_val); SEs.append(F_se)
        f_data[trt] = (gens, Fs, SEs)
        ax_a.errorbar(gens, Fs, yerr=SEs, color=trt_colors[trt], marker="o",
                      ms=4, lw=1.5, capsize=2, label=trt)

    slopes_4rep = {}
    for trt in ["B", "T", "M"]:
        g, f, _ = f_data[trt]
        g4 = [gi for gi, fi in zip(g, f) if gi not in [3, 4, 5]]
        f4 = [fi for gi, fi in zip(g, f) if gi not in [3, 4, 5]]
        sl, ic, _, _, _ = stats.linregress(g4, f4)
        slopes_4rep[trt] = (sl, ic)
    mean_ic = np.mean([v[1] for v in slopes_4rep.values()])
    x_d = np.linspace(0.5, 9.5, 50)
    ax_a.plot(x_d, mean_ic + x_d / (2 * 500), color=C_DRIFT, ls=":", lw=2,
              label="Drift (Ne=500)")

    ax_a.set_xlabel("Generation")
    ax_a.set_ylabel("Between-replicate F")
    ax_a.set_title("A. Between-replicate divergence", fontweight="bold",
                   loc="left", fontsize=9)
    ax_a.legend(fontsize=6, loc="upper left")
    ax_a.set_xlim(0.5, 9.5)
    ax_a.set_xticks(range(1, 10))

    ax_b = fig.add_subplot(gs[0, 1])

    for trt in ["B", "T", "M"]:
        het_by_gen = {}
        for gen in range(1, 10):
            H_vals = []
            for rep in range(1, 5):
                sname = f"{trt}{rep}G{gen:02d}"
                if sname not in samples:
                    continue
                idx = samples.index(sname)
                mask = (total[:, idx] >= 10) & np.isfinite(freq[:, idx])
                f_i = freq[mask, idx]
                d_i = total[mask, idx].astype(float)
                H_vals.append(np.mean(2 * f_i * (1 - f_i) * d_i / (d_i - 1)))
            if H_vals:
                het_by_gen[gen] = (np.mean(H_vals), np.std(H_vals) / np.sqrt(len(H_vals)))
                for h in H_vals:
                    ax_b.plot(gen, h, "o", color=trt_colors[trt], ms=2, alpha=0.15)

        gens = sorted(het_by_gen.keys())
        means = [het_by_gen[g][0] for g in gens]
        ses = [het_by_gen[g][1] for g in gens]
        ax_b.errorbar(gens, means, yerr=ses, color=trt_colors[trt], marker="o",
                      ms=4, lw=1.5, capsize=2, label=trt)

    ax_b.set_xlabel("Generation")
    ax_b.set_ylabel("Depth-corrected H")
    ax_b.set_title("B. Within-replicate diversity", fontweight="bold",
                   loc="left", fontsize=9)
    ax_b.legend(fontsize=6, loc="lower left")
    ax_b.set_xlim(0.5, 9.5)
    ax_b.set_xticks(range(1, 10))

    ax_c = fig.add_subplot(gs[0, 2])

    ne_true = sim_phi278["Ne_true"].values
    ne_between = sim_phi278["Ne_between_median"].values
    ne_between_lo = sim_phi278["Ne_between_25"].values
    ne_between_hi = sim_phi278["Ne_between_75"].values

    ax_c.fill_between(ne_true, ne_between_lo, ne_between_hi,
                      color=C_B, alpha=0.2)
    ax_c.plot(ne_true, ne_between, "o-", color=C_B, lw=2, ms=6,
              label="Between-rep Ne (simulated)")
    ax_c.plot([100, 1100], [100, 1100], "k--", lw=1, label="1:1 (unbiased)")
    ax_c.set_xlabel("True Ne (simulated)")
    ax_c.set_ylabel("Estimated Ne")
    ax_c.set_title("C. Between-rep estimator\n(unbiased under drift)",
                   fontweight="bold", loc="left", fontsize=9)
    ax_c.legend(fontsize=6)
    ax_c.set_xlim(100, 1100)
    ax_c.set_ylim(100, 1100)

    ax_d = fig.add_subplot(gs[1, 0])

    ne_temporal = sim_phi278["Ne_temporal_median"].values
    ne_temporal_lo = sim_phi278["Ne_temporal_25"].values
    ne_temporal_hi = sim_phi278["Ne_temporal_75"].values

    ax_d.fill_between(ne_true, ne_temporal_lo, ne_temporal_hi,
                      color=C_T, alpha=0.2)
    ax_d.plot(ne_true, ne_temporal, "o-", color=C_T, lw=2, ms=6,
              label="Temporal Ne (simulated)")
    ax_d.plot([100, 1100], [100, 1100], "k--", lw=1, label="1:1")
    ax_d.plot([100, 1100], [130, 1430], color=C_T, ls=":", lw=1,
              label="1.3× bias")

    ax_d.axhline(665, color=C_DRIFT, ls="-.", lw=1, alpha=0.5)
    ax_d.text(150, 680, "Observed temporal Ne = 665", fontsize=6, color=C_DRIFT)
    ax_d.annotate("", xy=(510, 665), xytext=(510, 620),
                  arrowprops=dict(arrowstyle="->", color=C_DRIFT, lw=1))
    ax_d.text(520, 590, "→ true Ne ≈ 510", fontsize=6, color=C_DRIFT)

    ax_d.set_xlabel("True Ne (simulated)")
    ax_d.set_ylabel("Estimated Ne")
    ax_d.set_title("D. Temporal estimator\n(1.3× upward bias under drift)",
                   fontweight="bold", loc="left", fontsize=9)
    ax_d.legend(fontsize=6)
    ax_d.set_xlim(100, 1100)
    ax_d.set_ylim(100, 1500)

    ax_e = fig.add_subplot(gs[1, 1])

    null_ne = 507
    null_lo = 492
    null_hi = 523

    observed_ne = {"B": 211, "T": 440, "M": 170}
    x_pos = np.arange(3)
    trt_list = ["B", "T", "M"]

    bars = ax_e.bar(x_pos, [observed_ne[t] for t in trt_list],
                    color=[C_B, C_T, C_M], alpha=0.8, edgecolor="black", lw=0.8)

    ax_e.axhspan(null_lo, null_hi, color=C_DRIFT, alpha=0.2, zorder=0)
    ax_e.axhline(null_ne, color=C_DRIFT, ls="--", lw=1.5,
                 label=f"Drift null (Ne=500)\n[{null_lo}, {null_hi}]")
    ax_e.axhline(500, color="black", ls=":", lw=1, alpha=0.4, label="Census Ne=500")

    ax_e.set_xticks(x_pos)
    ax_e.set_xticklabels(trt_list, fontsize=10)
    ax_e.set_ylabel("Between-replicate Ne")
    ax_e.set_title("E. Observed Ne vs. drift null",
                   fontweight="bold", loc="left", fontsize=9)
    ax_e.legend(fontsize=6, loc="upper right")

    for i, trt in enumerate(trt_list):
        ne = observed_ne[trt]
        if ne < null_lo:
            ax_e.text(i, ne - 20, f"Ne={ne}\n(< drift)", ha="center",
                      fontsize=6.5, fontweight="bold", color="black")
        else:
            ax_e.text(i, ne + 15, f"Ne={ne}\n(≈ drift)", ha="center",
                      fontsize=6.5, color=C_DRIFT)

    ax_f = fig.add_subplot(gs[1, 2])

    top4 = [c for c, _ in Counter(chroms).most_common(4)]
    dist_bins = [
        (0, 200), (200, 500), (500, 1000), (1000, 2000), (2000, 5000),
        (5000, 10000), (10000, 25000), (25000, 50000), (50000, 100000),
        (100000, 250000), (250000, 500000), (500000, 1000000),
        (1000000, 2500000), (2500000, 5000000), (5000000, 10000000),
    ]
    dist_mids = [(lo + hi) / 2 for lo, hi in dist_bins]

    for trt in ["B", "T", "M"]:
        deltas = []
        for rep in [1, 2, 3, 4]:
            s1, s9 = f"{trt}{rep}G01", f"{trt}{rep}G09"
            if s1 in samples and s9 in samples:
                i1, i9 = samples.index(s1), samples.index(s9)
                d = freq[:, i9] - freq[:, i1]
                bad = (total[:, i1] < 15) | (total[:, i9] < 15) | ~np.isfinite(d)
                d[bad] = np.nan
                deltas.append(d)
        mean_d = np.nanmean(deltas, axis=0)

        scaffold_rs = []
        for sc in top4:
            sc_mask = chroms == sc
            r_v = compute_autocorrelation(mean_d[sc_mask], positions[sc_mask],
                                           dist_bins, seed=42 + hash(sc) % 1000)
            ax_f.semilogx(dist_mids, r_v, color="#CCCCCC", lw=0.5, alpha=0.4, zorder=1)
            scaffold_rs.append(r_v)

        mean_r = np.nanmean(scaffold_rs, axis=0)
        ax_f.semilogx(dist_mids, mean_r, color=trt_colors[trt], marker="o",
                      ms=3, lw=1.5, label=trt, zorder=3)

    ax_f.axhline(0, color="grey", ls="--", lw=0.8, alpha=0.5)
    ax_f.set_xlabel("Physical distance (bp)")
    ax_f.set_ylabel("Autocorrelation of ΔAF")
    ax_f.set_title("F. Linkage of AF changes",
                   fontweight="bold", loc="left", fontsize=9)
    ax_f.legend(fontsize=6, loc="upper right")
    ax_f.set_xlim(100, 15000000)
    ax_f.set_ylim(bottom=-0.01)

    fig.suptitle(
        "Supplementary Figure X. Variance decomposition, drift estimation, and linkage structure",
        fontsize=11, fontweight="bold", y=0.98)

    fig.savefig(os.path.join(OUTDIR, "fig_supp_variance_final.png"),
                dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(os.path.join(OUTDIR, "fig_supp_variance_final.pdf"),
                bbox_inches="tight", facecolor="white")
    fig.savefig(os.path.join(OUTDIR, "fig_supp_variance_final.svg"),
                bbox_inches="tight", facecolor="white")

if __name__ == "__main__":
    main()
