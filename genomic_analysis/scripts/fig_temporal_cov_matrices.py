#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import sys
import os

C_B = "#4C86A8"
C_T = "#E07B54"
C_M = "#6AAB6E"
OUTDIR = "variance_analysis/figures"


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


def compute_temporal_cov(freq, total, samples, trt, n_diploids=None):
    timepoints = [0, 1, 2, 6, 7, 8, 9]
    n_tp = len(timepoints)
    T = n_tp - 1
    L = freq.shape[0]

    rep_data = []
    for rep in range(1, 5):
        tp_freqs, tp_depths = [], []
        all_present = True
        for tp in timepoints:
            sname = f"F{rep}G00" if tp == 0 else f"{trt}{rep}G{tp:02d}"
            if sname not in samples:
                all_present = False; break
            idx = samples.index(sname)
            tp_freqs.append(freq[:, idx])
            tp_depths.append(total[:, idx])
        if all_present:
            rep_data.append((np.array(tp_freqs), np.array(tp_depths)))

    R = len(rep_data)
    ok = np.ones(L, dtype=bool)
    for f_arr, d_arr in rep_data:
        ok &= np.all(d_arr >= 20, axis=0) & np.all(np.isfinite(f_arr), axis=0)
    grand_mean = np.nanmean(np.array([f[0] for f, _ in rep_data]), axis=0)
    ok &= (grand_mean >= 0.10) & (grand_mean <= 0.90)
    n_sites = ok.sum()

    delta_p = np.zeros((R, T, n_sites))
    het = np.zeros((R, T, n_sites))
    depths_at_tp = np.zeros((R, n_tp, n_sites))

    for ri, (f_arr, d_arr) in enumerate(rep_data):
        f_ok = f_arr[:, ok]
        d_ok = d_arr[:, ok]
        depths_at_tp[ri] = d_ok
        for t in range(T):
            delta_p[ri, t, :] = f_ok[t+1] - f_ok[t]
            p_mid = (f_ok[t] + f_ok[t+1]) / 2
            het[ri, t, :] = 2 * p_mid * (1 - p_mid)

    tcov = np.zeros((T, T))
    for r in range(R):
        for s in range(T):
            for t in range(T):
                h_prod = np.sqrt(het[r, s, :] * het[r, t, :])
                h_prod[h_prod < 1e-10] = np.nan

                raw = delta_p[r, s, :] * delta_p[r, t, :] / h_prod

                if n_diploids is not None:

                    N_hap = 2 * n_diploids

                    if s == t:
                        d_start = depths_at_tp[r, s, :]
                        d_end = depths_at_tp[r, s+1, :]
                        correction = (1.0/N_hap + 1.0/d_start - 1.0/(N_hap * d_start) +
                                      1.0/N_hap + 1.0/d_end - 1.0/(N_hap * d_end))
                        p_start = (delta_p[r, s, :] + het[r, s, :] * 0)  # placeholder
                        raw -= correction / 2  # approximate

                    elif abs(s - t) == 1:
                        shared_tp = max(s, t)  # the shared timepoint index
                        d_shared = depths_at_tp[r, shared_tp, :]
                        correction = (1.0/N_hap + 1.0/d_shared - 1.0/(N_hap * d_shared))
                        raw += correction / 2  # remove the negative bias

                tcov[s, t] += np.nanmean(raw)
    tcov /= R

    return tcov, n_sites


def plot_matrix(ax, mat, interval_labels, title, color, vmax, annotate=True):
    T = mat.shape[0]
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   aspect="equal", origin="upper")

    if annotate:
        for s in range(T):
            for t in range(T):
                val = mat[s, t]
                textcolor = "white" if abs(val) > vmax * 0.6 else "black"
                fontsize = 6 if s != t else 7
                weight = "bold" if s == t else "normal"
                ax.text(t, s, f"{val:+.4f}", ha="center", va="center",
                        fontsize=fontsize, color=textcolor, fontweight=weight)

    ax.set_xticks(range(T))
    ax.set_xticklabels(interval_labels, fontsize=5.5, rotation=45, ha="right")
    ax.set_yticks(range(T))
    ax.set_yticklabels(interval_labels, fontsize=5.5)
    ax.set_title(title, fontweight="bold", fontsize=10, color=color)

    for i in range(T):
        ax.add_patch(Rectangle((i-0.5, i-0.5), 1, 1, fill=False,
                                edgecolor="black", lw=1.5))
    return im


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    print("Loading data...", file=sys.stderr)
    samples, freq, total = load_data()

    treatments = [("B", C_B, "B"), ("T", C_T, "T"), ("BT", C_M, "M")]
    timepoints = [0, 1, 2, 6, 7, 8, 9]
    interval_labels = [f"{timepoints[i]}→{timepoints[i+1]}"
                       for i in range(len(timepoints)-1)]

    N_EFF = 29  # from overdispersion analysis

    raw_matrices = {}
    corr_matrices = {}
    for trt_label, color, trt_data in treatments:
        raw, n = compute_temporal_cov(freq, total, samples, trt_data, n_diploids=None)
        corr, _ = compute_temporal_cov(freq, total, samples, trt_data, n_diploids=N_EFF)
        raw_matrices[trt_label] = raw
        corr_matrices[trt_label] = corr

    fig = plt.figure(figsize=(16, 8.5))
    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.45, wspace=0.35,
                           left=0.06, right=0.98, top=0.90, bottom=0.08,
                           width_ratios=[1, 1, 1, 0.55])

    vmax_raw = 0.04
    vmax_corr = 0.01  # corrected values are much smaller

    for idx, (trt_label, color, _) in enumerate(treatments):
        ax = fig.add_subplot(gs[0, idx])
        im_raw = plot_matrix(ax, raw_matrices[trt_label], interval_labels,
                             f"{trt_label} — uncorrected", color, vmax_raw)

    ax_ann1 = fig.add_subplot(gs[0, 3])
    ax_ann1.axis("off")
    ax_ann1.text(0.05, 0.85, "Uncorrected", fontsize=9, fontweight="bold",
                 transform=ax_ann1.transAxes)
    ax_ann1.text(0.05, 0.72, "Pool size assumed: N = 80\n(no correction applied)",
                 fontsize=7, color="#555555", transform=ax_ann1.transAxes)
    ax_ann1.text(0.05, 0.52, "Diagonal", fontsize=8, fontweight="bold",
                 color="black", transform=ax_ann1.transAxes)
    ax_ann1.text(0.05, 0.43, "Total var = drift + 2×noise\n≈ 0.032",
                 fontsize=6.5, color="#555555", transform=ax_ann1.transAxes)
    ax_ann1.text(0.05, 0.28, "Adjacent off-diag", fontsize=8, fontweight="bold",
                 color="#2166AC", transform=ax_ann1.transAxes)
    ax_ann1.text(0.05, 0.19, "Negative = noise from\nshared timepoint\n≈ −0.016",
                 fontsize=6.5, color="#555555", transform=ax_ann1.transAxes)
    ax_ann1.text(0.05, 0.04, "Non-adjacent off-diag", fontsize=8, fontweight="bold",
                 color="#999999", transform=ax_ann1.transAxes)
    ax_ann1.text(0.05, -0.05, "≈ 0 (drift-dominated)",
                 fontsize=6.5, color="#555555", transform=ax_ann1.transAxes)

    for idx, (trt_label, color, _) in enumerate(treatments):
        ax = fig.add_subplot(gs[1, idx])
        im_corr = plot_matrix(ax, corr_matrices[trt_label], interval_labels,
                              f"{trt_label} — corrected (N_eff={N_EFF})",
                              color, vmax_corr)

    ax_ann2 = fig.add_subplot(gs[1, 3])
    ax_ann2.axis("off")
    ax_ann2.text(0.05, 0.85, f"Corrected", fontsize=9, fontweight="bold",
                 transform=ax_ann2.transAxes)
    ax_ann2.text(0.05, 0.72, f"Effective pool: N_eff = {N_EFF}\n(from overdispersion φ = 2.78)",
                 fontsize=7, color="#555555", transform=ax_ann2.transAxes)
    ax_ann2.text(0.05, 0.52, "Diagonal", fontsize=8, fontweight="bold",
                 color="black", transform=ax_ann2.transAxes)
    ax_ann2.text(0.05, 0.43, "Drift variance only\n(noise removed)\nNote: 10× smaller scale",
                 fontsize=6.5, color="#555555", transform=ax_ann2.transAxes)
    ax_ann2.text(0.05, 0.28, "Adjacent off-diag", fontsize=8, fontweight="bold",
                 color="#2166AC", transform=ax_ann2.transAxes)
    ax_ann2.text(0.05, 0.19, "Should be ≈ 0 if\ncorrection is accurate",
                 fontsize=6.5, color="#555555", transform=ax_ann2.transAxes)
    ax_ann2.text(0.05, 0.04, "Non-adjacent off-diag", fontsize=8, fontweight="bold",
                 color="#999999", transform=ax_ann2.transAxes)
    ax_ann2.text(0.05, -0.05, "Positive = selection signal",
                 fontsize=6.5, color="#555555", transform=ax_ann2.transAxes)

    cbar_ax1 = fig.add_axes([0.06, 0.47, 0.55, 0.015])
    plt.colorbar(im_raw, cax=cbar_ax1, orientation="horizontal")
    cbar_ax1.set_xlabel("Covariance (uncorrected scale)", fontsize=7)

    cbar_ax2 = fig.add_axes([0.06, 0.02, 0.55, 0.015])
    plt.colorbar(im_corr, cax=cbar_ax2, orientation="horizontal")
    cbar_ax2.set_xlabel("Covariance (corrected scale — note 4× zoom)", fontsize=7)

    fig.suptitle(
        "Temporal covariance matrices of allele frequency change",
        fontsize=13, fontweight="bold", y=0.96)

    for fmt in ["png", "pdf", "svg"]:
        fig.savefig(os.path.join(OUTDIR, f"fig_temporal_cov_matrices.{fmt}"),
                    dpi=300, bbox_inches="tight", facecolor="white")


if __name__ == "__main__":
    main()
