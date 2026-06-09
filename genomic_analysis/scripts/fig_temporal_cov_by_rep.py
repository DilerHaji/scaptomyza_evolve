#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import TwoSlopeNorm
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


def compute_per_rep_cov(freq, total, samples, trt):
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
            rep_data.append((rep, np.array(tp_freqs), np.array(tp_depths)))

    ok = np.ones(L, dtype=bool)
    for _, f_arr, d_arr in rep_data:
        ok &= np.all(d_arr >= 20, axis=0) & np.all(np.isfinite(f_arr), axis=0)
    grand_mean = np.nanmean(np.array([f[0] for _, f, _ in rep_data]), axis=0)
    ok &= (grand_mean >= 0.10) & (grand_mean <= 0.90)

    results = []
    for rep_num, f_arr, d_arr in rep_data:
        f_ok = f_arr[:, ok]
        delta_p = np.zeros((T, ok.sum()))
        het = np.zeros((T, ok.sum()))
        for t in range(T):
            delta_p[t, :] = f_ok[t+1] - f_ok[t]
            p_mid = (f_ok[t] + f_ok[t+1]) / 2
            het[t, :] = 2 * p_mid * (1 - p_mid)

        tcov = np.zeros((T, T))
        for s in range(T):
            for t in range(T):
                h_prod = np.sqrt(het[s, :] * het[t, :])
                h_prod[h_prod < 1e-10] = np.nan
                tcov[s, t] = np.nanmean(delta_p[s, :] * delta_p[t, :] / h_prod)
        results.append((rep_num, tcov))

    return results, ok.sum()


def plot_rep_matrix(ax, cov, interval_labels, title, show_ylabels=True):
    T = cov.shape[0]
    display = np.full((T, T), np.nan)
    grey_mask = np.zeros((T, T), dtype=bool)

    for s in range(T):
        for t in range(T):
            if s == t or abs(s - t) == 1:
                grey_mask[s, t] = True
            else:
                display[s, t] = cov[s, t]

    vmax_nonadj = 0.003
    im = ax.imshow(display, cmap="RdBu_r", vmin=-vmax_nonadj, vmax=vmax_nonadj,
                   aspect="equal", origin="upper")

    grey_display = np.full((T, T), np.nan)
    for s in range(T):
        for t in range(T):
            if grey_mask[s, t]:
                grey_display[s, t] = 0
    ax.imshow(grey_display, cmap="Greys", vmin=-1, vmax=1,
              aspect="equal", origin="upper", alpha=0.15)

    for s in range(T):
        for t in range(T):
            val = cov[s, t]
            if s == t:
                ax.text(t, s, f"{val:+.3f}", ha="center", va="center",
                        fontsize=5.5, color="black", fontweight="bold")
            elif abs(s - t) == 1:
                ax.text(t, s, f"{val:+.3f}", ha="center", va="center",
                        fontsize=5, color="#666666")
            else:
                textcolor = "black"
                if abs(val) > 0.002:
                    textcolor = "#B2182B" if val > 0 else "#2166AC"
                ax.text(t, s, f"{val:+.4f}", ha="center", va="center",
                        fontsize=5, color=textcolor, fontweight="bold")

    ax.set_xticks(range(T))
    ax.set_xticklabels(interval_labels, fontsize=5, rotation=45, ha="right")
    ax.set_yticks(range(T))
    ax.set_yticklabels(interval_labels if show_ylabels else [], fontsize=5)
    ax.set_title(title, fontsize=9, fontweight="bold")

    for i in range(T):
        ax.add_patch(Rectangle((i-0.5, i-0.5), 1, 1, fill=False,
                                edgecolor="black", lw=1.5))
    for i in range(T-1):
        ax.add_patch(Rectangle((i+0.5, i-0.5), 1, 1, fill=False,
                                edgecolor="#888888", lw=0.8, ls="--"))
        ax.add_patch(Rectangle((i-0.5, i+0.5), 1, 1, fill=False,
                                edgecolor="#888888", lw=0.8, ls="--"))

    return im


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    samples, freq, total = load_data()

    treatments = [("B", C_B, "B"), ("T", C_T, "T"), ("BT", C_M, "M")]
    timepoints = [0, 1, 2, 6, 7, 8, 9]
    interval_labels = [f"G{timepoints[i]:02d}→G{timepoints[i+1]:02d}"
                       for i in range(len(timepoints)-1)]
    T = len(interval_labels)

    for trt_label, color, trt_data in treatments:
        rep_covs, n_sites = compute_per_rep_cov(freq, total, samples, trt_data)
        mean_cov = np.mean([c for _, c in rep_covs], axis=0)

        fig, axes = plt.subplots(1, 5, figsize=(20, 4))
        fig.subplots_adjust(wspace=0.25, left=0.04, right=0.92, top=0.80, bottom=0.18)

        for idx, (rep_num, cov) in enumerate(rep_covs):
            diag = np.mean(np.diag(cov))
            adj = np.mean([cov[i, i+1] for i in range(T-1)])
            nonadj = [cov[s, t] for s in range(T) for t in range(s+2, T)]
            nonadj_mean = np.mean(nonadj)

            title = (f"Rep {rep_num}\n"
                     f"diag={diag:.4f}  adj={adj:+.4f}  non-adj={nonadj_mean:+.5f}")
            im = plot_rep_matrix(axes[idx], cov, interval_labels, title,
                                 show_ylabels=(idx == 0))

        diag = np.mean(np.diag(mean_cov))
        adj = np.mean([mean_cov[i, i+1] for i in range(T-1)])
        nonadj = [mean_cov[s, t] for s in range(T) for t in range(s+2, T)]
        nonadj_mean = np.mean(nonadj)

        title = (f"Mean\n"
                 f"diag={diag:.4f}  adj={adj:+.4f}  non-adj={nonadj_mean:+.5f}")
        im = plot_rep_matrix(axes[4], mean_cov, interval_labels, title,
                             show_ylabels=False)

        cbar_ax = fig.add_axes([0.93, 0.25, 0.015, 0.45])
        cb = plt.colorbar(im, cax=cbar_ax)
        cb.set_label("Non-adjacent covariance\n(selection signal)", fontsize=7)

        fig.text(0.93, 0.78, "Bold black = diagonal\n(total variance)",
                 fontsize=6, va="top")
        fig.text(0.93, 0.18, "Grey = adjacent\n(measurement noise)\n\n"
                 "Colored = non-adjacent\n(red + = selection\nblue − = overcorrection)",
                 fontsize=6, va="top")

        fig.suptitle(f"Treatment {trt_label} — temporal covariance by replicate "
                     f"({n_sites:,} sites)",
                     fontsize=12, fontweight="bold", color=color, y=0.95)

        for fmt in ["png", "pdf", "svg"]:
            fig.savefig(os.path.join(OUTDIR, f"fig_tcov_by_rep_{trt_label}.{fmt}"),
                        dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)


if __name__ == "__main__":
    main()
