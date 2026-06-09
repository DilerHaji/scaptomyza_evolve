#!/usr/bin/env python3

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

ROOT = Path(".")
OMEGA_FILE = ROOT / "baypass_wild/wild_omega_mat_omega.out"
BETAI_FILE = ROOT / "baypass_wild/wild_trt_summary_betai_reg.out"
XTX_FILE = ROOT / "baypass_wild/wild_omega_summary_pi_xtx.out"
POS_FILE = ROOT / "baypass_wild/wild_snp_positions.csv"
OUT_PNG = ROOT / "final_plots/wild/wild_baypass_twopanel.png"
OUT_SVG = ROOT / "final_plots/wild/wild_baypass_twopanel.svg"

POOL_ORDER = ["AVB", "AVT", "PSB", "PST", "RMB", "RMT"]
HOST_OF = {"AVB": "B", "AVT": "T", "PSB": "B", "PST": "T", "RMB": "B", "RMT": "T"}

# Colors matched to Figure 1 / lab-pool scheme:
C_B = "#499FFF"       # blue  — B host
C_T = "#EDB72D"       # gold  — T host
C_DECISIVE = "#CC3333"   # red  — Jeffreys-decisive threshold
C_GREY_DARK = "#555555"
C_GREY_LIGHT = "#BBBBBB"
BF_THRESHOLD = 10.0
DF_CHI2 = 6               # n_populations — chi^2 df for X^tX null

N_TOP_SCAFFOLDS = 10
N_TOP_SNPS = 5
N_DENSITY_SCAFFOLDS = 4
BG_SUBSAMPLE = 60_000
RNG_SEED = 20260417


def scaffold_short(name: str) -> str:
    m = re.search(r"(\d+)_HRSCAF", str(name))
    return m.group(1) if m else str(name)


def load_omega() -> np.ndarray:
    return np.loadtxt(OMEGA_FILE)


def load_scan() -> pd.DataFrame:
    betai = pd.read_csv(BETAI_FILE, sep=r"\s+")
    pos = pd.read_csv(POS_FILE)
    df = pos.merge(
        betai[["MRK", "BF(dB)", "Beta_is"]].rename(columns={"MRK": "mrk"}),
        on="mrk",
    )
    df = df.rename(columns={"BF(dB)": "bf_db", "Beta_is": "beta_i"})
    return df


def load_xtx() -> np.ndarray:
    df = pd.read_csv(XTX_FILE, sep=r"\s+")
    return df["XtXst"].values

def plot_omega(ax_heat: plt.Axes, omega: np.ndarray) -> None:
    sd = np.sqrt(np.diag(omega))
    cor = omega / np.outer(sd, sd)

    dist = 1 - np.abs(cor)
    np.fill_diagonal(dist, 0)
    Z = linkage(squareform(dist, checks=False), method="average")
    leaf_order = dendrogram(Z, no_plot=True)["leaves"]
    cor_reordered = cor[np.ix_(leaf_order, leaf_order)]
    labels = [POOL_ORDER[i] for i in leaf_order]

    im = ax_heat.imshow(cor_reordered, cmap=plt.cm.Greys, vmin=0.0, vmax=1.0,
                        aspect="equal")

    n = len(labels)
    for i in range(n):
        for j in range(n):
            val = cor_reordered[i, j]
            tc = "white" if val > 0.55 else "black"
            ax_heat.text(j, i, f"{val:.2f}", ha="center", va="center",
                         fontsize=7.5, color=tc)

    ax_heat.set_xticks(range(n))
    ax_heat.set_yticks(range(n))
    ax_heat.set_xticklabels(labels, rotation=0)
    ax_heat.set_yticklabels(labels)
    for tick, lab in zip(ax_heat.get_xticklabels(), labels):
        tick.set_color(C_B if HOST_OF[lab] == "B" else C_T)
        tick.set_fontweight("bold")
    for tick, lab in zip(ax_heat.get_yticklabels(), labels):
        tick.set_color(C_B if HOST_OF[lab] == "B" else C_T)
        tick.set_fontweight("bold")

    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    cbar = ax_heat.figure.colorbar(im, ax=ax_heat, shrink=0.72, pad=0.03)
    cbar.set_label(r"Corr from $\mathbf{\Omega}$", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax_heat.set_title("A", loc="left", fontsize=16, fontweight="bold", pad=6)


def plot_calibration(ax: plt.Axes, xtxst: np.ndarray) -> None:
    sorted_obs = np.sort(xtxst)
    n = len(sorted_obs)
    pp = (np.arange(1, n + 1) - 0.5) / n
    theoretical = stats.chi2.ppf(pp, df=DF_CHI2)

    qmax = 25.0
    qmin = -0.5

    ax.scatter(theoretical, sorted_obs, s=3, c=C_GREY_DARK,
               alpha=0.4, linewidths=0, rasterized=True)
    ax.plot([qmin, qmax], [qmin, qmax], color="black", linewidth=1.2,
            linestyle="--", zorder=3, label="null expectation (y = x)")

    q95 = stats.chi2.ppf(0.95, df=DF_CHI2)
    q99 = stats.chi2.ppf(0.99, df=DF_CHI2)
    for q, txt in [(q95, "5%"), (q99, "1%")]:
        ax.axvline(q, color=C_DECISIVE, linestyle=":", linewidth=0.7, alpha=0.7, zorder=2)
        ax.text(q + 0.2, qmax * 0.04, txt, color=C_DECISIVE, fontsize=7,
                ha="left", va="bottom")

    pct_above_q95 = 100 * np.mean(xtxst > q95)
    pct_above_q99 = 100 * np.mean(xtxst > q99)
    observed_median = float(np.median(xtxst))
    expected_median = stats.chi2.ppf(0.5, df=DF_CHI2)
    stats_text = (
        f"observed median = {observed_median:.2f}\n"
        f"$\\chi^2_{{{DF_CHI2}}}$ median = {expected_median:.2f}\n"
        f"{pct_above_q95:.1f}% above 5% tail\n"
        f"{pct_above_q99:.1f}% above 1% tail"
    )
    ax.text(
        0.03, 0.97, stats_text,
        transform=ax.transAxes, ha="left", va="top", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3",
                  facecolor="white", edgecolor=C_GREY_DARK, linewidth=0.6),
    )

    ax.set_xlabel(fr"$\chi^2_{{{DF_CHI2}}}$ quantile (expected)", fontsize=10)
    ax.set_ylabel(r"$X^tX_{st}$ quantile (observed)", fontsize=10)
    ax.set_xlim(qmin, qmax)
    ax.set_ylim(qmin, qmax)
    ax.tick_params(labelsize=8)
    ax.legend(loc="lower right", fontsize=8, frameon=False, handlelength=1.6)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title("B", loc="left", fontsize=16, fontweight="bold", pad=6)

    divider = make_axes_locatable(ax)
    ax_marg = divider.append_axes("right", size="18%", pad=0.08, sharey=ax)

    clipped = np.clip(xtxst, qmin, qmax)
    bins = np.linspace(qmin, qmax, 50)
    ax_marg.hist(clipped, bins=bins, orientation="horizontal",
                 color=C_GREY_LIGHT, edgecolor=C_GREY_DARK, linewidth=0.3,
                 alpha=0.9)
    ax_marg.axhline(q95, color=C_DECISIVE, linestyle=":", linewidth=0.7,
                    alpha=0.7, zorder=2)
    ax_marg.axhline(q99, color=C_DECISIVE, linestyle=":", linewidth=0.7,
                    alpha=0.7, zorder=2)
    ax_marg.set_xlabel("SNPs", fontsize=8)
    ax_marg.tick_params(axis="x", labelsize=7)
    ax_marg.tick_params(axis="y", left=False, labelleft=False)
    for spine in ("top", "right", "left"):
        ax_marg.spines[spine].set_visible(False)
    ax_marg.set_ylim(qmin, qmax)

def build_layout(df: pd.DataFrame, top_n: int):
    chrom_sizes = df.groupby("chrom")["pos"].max()
    top = chrom_sizes.sort_values(ascending=False).head(top_n).index.tolist()
    gap = max(int(chrom_sizes[top].max() * 0.03), 500_000)
    offsets, spans, midpoints = {}, [], []
    cur = 0
    for ch in top:
        L = int(chrom_sizes[ch])
        offsets[ch] = cur
        spans.append((ch, cur, cur + L))
        midpoints.append(cur + L / 2)
        cur += L + gap
    return top, offsets, spans, midpoints, cur - gap


def plot_manhattan(ax: plt.Axes, df: pd.DataFrame) -> None:
    top, offsets, spans, midpoints, total_w = build_layout(df, N_TOP_SCAFFOLDS)
    sub = df[df["chrom"].isin(top)].copy()
    sub["x"] = sub.apply(lambda r: offsets[r["chrom"]] + r["pos"], axis=1)
    sub["chrom_idx"] = sub["chrom"].map({c: i for i, c in enumerate(top)})

    bg = sub[sub["bf_db"] <= BF_THRESHOLD]
    rng = np.random.default_rng(RNG_SEED)
    if len(bg) > BG_SUBSAMPLE:
        bg = bg.iloc[rng.choice(len(bg), BG_SUBSAMPLE, replace=False)]
    bg_colors = np.where(bg["chrom_idx"].values % 2 == 0, C_GREY_DARK, C_GREY_LIGHT)
    ax.scatter(bg["x"].values, bg["bf_db"].values, s=2,
               c=bg_colors, alpha=0.35, linewidths=0, rasterized=True)

    decisive = sub[sub["bf_db"] > BF_THRESHOLD]
    dec_b = decisive[decisive["beta_i"] > 0]
    dec_t = decisive[decisive["beta_i"] < 0]

    ax.scatter(dec_b["x"].values, dec_b["bf_db"].values,
               s=24, c=C_B, marker="^", edgecolor="black", linewidth=0.35,
               zorder=5, rasterized=True, label=rf"ALT on B ($\beta_i > 0$), n = {len(dec_b)}")
    ax.scatter(dec_t["x"].values, dec_t["bf_db"].values,
               s=24, c=C_T, marker="v", edgecolor="black", linewidth=0.35,
               zorder=5, rasterized=True, label=rf"ALT on T ($\beta_i < 0$), n = {len(dec_t)}")

    ax.axhline(BF_THRESHOLD, color=C_DECISIVE, linestyle="--", linewidth=0.9,
               alpha=0.9, zorder=3)

    top_snps = sub.nlargest(N_TOP_SNPS, "bf_db")
    for _, r in top_snps.iterrows():
        ax.annotate(
            f"BF = {r['bf_db']:.1f}",
            xy=(r["x"], r["bf_db"]),
            xytext=(0, 10), textcoords="offset points",
            ha="center", va="bottom", fontsize=7.5, color="black",
            fontweight="bold",
            arrowprops=dict(arrowstyle="-", color="black", lw=0.4),
        )

    ax.set_xticks(midpoints)
    ax.set_xticklabels([scaffold_short(ch) for ch in top], fontsize=8, rotation=30, ha="right")
    ax.set_xlabel("Scaffold", fontsize=10)
    ax.set_xlim(-total_w * 0.005, total_w * 1.005)

    ymax = max(float(sub["bf_db"].max()) + 5, 35)
    ax.set_ylim(-25, ymax + 2)
    ax.set_ylabel("BF(dB)   [host-covariate model]", fontsize=10)
    ax.tick_params(axis="y", labelsize=8)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    legend_handles = [
        plt.Line2D([0], [0], marker="^", linestyle="none",
                   markerfacecolor=C_B, markeredgecolor="black", markersize=7,
                   label=rf"ALT on B ($\beta_i > 0$), n = {len(dec_b)}"),
        plt.Line2D([0], [0], marker="v", linestyle="none",
                   markerfacecolor=C_T, markeredgecolor="black", markersize=7,
                   label=rf"ALT on T ($\beta_i < 0$), n = {len(dec_t)}"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=8,
              frameon=False, handletextpad=0.4)

    ax.set_title("C", loc="left", fontsize=16, fontweight="bold", pad=6)


def main() -> None:
    plt.rcParams.update({
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.8,
    })

    omega = load_omega()
    df = load_scan()
    xtxst = load_xtx()


    fig = plt.figure(figsize=(15, 7.0))
    gs = GridSpec(
        nrows=2, ncols=2,
        width_ratios=[1.0, 2.6],
        height_ratios=[1.0, 1.0],
        wspace=0.22, hspace=0.35,
        figure=fig,
    )
    ax_omega = fig.add_subplot(gs[0, 0])
    ax_calib = fig.add_subplot(gs[1, 0])
    ax_manh = fig.add_subplot(gs[:, 1])

    plot_omega(ax_omega, omega)
    plot_calibration(ax_calib, xtxst)
    plot_manhattan(ax_manh, df)

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")

if __name__ == "__main__":
    main()
