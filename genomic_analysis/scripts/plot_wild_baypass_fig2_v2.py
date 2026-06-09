#!/usr/bin/env python3

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

ROOT = Path(".")
OMEGA_FILE = ROOT / "baypass_wild/wild_omega_mat_omega.out"
BETAI_FILE = ROOT / "baypass_wild/wild_trt_summary_betai_reg.out"
XTX_FILE = ROOT / "baypass_wild/wild_omega_summary_pi_xtx.out"
POS_FILE = ROOT / "baypass_wild/wild_snp_positions.csv"
OUT_PNG = ROOT / "final_plots/wild/wild_baypass_fig2_v2.png"
OUT_SVG = ROOT / "final_plots/wild/wild_baypass_fig2_v2.svg"

POOL_ORDER = ["AVB", "AVT", "PSB", "PST", "RMB", "RMT"]
HOST_OF = {"AVB": "B", "AVT": "T", "PSB": "B", "PST": "T", "RMB": "B", "RMT": "T"}

C_B = "#499FFF"              # B host label (heatmap ticks)
C_T = "#EDB72D"              # T host label
C_BRIGHT = "#00A5FF"         # bright blue for both-high SNPs
C_THRESH = "#CC3333"         # red dashed threshold lines
C_SCAFF_DARK = "#111111"     # near-black for odd-indexed scaffolds
C_SCAFF_LIGHT = "#888888"    # mid-grey for even-indexed scaffolds
C_GREY_MID = "#555555"

BF_THRESHOLD = 10.0
THIN_STEP = 16
DF_CHI2 = 6
XTX_TAIL_P = 0.10                  # one-tail p cutoff for "high X^tX" (softer than 0.05)
N_TOP_SCAFFOLDS = 10
N_TOP_SNPS_LABEL = 5
BG_SUBSAMPLE_FULL = 80_000        # for BF Manhattan (812,670 SNPs)
BG_SUBSAMPLE_THIN = 45_000        # for X^tX Manhattan (50,792 SNPs)
BG_SUBSAMPLE_SCATTER = 45_000
RNG_SEED = 20260418


def scaffold_short(name: str) -> str:
    m = re.search(r"(\d+)_HRSCAF", str(name))
    return m.group(1) if m else str(name)

def load_all():
    omega = np.loadtxt(OMEGA_FILE)

    pos = pd.read_csv(POS_FILE)
    betai = pd.read_csv(BETAI_FILE, sep=r"\s+").rename(
        columns={"MRK": "mrk", "BF(dB)": "bf_db", "Beta_is": "beta_i"}
    )
    bf_full = pos.merge(betai[["mrk", "bf_db", "beta_i"]], on="mrk")

    xtx = pd.read_csv(XTX_FILE, sep=r"\s+").rename(
        columns={"XtXst": "xtxst"}
    )
    pos_thin = pos.iloc[::THIN_STEP].reset_index(drop=True)
    if len(pos_thin) == len(xtx) + 1:
        pos_thin = pos_thin.iloc[:len(xtx)].copy()
    assert len(pos_thin) == len(xtx)
    xtx_thin = pd.concat(
        [pos_thin.reset_index(drop=True), xtx[["xtxst"]].reset_index(drop=True)],
        axis=1,
    )

    betai_thin = betai.iloc[::THIN_STEP].reset_index(drop=True)
    if len(betai_thin) == len(xtx_thin) + 1:
        betai_thin = betai_thin.iloc[:len(xtx_thin)].copy()
    xtx_thin["bf_db"] = betai_thin["bf_db"].values
    xtx_thin["beta_i"] = betai_thin["beta_i"].values

    q_xtx = stats.chi2.ppf(1 - XTX_TAIL_P, df=DF_CHI2)
    both_mask = (xtx_thin["bf_db"] > BF_THRESHOLD) & (xtx_thin["xtxst"] > q_xtx)
    both_positions = set(
        zip(xtx_thin.loc[both_mask, "chrom"], xtx_thin.loc[both_mask, "pos"])
    )

    return omega, bf_full, xtx_thin, both_positions


def build_layout(df: pd.DataFrame, top_n: int):
    chrom_sizes = df.groupby("chrom")["pos"].max()
    top = chrom_sizes.sort_values(ascending=False).head(top_n).index.tolist()
    gap = max(int(chrom_sizes[top].max() * 0.03), 500_000)
    offsets, midpoints = {}, []
    cur = 0
    for ch in top:
        L = int(chrom_sizes[ch])
        offsets[ch] = cur
        midpoints.append(cur + L / 2)
        cur += L + gap
    return top, offsets, midpoints, cur - gap


def plot_omega(ax: plt.Axes, omega: np.ndarray) -> None:
    sd = np.sqrt(np.diag(omega))
    cor = omega / np.outer(sd, sd)
    dist = 1 - np.abs(cor)
    np.fill_diagonal(dist, 0)
    Z = linkage(squareform(dist, checks=False), method="average")
    leaf_order = dendrogram(Z, no_plot=True)["leaves"]
    cor_reordered = cor[np.ix_(leaf_order, leaf_order)]
    labels = [POOL_ORDER[i] for i in leaf_order]

    im = ax.imshow(cor_reordered, cmap=plt.cm.Greys, vmin=0.0, vmax=1.0, aspect="equal")

    n = len(labels)
    for i in range(n):
        for j in range(n):
            val = cor_reordered[i, j]
            tc = "white" if val > 0.55 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=8, color=tc)

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=0)
    ax.set_yticklabels(labels)
    for tick, lab in zip(ax.get_xticklabels(), labels):
        tick.set_color(C_B if HOST_OF[lab] == "B" else C_T)
        tick.set_fontweight("bold")
    for tick, lab in zip(ax.get_yticklabels(), labels):
        tick.set_color(C_B if HOST_OF[lab] == "B" else C_T)
        tick.set_fontweight("bold")
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.72, pad=0.03)
    cbar.set_label(r"Corr from $\mathbf{\Omega}$", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax.set_title("A", loc="left", fontsize=18, fontweight="bold", pad=12)


def plot_bf_vs_xtx(ax: plt.Axes, df_thin: pd.DataFrame) -> None:
    q_xtx = stats.chi2.ppf(1 - XTX_TAIL_P, df=DF_CHI2)
    bf = df_thin["bf_db"].values
    xt = df_thin["xtxst"].values

    q_both = (bf > BF_THRESHOLD) & (xt > q_xtx)
    q_host_only = (bf > BF_THRESHOLD) & (xt <= q_xtx)
    q_diff_only = (bf <= BF_THRESHOLD) & (xt > q_xtx)
    q_neither = (bf <= BF_THRESHOLD) & (xt <= q_xtx)

    rng = np.random.default_rng(RNG_SEED)
    idx_bg = np.where(q_neither)[0]
    if len(idx_bg) > BG_SUBSAMPLE_SCATTER:
        idx_bg = rng.choice(idx_bg, BG_SUBSAMPLE_SCATTER, replace=False)
    ax.scatter(bf[idx_bg], xt[idx_bg], s=2, c=C_SCAFF_LIGHT,
               alpha=0.4, linewidths=0, rasterized=True)

    ax.scatter(bf[q_diff_only], xt[q_diff_only], s=7,
               c=C_SCAFF_DARK, alpha=0.6, linewidths=0, rasterized=True)
    ax.scatter(bf[q_host_only], xt[q_host_only], s=12,
               c=C_SCAFF_DARK, alpha=0.8, linewidths=0, rasterized=True)
    ax.scatter(bf[q_both], xt[q_both], s=26, c=C_BRIGHT,
               edgecolor="black", linewidth=0.5, zorder=6, rasterized=True)

    ax.axhline(q_xtx, color=C_THRESH, linestyle="--", linewidth=0.9, zorder=3)
    ax.axvline(BF_THRESHOLD, color=C_THRESH, linestyle="--", linewidth=0.9, zorder=3)

    xmin = float(np.quantile(bf, 0.001)) - 1
    xmax = float(np.quantile(bf, 0.9999)) + 2
    ymin = -0.5
    ymax = max(35, float(np.quantile(xt, 0.9999)) + 2)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    n_tot = len(df_thin)
    pct = lambda n: 100.0 * n / n_tot

    def _box(x, y, text, ha="left", va="top", edge=C_GREY_MID, weight="normal", color="black"):
        ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va,
                fontsize=7.5, color=color, fontweight=weight,
                bbox=dict(boxstyle="round,pad=0.22",
                          facecolor="white", edgecolor=edge, linewidth=0.6))

    _box(0.04, 0.97, f"X^tX only\nn = {q_diff_only.sum():,} ({pct(q_diff_only.sum()):.1f}%)")
    _box(0.96, 0.97,
         f"both\nn = {q_both.sum()} ({pct(q_both.sum()):.2f}%)",
         ha="right", edge=C_BRIGHT, weight="bold")
    _box(0.96, 0.05,
         f"BF only\nn = {q_host_only.sum()} ({pct(q_host_only.sum()):.2f}%)",
         ha="right", va="bottom", edge=C_THRESH)
    _box(0.04, 0.05, f"background\nn = {q_neither.sum():,} ({pct(q_neither.sum()):.1f}%)",
         va="bottom")

    rho, _ = stats.spearmanr(bf, xt)
    ax.text(0.04, 0.55, f"Spearman $\\rho$ = {rho:+.2f}",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.25",
                      facecolor="white", edgecolor=C_GREY_MID, linewidth=0.6))

    ax.set_xlabel("BF(dB)   [directional, host axis]", fontsize=10)
    ax.set_ylabel(r"$X^tX_{st}$   [undirected, overall differentiation]", fontsize=10)
    ax.tick_params(labelsize=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax.set_title("B", loc="left", fontsize=18, fontweight="bold", pad=12)


def plot_bf_manhattan(ax: plt.Axes, df: pd.DataFrame, both_positions,
                      layout) -> tuple[list, dict, list, int]:
    top, offsets, midpoints, total_w = layout
    sub = df[df["chrom"].isin(top)].copy()
    sub["x"] = sub.apply(lambda r: offsets[r["chrom"]] + r["pos"], axis=1)
    sub["chrom_idx"] = sub["chrom"].map({c: i for i, c in enumerate(top)})

    rng = np.random.default_rng(RNG_SEED)
    bg = sub[sub["bf_db"] <= BF_THRESHOLD]
    if len(bg) > BG_SUBSAMPLE_FULL:
        bg = bg.iloc[rng.choice(len(bg), BG_SUBSAMPLE_FULL, replace=False)]
    bg_colors = np.where(bg["chrom_idx"].values % 2 == 0, C_SCAFF_DARK, C_SCAFF_LIGHT)
    ax.scatter(bg["x"].values, bg["bf_db"].values, s=2, c=bg_colors,
               alpha=0.35, linewidths=0, rasterized=True)

    decisive = sub[sub["bf_db"] > BF_THRESHOLD].copy()
    decisive["is_both"] = [(c, p) in both_positions
                           for c, p in zip(decisive["chrom"], decisive["pos"])]
    dec_plain = decisive[~decisive["is_both"]]
    dec_both = decisive[decisive["is_both"]]

    plain_colors = np.where(dec_plain["chrom_idx"].values % 2 == 0,
                            C_SCAFF_DARK, C_SCAFF_LIGHT)
    ax.scatter(dec_plain["x"].values, dec_plain["bf_db"].values,
               s=14, c=plain_colors, edgecolor="black", linewidth=0.2,
               alpha=0.85, zorder=4, rasterized=True)
    ax.scatter(dec_both["x"].values, dec_both["bf_db"].values,
               s=28, c=C_BRIGHT, edgecolor="black", linewidth=0.4,
               zorder=6, rasterized=True,
               label=f"both-high (n = {len(dec_both)})")

    ax.axhline(BF_THRESHOLD, color=C_THRESH, linestyle="--", linewidth=0.9, zorder=3)
    ax.text(total_w * 0.999, BF_THRESHOLD + 1.0,
            f"BF = {int(BF_THRESHOLD)}",
            color=C_THRESH, fontsize=7.5, ha="right", va="bottom")

    top_snps = sub.nlargest(N_TOP_SNPS_LABEL, "bf_db")
    for _, r in top_snps.iterrows():
        ax.annotate(f"BF = {r['bf_db']:.1f}",
                    xy=(r["x"], r["bf_db"]),
                    xytext=(0, 10), textcoords="offset points",
                    ha="center", va="bottom", fontsize=7.5,
                    fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color="black", lw=0.4))

    ax.set_ylabel("BF(dB)  [host-covariate]", fontsize=10)
    ax.tick_params(axis="y", labelsize=8)
    ax.set_xlim(-total_w * 0.005, total_w * 1.005)
    ymax = max(float(sub["bf_db"].max()) + 5, 35)
    ax.set_ylim(-25, ymax)
    ax.set_xticks([])
    for spine in ("top", "right", "bottom"):
        ax.spines[spine].set_visible(False)

    if len(dec_both) > 0:
        ax.legend(loc="upper left", fontsize=8, frameon=False)

    ax.set_title("C", loc="left", fontsize=18, fontweight="bold", pad=6)

    return top, offsets, midpoints, total_w

def plot_xtx_manhattan(ax: plt.Axes, df_thin: pd.DataFrame, both_positions,
                       layout) -> None:
    top, offsets, midpoints, total_w = layout
    q_xtx = stats.chi2.ppf(1 - XTX_TAIL_P, df=DF_CHI2)
    sub = df_thin[df_thin["chrom"].isin(top)].copy()
    sub["x"] = sub.apply(lambda r: offsets[r["chrom"]] + r["pos"], axis=1)
    sub["chrom_idx"] = sub["chrom"].map({c: i for i, c in enumerate(top)})

    rng = np.random.default_rng(RNG_SEED)
    bg = sub[sub["xtxst"] <= q_xtx]
    if len(bg) > BG_SUBSAMPLE_THIN:
        bg = bg.iloc[rng.choice(len(bg), BG_SUBSAMPLE_THIN, replace=False)]
    bg_colors = np.where(bg["chrom_idx"].values % 2 == 0, C_SCAFF_DARK, C_SCAFF_LIGHT)
    ax.scatter(bg["x"].values, bg["xtxst"].values, s=2, c=bg_colors,
               alpha=0.4, linewidths=0, rasterized=True)

    outliers = sub[sub["xtxst"] > q_xtx].copy()
    outliers["is_both"] = [(c, p) in both_positions
                           for c, p in zip(outliers["chrom"], outliers["pos"])]
    out_plain = outliers[~outliers["is_both"]]
    out_both = outliers[outliers["is_both"]]

    plain_colors = np.where(out_plain["chrom_idx"].values % 2 == 0,
                            C_SCAFF_DARK, C_SCAFF_LIGHT)
    ax.scatter(out_plain["x"].values, out_plain["xtxst"].values,
               s=12, c=plain_colors, edgecolor="black", linewidth=0.2,
               alpha=0.85, zorder=4, rasterized=True)
    ax.scatter(out_both["x"].values, out_both["xtxst"].values,
               s=28, c=C_BRIGHT, edgecolor="black", linewidth=0.4,
               zorder=6, rasterized=True)

    ax.axhline(q_xtx, color=C_THRESH, linestyle="--", linewidth=0.9, zorder=3)
    ax.text(total_w * 0.999, q_xtx + 0.7,
            f"$\\chi^2_{{{DF_CHI2}}}$ {int(XTX_TAIL_P * 100)}% = {q_xtx:.2f}",
            color=C_THRESH, fontsize=7.5, ha="right", va="bottom")

    ax.set_xticks(midpoints)
    ax.set_xticklabels([scaffold_short(ch) for ch in top], fontsize=8, rotation=30, ha="right")
    ax.set_xlabel("Scaffold", fontsize=10)
    ax.set_xlim(-total_w * 0.005, total_w * 1.005)
    ymax = max(float(sub["xtxst"].max()) + 3, 35)
    ax.set_ylim(-0.5, ymax)
    ax.set_ylabel(r"$X^tX_{st}$  [overall differentiation]", fontsize=10)
    ax.tick_params(axis="y", labelsize=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax.set_title("D", loc="left", fontsize=18, fontweight="bold", pad=6)



def main() -> None:
    plt.rcParams.update({
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.8,
    })

    omega, bf_full, xtx_thin, both_positions = load_all()

    fig = plt.figure(figsize=(16, 8.5))
    gs = GridSpec(
        nrows=2, ncols=2,
        width_ratios=[1.0, 2.3],
        height_ratios=[1.0, 1.0],
        wspace=0.22, hspace=0.12,
        figure=fig,
    )
    ax_omega = fig.add_subplot(gs[0, 0])
    ax_scatter = fig.add_subplot(gs[1, 0])
    ax_bf = fig.add_subplot(gs[0, 1])
    ax_xtx = fig.add_subplot(gs[1, 1])

    plot_omega(ax_omega, omega)
    plot_bf_vs_xtx(ax_scatter, xtx_thin)

    layout = build_layout(bf_full, N_TOP_SCAFFOLDS)
    plot_bf_manhattan(ax_bf, bf_full, both_positions, layout)
    plot_xtx_manhattan(ax_xtx, xtx_thin, both_positions, layout)

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")

if __name__ == "__main__":
    main()
