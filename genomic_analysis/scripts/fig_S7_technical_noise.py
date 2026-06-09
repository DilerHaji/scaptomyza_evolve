#!/usr/bin/env python3

import os
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import LogNorm

IND_FOUNDER_DIR = "../founder_popstructure/af_comparison"
OUTDIR = "variance_analysis/section1_rigorous"

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",
})


def load_pool_and_ind():
    pool = pd.read_csv(f"{IND_FOUNDER_DIR}/pool_freq_table.tsv", sep="\t")
    pool["pool_mean"] = pool[["F1G00", "F2G00", "F3G00", "F4G00"]].mean(axis=1)
    pool["pool_var"] = pool[["F1G00", "F2G00", "F3G00", "F4G00"]].var(axis=1, ddof=1)

    ind = pd.read_csv(f"{IND_FOUNDER_DIR}/ind_pseudopool.mafs.gz",
                      sep="\t", compression="gzip")

    ind = ind.rename(columns={"chromo": "chrom", "position": "pos"})
    m = pool.merge(ind[["chrom", "pos", "major", "minor", "knownEM"]],
                   on=["chrom", "pos"], how="inner")

    m["alt"] = m["alt"].str.upper()
    m["major"] = m["major"].str.upper()
    m["minor"] = m["minor"].str.upper()
    m["ind_pseudo_af"] = np.where(m["minor"] == m["alt"], m["knownEM"],
                          np.where(m["major"] == m["alt"], 1 - m["knownEM"], np.nan))
    m = m.dropna(subset=["ind_pseudo_af"])
    return m


def panel_a_scatter(ax, df):
    ax.set_title("a  Individual vs pool-seq AF", loc="left", fontweight="bold")

    hb = ax.hexbin(df["ind_pseudo_af"], df["pool_mean"],
                   gridsize=60, cmap="viridis", mincnt=1,
                   norm=LogNorm())

    ax.plot([0, 1], [0, 1], "r--", lw=1, alpha=0.6, label="$y = x$")

    from scipy import stats as sp_stats
    slope, intercept, r, p, se = sp_stats.linregress(
        df["ind_pseudo_af"].values, df["pool_mean"].values)
    xfit = np.linspace(0, 1, 50)
    ax.plot(xfit, intercept + slope * xfit, "-", color="orange", lw=1.2,
            alpha=0.8, label=f"OLS ($R^2 = {r**2:.2f}$)")

    ax.set_xlabel("Individual pseudo-pool AF\n(192 flies, ~0.5$\\times$ each)")
    ax.set_ylabel("Pool-seq mean AF\n($F_1$–$F_4$)")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(frameon=False, loc="upper left", fontsize=7)
    cbar = plt.colorbar(hb, ax=ax, shrink=0.8)
    cbar.set_label("$n$ sites", fontsize=7)

def panel_b_variance(ax, df):
    ax.set_title("b  Between-founder-pool variance vs expected",
                 loc="left", fontweight="bold")
    N_pool = 80
    df = df.copy()
    df["expected_var"] = df["pool_mean"] * (1 - df["pool_mean"]) / (2 * N_pool)

    df["maf"] = np.minimum(df["pool_mean"], 1 - df["pool_mean"])

    bins = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
    labels = ["0–0.1", "0.1–0.2", "0.2–0.3", "0.3–0.4", "0.4–0.5"]
    df["maf_bin"] = pd.cut(df["maf"], bins=bins, labels=labels, include_lowest=True)

    obs = df.groupby("maf_bin", observed=False)["pool_var"].mean()
    exp = df.groupby("maf_bin", observed=False)["expected_var"].mean()
    counts = df.groupby("maf_bin", observed=False).size()
    inflation = obs / exp

    x = np.arange(len(labels))
    w = 0.35
    ax.bar(x - w/2, obs.values, w, label="Observed",
           color="coral", edgecolor="black", linewidth=0.5)
    ax.bar(x + w/2, exp.values, w, label="Expected binomial\n(N = 80)",
           color="steelblue", edgecolor="black", linewidth=0.5)

    for i, (lab, infl, n) in enumerate(zip(labels, inflation, counts)):
        if np.isfinite(infl):
            ax.text(i, obs.iloc[i] * 1.03, f"{infl:.1f}$\\times$",
                    ha="center", fontsize=7, color="firebrick")

    overall_obs = df["pool_var"].mean()
    overall_exp = df["expected_var"].mean()
    overall_inflation = overall_obs / overall_exp
    ax.text(0.98, 0.98,
            f"Genome-wide inflation: {overall_inflation:.2f}$\\times$\n"
            f"$\\Rightarrow$ $N_{{\\mathrm{{eff}}}} = {N_pool/overall_inflation:.0f}$ diploids",
            transform=ax.transAxes, va="top", ha="right",
            fontsize=8, color="firebrick",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="firebrick", alpha=0.9))

    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_xlabel("MAF bin (pool mean)")
    ax.set_ylabel("Allele frequency variance")
    ax.legend(frameon=False, loc="upper left", fontsize=7.5)

    return overall_inflation


def panel_c_summary(ax, inflation_factor):
    ax.set_title("c  $N_{\\mathrm{eff}}$ from three independent methods",
                 loc="left", fontweight="bold")
    N_pool = 80
    methods = [
        ("Founder pool\noverdispersion", 80 / inflation_factor, "#4EA2FF"),
        ("Temporal covariance\n($k = 1$ diagonal)", 24, "#EDB72F"),
        ("Pool technical\nreplicates", None, "#9BAB96"),  # TBD
    ]
    y_positions = np.arange(len(methods))[::-1]  # top to bottom
    for y, (label, val, color) in zip(y_positions, methods):
        if val is None:
            ax.errorbar([29], [y], xerr=[[5], [3]], fmt="o",
                        color=color, markersize=10, markeredgecolor="black",
                        markeredgewidth=0.5, capsize=4, alpha=0.3, zorder=5)
            ax.text(40, y, "(pending pool_variation\npipeline)",
                    fontsize=7, color="gray", style="italic", va="center")
        else:
            ax.scatter([val], [y], s=100, color=color, edgecolor="black",
                        linewidth=0.5, zorder=5)
            ax.text(val + 2, y, f"{val:.0f}", fontsize=9, va="center",
                    fontweight="bold")

    ax.axvline(N_pool, color="gray", ls="--", lw=1, alpha=0.6)
    ax.text(N_pool + 1, -0.55, f"$N_{{\\mathrm{{pool}}}} = {N_pool}$\n(nominal)",
            fontsize=7, color="gray", va="top")

    ax.axvspan(24, 32, alpha=0.1, color="gray")
    ax.text(28, -0.55, "range of estimates\n($N_{\\mathrm{eff}} = 29$ adopted)",
            fontsize=7, color="gray", va="top", ha="center")

    ax.set_yticks(y_positions)
    ax.set_yticklabels([m[0] for m in methods])
    ax.set_xlabel("$N_{\\mathrm{eff}}$  (diploid individuals)")
    ax.set_xlim(0, 90)
    ax.set_ylim(-1.2, len(methods) - 0.5)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    df = load_pool_and_ind()

    fig = plt.figure(figsize=(8, 4), dpi=150)
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.3)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    panel_a_scatter(ax_a, df)
    inflation = panel_b_variance(ax_b, df)

    plt.tight_layout()
    for ext in ["svg", "png"]:
        plt.savefig(os.path.join(OUTDIR, f"fig_S7_technical_noise.{ext}"),
                    dpi=300 if ext == "png" else None, bbox_inches="tight")

if __name__ == "__main__":
    main()
