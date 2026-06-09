#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

OUTDIR = "variance_analysis/section1_rigorous"
COLORS = {"B": "#4EA2FF", "T": "#EDB72F", "M": "#9BAB96"}
TREATMENTS = ["B", "T", "M"]
GENS_4REP = [1, 2, 6, 7, 8, 9]

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


def load_depths():
    samples = [l.strip() for l in open("variance_analysis/sample_list.txt") if l.strip()]
    sums = np.zeros(len(samples)); counts = np.zeros(len(samples))
    with open("variance_analysis/merged_depth.tsv") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            for i, x in enumerate(parts[2:]):
                if x and x != "." and x != "0":
                    try:
                        v = int(x); sums[i] += v; counts[i] += 1
                    except ValueError:
                        pass
    means = sums / np.maximum(counts, 1)
    return dict(zip(samples, means))


def panel_depth(ax, sd):
    ax.set_title("a  Mean sequencing depth per pool", loc="left", fontweight="bold")
    for trt in TREATMENTS:
        rep_depths = {g: [sd.get(f"{trt}{r}G{g:02d}") for r in range(1, 5)]
                      for g in GENS_4REP}
        rep_depths = {g: [v for v in vs if v is not None] for g, vs in rep_depths.items()}
        means = [np.mean(rep_depths[g]) for g in GENS_4REP]
        sems = [np.std(rep_depths[g], ddof=1) / np.sqrt(len(rep_depths[g]))
                for g in GENS_4REP]
        ax.errorbar(GENS_4REP, means, yerr=sems,
                    marker="o", ms=5, lw=1.5, capsize=2,
                    color=COLORS[trt], label=trt)
        for g in GENS_4REP:
            xs = np.full(len(rep_depths[g]), g) + np.random.default_rng(hash(trt+str(g)) % 2**31).uniform(-0.15, 0.15, len(rep_depths[g]))
            ax.scatter(xs, rep_depths[g], s=10, alpha=0.4, color=COLORS[trt],
                       edgecolors="none")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Mean depth")
    ax.set_xticks(GENS_4REP)
    ax.legend(title="Treatment", loc="upper right", frameon=False)
    ax.set_ylim(bottom=0)


def panel_F(ax, df, key, summary, title, show_intercept_label=False):
    ax.set_title(title, loc="left", fontweight="bold")
    for trt in TREATMENTS:
        sub = df[df["treatment"] == trt]
        xs = sub["generation"].values
        ys = sub[key].values
        ax.plot(xs, ys, "o", ms=5, color=COLORS[trt],
                markeredgecolor="white", markeredgewidth=0.5)
        row = summary[summary["treatment"] == trt].iloc[0]
        if "raw" in key:
            slope, intercept = row["slope_raw"], row["int_raw"]
            ne, lo, hi = row["Ne_raw_point"], row["Ne_raw_lo"], row["Ne_raw_hi"]
        else:
            slope, intercept = row["slope_bio"], row["int_bio"]
            ne, lo, hi = row["Ne_bio_point"], row["Ne_bio_lo"], row["Ne_bio_hi"]
        xfit = np.linspace(0, 10, 50)
        yfit = intercept + slope * xfit
        ax.plot(xfit, yfit, "-", lw=1.5, color=COLORS[trt], alpha=0.85)
    y_base = 0.97
    for i, trt in enumerate(TREATMENTS):
        row = summary[summary["treatment"] == trt].iloc[0]
        if "raw" in key:
            ne, lo, hi = row["Ne_raw_point"], row["Ne_raw_lo"], row["Ne_raw_hi"]
        else:
            ne, lo, hi = row["Ne_bio_point"], row["Ne_bio_lo"], row["Ne_bio_hi"]
        ax.text(0.02, y_base - i * 0.08,
                f"{trt}: N$_e$ = {ne:.0f} [{lo:.0f}, {hi:.0f}]",
                transform=ax.transAxes, fontsize=8,
                color=COLORS[trt], fontweight="bold", va="top")
    ax.axhline(0, color="gray", lw=0.5, ls="--", alpha=0.5)
    ax.set_xlabel("Generation")
    ax.set_ylabel("F  (between-rep divergence)")
    ax.set_xticks(GENS_4REP)
    ax.set_xlim(0, 10)


def main():
    sd = load_depths()
    df = pd.read_csv(os.path.join(OUTDIR, "F_noise_corrected_full.tsv"), sep="\t")
    summary = pd.read_csv(os.path.join(OUTDIR, "F_noise_corrected_bootstrap.tsv"), sep="\t")

    fig, axes = plt.subplots(3, 1, figsize=(5.5, 8.5), dpi=150)
    panel_depth(axes[0], sd)
    panel_F(axes[1], df, "F_raw", summary,
            "b  Raw F between replicates  (uncorrected)")
    y1 = axes[1].get_ylim()
    panel_F(axes[2], df, "F_bio", summary,
            "c  Noise-corrected F  (pool + sequencing sampling subtracted)")
    y2 = axes[2].get_ylim()
    ymin = min(y1[0], y2[0], -0.005)
    ymax = max(y1[1], y2[1])
    axes[1].set_ylim(ymin, ymax)
    axes[2].set_ylim(ymin, ymax)

    plt.tight_layout()
    svg_path = os.path.join(OUTDIR, "fig_noise_correction.svg")
    png_path = os.path.join(OUTDIR, "fig_noise_correction.png")
    plt.savefig(svg_path)
    plt.savefig(png_path, dpi=300)

if __name__ == "__main__":
    main()
