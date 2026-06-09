#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Patch

OUTDIR = "variance_analysis/section1_rigorous"
COLORS = {"B": "#4EA2FF", "T": "#EDB72F", "BT": "#9BAB96"}
N_EFF = 29
N_CENSUS = 500

f_reg = {
    "B":  (242, 238, 245),
    "T":  (347, 341, 355),
    "BT": (219, 216, 221),
}

poolseq = {
    "B":  (208, None, None),
    "T":  (211, None, None),
    "BT": (226, None, None),
}

manual = {
    "B":  (173, None, 450),
    "T":  (306, None, 800),
    "BT": (183, None, 450),
}

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 10,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",
})


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)

    treatments = ["B", "T", "BT"]
    method_offset = {"F-regression": +0.28,
                     "poolSeq (windowed)": 0.0,
                     "Method-of-moments": -0.28}
    method_marker = {"F-regression": "o",
                     "poolSeq (windowed)": "s",
                     "Method-of-moments": "^"}
    method_data = {
        "F-regression": f_reg,
        "poolSeq (windowed)": poolseq,
        "Method-of-moments": manual,
    }

    ax.axvspan(200, 350, alpha=0.10, color="gray", zorder=0)

    ax.axvline(N_EFF, ls="-", color="firebrick", lw=1.2, alpha=0.7, zorder=1)
    ax.axvline(N_CENSUS, ls="--", color="black", lw=1, alpha=0.7, zorder=1)

    y_positions = {t: i for i, t in enumerate(treatments[::-1])}

    for method, data in method_data.items():
        offset = method_offset[method]
        marker = method_marker[method]
        for trt in treatments:
            y = y_positions[trt] + offset
            point, lo, hi = data[trt]

            if lo is not None and np.isfinite(lo):
                ax.plot([lo, point], [y, y], "-", color=COLORS[trt],
                        lw=1.5, alpha=0.8, solid_capstyle="butt", zorder=4)
            if hi is not None and np.isfinite(hi):
                ax.plot([point, hi], [y, y], "-", color=COLORS[trt],
                        lw=1.5, alpha=0.8, solid_capstyle="butt", zorder=4)
                ax.plot([hi, hi], [y-0.06, y+0.06], "-",
                        color=COLORS[trt], lw=1, alpha=0.8, zorder=4)

            ax.scatter([point], [y], marker=marker, s=80,
                       color=COLORS[trt], edgecolor="black", linewidth=0.7,
                       zorder=5)

    ax.text(N_EFF, -1.0, f"$N_{{\\mathrm{{eff}}}} = {N_EFF}$",
            color="firebrick", fontsize=8, ha="center", va="top",
            fontweight="bold")
    ax.text(N_CENSUS, -1.0, f"census = {N_CENSUS}",
            color="black", fontsize=8, ha="center", va="top")
    ax.text(275, 2.7, "$N_e \\approx 200$–$350$\nconsensus range",
            color="gray", fontsize=7.5, ha="center", va="center",
            style="italic")

    for trt in treatments:
        y = y_positions[trt]
        ax.axhspan(y - 0.45, y + 0.45, color=COLORS[trt], alpha=0.08, zorder=0)

    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(list(y_positions.keys()), fontweight="bold")
    ax.set_xlim(0, 900)
    ax.set_ylim(-1.3, 3.0)
    ax.set_xlabel("Effective population size  $N_e$")

    legend_elems = [
        plt.Line2D([0], [0], marker="o", linestyle="", color="gray",
                   markeredgecolor="black", markersize=9,
                   label="$F$-regression (95% block bootstrap)"),
        plt.Line2D([0], [0], marker="s", linestyle="", color="gray",
                   markeredgecolor="black", markersize=9,
                   label="poolSeq windowed (point)"),
        plt.Line2D([0], [0], marker="^", linestyle="", color="gray",
                   markeredgecolor="black", markersize=9,
                   label="Method-of-moments ($\\to$ Chebyshev 95% UB)"),
    ]
    ax.legend(handles=legend_elems, loc="upper right", frameon=False,
              fontsize=8)

    plt.tight_layout()
    for ext in ["svg", "png"]:
        plt.savefig(os.path.join(OUTDIR, f"fig_SYY_ne_summary.{ext}"),
                    dpi=300 if ext == "png" else None, bbox_inches="tight")

if __name__ == "__main__":
    main()
