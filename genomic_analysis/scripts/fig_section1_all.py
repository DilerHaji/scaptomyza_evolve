#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

OUTDIR = "variance_analysis/section1_rigorous"
COLORS = {"B": "#4EA2FF", "T": "#EDB72F", "M": "#9BAB96"}

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


def fig_poolseq_sweep():
    df = pd.read_csv(os.path.join(OUTDIR, "poolseq_sweep_summary.tsv"), sep="\t")
    sub = df[(df["method"] == "P.planI") & (df["time_window"] == "G01-G09")]

    fig, ax = plt.subplots(figsize=(4.5, 3.5), dpi=150)
    pool_sizes = sorted(sub["pool_size"].unique())
    x = np.arange(len(pool_sizes))
    width = 0.22

    for i, trt in enumerate(["B", "T", "M"]):
        vals = []
        for ps in pool_sizes:
            row = sub[(sub["treatment"] == trt) & (sub["pool_size"] == ps)]
            vals.append(float(row["Ne_mean_of_medians"].values[0]) if len(row) else np.nan)
        ax.bar(x + (i - 1) * width, vals, width, label=trt, color=COLORS[trt],
               edgecolor="white", linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels([str(ps) for ps in pool_sizes])
    ax.set_xlabel("Pool size used in correction (diploids)")
    ax.set_ylabel("$N_e$ (poolSeq, Plan I, G01–G09)")
    ax.legend(title="Treatment", frameon=False)

    ax.annotate("", xy=(0.95, 0.85), xytext=(2.05, 0.85),
                xycoords="data", textcoords="data",
                arrowprops=dict(arrowstyle="->", color="gray", lw=1.5))
    ax.text(1.5, 0.88 * ax.get_ylim()[1], "effective pool\nsize correction",
            ha="center", va="bottom", fontsize=7, color="gray", style="italic")

    ax.axhline(250, ls="--", color="gray", alpha=0.4, lw=0.8)
    ax.text(0.02, 255, "$N_e$ = 250 (census)", fontsize=7, color="gray",
            transform=ax.get_yaxis_transform())

    plt.tight_layout()
    for ext in ["svg", "png"]:
        plt.savefig(os.path.join(OUTDIR, f"fig_poolseq_sweep.{ext}"),
                    dpi=300 if ext == "png" else None)
    plt.close()


def fig_simulation():
    sim_data = {
        "Drift baseline": {"Ne_F": 292, "Ne_pool": 246, "G": 0.000, "cc_d": 0.000, "cc_t": 0.000},
        "Sustained\n(s=0.20, 10%)": {"Ne_F": 390, "Ne_pool": 72, "G": 0.123, "cc_d": 0.035, "cc_t": 0.332},
        "Episodic\n(s=0.20, 5%/int)": {"Ne_F": 290, "Ne_pool": 167, "G": -0.001, "cc_d": 0.022, "cc_t": 0.091},
        "Clustered\n(10 reg, s=0.30)": {"Ne_F": 307, "Ne_pool_wnd": 232, "G": 0.050, "cc_d": 0.014, "cc_t": 0.143},
    }
    obs_data = {
        "B": {"Ne_F": 242, "Ne_pool": 230, "G": -0.08, "cc_d": 0.035, "cc_t": 0.108},
        "T": {"Ne_F": 347, "Ne_pool": 234, "G": 0.14, "cc_d": -0.010, "cc_t": 0.076},
        "M": {"Ne_F": 219, "Ne_pool": 251, "G": 0.16, "cc_d": 0.053, "cc_t": 0.022},
    }

    fig, axes = plt.subplots(1, 3, figsize=(9, 3.5), dpi=150)

    ax = axes[0]
    ax.set_title("a  $N_e$ estimates", loc="left", fontweight="bold")
    scenarios = list(sim_data.keys())
    y_sim = np.arange(len(scenarios))
    for yi, sc in enumerate(scenarios):
        d = sim_data[sc]
        ne_p = d.get("Ne_pool_wnd", d.get("Ne_pool", np.nan))
        ax.plot(d["Ne_F"], yi, "s", ms=7, color="steelblue", zorder=5)
        ax.plot(ne_p, yi, "o", ms=7, color="coral", zorder=5)
    y_obs_start = len(scenarios) + 0.5
    for oi, (trt, d) in enumerate(obs_data.items()):
        y = y_obs_start + oi
        ax.plot(d["Ne_F"], y, "s", ms=8, color=COLORS[trt], zorder=5,
                markeredgecolor="black", markeredgewidth=0.5)
        ax.plot(d["Ne_pool"], y, "o", ms=8, color=COLORS[trt], zorder=5,
                markeredgecolor="black", markeredgewidth=0.5)
    ax.set_yticks(list(range(len(scenarios))) + [y_obs_start + i for i in range(3)])
    ax.set_yticklabels(scenarios + ["Obs B", "Obs T", "Obs M"])
    ax.axvline(250, ls="--", color="gray", alpha=0.4, lw=0.8)
    ax.axhline(y_obs_start - 0.25, color="gray", alpha=0.3, lw=0.8)
    ax.set_xlabel("$N_e$")
    ax.plot([], [], "s", ms=6, color="gray", label="F-regression")
    ax.plot([], [], "o", ms=6, color="gray", label="poolSeq")
    ax.legend(frameon=False, fontsize=7, loc="lower right")

    ax = axes[1]
    ax.set_title("b  $G(k{\\geq}2)$", loc="left", fontweight="bold")
    for yi, sc in enumerate(scenarios):
        ax.barh(yi, sim_data[sc]["G"], height=0.6, color="lightsteelblue",
                edgecolor="steelblue", linewidth=0.5)
    for oi, (trt, d) in enumerate(obs_data.items()):
        y = y_obs_start + oi
        ax.barh(y, d["G"], height=0.6, color=COLORS[trt],
                edgecolor="black", linewidth=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.axhline(y_obs_start - 0.25, color="gray", alpha=0.3, lw=0.8)
    ax.set_yticks([])
    ax.set_xlabel("$G(k{\\geq}2)$")

    ax = axes[2]
    ax.set_title("c  Within-interval cc", loc="left", fontweight="bold")
    for yi, sc in enumerate(scenarios):
        ax.barh(yi, sim_data[sc]["cc_d"], height=0.6, color="lightsteelblue",
                edgecolor="steelblue", linewidth=0.5)
    for oi, (trt, d) in enumerate(obs_data.items()):
        y = y_obs_start + oi
        ax.barh(y, d["cc_d"], height=0.6, color=COLORS[trt],
                edgecolor="black", linewidth=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.axhline(y_obs_start - 0.25, color="gray", alpha=0.3, lw=0.8)
    ax.set_yticks([])
    ax.set_xlabel("Diagonal cc")

    plt.tight_layout()
    for ext in ["svg", "png"]:
        plt.savefig(os.path.join(OUTDIR, f"fig_simulation.{ext}"),
                    dpi=300 if ext == "png" else None)
    plt.close()


if __name__ == "__main__":
    fig_poolseq_sweep()
    fig_simulation()
