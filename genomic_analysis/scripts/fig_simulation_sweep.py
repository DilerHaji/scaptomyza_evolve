#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import TwoSlopeNorm

INFILE = "variance_analysis/section1_rigorous/simulation_sweep.tsv"
OUTDIR = "variance_analysis/section1_rigorous"
COLORS = {"B": "#4EA2FF", "T": "#EDB72F", "M": "#9BAB96"}

OBS = {
    "B": {"G_k2": -0.08, "cc_diag": 0.035, "cc_total": 0.108, "Ne_pool": 230},
    "T": {"G_k2":  0.14, "cc_diag":-0.010, "cc_total": 0.076, "Ne_pool": 234},
    "M": {"G_k2":  0.16, "cc_diag": 0.053, "cc_total": 0.022, "Ne_pool": 251},
}

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",
})


def main():
    df = pd.read_csv(INFILE, sep="\t")

    stats = [
        ("G_k2_median", "$G(k \\geq 2)$", "RdBu_r", True),
        ("cc_diag_median", "Within-interval cc", "RdBu_r", True),
        ("cc_total_median", "Total-trajectory cc", "RdBu_r", True),
        ("Ne_pool_median", "poolSeq $N_e$ (windowed)", "viridis", False),
    ]
    ld_vals = sorted(df["ld_size"].unique())
    s_vals = sorted(df["s"].unique())
    frac_vals = sorted(df["frac"].unique())

    drift_row = df[df["s"] == 0.0].iloc[0] if len(df[df["s"] == 0.0]) else None
    df_sel = df[df["s"] > 0].copy()
    s_sel = sorted(df_sel["s"].unique())

    n_stats = len(stats)
    n_ld = len(ld_vals)

    fig, axes = plt.subplots(n_stats, n_ld, figsize=(3.2 * n_ld, 2.8 * n_stats),
                              squeeze=False, dpi=150)

    for si, (stat_col, stat_label, cmap, diverging) in enumerate(stats):
        for li, ld in enumerate(ld_vals):
            ax = axes[si, li]

            mat = np.full((len(s_sel), len(frac_vals)), np.nan)
            for ri, s in enumerate(s_sel):
                for ci, frac in enumerate(frac_vals):
                    row = df_sel[(df_sel["s"] == s) &
                                 (df_sel["frac"] == frac) &
                                 (df_sel["ld_size"] == ld)]
                    if len(row) == 1:
                        mat[ri, ci] = row[stat_col].values[0]

            if diverging:
                vmax = max(abs(np.nanmin(mat)), abs(np.nanmax(mat)), 0.01)
                norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
            else:
                norm = None

            im = ax.imshow(mat, aspect="auto", cmap=cmap, norm=norm,
                           origin="lower")

            ax.set_xticks(range(len(frac_vals)))
            ax.set_xticklabels([f"{f:.0%}" for f in frac_vals])
            ax.set_yticks(range(len(s_sel)))
            ax.set_yticklabels([f"{s:.2f}" for s in s_sel])

            for ri in range(len(s_sel)):
                for ci in range(len(frac_vals)):
                    v = mat[ri, ci]
                    if np.isfinite(v):
                        fmt = f"{v:.0f}" if stat_col.startswith("Ne") else f"{v:.3f}"
                        color = "white" if abs(v) > 0.6 * np.nanmax(abs(mat)) else "black"
                        ax.text(ci, ri, fmt, ha="center", va="center",
                                fontsize=6, color=color)

            if si == 0:
                ax.set_title(f"LD = {ld}", fontweight="bold")
            if si == n_stats - 1:
                ax.set_xlabel("Fraction selected")
            if li == 0:
                ax.set_ylabel(f"{stat_label}\n$s$")

            plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02)

            if drift_row is not None and stat_col in drift_row.index:
                drift_val = drift_row[stat_col]
                if np.isfinite(drift_val):
                    fmt = f"{drift_val:.0f}" if stat_col.startswith("Ne") else f"{drift_val:.3f}"
                    ax.text(0.02, 0.98, f"drift: {fmt}", transform=ax.transAxes,
                            fontsize=6, va="top", ha="left", color="gray",
                            style="italic")

    plt.tight_layout()
    for ext in ["svg", "png"]:
        plt.savefig(f"{OUTDIR}/fig_simulation_sweep.{ext}",
                    dpi=300 if ext == "png" else None,
                    bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
