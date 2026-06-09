#!/usr/bin/env python3
from __future__ import annotations
import re
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(".")
MASTER_TSV = ROOT / "final_plots/wild/section2_candidate_master_v2.tsv"
GLMM_PER_WINDOW_TSV = ROOT / "final_plots/wild/section2_glmm_per_window.tsv"
OUT_BASE = ROOT / "final_plots/wild/section2_fig3c_pointcloud_diversity"

WIN = 200_000
SCAFF_MIN_WINDOWS = 50
TOP_FRAC = 0.05
PEAK_CHROM = "chr_ScDA7r2_439_HRSCAF_779"
PEAK_START = 2_800_000

GREY = "#bdbdbd"
RED = "#C84A45"

ROWS = [
    ("θπ slope", [
        ("B_thetaPi_slope",       "B",   r"$\Delta\theta_\pi$/gen"),
        ("T_thetaPi_slope",       "T",   r"$\Delta\theta_\pi$/gen"),
        ("BT_thetaPi_slope_mean", "B+T", r"$\Delta\theta_\pi$/gen"),
    ]),
    ("θW slope", [
        ("B_thetaW_slope",        "B",   r"$\Delta\theta_W$/gen"),
        ("T_thetaW_slope",        "T",   r"$\Delta\theta_W$/gen"),
        ("BT_thetaW_slope_mean",  "B+T", r"$\Delta\theta_W$/gen"),
    ]),
]
DIRECTION = "low"  # diversity loss → negative slope under selection


def main() -> None:
    plt.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.6,
    })

    df = pd.read_csv(MASTER_TSV, sep="\t")
    glmm = pd.read_csv(GLMM_PER_WINDOW_TSV, sep="\t")
    df = df.merge(glmm[["chrom", "start", "glmm_lrt"]],
                  on=["chrom", "start"], how="left")

    chrom_sizes = df.groupby("chrom").size().sort_values(ascending=False)
    keep = chrom_sizes[chrom_sizes >= SCAFF_MIN_WINDOWS].index.tolist()
    df_p = df[df["chrom"].isin(keep)].copy().reset_index(drop=True)
    df_p["BT_thetaPi_slope_mean"] = df_p[["B_thetaPi_slope", "T_thetaPi_slope"]].mean(axis=1)
    df_p["BT_thetaW_slope_mean"]  = df_p[["B_thetaW_slope",  "T_thetaW_slope"]].mean(axis=1)

    glmm_thr = df_p["glmm_lrt"].quantile(1 - TOP_FRAC)
    df_p["IND_glmm_lrt"] = (df_p["glmm_lrt"] >= glmm_thr).fillna(False).astype(int)
    IND_COLS = [
        "IND_HV_blocks_B", "IND_HV_blocks_T", "IND_HV_blocks_M",
        "IND_cov_BT_neg", "IND_slope_div", "IND_permFST",
        "IND_glmm_lrt", "IND_wild_C2",
    ]
    df_p["votes_v3"] = df_p[IND_COLS].fillna(0).astype(int).sum(axis=1)

    n_rows = len(ROWS)
    n_sig_per_row = 3
    fig, axes = plt.subplots(n_rows, n_sig_per_row * 2,
                             figsize=(n_sig_per_row * 2 * 1.55, n_rows * 2.3))
    axes = np.array(axes).reshape(n_rows, n_sig_per_row * 2)

    rng = np.random.default_rng(0)
    peak = df_p[(df_p["chrom"] == PEAK_CHROM) & (df_p["start"] == PEAK_START)]

    for r, (row_label, signals) in enumerate(ROWS):
        for s_i, (col, sublabel, ylab) in enumerate(signals):
            ax_raw = axes[r, s_i * 2]
            ax_dev = axes[r, s_i * 2 + 1]

            sub = df_p[df_p[col].notna()].copy()
            null_vals = sub.loc[sub["votes_v3"] == 0, col].dropna().values
            base = null_vals.mean() if len(null_vals) else np.nan
            peak_y = float(peak[col].iloc[0]) if (len(peak) and peak[col].notna().any()) else np.nan

            if not np.isnan(peak_y) and len(null_vals):
                if DIRECTION == "high":
                    n_extreme = int((null_vals >= peak_y).sum())
                else:
                    n_extreme = int((null_vals <= peak_y).sum())
                p_emp = (n_extreme + 1) / (len(null_vals) + 1)
                p_str = f"p={p_emp:.3f}" if p_emp >= 1e-3 else "p<0.001"
                stat_str = f"{p_str} (n_null={len(null_vals)})"
            else:
                stat_str = "p=NA"

            x_jit = 1 + rng.uniform(-0.18, 0.18, len(sub))
            ax_raw.scatter(x_jit, sub[col], s=3, color=GREY, alpha=0.55,
                           edgecolor="none", zorder=1)
            if not np.isnan(peak_y):
                ax_raw.scatter([1], [peak_y], s=22, color=RED,
                               edgecolor="black", linewidth=0.4, zorder=5)
            if not np.isnan(base):
                ax_raw.axhline(base, color="#555555", linestyle="--",
                               linewidth=0.5, alpha=0.7)

            dev = sub[col] - base
            ax_dev.scatter(x_jit, dev, s=3, color=GREY, alpha=0.55,
                           edgecolor="none", zorder=1)
            if not np.isnan(peak_y) and not np.isnan(base):
                ax_dev.scatter([1], [peak_y - base], s=22, color=RED,
                               edgecolor="black", linewidth=0.4, zorder=5)
            ax_dev.axhline(0, color="#555555", linestyle="--",
                           linewidth=0.5, alpha=0.7)

            for ax in (ax_raw, ax_dev):
                ax.set_xlim(0.4, 1.6)
                ax.set_xticks([])
                for sp in ("top", "right"):
                    ax.spines[sp].set_visible(False)
                ax.tick_params(axis="both", labelsize=5.5, length=2)

            if s_i == 0:
                ax_raw.set_ylabel(f"{row_label}\n{ylab}", fontsize=7)
            else:
                ax_raw.set_ylabel(ylab, fontsize=6)
            ax_raw.set_title(f"{sublabel} · raw\n{stat_str}",
                             fontsize=6.0, loc="left", pad=2)
            ax_dev.set_title(f"{sublabel} · Δ vs vote=0",
                             fontsize=6.0, loc="left", pad=2)

    fig.tight_layout()
    fig.savefig(f"{OUT_BASE}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.svg", bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.pdf", bbox_inches="tight")
    svg = Path(f"{OUT_BASE}.svg")
    txt = svg.read_text()
    txt = re.sub(r"<clipPath[^>]*>.*?</clipPath>", "", txt, flags=re.DOTALL)
    txt = re.sub(r'\s*clip-path="url\([^)]+\)"', "", txt)
    svg.write_text(txt)

if __name__ == "__main__":
    main()
