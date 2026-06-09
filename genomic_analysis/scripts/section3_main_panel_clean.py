#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import re
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from section3_x_enrichment_analysis import load_q_pred_with_filters

ROOT = Path(".")
F100 = ROOT / "grenfst/diversity_attrition/wild_pi_100000_n100_diversity.csv"
OUT_BASE = ROOT / "final_plots/wild/section3_fig4_main_clean"

WILD = ["AVB", "AVT", "PSB", "PST", "RMB", "RMT"]
CHROM_439 = "chr_ScDA7r2_439_HRSCAF_779"
LD_S, LD_E = 2_640_000, 3_610_000


def strip_svg(svg_path: Path):
    txt = svg_path.read_text()
    txt = re.sub(r"<clipPath\b[^>]*>.*?</clipPath>", "", txt, flags=re.DOTALL)
    txt = re.sub(r"\s+clip-path\s*=\s*['\"]url\(#[^)]+\)['\"]", "", txt)
    txt = txt.replace('<svg ',
        '<svg xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" ', 1)
    txt = re.sub(r'(<g\s+id="([^"]+)")',
                 r'\1 inkscape:label="\2" inkscape:groupmode="layer"', txt)
    svg_path.write_text(txt)


def assign_local_D(df, w):
    out = df.copy()
    out["local_D"] = np.nan
    for chrom, sub in out.groupby("chrom"):
        ww = w[w["chrom"] == chrom].sort_values("start").reset_index(drop=True)
        if ww.empty:
            continue
        starts = ww["start"].values
        ends = ww["end"].values
        positions = sub["pos"].values.astype(np.int64)
        idx = np.clip(np.searchsorted(starts, positions, side="right") - 1,
                       0, len(starts) - 1)
        in_w = (starts[idx] <= positions) & (positions <= ends[idx])
        D_vals = ww["D"].values[idx]
        D_vals = np.where(in_w, D_vals, np.nan)
        out.loc[sub.index, "local_D"] = D_vals
    return out


def main():
    plt.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.7,
    })

    df, _ = load_q_pred_with_filters()
    df["fg_pass"] = df["bg_pass"] & (df["within_host_SD"] <= 0.07)

    print("loading wild diversity for local D...")
    w = pd.read_csv(F100)
    w = w[w["total.passed"] >= 30].copy()
    for arm, suf in [("pi","theta_pi"),("thW","theta_watterson"),
                      ("D","tajimas_d")]:
        w[arm] = w[[f"{p}.1.{suf}" for p in WILD]].mean(axis=1)

    df = assign_local_D(df, w)

    bg = df[df["bg_pass"] & ~df["fg_pass"]]
    fg = df[df["fg_pass"] & df["local_D"].notna()]
    fg_chr439_LD = fg[(fg["chrom"] == CHROM_439)
                      & (fg["pos"].between(LD_S, LD_E))]

    fig, ax = plt.subplots(figsize=(6.5, 6.0))

    ax.scatter(bg["qstar_pred"], bg["wild_AF_pol"],
                s=3, color="#cccccc", alpha=0.35, edgecolors="none",
                rasterized=True, zorder=1)

    norm = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
    sc = ax.scatter(fg["qstar_pred"], fg["wild_AF_pol"],
                     c=fg["local_D"], cmap="RdBu_r", norm=norm,
                     s=14, edgecolors="none", alpha=0.85, zorder=3)

    ax.scatter(fg_chr439_LD["qstar_pred"], fg_chr439_LD["wild_AF_pol"],
                marker="D", s=70, c=fg_chr439_LD["local_D"], cmap="RdBu_r",
                norm=norm, edgecolor="black", linewidth=0.9, zorder=10)

    ax.plot([0, 1], [0, 1], color="black", linewidth=1.0,
             linestyle="--", alpha=0.55, zorder=2)

    x = fg["qstar_pred"].values
    y = fg["wild_AF_pol"].values
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() >= 3:
        b, a = np.polyfit(x[keep], y[keep], 1)
        xs = np.linspace(0.05, 0.95, 50)
        ax.plot(xs, a + b * xs, color="#222222", linewidth=2.0, zorder=8)

    # Stats
    rho, p = spearmanr(x[keep], y[keep])
    n_fg = int(keep.sum())
    n_LD = len(fg_chr439_LD)
    ax.text(0.03, 0.97,
             f"foreground (n={n_fg:,}):  ρ = {rho:+.3f}   p = {p:.1e}\n"
             f"chr_439 LD-block subset (♦):  n = {n_LD}",
             transform=ax.transAxes, ha="left", va="top",
             fontsize=8.5, family="monospace",
             bbox=dict(facecolor="white", alpha=0.95, edgecolor="#aaa",
                        linewidth=0.4, boxstyle="round,pad=0.4"))

    cbar = plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Local wild Tajima's D", fontsize=8.5)
    cbar.ax.tick_params(labelsize=7.5)

    handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                    markerfacecolor="#cccccc", markeredgecolor="none",
                    markersize=4, label=f"background  (n={len(bg):,})"),
        plt.Line2D([0], [0], marker="o", color="w",
                    markerfacecolor="#888888", markeredgecolor="none",
                    markersize=5, label=f"foreground  (n={n_fg:,})"),
        plt.Line2D([0], [0], marker="D", color="w",
                    markerfacecolor="#888888", markeredgecolor="black",
                    markersize=7, label=f"chr_439 LD block  (n={n_LD})"),
        plt.Line2D([0], [0], color="black", linewidth=2.0, label="OLS fit"),
        plt.Line2D([0], [0], color="black", linewidth=1.0, linestyle="--",
                    alpha=0.6, label="y = x"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=8,
              frameon=False, bbox_to_anchor=(1.0, 0.0))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"Lab Dempster $q^{*}_{\rm pred}$", fontsize=10)
    ax.set_ylabel("Wild allele frequency (polarized)", fontsize=10)
    ax.set_title(
        "Lab-derived equilibrium predicts wild allele frequency",
        fontsize=11, pad=6)
    ax.tick_params(labelsize=9)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)

    fig.tight_layout()
    fig.savefig(f"{OUT_BASE}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.svg", bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.pdf", bbox_inches="tight")
    strip_svg(Path(f"{OUT_BASE}.svg"))


if __name__ == "__main__":
    main()
