#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"

PALETTE = {"B": "#96CDFF", "F": "#B1B1B1", "M": "#901442", "T": "#F3C43C"}
TREAT_LABEL = {"B": "Barbarea", "T": "Turitus", "M": "B+T (Mixed)", "F": "Founders"}


def classify(name):
    m = re.match(r"^([BTM])([1-4])G(\d{2})$", name)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    m = re.match(r"^F([1-4])(G00)?$", name)
    if m:
        return "F", int(m.group(1)), 0
    return None, None, None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alt-dir", required=True,
                    help="matched_null_*_alt_direction.tsv (has both direction "
                         "flavors per sample for all generations)")
    ap.add_argument("--scores",  required=True,
                    help="levene_scores_81_scores.tsv (has per-rep trajectory, "
                         "used here only for gen metadata)")
    ap.add_argument("--null",    required=True)
    ap.add_argument("--het-asc", required=True)
    ap.add_argument("--het-mat", required=True)
    ap.add_argument("--out",     required=True)
    args = ap.parse_args()

    alt = pd.read_csv(args.alt_dir, sep="\t")
    null = pd.read_csv(args.null, sep="\t")
    het_asc = pd.read_csv(args.het_asc, sep="\t")
    het_mat = pd.read_csv(args.het_mat, sep="\t")

    alt["score"] = alt["score_g10only_direction"]

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    ax = axes[0, 0]
    for treat, color in PALETTE.items():
        sub = alt[alt["treat"] == treat]
        if len(sub) == 0:
            continue
        first = True
        for rp, g in sub.groupby("rp"):
            g = g.sort_values("gen")
            ax.plot(g["gen"], g["score"], color=color, marker="o",
                    markersize=4, linewidth=1.0, alpha=0.75,
                    label=TREAT_LABEL[treat] if first else None)
            first = False
    ax.axhline(0, color="#888", linewidth=0.6, linestyle="--")
    ax.set_xlabel("Generation", fontsize=11)
    ax.set_ylabel("B-T projection score\n(G10-only direction vector)", fontsize=10)
    ax.set_title("A - Projection score trajectories", fontsize=11, loc="left")
    ax.legend(fontsize=9, frameon=False, loc="upper left")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)

    ax = axes[0, 1]
    positions = {"F": 0, "B": 1, "M": 2, "T": 3}
    for treat, pos in positions.items():
        if treat == "F":
            sub = alt[alt["treat"] == treat]
        else:
            sub = alt[(alt["treat"] == treat) & (alt["gen"] == 10)]
        if len(sub) == 0:
            continue
        jitter = 0.08 * (pd.Series(range(len(sub))).mod(3) - 1).values
        ax.scatter(np.full(len(sub), pos) + jitter, sub["score"], s=55,
                   c=PALETTE[treat], edgecolors="black", linewidths=0.5,
                   alpha=0.85, label=TREAT_LABEL[treat])
    ax.axhline(0, color="#888", linewidth=0.6, linestyle="--")
    ax.set_xticks(list(positions.values()))
    ax.set_xticklabels([TREAT_LABEL[t] for t in positions.keys()],
                       rotation=15, ha="right")
    ax.set_ylabel("Endpoint score\n(G10; F at G0)", fontsize=10)
    ax.set_title("B - Endpoint projection by treatment", fontsize=11, loc="left")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)

    ax = axes[1, 0]
    m_cols = [f"M{rp}G10" for rp in [1, 2, 3, 4] if f"M{rp}G10" in null.columns]
    null_std = null[m_cols].std(axis=1, ddof=1)

    obs_m_g10 = alt[(alt["treat"] == "M") & (alt["gen"] == 10)]
    obs_std = float(obs_m_g10["score_delta_direction"].std(ddof=1))

    ax.hist(null_std, bins=40, color="#CCCCCC", edgecolor="white",
            label="Matched-null draws (n=1000)")
    ax.axvline(obs_std, color="#901442", linewidth=2.5,
               label=f"Observed (ascertained 81)\nstd = {obs_std:.5f}")
    p_low = float((null_std <= obs_std).mean())
    ax.text(0.98, 0.95,
            f"p(null ≤ obs) = {p_low:.3f}\nnull median = {null_std.median():.5f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=10)
    ax.set_xlabel("Across-replicate std of M G10 score", fontsize=10)
    ax.set_ylabel("# null draws", fontsize=10)
    ax.set_title("C - M replicate convergence at ascertained windows",
                 fontsize=11, loc="left")
    ax.legend(fontsize=9, frameon=False, loc="upper left")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)

    ax = axes[1, 1]
    r_asc = het_asc["ratio_M_F"].replace([np.inf, -np.inf], np.nan).dropna()
    r_mat = het_mat["ratio_M_F"].replace([np.inf, -np.inf], np.nan).dropna()
    bins = np.linspace(0, 2.0, 80)
    ax.hist(r_asc, bins=bins, density=True, alpha=0.55, color="#901442",
            label=f"Ascertained (n={len(r_asc):,}, median={r_asc.median():.3f})")
    ax.hist(r_mat, bins=bins, density=True, alpha=0.40, color="#888888",
            label=f"Matched-AF genome (n={len(r_mat):,}, median={r_mat.median():.3f})")
    ax.axvline(1.0, color="#444", linewidth=0.6, linestyle="--")
    ax.set_xlabel("H_M_G10 / H_F  (per-SNP heterozygosity ratio)", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_title("D - M heterozygosity preservation vs founders",
                 fontsize=11, loc="left")
    ax.legend(fontsize=9, frameon=False, loc="upper right")
    ax.set_xlim(0, 2.0)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)

    plt.tight_layout()
    out_base = re.sub(r"\.(png|svg|pdf)$", "", args.out, flags=re.IGNORECASE)
    Path(out_base).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_base}.png", dpi=180, bbox_inches="tight")
    fig.savefig(f"{out_base}.svg",            bbox_inches="tight")

if __name__ == "__main__":
    main()
