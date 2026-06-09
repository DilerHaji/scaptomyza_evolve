#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(".")

WILD_POS = ROOT / "baypass_wild/wild_snp_positions.csv"
WILD_BF = ROOT / "baypass_wild/wild_trt_summary_betai_reg.out"
BLOCKS = ROOT / "hv_results_cluster/block_dynamics/block_dynamics_summary.tsv"
PBS = ROOT / "variance_analysis/cvtkpy_final/pbs_vs_bt_aligned.tsv"

OUT_PNG = ROOT / "final_plots/wild/wild_lab_convergence_fig3.png"
OUT_SVG = ROOT / "final_plots/wild/wild_lab_convergence_fig3.svg"

SCAFFS = [
    "chr_ScDA7r2_110_HRSCAF_295",
    "chr_ScDA7r2_126_HRSCAF_325",
    "chr_ScDA7r2_439_HRSCAF_779",
    "chr_ScDA7r2_597_HRSCAF_953",
]

WIN = 500_000
STEP = 100_000
MIN_SNPS_PER_WIN = 10
BF_TOP_PCTILE = 0.95           # top-5%
WIN_TOP_PCTILE = 0.95          # top-5% of windows by frac_top5pct (used for Fisher)
SYM_CUTOFF = 0.5

TRT_COLOR = {
    "B": "#499FFF",
    "T": "#EDB72D",
    "B+T": "#9BAB96",
}
TRT_DISPLAY = {"B": "B", "T": "T", "M": "B+T"}

C_NO_BLOCK = "#BBBBBB"
C_LOW_SYM = "#F4BFBF"
C_HIGH_SYM = "#CC3333"
C_NULL = "#7F7F7F"


def load_wild_bf() -> pd.DataFrame:
    pos = pd.read_csv(WILD_POS)
    bf = pd.read_csv(WILD_BF, sep=r"\s+")
    df = pos.merge(
        bf[["MRK", "BF(dB)"]].rename(columns={"MRK": "mrk", "BF(dB)": "bf_db"}),
        on="mrk",
    )
    df = df[df["chrom"].isin(SCAFFS)].copy()
    thr = df["bf_db"].quantile(BF_TOP_PCTILE)
    df["is_top5pct"] = df["bf_db"] >= thr
    return df


def load_blocks_with_pbs() -> pd.DataFrame:
    dyn = pd.read_csv(BLOCKS, sep="\t")
    pbs = pd.read_csv(PBS, sep="\t").rename(columns={"chrom": "chr"})
    rz_b, rz_t = [], []
    for _, b in dyn.iterrows():
        sub = pbs[(pbs["chr"] == b["chr"]) &
                  (pbs["start"] <= b["end"]) &
                  (pbs["end"] >= b["start"])]
        if len(sub):
            rz_b.append(sub["rz_pbs_B"].mean(skipna=True))
            rz_t.append(sub["rz_pbs_T"].mean(skipna=True))
        else:
            rz_b.append(np.nan)
            rz_t.append(np.nan)
    dyn["rz_pbs_B"] = rz_b
    dyn["rz_pbs_T"] = rz_t
    dyn = dyn.dropna(subset=["rz_pbs_B", "rz_pbs_T"]).copy()

    combined_z = (dyn["rz_pbs_B"] + dyn["rz_pbs_T"]) / 2
    denom = (dyn["rz_pbs_B"] + dyn["rz_pbs_T"]).replace(0, np.nan)
    direction = np.where(combined_z > 0,
                         (dyn["rz_pbs_T"] - dyn["rz_pbs_B"]) / denom, 0)
    direction = np.where(np.isnan(direction), 0, direction)
    sym_score = np.where(combined_z > 0,
                         combined_z * (1 - np.minimum(np.abs(direction), 1)), 0)
    dyn["combined_z"] = combined_z
    dyn["sym_score"] = sym_score
    dyn["treatment_disp"] = dyn["treatment"].map(TRT_DISPLAY)
    return dyn[dyn["treatment_disp"].isin(TRT_COLOR.keys())].copy()


def build_windows(wild: pd.DataFrame, blocks: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for ch in SCAFFS:
        wc = wild[wild["chrom"] == ch]
        if wc.empty:
            continue
        max_pos = int(wc["pos"].max())
        starts = np.arange(0, max_pos, STEP, dtype=int)
        for s in starts:
            e = s + WIN
            sub = wc[(wc["pos"] >= s) & (wc["pos"] < e)]
            if len(sub) < MIN_SNPS_PER_WIN:
                continue
            row = {
                "chr": ch, "start": s, "end": e,
                "n_snps": len(sub),
                "frac_top5pct": sub["is_top5pct"].mean(),
            }
            for trt in TRT_COLOR:
                ov = blocks[
                    (blocks["chr"] == ch) &
                    (blocks["treatment_disp"] == trt) &
                    (blocks["start"] <= e) &
                    (blocks["end"] >= s)
                ]
                if len(ov):
                    row[f"max_sym_{trt}"] = float(ov["sym_score"].max())
                    row[f"any_block_{trt}"] = True
                else:
                    row[f"max_sym_{trt}"] = 0.0
                    row[f"any_block_{trt}"] = False
            rows.append(row)
    return pd.DataFrame(rows)


def scatter_panel(ax, wins: pd.DataFrame, trt: str,
                  max_sym_overall: float) -> tuple[float, float]:
    x = wins[f"max_sym_{trt}"].values
    y = wins["frac_top5pct"].values * 100.0

    cat = np.where(x > SYM_CUTOFF, C_HIGH_SYM,
                   np.where(x > 0, C_LOW_SYM, C_NO_BLOCK))
    ax.scatter(x, y, s=14, c=cat, edgecolor="black",
               linewidth=0.15, alpha=0.8, zorder=2)

    valid = ~np.isnan(x) & ~np.isnan(y)
    rho, pval = stats.spearmanr(x[valid], y[valid])

    has_block = wins[f"any_block_{trt}"].values
    if has_block.sum() > 10:
        m, b = np.polyfit(x[has_block], y[has_block], 1)
        xs = np.linspace(0, max_sym_overall, 100)
        ax.plot(xs, m * xs + b, color=TRT_COLOR[trt],
                linewidth=1.8, zorder=3)

    ax.axhline(5.0, color=C_NULL, linestyle=":", linewidth=0.8, zorder=1)
    ax.set_xlim(-0.04, max_sym_overall + 0.05)
    ax.set_ylim(-0.8, max(y.max() * 1.15, 12))
    ax.set_xlabel("Max lab-block symmetry score\n(per 500-kb window)", fontsize=9)
    ax.set_ylabel("% of SNPs per window\nin wild top-5% BF(dB)", fontsize=9)
    ax.tick_params(labelsize=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax.text(
        0.03, 0.97,
        f"Spearman $\\rho$ = {rho:+.2f}\n"
        f"p = {pval:.1e}",
        transform=ax.transAxes, ha="left", va="top",
        fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.25",
                  facecolor="white", edgecolor=TRT_COLOR[trt], linewidth=0.8),
    )

    ax.set_title(trt, loc="left", fontsize=12, fontweight="bold",
                 color=TRT_COLOR[trt])

    return float(rho), float(pval)


def enrichment_panel(ax, wins: pd.DataFrame, trt: str) -> tuple[float, float, int, int, int]:
    x = wins[f"max_sym_{trt}"].values
    cat = np.where(x > SYM_CUTOFF, "sym>0.5",
                   np.where(x > 0, "0<sym≤0.5", "no block"))

    thr = wins["frac_top5pct"].quantile(WIN_TOP_PCTILE)
    is_top_wild = wins["frac_top5pct"].values >= thr

    categories = ["no block", "0<sym≤0.5", "sym>0.5"]
    colors = [C_NO_BLOCK, C_LOW_SYM, C_HIGH_SYM]
    n_per = []
    pct_per = []
    for c in categories:
        mask = cat == c
        n = mask.sum()
        pct = is_top_wild[mask].mean() * 100 if n else 0
        n_per.append(n)
        pct_per.append(pct)

    xpos = np.arange(len(categories))
    bars = ax.bar(xpos, pct_per, color=colors, edgecolor="black", linewidth=0.6)
    for i, (n, p) in enumerate(zip(n_per, pct_per)):
        ax.text(i, p + max(pct_per) * 0.04,
                f"n = {n}\n{p:.1f}%",
                ha="center", va="bottom", fontsize=8)
    ax.axhline(5.0, color=C_NULL, linestyle=":", linewidth=0.8, zorder=1)
    ax.set_xticks(xpos)
    ax.set_xticklabels(categories, fontsize=8.5)
    ax.set_ylabel("% of windows in top-5%\nby wild BF density", fontsize=9)
    ax.set_ylim(0, max(pct_per) * 1.35 + 1)
    ax.tick_params(axis="y", labelsize=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    a = int(((cat == "sym>0.5") & is_top_wild).sum())
    b = int(((cat == "sym>0.5") & ~is_top_wild).sum())
    c = int(((cat != "sym>0.5") & is_top_wild).sum())
    d = int(((cat != "sym>0.5") & ~is_top_wild).sum())
    table = np.array([[a, b], [c, d]])
    if a > 0 and b > 0 and c > 0 and d > 0:
        odds, p = stats.fisher_exact(table, alternative="greater")
    else:
        odds, p = (a * d) / max(b * c, 1), 1.0

    ax.text(
        0.97, 0.97,
        f"OR = {odds:.1f}\np = {p:.1e}",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.25",
                  facecolor="white", edgecolor=TRT_COLOR[trt], linewidth=0.8),
    )

    return odds, p, n_per[2], n_per[1], n_per[0]


def main() -> None:
    plt.rcParams.update({
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.8,
    })

    wild = load_wild_bf()

    blocks = load_blocks_with_pbs()

    wins = build_windows(wild, blocks)

    max_sym = max(blocks["sym_score"].max(), 0.5)

    fig, axes = plt.subplots(
        2, 3, figsize=(12.5, 7.2),
        gridspec_kw={"height_ratios": [1.15, 1.0], "hspace": 0.45, "wspace": 0.32},
    )

    for j, trt in enumerate(["B", "T", "B+T"]):
        rho, pval = scatter_panel(axes[0, j], wins, trt, max_sym)

    for j, trt in enumerate(["B", "T", "B+T"]):
        odds, p, n_high, n_low, n_none = enrichment_panel(axes[1, j], wins, trt)

    axes[0, 0].text(-0.17, 1.08, "A", transform=axes[0, 0].transAxes,
                    fontsize=16, fontweight="bold", ha="left", va="top")
    axes[1, 0].text(-0.17, 1.12, "B", transform=axes[1, 0].transAxes,
                    fontsize=16, fontweight="bold", ha="left", va="top")

    fig.suptitle(
        "Wild host-associated BF overlaps lab selection signal, by treatment",
        fontsize=12, y=0.995, x=0.05, ha="left", fontweight="bold",
    )

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")

if __name__ == "__main__":
    main()
