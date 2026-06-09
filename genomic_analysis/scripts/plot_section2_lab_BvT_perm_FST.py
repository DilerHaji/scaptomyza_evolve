#!/usr/bin/env python3

from __future__ import annotations
from pathlib import Path
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(".")
AF_MATRIX = ROOT / "final_plots/wild/af_matrix_22pools.csv"
OUT_DIR = ROOT / "final_plots/wild"
OUT_BASE = OUT_DIR / "section2_lab_BvT_permFST_manhattan"
OUT_WIN_TSV = OUT_DIR / "section2_lab_BvT_permFST_windows.tsv"
OUT_TOP_TSV = OUT_DIR / "section2_lab_BvT_permFST_top.tsv"

POOL_COLS = ["B1G10", "B2G10", "B3G10", "B4G10",
              "T1G10", "T2G10", "T3G10", "T4G10"]
TRUE_SPLIT = frozenset([0, 1, 2, 3])   # B indices in POOL_COLS

N_EFF = 58    # pool-seq-corrected effective chromosomes per pool (supplement)
WINDOW_BP = 200_000
MIN_SNPS_PER_WINDOW = 20
TOP_N_TABLE = 40


def hudson_numden(p1: np.ndarray, p2: np.ndarray, n1=N_EFF, n2=N_EFF):
    num = (p1 - p2) ** 2 - (p1 * (1 - p1)) / (n1 - 1) - (p2 * (1 - p2)) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return num, den


def main() -> None:
    plt.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.8,
    })

    af = pd.read_csv(AF_MATRIX, usecols=["chrom_pos", *POOL_COLS])
    af[["chrom", "pos"]] = af["chrom_pos"].str.split(":", expand=True)
    af["pos"] = af["pos"].astype(int)
    af["bin_start"] = (af["pos"] // WINDOW_BP) * WINDOW_BP

    P = af[POOL_COLS].to_numpy(dtype=float)
    n_snps = P.shape[0]

    pairs = list(combinations(range(8), 2))

    pair_num = np.empty((len(pairs), n_snps))
    pair_den = np.empty((len(pairs), n_snps))
    for k, (i, j) in enumerate(pairs):
        num, den = hudson_numden(P[:, i], P[:, j])
        pair_num[k] = num
        pair_den[k] = den

    keys = af[["chrom", "bin_start"]].to_numpy()
    agg = pd.DataFrame({"chrom": af["chrom"], "bin_start": af["bin_start"]})
    for k in range(len(pairs)):
        agg[f"num_{k}"] = pair_num[k]
        agg[f"den_{k}"] = pair_den[k]
    agg_sum = agg.groupby(["chrom", "bin_start"], as_index=False).agg(
        {**{f"num_{k}": "sum" for k in range(len(pairs))},
         **{f"den_{k}": "sum" for k in range(len(pairs))},
        }
    )
    agg_sum["n_snps"] = af.groupby(["chrom", "bin_start"]).size().values
    agg_sum = agg_sum[agg_sum["n_snps"] >= MIN_SNPS_PER_WINDOW].reset_index(drop=True)

    num_w = agg_sum[[f"num_{k}" for k in range(len(pairs))]].to_numpy()
    den_w = agg_sum[[f"den_{k}" for k in range(len(pairs))]].to_numpy()

    splits = list(combinations(range(8), 4))
    cross_mask = np.zeros((len(splits), len(pairs)), dtype=bool)
    for s_idx, group_a in enumerate(splits):
        group_a_set = set(group_a)
        for k, (i, j) in enumerate(pairs):
            if (i in group_a_set) ^ (j in group_a_set):
                cross_mask[s_idx, k] = True
    assert cross_mask.sum(axis=1).max() == 16
    assert cross_mask.sum(axis=1).min() == 16  # each 4/4 split has exactly 16 cross-pairs


    sum_cross_num = num_w @ cross_mask.T
    sum_cross_den = den_w @ cross_mask.T
    with np.errstate(invalid="ignore", divide="ignore"):
        fst_splits = sum_cross_num / np.maximum(sum_cross_den, 1e-12)

    true_split_idx = splits.index(tuple(sorted(TRUE_SPLIT)))
    obs_fst = fst_splits[:, true_split_idx]

    ranks = (fst_splits >= obs_fst[:, None]).sum(axis=1)
    emp_p = ranks / fst_splits.shape[1]

    agg_sum["fst_BvT"] = obs_fst
    agg_sum["fst_median_perm"] = np.nanmedian(fst_splits, axis=1)
    agg_sum["fst_ratio"] = obs_fst / np.maximum(agg_sum["fst_median_perm"].to_numpy(), 1e-9)
    agg_sum["emp_p"] = emp_p
    agg_sum["neg_log10_p"] = -np.log10(np.maximum(emp_p, 1.0 / fst_splits.shape[1]))


    top = agg_sum.sort_values("fst_BvT", ascending=False).head(TOP_N_TABLE)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    agg_sum[["chrom", "bin_start", "n_snps", "fst_BvT",
              "fst_median_perm", "fst_ratio", "emp_p", "neg_log10_p"]] \
        .to_csv(OUT_WIN_TSV, sep="\t", index=False)
    top[["chrom", "bin_start", "n_snps", "fst_BvT",
          "fst_median_perm", "fst_ratio", "emp_p", "neg_log10_p"]] \
        .to_csv(OUT_TOP_TSV, sep="\t", index=False)

    for _, r in top.head(10).iterrows():
        sc = r["chrom"].split("_")[2] if "_" in r["chrom"] else r["chrom"]







    keep = (agg_sum.groupby("chrom").size()
                     .loc[lambda s: s >= 10].index.tolist())
    keep = sorted(keep, key=lambda s: -agg_sum.loc[agg_sum["chrom"] == s,
                                                     "bin_start"].max())
    win_p = agg_sum[agg_sum["chrom"].isin(keep)].copy()

    offsets = {}
    cum = 0
    gap = 5_000_000
    for s in keep:
        offsets[s] = cum
        cum += agg_sum.loc[agg_sum["chrom"] == s, "bin_start"].max() + WINDOW_BP + gap
    win_p["x"] = win_p.apply(lambda r: offsets[r["chrom"]] + r["bin_start"], axis=1)

    fig, axes = plt.subplots(2, 1, figsize=(13, 6.5), sharex=True,
                              gridspec_kw=dict(hspace=0.10))

    ax1 = axes[0]
    pal = ["#888888", "#444444"]
    for i, s in enumerate(keep):
        sub = win_p[win_p["chrom"] == s]
        ax1.scatter(sub["x"], sub["fst_BvT"],
                    s=8, color=pal[i % 2], linewidths=0,
                    alpha=0.55, zorder=1)

    sig = win_p[win_p["emp_p"] <= 5.0 / 70]
    ax1.scatter(sig["x"], sig["fst_BvT"], s=16, color="#C84A45",
                edgecolor="black", linewidth=0.3, zorder=3,
                label=f"emp p <= 5/70 (n = {len(sig)})")

    top_label = win_p.sort_values("fst_BvT", ascending=False).head(8)
    for _, r in top_label.iterrows():
        sc = r["chrom"].split("_")[2] if "_" in r["chrom"] else r["chrom"]
        ax1.annotate(f"{sc}:{r['bin_start']/1e6:.1f}M",
                      xy=(r["x"], r["fst_BvT"]),
                      xytext=(0, 6), textcoords="offset points",
                      fontsize=7.5, ha="center", color="#C84A45",
                      fontweight="bold")

    ax1.set_ylabel("Hudson FST(B, T) at G10\n(per 200 kb window)", fontsize=9.5)
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)
    ax1.tick_params(axis="y", labelsize=9)
    ax1.legend(loc="upper right", frameon=False, fontsize=8)
    ax1.set_title("Section 2 — lab B vs T genome scan (G10, permutation null; 70 splits of 8 pools)",
                  loc="left", fontsize=11, fontweight="bold", pad=6)

    ax2 = axes[1]
    for i, s in enumerate(keep):
        sub = win_p[win_p["chrom"] == s]
        ax2.scatter(sub["x"], sub["neg_log10_p"],
                    s=8, color=pal[i % 2], linewidths=0,
                    alpha=0.55, zorder=1)
    ax2.scatter(sig["x"], sig["neg_log10_p"], s=16, color="#C84A45",
                edgecolor="black", linewidth=0.3, zorder=3)
    ax2.axhline(-np.log10(5.0 / 70), color="#C84A45", linestyle="--",
                 linewidth=0.8, alpha=0.6, zorder=0)

    tick_x = []
    for s in keep:
        sub_x = win_p.loc[win_p["chrom"] == s, "x"]
        tick_x.append(sub_x.mean())
    short_labels = [s.split("_")[2] if "_" in s else s for s in keep]
    ax2.set_xticks(tick_x)
    ax2.set_xticklabels(short_labels, fontsize=7)
    ax2.set_xlabel("Genome (scaffolds, ordered by length)", fontsize=10)
    ax2.set_ylabel("-log10(empirical p) from\n70-split permutation", fontsize=9.5)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)
    ax2.tick_params(axis="y", labelsize=9)

    import re
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
