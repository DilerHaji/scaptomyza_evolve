#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42


PALETTE_BAND = ["#bdbdbd", "#2f2f2f"]         # alternating scaffold colors
HL_ANTAG     = "#D76161"                      # negative cov_BT highlight
HL_M_BIMODAL = "#1f78b4"                      # high var_M highlight


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--windows", required=True)
    ap.add_argument("--out",     required=True,
                    help="Output prefix. .png + .svg written.")
    ap.add_argument("--top-n",   type=int, default=15,
                    help="Annotate top-N composite-score windows (default 15)")
    ap.add_argument("--min-scaffolds", type=int, default=50,
                    help="Only show scaffolds with at least this many windows (default 50)")
    ap.add_argument("--min-snps",      type=int, default=0,
                    help="Only consider windows with at least this many SNPs (default 0 = no filter)")
    return ap.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.windows, sep="\t")

    scaffold_counts = df["chrom"].value_counts()
    keep = scaffold_counts[scaffold_counts >= args.min_scaffolds].index.tolist()
    df = df[df["chrom"].isin(keep)].copy()

    if args.min_snps > 0:
        before = len(df)
        df = df[df["n_snps"] >= args.min_snps].reset_index(drop=True)
    if len(df) == 0:
        sys.exit("no windows left after scaffold filter")

    scaffold_order = scaffold_counts.loc[keep].sort_values(ascending=False).index.tolist()
    offsets = {}
    cum = 0
    gap = 2_000_000
    scaffold_midpoints = []
    for i, sc in enumerate(scaffold_order):
        sub = df[df["chrom"] == sc]
        offsets[sc] = cum
        sc_max = sub["start"].max() + 200_000
        scaffold_midpoints.append((sc, cum + sc_max / 2))
        cum += sc_max + gap

    df["x"] = df.apply(lambda r: offsets[r["chrom"]] + r["start"], axis=1)
    df["band"] = df["chrom"].apply(lambda c: scaffold_order.index(c) % 2)

    if ("z_perm" in df.columns and df["z_perm"].notna().any()
            and "ratio_M_over_F" in df.columns and "n_genes" in df.columns):
        gene_mask = df["n_genes"] >= 10
        df["levene_score"] = -df["z_perm"] * df["ratio_M_over_F"]
        df.loc[~gene_mask, "levene_score"] = np.nan  # hide gene-poor
        score_label = "−z_perm × (var_M / var_F), gene-dense only"
    elif "z_cov_BT" in df.columns and "ratio_M_over_F" in df.columns:
        df["levene_score"] = -df["z_cov_BT"] * df["ratio_M_over_F"]
        score_label = "−z_cov_BT × (var_M / var_F)"
    elif "ratio_M_over_F" in df.columns and df["ratio_M_over_F"].notna().any():
        df["levene_score"] = -df["cov_BT"] * df["ratio_M_over_F"]
        score_label = "−cov_BT × (var_M / var_F)"
    else:
        df["levene_score"] = -df["cov_BT"] * df["across_rep_var_M"]
        score_label = "−cov_BT × var_M"
    top_hits = df.nlargest(args.top_n, "levene_score")

    def short_sc(s):
        m = re.search(r"_(\d+)_HRSCAF", s)
        return f"sc{m.group(1)}" if m else s

    has_ratio = "ratio_M_over_F" in df.columns and df["ratio_M_over_F"].notna().any()
    has_genes = "n_genes" in df.columns

    extra_rows = [True, True, has_genes, True]  # always include SNP density at bottom
    n_rows = sum(extra_rows)
    fig, axes = plt.subplots(n_rows, 1, figsize=(16, 3.0 * n_rows), sharex=True)
    axA = axes[0]
    axB = axes[1]
    axC = axes[2] if has_genes else None
    axD = axes[-1]

    has_zperm = "z_perm" in df.columns and df["z_perm"].notna().any()
    if has_zperm:
        ycol_A = "z_perm"
        ylab_A = "−z_perm per window\n(up = antagonistic beyond matched-n permutation null)"
        thresh_line = 3  # z_perm < -3 roughly matches p_perm_neg < 0.001
    elif "z_cov_BT" in df.columns:
        ycol_A = "z_cov_BT"
        ylab_A = "−z_cov_BT per window\n(up = antagonistic, |z|>2 statistically resolvable)"
        thresh_line = 2
    else:
        ycol_A = "cov_BT"
        ylab_A = "−cov(ΔAF_B, ΔAF_T) per window\n(up = antagonistic)"
        thresh_line = None
    use_z = ycol_A.startswith("z_")
    for band_id in (0, 1):
        sub = df[df["band"] == band_id]
        axA.scatter(sub["x"], -sub[ycol_A], s=5, c=PALETTE_BAND[band_id],
                    alpha=0.7, edgecolors="none", rasterized=True)
    if use_z:
        antag_mask = df[ycol_A] < -thresh_line
    else:
        antag_mask = df["cov_BT"] < -df["cov_BT"].abs().quantile(0.95)
    antag = df[antag_mask]
    axA.scatter(antag["x"], -antag[ycol_A], s=18, facecolors="none",
                edgecolors=HL_ANTAG, linewidths=0.8,
                label=f"{ycol_A} < -{thresh_line}" if use_z else f"{ycol_A} bottom-5%")

    if "n_genes" in df.columns and ycol_A == "z_perm":
        cand_mask = (df["z_perm"] < -thresh_line) & (df["n_genes"] >= 10)
        cand = df[cand_mask]
        axA.scatter(cand["x"], -cand[ycol_A], s=32, facecolors="#FFD700",
                    edgecolors="black", linewidths=0.6,
                    label=f"candidate (z_perm<−{thresh_line} AND n_genes≥10, n={len(cand)})",
                    zorder=5)
    axA.axhline(0, color="#888", linewidth=0.5, linestyle=":")
    if thresh_line is not None:
        axA.axhline(thresh_line, color=HL_ANTAG, linewidth=0.6, linestyle="--")
    axA.set_ylabel(ylab_A, fontsize=10)
    axA.legend(loc="upper right", fontsize=8, frameon=False)
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)

    if "z_perm_bimodal" in df.columns and df["z_perm_bimodal"].notna().any():
        ycol_B = "z_perm_bimodal"
        ylabel_B = "z_perm_bimodal per window\n(up = var_M/var_F beyond matched-n permutation null)"
    elif "z_log_ratio_M_over_F" in df.columns and df["z_log_ratio_M_over_F"].notna().any():
        ycol_B = "z_log_ratio_M_over_F"
        ylabel_B = "z(log var_M / var_F)\n(up = M reps disagree BEYOND founder pool-seq noise, SE-corrected)"
    elif has_ratio:
        ycol_B = "ratio_M_over_F"
        ylabel_B = "Across-rep var(ΔAF_M) / var(AF_F)\n(up = M reps disagree BEYOND founder pool-seq noise)"
    else:
        ycol_B = "across_rep_var_M"
        ylabel_B = "Mean across-rep var(ΔAF_M) per window\n(up = M reps disagree)"
    for band_id in (0, 1):
        sub = df[df["band"] == band_id]
        axB.scatter(sub["x"], sub[ycol_B], s=5, c=PALETTE_BAND[band_id],
                    alpha=0.7, edgecolors="none", rasterized=True)
    bimodal = df[df[ycol_B] > df[ycol_B].quantile(0.95)]
    axB.scatter(bimodal["x"], bimodal[ycol_B], s=18,
                facecolors="none", edgecolors=HL_M_BIMODAL, linewidths=0.8,
                label=f"{ycol_B} top-5%")
    axB.set_ylabel(ylabel_B, fontsize=10)
    axB.legend(loc="upper right", fontsize=8, frameon=False)
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)

    if axC is not None:
        for band_id in (0, 1):
            sub = df[df["band"] == band_id]
            axC.bar(sub["x"], sub["n_genes"], width=150_000, color=PALETTE_BAND[band_id],
                    alpha=0.7, edgecolor="none")
        axC.set_ylabel("n genes / 200 kb window", fontsize=10)
        axC.spines["top"].set_visible(False)
        axC.spines["right"].set_visible(False)

    ycol_B_ann = ycol_B
    for _, r in top_hits.iterrows():
        label = f"{short_sc(r['chrom'])}:{int(r['start']/1e6)}Mb"
        y_a = -r[ycol_A]
        axA.annotate(label, xy=(r["x"], y_a),
                     xytext=(4, 4), textcoords="offset points",
                     fontsize=6, color="#7F0000", alpha=0.95)
        axB.annotate(label, xy=(r["x"], r[ycol_B_ann]),
                     xytext=(4, 4), textcoords="offset points",
                     fontsize=6, color="#0A3D91", alpha=0.95)

    for band_id in (0, 1):
        sub = df[df["band"] == band_id]
        axD.bar(sub["x"], sub["n_snps"], width=150_000, color=PALETTE_BAND[band_id],
                alpha=0.7, edgecolor="none")
    axD.axhline(df["n_snps"].median(), color="#444", linewidth=0.6, linestyle="--",
                label=f"genome median = {int(df['n_snps'].median())}")
    axD.set_ylabel("n SNPs / 200 kb window", fontsize=10)
    axD.legend(loc="upper right", fontsize=8, frameon=False)
    axD.spines["top"].set_visible(False)
    axD.spines["right"].set_visible(False)

    axD.set_xticks([m for _, m in scaffold_midpoints])
    axD.set_xticklabels([short_sc(s) for s, _ in scaffold_midpoints],
                         rotation=60, ha="right", fontsize=7)
    axD.set_xlabel("Scaffold (ordered by size)", fontsize=10)

    plt.tight_layout()

    out_base = re.sub(r"\.(png|svg|pdf)$", "", args.out, flags=re.IGNORECASE)
    Path(out_base).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_base}.png", dpi=180, bbox_inches="tight")
    fig.savefig(f"{out_base}.svg",            bbox_inches="tight")

    base_cols = ["chrom", "start", "end", "n_snps", "cov_BT",
                 "across_rep_var_M", "across_rep_var_F"]
    if "ratio_M_over_F" in top_hits.columns:
        base_cols.append("ratio_M_over_F")
    if "n_genes" in top_hits.columns:
        base_cols.append("n_genes")
    base_cols.append("levene_score")

if __name__ == "__main__":
    main()
