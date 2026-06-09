#!/usr/bin/env python3

from __future__ import annotations
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(".")
QSTAR_TSV = ROOT / "final_plots/wild/section2_wild_qstar_correlation_gw_top50_perrep.tsv"
BAYPASS_C2 = ROOT / "baypass_wild/wild_contrast_summary_contrast.out"
BAYPASS_POS = ROOT / "baypass_wild/wild_snp_positions.csv"
AF_MATRIX = ROOT / "final_plots/wild/af_matrix_22pools.csv"
SNP_MULLER = ROOT / "final_plots/wild/sfla_v2_snp_muller.tsv"
OUT_BASE = ROOT / "final_plots/wild/section3_fig4_three_panels_v2"

WILD_B = ["AVB", "PSB", "RMB"]
WILD_T = ["AVT", "PST", "RMT"]
LAB_B_G10 = ["B1G10", "B2G10", "B3G10", "B4G10"]
LAB_T_G10 = ["T1G10", "T2G10", "T3G10", "T4G10"]
CHROM_439 = "chr_ScDA7r2_439_HRSCAF_779"
SIG_S, SIG_E = 2_800_000, 3_000_000

BLOCK_BP = 1_000_000
N_BOOT = 2000

MULLER_LABEL = {"A": "A (X)", "B": "B (2L)", "C": "C (2R)",
                "D": "D (3L)", "E": "E (3R)"}
MULLER_COLOR = {"A": "#d62728", "B": "#3d6cb5", "C": "#33a02c",
                "D": "#9467bd", "E": "#e6ab02"}
COL_BG = "#bbbbbb"
COL_SIG = "#762a83"


def strip_svg(svg_path: Path):
    txt = svg_path.read_text()
    txt = re.sub(r"<clipPath\b[^>]*>.*?</clipPath>", "", txt, flags=re.DOTALL)
    txt = re.sub(r"\s+clip-path\s*=\s*['\"]url\(#[^)]+\)['\"]", "", txt)
    txt = txt.replace('<svg ',
        '<svg xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" ', 1)
    txt = re.sub(r'(<g\s+id="([^"]+)")',
                 r'\1 inkscape:label="\2" inkscape:groupmode="layer"', txt)
    svg_path.write_text(txt)


def block_bootstrap_ci(df, x_col, y_col, n_boot=N_BOOT, block_bp=BLOCK_BP,
                         rng=None):
    if rng is None:
        rng = np.random.default_rng(42)
    df = df.dropna(subset=[x_col, y_col]).copy()
    if len(df) < 5:
        return np.nan, np.nan, np.nan
    df["block"] = (df["chrom"].astype(str) + "_"
                    + (df["pos"] // block_bp).astype(str))
    blocks = df["block"].unique()
    n_blocks = len(blocks)
    obs = float(spearmanr(df[x_col], df[y_col]).correlation)
    if n_blocks < 5:
        rs = np.empty(n_boot)
        n = len(df)
        for i in range(n_boot):
            idx = rng.integers(0, n, n)
            rs[i] = spearmanr(df[x_col].values[idx],
                                df[y_col].values[idx]).correlation
        return obs, float(np.percentile(rs, 2.5)), float(np.percentile(rs, 97.5))
    rs = []
    for _ in range(n_boot):
        sampled_blocks = rng.choice(blocks, size=n_blocks, replace=True)
        idx = df["block"].isin(sampled_blocks)
        sub = df[idx]
        if len(sub) >= 5:
            rs.append(spearmanr(sub[x_col], sub[y_col]).correlation)
    rs = np.array(rs)
    return obs, float(np.percentile(rs, 2.5)), float(np.percentile(rs, 97.5))


def main():
    plt.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.7,
    })

    df = pd.read_csv(QSTAR_TSV, sep="\t",
                      usecols=["chrom", "pos", "lrt", "sB", "sT",
                                "sB_SE", "sT_SE"])
    df["lab_signed"] = np.sign(df["sB"]) * np.sqrt(df["lrt"].clip(lower=0))
    df["zB"] = df["sB"].abs() / df["sB_SE"].clip(lower=1e-9)
    df["zT"] = df["sT"].abs() / df["sT_SE"].clip(lower=1e-9)
    INPUT_TOP_PCT = 0.50
    q95 = max(0.0, (0.95 - INPUT_TOP_PCT) / (1 - INPUT_TOP_PCT))
    thr5 = df["lrt"].quantile(q95)
    df["bg_pass"] = (df["lrt"] >= thr5) & (df["zB"] >= 2) & (df["zT"] >= 2)

    bp = pd.read_csv(BAYPASS_C2, sep=r"\s+",
                      names=["contrast", "mrk", "C2", "log10p"], skiprows=1)
    pos = pd.read_csv(BAYPASS_POS)
    bp = bp.merge(pos, on="mrk", how="inner")
    af = pd.read_csv(AF_MATRIX,
                       usecols=["chrom_pos"] + WILD_B + WILD_T
                                + LAB_B_G10 + LAB_T_G10)
    af[["chrom", "pos"]] = af["chrom_pos"].str.split(":", n=1, expand=True)
    af["pos"] = af["pos"].astype(int)
    af["wild_dp"] = af[WILD_B].mean(axis=1) - af[WILD_T].mean(axis=1)
    af["lab_dp"] = af[LAB_B_G10].mean(axis=1) - af[LAB_T_G10].mean(axis=1)
    bp = bp.merge(af[["chrom", "pos", "wild_dp", "lab_dp"]],
                   on=["chrom", "pos"], how="left")
    bp["wild_signed"] = np.sign(bp["wild_dp"]) * np.sqrt(bp["C2"].clip(lower=0))

    sm = pd.read_csv(SNP_MULLER, sep="\t")[["chrom", "pos", "muller"]]

    df = df.merge(bp[["chrom", "pos", "wild_signed", "wild_dp", "lab_dp"]],
                   on=["chrom", "pos"], how="inner")
    df = df.merge(sm, on=["chrom", "pos"], how="left")
    bg = df[df["bg_pass"]].copy()

    rng = np.random.default_rng(42)
    rows = []
    scaff_counts = bg.groupby("chrom").size()
    big_scaffs = scaff_counts[scaff_counts >= 30].sort_values(ascending=False)
    scaff_muller = {}
    for s in big_scaffs.index:
        sub = bg[bg["chrom"] == s]
        if sub["muller"].notna().any():
            scaff_muller[s] = sub["muller"].mode().iloc[0]

    for s in big_scaffs.index:
        sub = bg[bg["chrom"] == s]
        rho, lo, hi = block_bootstrap_ci(sub, "lab_signed", "wild_signed",
                                            rng=rng)
        m = scaff_muller.get(s, "")
        m_lab = MULLER_LABEL.get(m, m or "?")
        scaff_short = "chr_" + s.split("_")[2]
        if s == "chr_ScDA7r2_597_HRSCAF_953":
            m_lab = "C + D fusion"
        rows.append({"label": scaff_short, "muller_lab": m_lab,
                      "n": len(sub), "rho": rho, "ci_lo": lo, "ci_hi": hi,
                      "kind": "scaffold", "muller": m})

    for m in ["A", "B", "C", "D", "E"]:
        sub = bg[bg["muller"] == m]
        if len(sub) < 30:
            continue
        rho, lo, hi = block_bootstrap_ci(sub, "lab_signed", "wild_signed",
                                            rng=rng)
        rows.append({"label": f"All Muller {m}", "muller_lab": MULLER_LABEL[m],
                      "n": len(sub), "rho": rho, "ci_lo": lo, "ci_hi": hi,
                      "kind": "muller", "muller": m})

    A_df = pd.DataFrame(rows)
    A_df["muller_order"] = A_df["muller"].map(
        {m: i for i, m in enumerate(["A", "B", "C", "D", "E"])}).fillna(99)
    A_df["kind_order"] = (A_df["kind"] == "muller").astype(int)
    A_df = A_df.sort_values(["muller_order", "kind_order", "n"],
                              ascending=[True, True, False]).reset_index(drop=True)

    A_subset = bg[bg["muller"] == "A"].copy()
    sig_block_mask = ((A_subset["chrom"] == CHROM_439)
                       & (A_subset["pos"].between(SIG_S, SIG_E)))
    A_subset["is_signal"] = sig_block_mask

    A_subset["lab_abs"] = A_subset["lab_signed"].abs()
    edges = np.quantile(A_subset["lab_abs"], [0, 0.25, 0.50, 0.75, 1.00])
    edges[0] = -np.inf; edges[-1] = np.inf
    bin_labels = ["lo (Q1)", "Q2", "Q3", "hi (Q4)"]
    A_subset["lab_bin"] = pd.cut(A_subset["lab_abs"], bins=edges,
                                   labels=bin_labels, include_lowest=True)

    chr_439_signal_dp = bg[(bg["chrom"] == CHROM_439)
                            & (bg["pos"].between(SIG_S, SIG_E))].copy()
    bg_other = bg[~bg.index.isin(chr_439_signal_dp.index)].copy()
    rho_all_dp, lo_all_dp, hi_all_dp = block_bootstrap_ci(
        bg_other, "lab_dp", "wild_dp", rng=rng)
    rho_sig_dp, lo_sig_dp, hi_sig_dp = block_bootstrap_ci(
        chr_439_signal_dp, "lab_dp", "wild_dp", rng=rng)

    fig = plt.figure(figsize=(15.0, 5.5))
    gs = fig.add_gridspec(1, 3, wspace=0.45, width_ratios=[1.4, 1.0, 1.0])

    axA = fig.add_subplot(gs[0])
    axA.axvline(0, color="black", lw=0.5, alpha=0.6)
    y = np.arange(len(A_df))
    for i, (_, r) in enumerate(A_df.iterrows()):
        if not np.isfinite(r["rho"]):
            continue
        c = MULLER_COLOR.get(r["muller"], "#888")

        if r["kind"] == "scaffold":
            face = "white"; edge = c; size = 8; lw = 1.5
        else:
            face = c; edge = "black"; size = 10; lw = 0.6
        axA.errorbar(r["rho"], i,
                      xerr=[[r["rho"] - r["ci_lo"]],
                             [r["ci_hi"] - r["rho"]]],
                      fmt="o", markerfacecolor=face, markeredgecolor=edge,
                      markersize=size, markeredgewidth=lw, ecolor=c,
                      elinewidth=1.2, capsize=3, zorder=5)
    ylabels = [f"{r['label']}  [{r['muller_lab']}]  n={r['n']:,}"
               for _, r in A_df.iterrows()]
    axA.set_yticks(y)
    axA.set_yticklabels(ylabels, fontsize=7.5)
    axA.invert_yaxis()
    axA.set_xlabel(r"$\rho$  (signed √LRT  vs  signed √C₂)  [block bootstrap]",
                    fontsize=8.5)
    axA.set_title("A. Per-scaffold and per-Muller correspondence",
                   fontsize=10, pad=4)
    axA.tick_params(labelsize=7.5)


    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="white",
                markeredgecolor="#888", markeredgewidth=1.5, markersize=7,
                label="single scaffold"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#888",
                markeredgecolor="black", markersize=8,
                label="all scaffolds in Muller"),
    ]
    axA.legend(handles=handles, fontsize=7, frameon=False, loc="lower right")
    for sp in ("top", "right"): axA.spines[sp].set_visible(False)


    axB = fig.add_subplot(gs[1])
    axB.axhline(0, color="black", lw=0.5, alpha=0.5)


    rng_jit = np.random.default_rng(7)
    bin_centers = np.arange(len(bin_labels))
    for i, b in enumerate(bin_labels):
        sub_bg = A_subset[(A_subset["lab_bin"] == b) & (~A_subset["is_signal"])]
        sub_sig = A_subset[(A_subset["lab_bin"] == b) & (A_subset["is_signal"])]
        n_bg = len(sub_bg)
        n_sig = len(sub_sig)
        if n_bg > 0:
            jit = rng_jit.uniform(-0.30, 0.30, size=n_bg)
            axB.scatter(np.full(n_bg, i) + jit, sub_bg["wild_signed"],
                         s=14, color=COL_BG, alpha=0.55, edgecolors="none",
                         zorder=2)
        if n_sig > 0:
            jit = rng_jit.uniform(-0.20, 0.20, size=n_sig)
            axB.scatter(np.full(n_sig, i) + jit, sub_sig["wild_signed"],
                         s=70, facecolor=COL_SIG, edgecolor="black",
                         linewidth=0.8, zorder=10)

        ymax = A_subset["wild_signed"].max()
        axB.text(i, ymax * 1.05, f"X-A: {n_bg}\nsig: {n_sig}",
                  ha="center", va="bottom", fontsize=7,
                  family="monospace")

    axB.set_xticks(bin_centers)
    axB.set_xticklabels(bin_labels, fontsize=8)
    axB.set_xlabel("Lab |signed √LRT| quartile  (within Muller A)",
                    fontsize=9)
    axB.set_ylabel(r"Wild  sign($\Delta p$) × $\sqrt{C_2}$", fontsize=9)
    axB.set_title("B. Muller A by lab significance — signal block on top",
                   fontsize=10, pad=4)
    axB.tick_params(labelsize=8)

    handles_B = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_BG,
                markeredgecolor="none", markersize=7,
                label=f"Muller A (excl signal)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_SIG,
                markeredgecolor="black", markersize=9,
                label=f"chr_439 signal block"),
    ]
    axB.legend(handles=handles_B, fontsize=7, frameon=False, loc="lower right")
    for sp in ("top", "right"): axB.spines[sp].set_visible(False)


    axC = fig.add_subplot(gs[2])
    axC.axhline(0, color="black", lw=0.5, alpha=0.5)
    axC.axvline(0, color="black", lw=0.5, alpha=0.5)
    axC.plot([-0.5, 0.5], [-0.5, 0.5], color="black", lw=0.7, ls="--",
              alpha=0.4)
    axC.scatter(bg_other["lab_dp"], bg_other["wild_dp"],
                 s=4, color="#cccccc", alpha=0.30, edgecolors="none",
                 rasterized=True, zorder=2,
                 label=f"genome-wide bg (n={len(bg_other):,})")
    axC.scatter(chr_439_signal_dp["lab_dp"], chr_439_signal_dp["wild_dp"],
                 marker="o", s=85, facecolor=COL_SIG, edgecolor="black",
                 linewidth=0.9, zorder=10,
                 label=f"chr_439 signal block (n={len(chr_439_signal_dp)})")
    if len(chr_439_signal_dp) >= 5:
        b, a = np.polyfit(chr_439_signal_dp["lab_dp"],
                            chr_439_signal_dp["wild_dp"], 1)
        xs = np.linspace(chr_439_signal_dp["lab_dp"].min(),
                          chr_439_signal_dp["lab_dp"].max(), 50)
        axC.plot(xs, a + b * xs, color=COL_SIG, lw=2.0, zorder=8)

    txt = (f"genome-wide bg: ρ={rho_all_dp:+.3f}  "
           f"CI=[{lo_all_dp:+.2f},{hi_all_dp:+.2f}]\n"
           f"signal block:   ρ={rho_sig_dp:+.3f}  "
           f"CI=[{lo_sig_dp:+.2f},{hi_sig_dp:+.2f}]")
    axC.text(0.04, 0.97, txt, transform=axC.transAxes,
              ha="left", va="top", fontsize=7.5, family="monospace",
              bbox=dict(facecolor="white", alpha=0.95, edgecolor="#aaa",
                         linewidth=0.4, boxstyle="round,pad=0.3"))
    axC.set_xlim(-0.45, 0.45); axC.set_ylim(-0.45, 0.45)
    axC.set_xlabel(r"Lab $\Delta p$  (B G10 − T G10)", fontsize=9)
    axC.set_ylabel(r"Wild $\Delta p$  (B − T)", fontsize=9)
    axC.set_title("C. Lab Δp vs wild Δp: signal block highlighted",
                   fontsize=10, pad=4)
    axC.tick_params(labelsize=8)
    axC.legend(fontsize=7, frameon=False, loc="lower right")
    for sp in ("top", "right"): axC.spines[sp].set_visible(False)

    fig.suptitle(
        "Lab vs wild antagonism: chr_439 signal block as the focal locus  (block-bootstrap CIs)",
        fontsize=11, y=0.995)
    fig.tight_layout()
    fig.savefig(f"{OUT_BASE}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.svg", bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.pdf", bbox_inches="tight")
    strip_svg(Path(f"{OUT_BASE}.svg"))


if __name__ == "__main__":
    main()
