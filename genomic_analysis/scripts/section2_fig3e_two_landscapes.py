#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(".")
S_TSV = ROOT / "final_plots/wild/section2_fig3e_per_snp_s_estimates.tsv"
OUT_DIR = ROOT / "final_plots/wild"
C_HABITAT = 0.5

W_AA_B_P, W_Aa_B_P, W_aa_B_P = 1.000, 0.913, 0.612
W_AA_T_P, W_Aa_T_P, W_aa_T_P = 0.611, 0.806, 0.902
sB_P = (W_AA_B_P - W_aa_B_P) / W_AA_B_P
hB_P = (W_Aa_B_P - W_aa_B_P) / (W_AA_B_P - W_aa_B_P)
sT_P = (W_aa_T_P - W_AA_T_P) / W_aa_T_P
hT_P = (W_Aa_T_P - W_AA_T_P) / (W_aa_T_P - W_AA_T_P)
qstar_P = (C_HABITAT * sB_P * hB_P) / (C_HABITAT * sB_P * hB_P
                                          + (1 - C_HABITAT) * sT_P * hT_P)

def landscape_curves(W_AA_B, W_Aa_B, W_aa_B, W_AA_T, W_Aa_T, W_aa_T, q=None):
    if q is None:
        q = np.linspace(0, 1, 200)
    Wb = q*q*W_AA_B + 2*q*(1-q)*W_Aa_B + (1-q)*(1-q)*W_aa_B
    Wt = q*q*W_AA_T + 2*q*(1-q)*W_Aa_T + (1-q)*(1-q)*W_aa_T
    Wm = C_HABITAT * Wb + (1 - C_HABITAT) * Wt
    return q, Wb, Wt, Wm


def qstar_dempster(sB, hB, sT, hT, c=C_HABITAT):
    return (c * sB * hB) / (c * sB * hB + (1 - c) * sT * hT)


def plot_landscape(ax, W_AA_B, W_Aa_B, W_aa_B,
                    W_AA_T, W_Aa_T, W_aa_T,
                    qstar, title_prefix, sB, sT, hB, hT,
                    extra_text=""):
    q, Wb, Wt, Wm = landscape_curves(W_AA_B, W_Aa_B, W_aa_B,
                                        W_AA_T, W_Aa_T, W_aa_T)

    norm = Wm.max()
    Wb_n = Wb / norm
    Wt_n = Wt / norm
    Wm_n = Wm / norm

    ax.plot(q, Wb_n, color="#a04848", linewidth=1.4, label="W on B host")
    ax.plot(q, Wt_n, color="#3d6cb5", linewidth=1.4, label="W on T host")
    ax.plot(q, Wm_n, color="#1a1a1a", linewidth=2.2, label="Marginal (c=0.5)")
    ax.axvline(qstar, color="#1a5e1a", linestyle="--", linewidth=1.2,
                zorder=4, label=f"q* = {qstar:.3f}")
    ax.scatter([qstar], [Wm_n.max()], s=120, color="#1a5e1a",
                marker="v", edgecolor="black", linewidth=0.7, zorder=5)

    for q_pt, label in zip([0.0, 0.5, 1.0], ["aa", "Aa", "AA"]):
        ax.scatter([q_pt], [(W_aa_B if label=='aa'
                              else W_Aa_B if label=='Aa'
                              else W_AA_B) / norm],
                    s=22, color="#a04848", edgecolor="white", linewidth=0.5,
                    zorder=6)
        ax.scatter([q_pt], [(W_aa_T if label=='aa'
                              else W_Aa_T if label=='Aa'
                              else W_AA_T) / norm],
                    s=22, color="#3d6cb5", edgecolor="white", linewidth=0.5,
                    zorder=6)
        ax.text(q_pt, -0.04, label, transform=ax.get_xaxis_transform(),
                 ha="center", va="top", fontsize=7, color="#444")

    ax.set_xlim(0, 1)
    ax.set_ylim(min(Wb_n.min(), Wt_n.min(), Wm_n.min()) * 0.96, 1.06)
    ax.set_xlabel("allele frequency q   (A = B-favored)", fontsize=8)
    ax.set_ylabel("normalized mean fitness", fontsize=8)
    ax.set_title(title_prefix, fontsize=10)
    ax.tick_params(labelsize=7)
    ax.legend(loc="lower center", fontsize=6.5, frameon=False, ncol=2)
    txt = (f"s_B = {sB:+.3f}   h_B = {hB:.2f}\n"
            f"s_T = {sT:+.3f}   h_T = {hT:.2f}")
    if extra_text:
        txt += "\n" + extra_text
    ax.text(0.04, 0.96, txt, transform=ax.transAxes,
             ha="left", va="top", fontsize=6.5, family="monospace",
             color="#222", bbox=dict(facecolor="white", alpha=0.7,
                                       edgecolor="#bbbbbb", linewidth=0.4,
                                       boxstyle="round,pad=0.25"))


def bootstrap_q_star(df_subset, n_boot=2000, seed=42):
    rng = np.random.default_rng(seed)
    n = len(df_subset)
    sB_vals = df_subset["sB"].values
    sT_vals = df_subset["sT"].values
    out = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        med_sB = np.median(sB_vals[idx])
        med_sT = np.median(sT_vals[idx])
        if med_sB > 0 and med_sT < 0:
            out[i] = qstar_dempster(med_sB, 0.5, abs(med_sT), 0.5)
        else:
            out[i] = np.nan
    return out


def main():
    plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42,
                          "font.family": "sans-serif",
                          "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
                          "axes.linewidth": 0.6})

    df = pd.read_csv(S_TSV, sep="\t")
    top = df[df["in_top"]]
    bg  = df[~df["in_top"]]
    sB_G = float(top["sB"].median())   # already polarized so s_B > 0 globally
    sT_G = float(top["sT"].median())

    qstar_top_boot = bootstrap_q_star(top, n_boot=2000, seed=42)
    valid = ~np.isnan(qstar_top_boot)
    qstar_top_med = float(np.nanmedian(qstar_top_boot))
    qstar_top_lo  = float(np.nanquantile(qstar_top_boot, 0.025))
    qstar_top_hi  = float(np.nanquantile(qstar_top_boot, 0.975))

    n_top = int(top.shape[0])
    rng = np.random.default_rng(0)
    null_q = []
    n_null = 2000
    for i in range(n_null):
        sub = bg.iloc[rng.integers(0, len(bg), size=n_top)]
        msB = sub["sB"].median()
        msT = sub["sT"].median()
        if msB > 0 and msT < 0:
            null_q.append(qstar_dempster(msB, 0.5, abs(msT), 0.5))
    null_q = np.array(null_q)

    p_one_tailed = float(np.mean(null_q <= qstar_top_med))

    hB_G, hT_G = 0.5, 0.5

    W_aa_B_G = 1.0
    W_Aa_B_G = 1.0 + sB_G * hB_G
    W_AA_B_G = 1.0 + sB_G

    sT_pos = abs(sT_G)
    W_AA_T_G = 1.0
    W_Aa_T_G = 1.0 + sT_pos * hT_G
    W_aa_T_G = 1.0 + sT_pos

    qstar_G = qstar_dempster(sB_G, hB_G, sT_pos, hT_G)

    fig = plt.figure(figsize=(11.5, 3.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.85], wspace=0.35)
    ax_p = fig.add_subplot(gs[0, 0])
    ax_g = fig.add_subplot(gs[0, 1])
    ax_n = fig.add_subplot(gs[0, 2])

    plot_landscape(ax_p, W_AA_B_P, W_Aa_B_P, W_aa_B_P,
                    W_AA_T_P, W_Aa_T_P, W_aa_T_P, qstar_P,
                    "Phenotypic perspective\n(F3 reciprocal-performance, Panel A)",
                    sB_P, -sT_P, hB_P, hT_P)

    extra = (f"derived from {len(top)} top-LRT SNPs\n"
              f"in chr_439 LD block (970 kb)")
    plot_landscape(ax_g, W_AA_B_G, W_Aa_B_G, W_aa_B_G,
                    W_AA_T_G, W_Aa_T_G, W_aa_T_G, qstar_G,
                    "Genomic perspective\n(chr_439 signal-region AF dynamics)",
                    sB_G, sT_G, hB_G, hT_G,
                    extra_text=extra)

    bins = np.linspace(0, 1, 41)
    ax_n.hist(null_q, bins=bins, color="#a8a8a8", alpha=0.7,
                edgecolor="none",
                label=f"Random matched samples\n(n={len(null_q)} draws)")

    ax_n.axvline(qstar_P, color="#a04848", linestyle="-", linewidth=1.4,
                  alpha=0.8, label=f"q*_pheno = {qstar_P:.3f}")

    ax_n.axvspan(qstar_top_lo, qstar_top_hi, color="#1a5e1a",
                  alpha=0.2, zorder=2)
    ax_n.axvline(qstar_top_med, color="#1a5e1a", linestyle="--", linewidth=1.4,
                  alpha=0.95,
                  label=f"q*_geno (signal) = {qstar_top_med:.3f}\n"
                         f"  95% boot CI [{qstar_top_lo:.3f}, {qstar_top_hi:.3f}]")
    ax_n.set_xlabel("q*  (Dempster equilibrium)", fontsize=8)
    ax_n.set_ylabel("count", fontsize=8)
    ax_n.set_xlim(0, 1)
    ax_n.set_title("Null vs observed q* (genomic)", fontsize=10)
    ax_n.tick_params(labelsize=7)
    ax_n.legend(loc="upper left", fontsize=6.5, frameon=False)
    ax_n.text(0.97, 0.97,
                f"empirical p\n(null ≤ observed) = {p_one_tailed:.3f}",
                transform=ax_n.transAxes, ha="right", va="top",
                fontsize=7, color="#222",
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="#bbb",
                          linewidth=0.4, boxstyle="round,pad=0.25"))

    fig.suptitle(f"Two Dempster perspectives:  q*_pheno = {qstar_P:.2f}, "
                  f"q*_geno = {qstar_top_med:.2f}  (matched-region p = {p_one_tailed:.3f})",
                  fontsize=11)
    fig.tight_layout()
    save(fig, "section2_fig3e_two_landscapes")


def save(fig, base):
    out = OUT_DIR / base
    fig.savefig(f"{out}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out}.svg", bbox_inches="tight")
    fig.savefig(f"{out}.pdf", bbox_inches="tight")
    import re
    txt = Path(f"{out}.svg").read_text()
    txt = re.sub(r"<clipPath[^>]*>.*?</clipPath>", "", txt, flags=re.DOTALL)
    txt = re.sub(r'\s*clip-path="url\([^)]+\)"', "", txt)
    Path(f"{out}.svg").write_text(txt)


if __name__ == "__main__":
    main()
