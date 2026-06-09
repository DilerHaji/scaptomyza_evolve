#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(".")
LRT_ITER_TSV = ROOT / "final_plots/wild/section2_fig3e_qstar_iteration.tsv"  # LRT-based
LYNCH_PANEL_TSV = ROOT / "final_plots/wild/section2_fig3e_lynch_panel.tsv"   # Lynch G10
OUT_DIR = ROOT / "final_plots/wild"

C = 0.5
H = 0.5

W_AA_B_P, W_Aa_B_P, W_aa_B_P = 1.000, 0.913, 0.612
W_AA_T_P, W_Aa_T_P, W_aa_T_P = 0.611, 0.806, 0.902
sB_P = (W_AA_B_P - W_aa_B_P) / W_AA_B_P
sT_P = (W_aa_T_P - W_AA_T_P) / W_aa_T_P
qstar_P_h05 = sB_P / (sB_P + sT_P)


def landscape_curves(W_AA_B, W_Aa_B, W_aa_B, W_AA_T, W_Aa_T, W_aa_T,
                      c=C, q=None):
    if q is None:
        q = np.linspace(0, 1, 400)
    Wb = q*q*W_AA_B + 2*q*(1-q)*W_Aa_B + (1-q)*(1-q)*W_aa_B
    Wt = q*q*W_AA_T + 2*q*(1-q)*W_Aa_T + (1-q)*(1-q)*W_aa_T
    Wm = c * Wb + (1 - c) * Wt
    return q, Wb, Wt, Wm


def plot_minimal_landscape(ax, Wgt, qstar, label):
    q, Wb, Wt, Wm = landscape_curves(*Wgt)
    norm = Wm.max()
    Wb_n, Wt_n, Wm_n = Wb / norm, Wt / norm, Wm / norm

    ax.plot(q, Wb_n, color="#a04848", linewidth=1.1)
    ax.plot(q, Wt_n, color="#3d6cb5", linewidth=1.1)
    ax.plot(q, Wm_n, color="#1a1a1a", linewidth=1.7)
    ax.axvline(qstar, color="#1a5e1a", linestyle="--", linewidth=1.0, zorder=4)

    q_grid = q
    idx = int(np.argmin(np.abs(q_grid - qstar)))
    y_at_qstar = Wm_n[idx]
    ax.scatter([qstar], [y_at_qstar], s=45, color="#1a5e1a",
                marker="v", edgecolor="black", linewidth=0.5, zorder=6)

    for q_pt, lab in zip([0.0, 1.0], ["aa", "AA"]):
        ax.text(q_pt, -0.06, lab, transform=ax.get_xaxis_transform(),
                 ha="center", va="top", fontsize=6.5, color="#444")

    ax.set_xlim(0, 1)
    y_min = min(Wb_n.min(), Wt_n.min(), Wm_n.min()) * 0.97
    ax.set_ylim(y_min, 1.04)
    ax.set_title(label, fontsize=8, pad=2)
    ax.tick_params(labelsize=6, length=2, pad=1)
    ax.set_xticks([])
    ax.set_yticks(np.arange(0.75, 1.05, 0.10))


def main():
    plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42,
                          "font.family": "sans-serif",
                          "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
                          "axes.linewidth": 0.5})

    L = pd.read_csv(LYNCH_PANEL_TSV, sep="\t")
    g10 = L[L["g_max"] == 10].iloc[0]
    sB_G_lynch = float(g10["med_sB"])
    sT_G_lynch = float(g10["med_sT"])
    qstar_G_lynch = float(g10["qstar_h05"])

    Wgt_G = (1 + sB_G_lynch, 1 + sB_G_lynch * H, 1.0,
             1.0, 1 + abs(sT_G_lynch) * H, 1 + abs(sT_G_lynch))
    Wgt_P = (W_AA_B_P, W_Aa_B_P, W_aa_B_P,
             W_AA_T_P, W_Aa_T_P, W_aa_T_P)

    R = pd.read_csv(LRT_ITER_TSV, sep="\t")
    R_all = R[R["g_max"] >= 2].reset_index(drop=True)

    fig = plt.figure(figsize=(4.4, 4.6))
    gs = fig.add_gridspec(2, 2,
                            height_ratios=[1.0, 1.4],
                            width_ratios=[1, 1],
                            hspace=0.50, wspace=0.30)
    ax_p = fig.add_subplot(gs[0, 0])
    ax_g = fig.add_subplot(gs[0, 1])
    ax_t = fig.add_subplot(gs[1, :])

    plot_minimal_landscape(ax_p, Wgt_P, qstar_P_h05,
                            f"Phenotypic   q*={qstar_P_h05:.2f}")
    plot_minimal_landscape(ax_g, Wgt_G, qstar_G_lynch,
                            f"Genomic (Lynch)   q*={qstar_G_lynch:.2f}")
    ax_p.set_ylabel("normalized\nmean fitness", fontsize=6.5, labelpad=1)

    ax_t.axhline(qstar_P_h05, color="#a04848", linewidth=1.0, alpha=0.85,
                  zorder=2, label=f"q*_pheno (h=0.5) = {qstar_P_h05:.3f}")
    ax_t.axhline(0.5, color="#888", linestyle=":", linewidth=0.4,
                  alpha=0.6, zorder=1)

    R_all = R_all.sort_values("g_max").reset_index(drop=True)
    Y_ND = 0.95
    R_valid = R_all[R_all["qstar"].notna()].reset_index(drop=True)
    R_undef = R_all[R_all["qstar"].isna()].reset_index(drop=True)

    y_pos = [r["qstar"] if pd.notna(r["qstar"]) else Y_ND
              for _, r in R_all.iterrows()]
    is_valid = [pd.notna(r["qstar"]) for _, r in R_all.iterrows()]
    for i in range(len(R_all) - 1):
        x0, x1 = R_all["g_max"].iloc[i], R_all["g_max"].iloc[i + 1]
        y0, y1 = y_pos[i], y_pos[i + 1]
        both = is_valid[i] and is_valid[i + 1]
        ax_t.plot([x0, x1], [y0, y1],
                   linestyle="-" if both else ":",
                   linewidth=1.3 if both else 0.7,
                   color="#1a5e1a",
                   alpha=1.0 if both else 0.55, zorder=3)

    ax_t.errorbar(R_valid["g_max"], R_valid["qstar"],
                   yerr=[R_valid["qstar"] - R_valid["q_lo"],
                          R_valid["q_hi"] - R_valid["qstar"]],
                   fmt="o", color="#1a5e1a", markersize=5,
                   capsize=2.5, linewidth=0, elinewidth=1.2,
                   label="q*_geno (top-LRT, h=0.5)", zorder=5)
    for _, r in R_valid.iterrows():
        ax_t.text(r["g_max"], r["qstar"] + 0.05, f"{r['qstar']:.2f}",
                   ha="center", fontsize=6, color="#1a5e1a")

    for _, r in R_undef.iterrows():
        ax_t.scatter([r["g_max"]], [Y_ND], s=40, facecolors="none",
                      edgecolors="#1a5e1a", linewidth=1.0, zorder=5,
                      marker="o")
        ax_t.text(r["g_max"], Y_ND + 0.03, "n.d.", ha="center",
                   va="bottom", fontsize=6, color="#1a5e1a", style="italic")

    ax_t.set_xlabel("included generations  (G1 to G_max)", fontsize=7)
    ax_t.set_ylabel("q*", fontsize=7)
    ax_t.set_xticks(R_all["g_max"])
    ax_t.set_xticklabels([f"G{g}" for g in R_all["g_max"]], fontsize=6)
    ax_t.set_ylim(0.4, 1.05)
    ax_t.tick_params(labelsize=6, length=2, pad=1)
    ax_t.legend(loc="lower left", fontsize=6, frameon=False)

    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.10, top=0.94)
    save(fig, "section2_fig3e_combined_panel")


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
