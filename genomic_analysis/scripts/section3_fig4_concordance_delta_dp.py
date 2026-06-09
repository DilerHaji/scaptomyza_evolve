#!/usr/bin/env python3

from __future__ import annotations
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(".")
QSTAR_TSV = ROOT / "final_plots/wild/section2_wild_qstar_correlation_gw_top50_perrep.tsv"
QSTAR_GLM_TSV = ROOT / "final_plots/wild/qstar_binomial_glm.tsv"
BAYPASS_BF = ROOT / "baypass_wild/wild_trt_summary_betai_reg.out"
BAYPASS_POS = ROOT / "baypass_wild/wild_snp_positions.csv"
AF_MATRIX = ROOT / "final_plots/wild/af_matrix_22pools.csv"
SNP_MULLER = ROOT / "final_plots/wild/sfla_v2_snp_muller.tsv"
OUT_BASE = ROOT / "final_plots/wild/section3_fig4_concordance_delta_dp"

WILD_B = ["AVB", "PSB", "RMB"]
WILD_T = ["AVT", "PST", "RMT"]
LAB_B_G10 = ["B1G10", "B2G10", "B3G10", "B4G10"]
LAB_T_G10 = ["T1G10", "T2G10", "T3G10", "T4G10"]
FOUNDER = ["F1G00", "F2G00", "F3G00", "F4G00"]

CHROM_439 = "chr_ScDA7r2_439_HRSCAF_779"
LD_S, LD_E = 2_640_000, 3_610_000
SIG_S, SIG_E = 2_800_000, 3_000_000

MULLER_ORDER = ["A", "B", "C", "D", "E"]
MULLER_LABEL = {"A": "A (X)", "B": "B (2L)", "C": "C (2R)",
                "D": "D (3L)", "E": "E (3R)"}
MULLER_COLOR = {"A": "#d62728", "B": "#3d6cb5", "C": "#33a02c",
                "D": "#9467bd", "E": "#e6ab02"}
COL_BLUE_LIGHT = "#5da7d6"   # signal-block-matched
COL_BLUE_DK = "#0d3a66"       # LD-block-matched
COL_PURPLE = "#762a83"

N_BOOT_MATCHED = 1000
AF_TOL = 0.05  # founder MAF tolerance for AF-matched sampling
SEED = 42


def strip_svg(svg_path: Path):
    txt = svg_path.read_text()
    txt = re.sub(r"<clipPath\b[^>]*>.*?</clipPath>", "", txt, flags=re.DOTALL)
    txt = re.sub(r"\s+clip-path\s*=\s*['\"]url\(#[^)]+\)['\"]", "", txt)
    txt = txt.replace('<svg ',
        '<svg xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" ', 1)
    txt = re.sub(r'(<g\s+id="([^"]+)")',
                 r'\1 inkscape:label="\2" inkscape:groupmode="layer"', txt)
    svg_path.write_text(txt)


def load_data():
    df = pd.read_csv(QSTAR_TSV, sep="\t",
                      usecols=["chrom", "pos", "lrt", "sB", "sT",
                                "sB_SE", "sT_SE", "qstar_pred", "qstar_pred_SE",
                                "wild_AF_pol", "antag"])
    df["zB"] = df["sB"].abs() / df["sB_SE"].clip(lower=1e-9)
    df["zT"] = df["sT"].abs() / df["sT_SE"].clip(lower=1e-9)
    INPUT_TOP_PCT = 0.50
    q95 = max(0.0, (0.95 - INPUT_TOP_PCT) / (1 - INPUT_TOP_PCT))
    thr5 = df["lrt"].quantile(q95)
    df["bg_pass"] = (df["lrt"] >= thr5) & (df["zB"] >= 2) & (df["zT"] >= 2)

    bp = pd.read_csv(BAYPASS_BF, sep=r"\s+",
                      names=["covariable", "mrk", "M_Pearson", "SD_Pearson",
                              "M_Spearman", "SD_Spearman", "BF_dB",
                              "Beta_is", "SD_Beta_is", "eBPis"], skiprows=1)
    pos = pd.read_csv(BAYPASS_POS)
    bp = bp.merge(pos, on="mrk", how="inner")

    af = pd.read_csv(AF_MATRIX,
                       usecols=["chrom_pos"] + WILD_B + WILD_T
                                + LAB_B_G10 + LAB_T_G10 + FOUNDER)
    af[["chrom", "pos"]] = af["chrom_pos"].str.split(":", n=1, expand=True)
    af["pos"] = af["pos"].astype(int)
    af["wild_dp"] = af[WILD_B].mean(axis=1) - af[WILD_T].mean(axis=1)
    af["lab_dp"] = af[LAB_B_G10].mean(axis=1) - af[LAB_T_G10].mean(axis=1)
    af["wild_AF_mean"] = af[WILD_B + WILD_T].mean(axis=1)
    af["founder_AF"] = af[FOUNDER].mean(axis=1)
    af["founder_MAF"] = np.minimum(af["founder_AF"], 1 - af["founder_AF"])

    Bp = af[WILD_B].to_numpy(dtype=float)
    Tp = af[WILD_T].to_numpy(dtype=float)
    af["wild_within_host_SD"] = 0.5 * (np.nanstd(Bp, axis=1, ddof=1)
                                         + np.nanstd(Tp, axis=1, ddof=1))

    sm = pd.read_csv(SNP_MULLER, sep="\t")[["chrom", "pos", "muller"]]

    bg = (df[df["bg_pass"]]
            .merge(bp[["chrom", "pos", "BF_dB", "Beta_is"]],
                    on=["chrom", "pos"], how="inner")
            .merge(af[["chrom", "pos", "founder_MAF",
                         "wild_dp", "lab_dp", "wild_within_host_SD",
                         "wild_AF_mean"]],
                    on=["chrom", "pos"], how="inner")
            .merge(sm, on=["chrom", "pos"], how="left"))

    bg["wild_AF_pol"] = np.where(bg["sB"] >= 0,
                                  bg["wild_AF_mean"],
                                  1 - bg["wild_AF_mean"])

    bg["concord_sign"] = (np.sign(bg["lab_dp"]) == np.sign(bg["wild_dp"])).astype(int)
    return bg


def af_matched_bootstrap(target_mafs, candidate_pool, target_n, n_boot,
                           rng, af_tol=AF_TOL):
    cand_mafs = candidate_pool["founder_MAF"].values
    cand_signs = candidate_pool["concord_sign"].values
    out = []
    for _ in range(n_boot):
        chosen_signs = []
        for tm in target_mafs:
            mask = np.abs(cand_mafs - tm) <= af_tol
            idx_pool = np.where(mask)[0]
            if len(idx_pool) == 0:
                # widen tolerance
                idx_pool = np.where(np.abs(cand_mafs - tm) <= af_tol * 3)[0]
            if len(idx_pool) == 0:
                continue
            idx = rng.choice(idx_pool, 1)[0]
            chosen_signs.append(cand_signs[idx])
        if len(chosen_signs) >= target_n // 2:  # require at least half matched
            out.append(np.mean(chosen_signs))
    return np.asarray(out)


def size_matched_bootstrap(candidate_pool, target_n, n_boot, rng):
    signs = candidate_pool["concord_sign"].values
    n = len(signs)
    if n < target_n:
        return np.array([])
    out = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.choice(n, size=target_n, replace=True)
        out[i] = signs[idx].mean()
    return out


def main():
    plt.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.7,
    })

    bg = load_data()

    A_data = bg[bg["muller"] == "A"]
    sig_block = A_data[(A_data["chrom"] == CHROM_439)
                        & (A_data["pos"].between(SIG_S, SIG_E))]
    LD_block = A_data[(A_data["chrom"] == CHROM_439)
                       & (A_data["pos"].between(LD_S, LD_E))]

    sig_target_mafs = sig_block["founder_MAF"].values
    LD_target_mafs = LD_block["founder_MAF"].values
    target_n_sig = len(sig_block)
    target_n_LD = len(LD_block)

    rng = np.random.default_rng(SEED)
    rows = []
    for m in MULLER_ORDER:
        sub = bg[bg["muller"] == m]
        if len(sub) < 30:
            continue

        v = sub["concord_sign"].dropna().to_numpy()
        p = float(np.mean(v))
        se = float(np.sqrt(p * (1 - p) / len(v)))

        pool = sub[~((sub["chrom"] == CHROM_439)
                      & (sub["pos"].between(LD_S, LD_E)))]

        sig_boots_af = af_matched_bootstrap(
            sig_target_mafs, pool, target_n_sig, N_BOOT_MATCHED, rng)
        LD_boots_af = af_matched_bootstrap(
            LD_target_mafs, pool, target_n_LD, N_BOOT_MATCHED, rng)

        sig_boots_size = size_matched_bootstrap(
            pool, target_n_sig, N_BOOT_MATCHED, rng)
        LD_boots_size = size_matched_bootstrap(
            pool, target_n_LD, N_BOOT_MATCHED, rng)

        rows.append({
            "muller": m, "n": len(v),
            "baseline_p": p, "baseline_se": se,
            "sig_match_med": float(np.median(sig_boots_af)),
            "sig_match_lo": float(np.percentile(sig_boots_af, 2.5)),
            "sig_match_hi": float(np.percentile(sig_boots_af, 97.5)),
            "LD_match_med": float(np.median(LD_boots_af)),
            "LD_match_lo": float(np.percentile(LD_boots_af, 2.5)),
            "LD_match_hi": float(np.percentile(LD_boots_af, 97.5)),
            "sig_size_cloud": sig_boots_size,
            "LD_size_cloud": LD_boots_size,
        })

    sig_v = sig_block["concord_sign"].dropna().to_numpy()
    LD_v = LD_block["concord_sign"].dropna().to_numpy()
    sig_p = float(np.mean(sig_v))
    sig_se = float(np.sqrt(sig_p * (1 - sig_p) / len(sig_v)))
    LD_p = float(np.mean(LD_v))
    LD_se = float(np.sqrt(LD_p * (1 - LD_p) / len(LD_v)))

    bg_other = bg[~((bg["chrom"] == CHROM_439)
                     & (bg["pos"].between(SIG_S, SIG_E)))]
    sig_dp = bg[(bg["chrom"] == CHROM_439)
                 & (bg["pos"].between(SIG_S, SIG_E))]
    from scipy.stats import spearmanr
    bg_other_clean = bg_other[["lab_dp", "wild_dp"]].dropna()
    rho_other, p_other = spearmanr(bg_other_clean["lab_dp"],
                                       bg_other_clean["wild_dp"])
    rho_sig, p_sig = spearmanr(sig_dp["lab_dp"], sig_dp["wild_dp"])

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.5, 5.5),
                                            gridspec_kw={"width_ratios": [1.25, 1.0],
                                                         "wspace": 0.30})


    axA.axvline(0.5, color="black", lw=0.7, ls="--", alpha=0.6,
                 label="50% (no preference)")
    rng_jit = np.random.default_rng(11)
    for i, r in enumerate(rows):
        c = MULLER_COLOR[r["muller"]]

        sig_cloud = r["sig_size_cloud"]
        LD_cloud = r["LD_size_cloud"]

        sig_show = rng_jit.choice(sig_cloud, size=min(400, len(sig_cloud)),
                                     replace=False) if len(sig_cloud) else np.array([])
        LD_show = rng_jit.choice(LD_cloud, size=min(400, len(LD_cloud)),
                                    replace=False) if len(LD_cloud) else np.array([])

        sig_xjit = 0.5 / target_n_sig
        LD_xjit = 0.5 / target_n_LD
        if len(sig_show):
            jit_y = rng_jit.uniform(0.08, 0.30, size=len(sig_show))
            jit_x = rng_jit.uniform(-sig_xjit, sig_xjit, size=len(sig_show))
            axA.scatter(sig_show + jit_x,
                         np.full(len(sig_show), i) + jit_y,
                         s=3, color=COL_BLUE_LIGHT, alpha=0.20,
                         edgecolors="none", rasterized=True, zorder=2)
        if len(LD_show):
            jit_y = rng_jit.uniform(-0.30, -0.08, size=len(LD_show))
            jit_x = rng_jit.uniform(-LD_xjit, LD_xjit, size=len(LD_show))
            axA.scatter(LD_show + jit_x,
                         np.full(len(LD_show), i) + jit_y,
                         s=3, color=COL_BLUE_DK, alpha=0.20,
                         edgecolors="none", rasterized=True, zorder=2)

        axA.errorbar(r["baseline_p"], i, xerr=r["baseline_se"], fmt="o",
                      color="black", markeredgecolor="black",
                      markersize=10, ecolor="black", elinewidth=1.4,
                      capsize=4, markeredgewidth=0.7, zorder=10)

        m, lo, hi = r["sig_match_med"], r["sig_match_lo"], r["sig_match_hi"]
        axA.errorbar(m, i + 0.18, xerr=[[m - lo], [hi - m]], fmt="o",
                      color=COL_BLUE_LIGHT, markeredgecolor="black",
                      markersize=9, ecolor=COL_BLUE_LIGHT, elinewidth=1.2,
                      capsize=3.5, markeredgewidth=0.6, zorder=9)

        m, lo, hi = r["LD_match_med"], r["LD_match_lo"], r["LD_match_hi"]
        axA.errorbar(m, i - 0.18, xerr=[[m - lo], [hi - m]], fmt="o",
                      color=COL_BLUE_DK, markeredgecolor="black",
                      markersize=9, ecolor=COL_BLUE_DK, elinewidth=1.2,
                      capsize=3.5, markeredgewidth=0.6, zorder=9)

        if r["muller"] == "A":
            axA.errorbar(sig_p, i + 0.32, xerr=sig_se, fmt="*",
                          color=COL_BLUE_LIGHT, markeredgecolor="black",
                          markersize=18, ecolor=COL_BLUE_LIGHT, elinewidth=1.2,
                          capsize=3.5, markeredgewidth=1.0, zorder=11,
                          label=f"actual signal block (n={len(sig_v)})")
            axA.errorbar(LD_p, i - 0.32, xerr=LD_se, fmt="*",
                          color=COL_BLUE_DK, markeredgecolor="black",
                          markersize=18, ecolor=COL_BLUE_DK, elinewidth=1.2,
                          capsize=3.5, markeredgewidth=1.0, zorder=11,
                          label=f"actual LD block (n={len(LD_v)})")
    axA.set_yticks(range(len(rows)))
    axA.set_yticklabels(
        [f"{MULLER_LABEL[r['muller']]}\nn={r['n']:,}" for r in rows],
        fontsize=8.5)
    axA.invert_yaxis()
    axA.set_xlabel(r"Fraction with sign(lab $\Delta p$) = sign(wild $\Delta p$)",
                    fontsize=9.5)
    axA.set_title("A. Concordant-fraction per Muller — "
                   "AF/n-matched controls",
                   fontsize=10, pad=4)
    axA.tick_params(labelsize=8)

    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="black",
                markeredgecolor="black", markersize=8,
                label="full-Muller baseline (binomial SE)"),
        Line2D([0], [0], marker="o", color="w",
                markerfacecolor=COL_BLUE_LIGHT, markeredgecolor="black",
                markersize=8, label="signal-block-matched (bootstrap CI)"),
        Line2D([0], [0], marker="o", color="w",
                markerfacecolor=COL_BLUE_DK, markeredgecolor="black",
                markersize=8, label="LD-block-matched (bootstrap CI)"),
        Line2D([0], [0], marker="*", color="w",
                markerfacecolor=COL_BLUE_LIGHT, markeredgecolor="black",
                markersize=14, label="actual chr_439 signal block"),
        Line2D([0], [0], marker="*", color="w",
                markerfacecolor=COL_BLUE_DK, markeredgecolor="black",
                markersize=14, label="actual chr_439 LD block"),
    ]
    axA.legend(handles=handles, fontsize=7, frameon=False,
                 loc="lower right")
    for sp in ("top", "right"): axA.spines[sp].set_visible(False)


    axB.axhline(0, color="black", lw=0.5, alpha=0.5)
    axB.axvline(0, color="black", lw=0.5, alpha=0.5)
    axB.plot([-0.5, 0.5], [-0.5, 0.5], color="black", lw=0.7, ls="--",
              alpha=0.4)
    axB.scatter(bg_other["lab_dp"], bg_other["wild_dp"],
                 s=4, color="#cccccc", alpha=0.30, edgecolors="none",
                 rasterized=True, zorder=2,
                 label=f"genome-wide bg (n={len(bg_other):,})")
    axB.scatter(sig_dp["lab_dp"], sig_dp["wild_dp"],
                 s=85, facecolor=COL_PURPLE, edgecolor="black",
                 linewidth=0.9, zorder=10,
                 label=f"chr_439 signal block (n={len(sig_dp)})")
    if len(sig_dp) >= 5:
        b, a = np.polyfit(sig_dp["lab_dp"].dropna(),
                            sig_dp["wild_dp"].dropna(), 1)
        xs = np.linspace(sig_dp["lab_dp"].min(),
                          sig_dp["lab_dp"].max(), 50)
        axB.plot(xs, a + b * xs, color=COL_PURPLE, lw=2.0, zorder=8)
    txt = (f"genome-wide:    ρ={rho_other:+.3f}  p={p_other:.1e}\n"
           f"signal block:   ρ={rho_sig:+.3f}  p={p_sig:.1e}")
    axB.text(0.04, 0.97, txt, transform=axB.transAxes,
              ha="left", va="top", fontsize=8, family="monospace",
              bbox=dict(facecolor="white", alpha=0.95, edgecolor="#aaa",
                         linewidth=0.4, boxstyle="round,pad=0.3"))
    axB.set_xlim(-0.45, 0.45); axB.set_ylim(-0.45, 0.45)
    axB.set_xlabel(r"Lab $\Delta p$  (B G10 − T G10)", fontsize=9.5)
    axB.set_ylabel(r"Wild $\Delta p$  (B − T)", fontsize=9.5)
    axB.set_title("B. Lab Δp vs wild Δp — chr_439 signal block",
                   fontsize=10, pad=4)
    axB.tick_params(labelsize=8.5)
    axB.legend(fontsize=8, frameon=False, loc="lower right")
    for sp in ("top", "right"): axB.spines[sp].set_visible(False)

    fig.suptitle(
        "Lab vs wild antagonism at chr_439 signal block",
        fontsize=11, y=0.995)
    fig.tight_layout()
    fig.savefig(f"{OUT_BASE}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.svg", bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.pdf", bbox_inches="tight")
    strip_svg(Path(f"{OUT_BASE}.svg"))

if __name__ == "__main__":
    main()
