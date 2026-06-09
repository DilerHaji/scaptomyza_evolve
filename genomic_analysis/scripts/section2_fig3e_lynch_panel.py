#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(".")
LYNCH_DIR = ROOT / "lynch_s"
WB_TSV = ROOT / "final_plots/wild/section2_fig3e_lynch_chr439_wellbehaved.tsv"
OUT_DIR = ROOT / "final_plots/wild"

C = 0.5
N_INTERVALS_WINDOW = 3   # rolling window of 3 last intervals
G_MAX_VALUES = [6, 7, 8, 9, 10]

W_AA_B_P, W_Aa_B_P, W_aa_B_P = 1.000, 0.913, 0.612
W_AA_T_P, W_Aa_T_P, W_aa_T_P = 0.611, 0.806, 0.902
sB_P = (W_AA_B_P - W_aa_B_P) / W_AA_B_P
hB_P = (W_Aa_B_P - W_aa_B_P) / (W_AA_B_P - W_aa_B_P)
sT_P = (W_aa_T_P - W_AA_T_P) / W_aa_T_P
hT_P = (W_Aa_T_P - W_AA_T_P) / (W_aa_T_P - W_AA_T_P)
qstar_P = (C * sB_P * hB_P) / (C * sB_P * hB_P + (1 - C) * sT_P * hT_P)
qstar_P_h05 = sB_P / (sB_P + sT_P)


def landscape_curves(W_AA_B, W_Aa_B, W_aa_B, W_AA_T, W_Aa_T, W_aa_T,
                      c=C, q=None):
    if q is None:
        q = np.linspace(0, 1, 400)
    Wb = q*q*W_AA_B + 2*q*(1-q)*W_Aa_B + (1-q)*(1-q)*W_aa_B
    Wt = q*q*W_AA_T + 2*q*(1-q)*W_Aa_T + (1-q)*(1-q)*W_aa_T
    Wm = c * Wb + (1 - c) * Wt
    return q, Wb, Wt, Wm


def plot_minimal_landscape(ax, Wgt, qstar, qstar_add, label):
    q, Wb, Wt, Wm = landscape_curves(*Wgt)
    norm = Wm.max()
    Wb_n, Wt_n, Wm_n = Wb / norm, Wt / norm, Wm / norm

    ax.plot(q, Wb_n, color="#a04848", linewidth=1.1)
    ax.plot(q, Wt_n, color="#3d6cb5", linewidth=1.1)
    ax.plot(q, Wm_n, color="#1a1a1a", linewidth=1.7)
    ax.axvline(qstar, color="#1a5e1a", linestyle="--", linewidth=1.0, zorder=4)
    ax.scatter([qstar], [Wm_n.max()], s=45, color="#1a5e1a",
                marker="v", edgecolor="black", linewidth=0.5, zorder=6)
    ax.axvline(qstar_add, color="#7fb185", linestyle=":", linewidth=0.9,
                alpha=0.85, zorder=3)
    ax.scatter([qstar_add], [Wm_n.max()], s=30, color="#7fb185",
                marker="v", edgecolor="#444", linewidth=0.4,
                alpha=0.85, zorder=5)

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


def qstar_dempster(sB, sT_neg, hB=hB_P, hT=hT_P, c=C):
    if not (np.isfinite(sB) and np.isfinite(sT_neg) and sB > 0 and sT_neg < 0):
        return np.nan
    return (c * sB * hB) / (c * sB * hB + (1 - c) * abs(sT_neg) * hT)


def qstar_h05(sB, sT_neg, c=C):
    if not (np.isfinite(sB) and np.isfinite(sT_neg) and sB > 0 and sT_neg < 0):
        return np.nan
    return c * sB / (c * sB + (1 - c) * abs(sT_neg))


def lynch_columns_for_window(treatment_letter, rep, g_max,
                                n_intervals=N_INTERVALS_WINDOW):
    cols = []
    for g_end in range(g_max, g_max - n_intervals, -1):
        g_start = g_end - 1
        if g_start < 1:
            continue
        col = (f"s_{treatment_letter}{rep}G{g_start:02d}_"
                f"{treatment_letter}{rep}G{g_end:02d}")
        cols.append(col)
    return cols


def load_wb_lynch_for_window(wb_pos_set, treatment_letter, g_max):
    dfs = []
    for rep in range(1, 5):
        f = LYNCH_DIR / f"d{treatment_letter}{rep}_per_interval_s.csv"
        with open(f) as fh:
            header = fh.readline().rstrip().split(",")
        existing_s_cols = [c for c in header if c.startswith("s_")]
        keep_cols = []
        for c in existing_s_cols:
            try:
                parts = c.split("_")
                g_start = int(parts[1][3:])
                g_end = int(parts[2][3:])
            except (ValueError, IndexError):
                continue
            window_lo = g_max - N_INTERVALS_WINDOW
            if g_end <= g_max and g_start >= window_lo:
                keep_cols.append(c)
        cols = ["CHROM", "POS"] + keep_cols
        if not keep_cols:
            continue
        d = pd.read_csv(f, usecols=cols, low_memory=False)
        d = d[d["POS"].isin(wb_pos_set)].copy()
        d[f"sR{rep}"] = d[keep_cols].mean(axis=1, skipna=True)
        d = d[["CHROM", "POS", f"sR{rep}"]]
        dfs.append(d)
    if not dfs:
        return pd.DataFrame()
    out = dfs[0]
    for d in dfs[1:]:
        out = out.merge(d, on=["CHROM", "POS"], how="outer")
    rep_cols = [c for c in out.columns if c.startswith("sR")]
    out["s_treatment"] = out[rep_cols].mean(axis=1, skipna=True)
    return out


def main():
    plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42,
                          "font.family": "sans-serif",
                          "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
                          "axes.linewidth": 0.5})

    wb = pd.read_csv(WB_TSV, sep="\t")
    wb_pos = set(wb["POS"].astype(int).tolist())

    rows = []
    g10_per_snp = None  # to save the G10 per-SNP data for landscape derivation
    for g_max in G_MAX_VALUES:
        df_B = load_wb_lynch_for_window(wb_pos, "B", g_max)
        df_T = load_wb_lynch_for_window(wb_pos, "T", g_max)
        df = df_B.merge(df_T[["CHROM", "POS", "s_treatment"]],
                          on=["CHROM", "POS"], suffixes=("_B", "_T"))
        df = df.rename(columns={"s_treatment_B": "sB_pre",
                                  "s_treatment_T": "sT_pre"})
        # Polarize so sB > 0 per SNP
        flip = df["sB_pre"] < 0
        df["sB"] = df["sB_pre"]
        df["sT"] = df["sT_pre"]
        df.loc[flip, "sB"] = -df.loc[flip, "sB_pre"]
        df.loc[flip, "sT"] = -df.loc[flip, "sT_pre"]
        if g_max == 10:
            g10_per_snp = df.copy()

        n_valid = (df["sB"].notna() & df["sT"].notna()).sum()
        med_sB = df["sB"].median()
        med_sT = df["sT"].median()
        q_phenoH = qstar_dempster(med_sB, med_sT)
        q_h05 = qstar_h05(med_sB, med_sT)

        sub = df.dropna(subset=["sB", "sT"])
        n = len(sub)
        sB_v, sT_v = sub["sB"].values, sub["sT"].values
        boot_pH = np.empty(2000); boot_05 = np.empty(2000)
        rng = np.random.default_rng(g_max * 10)
        for i in range(2000):
            idx = rng.integers(0, n, size=n)
            boot_pH[i] = qstar_dempster(np.median(sB_v[idx]),
                                          np.median(sT_v[idx]))
            boot_05[i] = qstar_h05(np.median(sB_v[idx]),
                                     np.median(sT_v[idx]))
        valid_pH = ~np.isnan(boot_pH)
        q_pH_lo = float(np.nanquantile(boot_pH, 0.025)) if valid_pH.any() else np.nan
        q_pH_hi = float(np.nanquantile(boot_pH, 0.975)) if valid_pH.any() else np.nan
        valid_05 = ~np.isnan(boot_05)
        q_05_lo = float(np.nanquantile(boot_05, 0.025)) if valid_05.any() else np.nan
        q_05_hi = float(np.nanquantile(boot_05, 0.975)) if valid_05.any() else np.nan

        rows.append({"g_max": g_max, "n_snps": n_valid,
                      "med_sB": med_sB, "med_sT": med_sT,
                      "qstar_phenoH": q_phenoH, "qstar_h05": q_h05,
                      "qP_lo": q_pH_lo, "qP_hi": q_pH_hi,
                      "q05_lo": q_05_lo, "q05_hi": q_05_hi,
                      "valid_frac_phenoH": valid_pH.mean()})

    R = pd.DataFrame(rows)
    R.to_csv(OUT_DIR / "section2_fig3e_lynch_panel.tsv",
              sep="\t", index=False)

    g10_row = R[R["g_max"] == 10].iloc[0]
    sB_G, sT_G = float(g10_row["med_sB"]), float(g10_row["med_sT"])
    qstar_G = float(g10_row["qstar_phenoH"])
    qstar_G_add = float(g10_row["qstar_h05"])
    Wgt_G = (1 + sB_G, 1 + sB_G * hB_P, 1.0,
             1.0, 1 + abs(sT_G) * hT_P, 1 + abs(sT_G))
    Wgt_P = (W_AA_B_P, W_Aa_B_P, W_aa_B_P,
             W_AA_T_P, W_Aa_T_P, W_aa_T_P)

    fig = plt.figure(figsize=(4.4, 4.6))
    gs = fig.add_gridspec(2, 2,
                            height_ratios=[1.0, 1.4],
                            width_ratios=[1, 1],
                            hspace=0.50, wspace=0.30)
    ax_p = fig.add_subplot(gs[0, 0])
    ax_g = fig.add_subplot(gs[0, 1])
    ax_t = fig.add_subplot(gs[1, :])

    plot_minimal_landscape(ax_p, Wgt_P, qstar_P, qstar_P_h05,
                            f"Phenotypic   q*={qstar_P:.2f} | {qstar_P_h05:.2f}")
    plot_minimal_landscape(ax_g, Wgt_G, qstar_G, qstar_G_add,
                            f"Genomic   q*={qstar_G:.2f} | {qstar_G_add:.2f}")
    ax_p.set_ylabel("normalized\nmean fitness", fontsize=6.5, labelpad=1)

    ax_t.axhline(qstar_P, color="#a04848", linewidth=1.0, alpha=0.85,
                  zorder=2, label=f"q*_pheno (pheno h)={qstar_P:.3f}")
    ax_t.axhline(qstar_P_h05, color="#d09a9a", linestyle="--", linewidth=0.9,
                  alpha=0.85, zorder=2,
                  label=f"q*_pheno (h=0.5)={qstar_P_h05:.3f}")
    ax_t.axhline(0.5, color="#888", linestyle=":", linewidth=0.4, alpha=0.6, zorder=1)

    R_all = R.sort_values("g_max").reset_index(drop=True)
    Y_ND = 0.95
    def draw_traj(qcol, lo_col, hi_col, color, label, marker, linestyle,
                    alpha, label_offset):
        y_pos = [r[qcol] if pd.notna(r[qcol]) else Y_ND for _, r in R_all.iterrows()]
        is_valid = [pd.notna(r[qcol]) for _, r in R_all.iterrows()]
        for i in range(len(R_all) - 1):
            x0, x1 = R_all["g_max"].iloc[i], R_all["g_max"].iloc[i+1]
            y0, y1 = y_pos[i], y_pos[i+1]
            both = is_valid[i] and is_valid[i+1]
            ax_t.plot([x0, x1], [y0, y1],
                       linestyle=linestyle if both else ":",
                       linewidth=1.2 if both else 0.7, color=color,
                       alpha=alpha if both else alpha * 0.5, zorder=3)
        valid = R_all[R_all[qcol].notna()]
        ax_t.errorbar(valid["g_max"], valid[qcol],
                       yerr=[valid[qcol] - valid[lo_col],
                              valid[hi_col] - valid[qcol]],
                       fmt=marker, color=color, markersize=4.5,
                       capsize=2, linewidth=0, elinewidth=1.0, alpha=alpha,
                       label=label, zorder=5)
        for _, r in valid.iterrows():
            dy = 0.05 if label_offset > 0 else -0.07
            ax_t.text(r["g_max"], r[qcol] + dy, f"{r[qcol]:.2f}",
                       ha="center", fontsize=5.5, color=color, alpha=alpha)
        for i, v in enumerate(is_valid):
            if not v:
                gm = R_all["g_max"].iloc[i]
                ax_t.scatter([gm], [Y_ND], s=40, facecolors="none",
                              edgecolors=color, linewidth=1.0, zorder=5)
                ax_t.text(gm, Y_ND + 0.03, "n.d.", ha="center", va="bottom",
                           fontsize=5.5, color=color, style="italic")

    draw_traj("qstar_phenoH", "qP_lo", "qP_hi", "#1a5e1a",
                "q* (pheno h: 0.78, 0.67)", "o", "-", 1.0, +1)
    draw_traj("qstar_h05", "q05_lo", "q05_hi", "#7fb185",
                "q* (additive h=0.5)", "s", "--", 0.85, -1)

    ax_t.set_xlabel(
        f"included generations  (G_max, last {N_INTERVALS_WINDOW} intervals)",
        fontsize=7)
    ax_t.set_ylabel("q*", fontsize=7)
    ax_t.set_xticks(R_all["g_max"])
    ax_t.set_xticklabels([f"G{g}" for g in R_all["g_max"]], fontsize=6)
    ax_t.set_ylim(0.0, 1.05)
    ax_t.tick_params(labelsize=6, length=2, pad=1)
    ax_t.legend(loc="lower left", fontsize=5.5, frameon=False)
    ax_t.set_title(f"Lynch q*_geno (n={int(g10_row['n_snps'])} well-behaved chr_439 SNPs)",
                    fontsize=8)

    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.10, top=0.94)
    save(fig, "section2_fig3e_lynch_panel")


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
