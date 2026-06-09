#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(".")
DIVERSITY_TSV = ROOT / "grenfst/diversity/trajectory_pi_200000_nonoverlap.csv"
QSTAR_TSV = ROOT / "final_plots/wild/section2_wild_qstar_correlation_gw_top50.tsv"
OUT_BASE = ROOT / "final_plots/wild/section2_wild_qstar_diversity_filter_v3"

INPUT_TOP_PCT = 0.50

CHROM_439 = "chr_ScDA7r2_439_HRSCAF_779"
LD_S, LD_E = 2_640_000, 3_610_000


def main():
    plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42,
                          "font.family": "sans-serif",
                          "font.sans-serif": ["Helvetica","Arial","DejaVu Sans"],
                          "axes.linewidth": 0.7})

    div = pd.read_csv(DIVERSITY_TSV)
    div["start"] = div["start"].astype(int)

    for arm in ["B", "T", "M"]:
        rep_cols = [f"{arm}{r}G10.theta_watterson" for r in range(1, 5)
                    if f"{arm}{r}G10.theta_watterson" in div.columns]
        div[f"{arm}_SD_G10"] = div[rep_cols].std(axis=1, ddof=1)
        div[f"{arm}_mean_G10"] = div[rep_cols].mean(axis=1)
        ps_cols = [f"{arm}{r}G10.passed" for r in range(1, 5)
                   if f"{arm}{r}G10.passed" in div.columns]
        div[f"{arm}_SD_npass_G10"] = div[ps_cols].std(axis=1, ddof=1)

    div = div.dropna(subset=["B_SD_G10", "T_SD_G10", "M_SD_G10"]).copy()

    bSD = div['B_SD_G10'].mean(); tSD = div['T_SD_G10'].mean()
    mSD = div['M_SD_G10'].mean()
    bSDn = div['B_SD_npass_G10'].mean(); tSDn = div['T_SD_npass_G10'].mean()
    mSDn = div['M_SD_npass_G10'].mean()
    chr110 = div[div["chrom"] == "chr_ScDA7r2_110_HRSCAF_368"]
    div["M_var_excess"] = (div["M_SD_G10"]
                              - 0.5 * (div["B_SD_G10"] + div["T_SD_G10"]))
    div["M_var_excess_npass"] = (div["M_SD_npass_G10"]
                                    - 0.5 * (div["B_SD_npass_G10"]
                                              + div["T_SD_npass_G10"]))
    n_pos = int((div["M_var_excess"] > 0).sum())
    chr439_div = div[(div["chrom"] == CHROM_439) &
                       (div["start"] >= LD_S) &
                       (div["end"] <= LD_E)]
    pct_chr439 = (
        (div["M_var_excess"].rank(pct=True) >
         (1 - len(chr439_div[chr439_div["M_var_excess"] >
                                 chr439_div["M_var_excess"].median()]) / len(chr439_div)))
    )
    chr439_pct_rank = chr439_div["M_var_excess"].apply(
        lambda v: (div["M_var_excess"] <= v).mean()).mean()

    df = pd.read_csv(QSTAR_TSV, sep="\t")
    df = df[df["antag"] & df["qstar_pred"].notna() &
             df["wild_AF_pol"].notna()].copy()

    WILD_POOLS_ALL = ["AVB", "PSB", "RMB", "AVT", "PST", "RMT"]
    af = pd.read_csv(ROOT / "final_plots/wild/af_matrix_22pools.csv",
                       usecols=["chrom_pos"] + WILD_POOLS_ALL)
    af[["chrom", "pos"]] = af["chrom_pos"].str.split(":", expand=True)
    af["pos"] = af["pos"].astype(int)
    af = af.drop(columns=["chrom_pos"])
    df = df.merge(af, on=["chrom", "pos"], how="left")
    pool_arr = df[WILD_POOLS_ALL].values
    pol = df["polarized"].values if "polarized" in df.columns else \
            (df["wild_AF_pol"].values != df["wild_AF_orig"].values)
    pol_pool = np.where(pol[:, None], 1 - pool_arr, pool_arr)
    df["wild_AF_SD"] = np.nanstd(pol_pool, axis=1, ddof=1)
    n_valid = np.sum(~np.isnan(pol_pool), axis=1)
    df["wild_AF_SE"] = df["wild_AF_SD"] / np.sqrt(np.clip(n_valid, 1, None))

    snp_to_div = {}
    for chrom, sub_div in div.groupby("chrom"):
        sub = sub_div.sort_values("start").reset_index(drop=True)
        starts = sub["start"].values
        ends = sub["end"].values
        snps_chrom = df[df["chrom"] == chrom]["pos"].values
        for pos in snps_chrom:
            idx = np.searchsorted(starts, pos, side="right") - 1
            if 0 <= idx < len(starts) and starts[idx] <= pos <= ends[idx]:
                snp_to_div[(chrom, int(pos))] = int(starts[idx])
            else:
                snp_to_div[(chrom, int(pos))] = -1
    df["div_window_start"] = [snp_to_div.get((c, int(p)), -1)
                                for c, p in zip(df["chrom"], df["pos"])]
    n_mapped = (df["div_window_start"] >= 0).sum()

    div_lookup = {(r["chrom"], int(r["start"])): float(r["M_var_excess"])
                  for _, r in div.iterrows()}
    div_lookup_npass = {(r["chrom"], int(r["start"])):
                          float(r["M_var_excess_npass"])
                          for _, r in div.iterrows()}
    df["M_var_excess"] = [div_lookup.get((c, int(b)), np.nan)
                          for c, b in zip(df["chrom"], df["div_window_start"])]
    df["M_var_excess_npass"] = [div_lookup_npass.get((c, int(b)), np.nan)
                                  for c, b in zip(df["chrom"],
                                                  df["div_window_start"])]

    results = []
    lrt_thresholds = [
        ("top 50%", 0.50),
        ("top 5%", 0.95),
        ("top 0.1%", 0.999),
    ]
    div_filters = [
        ("All windows", None),
        ("M_var_excess > 0 (M more idiosyncratic)", "pos"),
        ("M_var_excess ≥ Q50", ("q", 0.50)),
        ("M_var_excess ≥ Q75", ("q", 0.75)),
    ]
    for lrt_label, lrt_q in lrt_thresholds:
        q_local = max(0.0, (lrt_q - INPUT_TOP_PCT) / (1 - INPUT_TOP_PCT))
        thr = df["lrt"].quantile(q_local)
        lrt_set = df[df["lrt"] >= thr]
        for div_label, div_key in div_filters:
            if div_key is None:
                sub = lrt_set
            elif div_key == "pos":
                sub = lrt_set[lrt_set["M_var_excess"] > 0]
            elif div_key == "neg":
                sub = lrt_set[lrt_set["M_var_excess"] < 0]
            elif div_key[0] == "q":
                pct = div_key[1]
                thrM = div["M_var_excess"].quantile(pct)
                sub = lrt_set[lrt_set["M_var_excess"] >= thrM]
            sub = sub.dropna(subset=["M_var_excess"])
            if len(sub) < 10:
                continue
            rs, ps = spearmanr(sub["qstar_pred"], sub["wild_AF_pol"])
            slope_, _ = np.polyfit(sub["qstar_pred"],
                                       sub["wild_AF_pol"], 1)
            results.append({"lrt": lrt_label, "div_filter": div_label,
                              "n": len(sub), "spearman": float(rs),
                              "spearman_p": float(ps), "slope": float(slope_)})
    R = pd.DataFrame(results)
    R.to_csv(f"{OUT_BASE}.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.0))

    ax = axes[0]
    div_labels = [l for l, _ in div_filters]
    div_colors = ["#1a1a1a", "#1a5e1a", "#3d6cb5", "#E07B2F"]
    div_markers = ["o", "D", "^", "s"]
    lrt_label_to_idx = {t[0]: i for i, t in enumerate(lrt_thresholds)}
    for div_label, color, marker in zip(div_labels, div_colors, div_markers):
        sub = R[R["div_filter"] == div_label].copy()
        if len(sub) == 0:
            continue
        xs = [lrt_label_to_idx[l] for l in sub["lrt"].values]
        ax.plot(xs, sub["spearman"].values,
                  marker=marker, color=color, linewidth=1.4, markersize=7,
                  markeredgecolor="black", markeredgewidth=0.4,
                  label=div_label)
        for i, (_, r) in zip(xs, sub.iterrows()):
            stars = ("***" if r["spearman_p"] < 0.001 else
                      "**" if r["spearman_p"] < 0.01 else
                      "*" if r["spearman_p"] < 0.05 else "")
            if stars:
                ax.text(i, r["spearman"] + 0.005, stars,
                          ha="center", fontsize=7, color=color,
                          fontweight="bold")
    ax.axhline(0, color="#888", linestyle=":", linewidth=0.5)
    ax.set_xticks(range(len(lrt_thresholds)))
    ax.set_xticklabels([t[0] for t in lrt_thresholds],
                          fontsize=8, rotation=15, ha="right")
    ax.set_ylabel("Spearman ρ\n(q*_pred vs wild_AF_polarized)", fontsize=9)
    ax.set_xlabel("LRT stringency", fontsize=9)
    ax.set_title("Fig-2-aligned filter: between-rep SD excess in M",
                  fontsize=10)
    ax.legend(loc="upper left", fontsize=6.5, frameon=False, ncol=1)
    ax.tick_params(labelsize=8)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)

    ax = axes[1]
    sd_avg_BT = 0.5 * (div["B_SD_G10"] + div["T_SD_G10"])
    sc = ax.scatter(sd_avg_BT, div["M_SD_G10"],
                      c=div["M_var_excess"], cmap="RdBu_r",
                      vmin=-div["M_var_excess"].abs().quantile(0.95),
                      vmax=div["M_var_excess"].abs().quantile(0.95),
                      s=10, alpha=0.6, edgecolor="none")
    cbar = plt.colorbar(sc, ax=ax, label="M_var_excess")
    lo = min(sd_avg_BT.min(), div["M_SD_G10"].min())
    hi = max(sd_avg_BT.max(), div["M_SD_G10"].max())
    ax.plot([lo, hi], [lo, hi], color="#333", linestyle=":",
              linewidth=0.7, label="identity (SD_M = SD_avg(B,T))")
    chr439_avg = 0.5 * (chr439_div["B_SD_G10"] + chr439_div["T_SD_G10"])
    ax.scatter(chr439_avg, chr439_div["M_SD_G10"],
                s=80, color="none", edgecolor="black", linewidth=1.2,
                zorder=5, label=f"chr_439 LD (n={len(chr439_div)})")
    ax.set_xlabel("Mean of (SD across B reps, SD across T reps)\n"
                    "θ_W at G10", fontsize=9)
    ax.set_ylabel("SD across M reps, θ_W at G10", fontsize=9)
    ax.set_title("Per-window between-rep variability\n"
                  "Above identity = M more variable across reps "
                  "than B,T (Fig 2C signature)", fontsize=10)
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    ax.tick_params(labelsize=8)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)

    ax = axes[2]
    q_local_top5 = max(0.0, (0.95 - INPUT_TOP_PCT) / (1 - INPUT_TOP_PCT))
    thr_top5 = df["lrt"].quantile(q_local_top5)
    sc_df = df[(df["lrt"] >= thr_top5)
                & df["qstar_pred"].notna()
                & df["wild_AF_pol"].notna()
                & df["M_var_excess"].notna()].copy()
    ax.scatter(sc_df["qstar_pred"], sc_df["wild_AF_pol"],
                 s=5, color="#D9D9D9", alpha=0.5, edgecolor="none",
                 rasterized=True)

    q_local_top01 = max(0.0, (0.999 - INPUT_TOP_PCT) / (1 - INPUT_TOP_PCT))
    thr_top01 = df["lrt"].quantile(q_local_top01)
    thrM_q75_pre = div["M_var_excess"].quantile(0.75)
    hl = sc_df[(sc_df["lrt"] >= thr_top01)
                & (sc_df["M_var_excess"] >= thrM_q75_pre)]
    ax.errorbar(hl["qstar_pred"], hl["wild_AF_pol"],
                  yerr=hl["wild_AF_SE"].values,
                  fmt="none", ecolor="#E07B2F", elinewidth=0.7,
                  capsize=2, capthick=0.6, alpha=0.85, zorder=4)
    ax.scatter(hl["qstar_pred"], hl["wild_AF_pol"],
                 s=28, color="#E07B2F", alpha=0.95,
                 edgecolor="black", linewidth=0.5, zorder=5,
                 label=f"top 0.1% LRT × Q75 (n={len(hl)})")

    ax.plot([0, 1], [0, 1], color="#333", linestyle=":",
              linewidth=0.7, label="y = x")
    x_fit = np.linspace(0, 1, 100)

    fits = []
    fits.append(("All top 5%", sc_df, "#1a1a1a", "--"))
    fits.append(("M_var_excess > 0",
                  sc_df[sc_df["M_var_excess"] > 0], "#1a5e1a", "-"))
    thrM_q50 = div["M_var_excess"].quantile(0.50)
    fits.append(("M_var_excess ≥ Q50",
                  sc_df[sc_df["M_var_excess"] >= thrM_q50], "#3d6cb5", "-"))
    thrM_q75 = div["M_var_excess"].quantile(0.75)
    fits.append(("M_var_excess ≥ Q75",
                  sc_df[sc_df["M_var_excess"] >= thrM_q75], "#E07B2F", "-"))
    for fit_label, sub, color, ls in fits:
        s_, i_ = np.polyfit(sub["qstar_pred"], sub["wild_AF_pol"], 1)
        rs_, ps_ = spearmanr(sub["qstar_pred"], sub["wild_AF_pol"])
        ax.plot(x_fit, s_ * x_fit + i_, color=color, linewidth=1.4,
                  linestyle=ls,
                  label=f"{fit_label}  ρ={rs_:+.3f}  slope={s_:+.2f}")
    ax.axhline(0.5, color="#888", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.axvline(0.5, color="#888", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.set_xlim(0.2, 0.8); ax.set_ylim(0.2, 0.8)
    ax.set_xlabel("q*_predicted from lab Dempster\n(per-SNP, additive h=0.5)",
                   fontsize=9)
    ax.set_ylabel("Wild allele frequency (polarized to lab B-favored)",
                   fontsize=9)
    ax.set_title(f"Top 5% LRT (n={len(sc_df):,}); color = M_var_excess",
                  fontsize=10)
    ax.legend(loc="lower right", fontsize=7.5, frameon=False)
    ax.tick_params(labelsize=8)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)

    fig.suptitle("Filter top-LRT SNPs by Fig-2 between-rep SD excess "
                  "(M idiosyncratic retention)", fontsize=11)
    fig.subplots_adjust(left=0.05, right=0.985, bottom=0.13, top=0.86,
                          wspace=0.32)
    fig.savefig(f"{OUT_BASE}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.svg", bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
