#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr, pearsonr

ROOT = Path(".")
GLM_CSV = ROOT / "glm_lrt_gw_final/glmV1full.csv"
AF_MATRIX = ROOT / "final_plots/wild/af_matrix_22pools.csv"
AD_TSV = ROOT / "variance_analysis/merged_ad.tsv"
SAMPLES = ROOT / "variance_analysis/sample_list.txt"
OUT_BASE = ROOT / "final_plots/wild/section2_wild_qstar_correlation_gw_top50_perrep"

WILD_B = ["AVB","PSB","RMB"]; WILD_T = ["AVT","PST","RMT"]
FOUNDERS = ["F1G00","F2G00","F3G00","F4G00"]
CHROM_439 = "chr_ScDA7r2_439_HRSCAF_779"
LD_S, LD_E = 2_640_000, 3_610_000

GENS = list(range(1, 11))
MIN_DEPTH = 30
MIN_AF, MAX_AF = 0.05, 0.95
TOP_PCT = 0.50 


def parse_samples():
    return {n: i + 4 for i, n in enumerate(open(SAMPLES).read().splitlines())}


def parse_ad(s):
    if "," not in s: return None, None
    try: r, a = (int(x) for x in s.split(","))
    except (ValueError, AttributeError): return None, None
    return r, a


def logit(p, eps=1e-3):
    p = np.clip(p, eps, 1 - eps)
    return np.log(p / (1 - p))


def slope(af, gens):
    if len(af) < 4: return np.nan
    x = np.array(gens, dtype=float)
    y = logit(np.array(af))
    if not np.all(np.isfinite(y)): return np.nan
    if np.std(x) == 0: return np.nan
    s, *_ = stats.linregress(x, y)
    return s


def main():
    plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42,
                          "font.family": "sans-serif",
                          "font.sans-serif": ["Helvetica","Arial","DejaVu Sans"],
                          "axes.linewidth": 0.7})

    cached = Path(f"{OUT_BASE}.tsv")
    if cached.exists():
        df = pd.read_csv(cached, sep="\t")
    else:
        g = pd.read_csv(GLM_CSV, usecols=["chrom","pos","LRT_chisq",
                                           "converged","error"])
        g = g[(g.converged == True) & (g.error == "OK") &
              g.LRT_chisq.notna()]
        thr = g.LRT_chisq.quantile(TOP_PCT)
        g_top = g[g.LRT_chisq >= thr][["chrom","pos","LRT_chisq"]]
        target_set = set(zip(g_top["chrom"], g_top["pos"].astype(int)))
        chrom_pos = {}
        for chrom, pos in target_set:
            chrom_pos.setdefault(chrom, set()).add(pos)
        lrt_lookup = dict(zip(zip(g_top["chrom"], g_top["pos"].astype(int)),
                                g_top["LRT_chisq"]))

        cols = parse_samples()
        F_idx = [cols[f"F{r}G00"] for r in range(1, 5) if f"F{r}G00" in cols]
        B_per_rep = {r: {g: cols[f"B{r}G{g:02d}"] for g in GENS
                          if f"B{r}G{g:02d}" in cols} for r in range(1, 5)}
        T_per_rep = {r: {g: cols[f"T{r}G{g:02d}"] for g in GENS
                          if f"T{r}G{g:02d}" in cols} for r in range(1, 5)}

        def af_at_cell(f, ci):
            ref, alt = parse_ad(f[ci])
            if ref is not None and ref + alt >= MIN_DEPTH:
                return alt / (ref + alt)
            return np.nan

        def af_at_idx_list(f, idx_list):
            vals = [af_at_cell(f, ci) for ci in idx_list]
            vals = [v for v in vals if not np.isnan(v)]
            return float(np.mean(vals)) if len(vals) >= 2 else np.nan

        def per_rep_slopes(f, per_rep_map):
            out = []
            for r, gen_to_col in per_rep_map.items():
                pts = [(g, af_at_cell(f, ci))
                       for g, ci in gen_to_col.items()]
                pts = [(g, a) for g, a in pts if not np.isnan(a)]
                if len(pts) < 4:
                    continue
                gs, afs = zip(*pts)
                s = slope(list(afs), list(gs))
                if np.isfinite(s):
                    out.append(s)
            return out

        F_all = [cols[f"F{r}G00"] for r in range(1, 5) if f"F{r}G00" in cols]
        rows_out = []
        n_seen = 0
        with open(AD_TSV) as fh:
            for line in fh:
                f = line.rstrip().split("\t")
                chrom, pos = f[0], int(f[1])
                ps = chrom_pos.get(chrom)
                if ps is None or pos not in ps:
                    continue
                af_F = af_at_idx_list(f, F_all)
                if np.isnan(af_F) or not (MIN_AF <= af_F <= MAX_AF):
                    continue
                slopes_B = per_rep_slopes(f, B_per_rep)
                slopes_T = per_rep_slopes(f, T_per_rep)
                if len(slopes_B) < 2 or len(slopes_T) < 2:
                    continue
                sB = float(np.mean(slopes_B))
                sT = float(np.mean(slopes_T))
                sB_SE = float(np.std(slopes_B, ddof=1) / np.sqrt(len(slopes_B)))
                sT_SE = float(np.std(slopes_T, ddof=1) / np.sqrt(len(slopes_T)))
                polarized = False
                if sB < 0:
                    sB, sT = -sB, -sT
                    polarized = True
                af_F_pol = af_F if not polarized else (1 - af_F)
                rows_out.append({
                    "chrom": chrom, "pos": pos,
                    "lrt": lrt_lookup[(chrom, pos)],
                    "founder_AF_orig": af_F,
                    "founder_AF_pol": af_F_pol,
                    "polarized": polarized,
                    "sB": sB, "sT": sT,
                    "sB_SE": sB_SE, "sT_SE": sT_SE,
                    "n_reps_B": len(slopes_B),
                    "n_reps_T": len(slopes_T),
                })
                n_seen += 1
                if n_seen % 5000 == 0:
                    pass
        df_lab = pd.DataFrame(rows_out)

        cols_af = ["chrom_pos"] + WILD_B + WILD_T
        af = pd.read_csv(AF_MATRIX, usecols=cols_af)
        af[["chrom","pos"]] = af["chrom_pos"].str.split(":", expand=True)
        af["pos"] = af["pos"].astype(int)
        af["wild_AF_orig"] = af[WILD_B + WILD_T].mean(axis=1)
        df = df_lab.merge(af[["chrom","pos","wild_AF_orig"]],
                            on=["chrom","pos"], how="inner")

        df["wild_AF_pol"] = np.where(df["polarized"],
                                       1 - df["wild_AF_orig"],
                                       df["wild_AF_orig"])

        df["antag"] = (df["sB"] > 0) & (df["sT"] < 0)
        df["qstar_pred"] = np.where(df["antag"],
                                       df["sB"] / (df["sB"] + np.abs(df["sT"])),
                                       np.nan)

        denom4 = (df["sB"] - df["sT"]).pow(4)
        var_q = (df["sT"].pow(2) * df["sB_SE"].pow(2)
                  + df["sB"].pow(2) * df["sT_SE"].pow(2)) / denom4
        df["qstar_pred_SE"] = np.where(df["antag"], np.sqrt(var_q), np.nan)
        df["in_chr439_LD"] = ((df["chrom"] == CHROM_439) &
                                (df["pos"] >= LD_S) &
                                (df["pos"] <= LD_E))
        df.to_csv(cached, sep="\t", index=False)

    df_antag = df[df["antag"] & df["qstar_pred"].notna() &
                   df["wild_AF_pol"].notna()].copy()

    sub_thresholds = [
        ("top 5% LRT", 0.95),
        ("top 1% LRT", 0.99),
        ("top 0.5% LRT", 0.995),
        ("top 0.1% LRT", 0.999),
    ]
    rows_results = []
    for label, q in sub_thresholds:
        thr = df_antag["lrt"].quantile((q - 0.95) / (1 - 0.95))
        sub = df_antag[df_antag["lrt"] >= thr]
        for set_label, mask in [
            (f"{label} - all genome-wide", np.ones(len(sub), dtype=bool)),
            (f"{label} - chr_439 LD block", sub["in_chr439_LD"].values),
            (f"{label} - rest of genome",  ~sub["in_chr439_LD"].values),
        ]:
            ssub = sub[mask]
            if len(ssub) < 10: continue
            rs, ps = spearmanr(ssub["qstar_pred"], ssub["wild_AF_pol"])
            rp, pp = pearsonr(ssub["qstar_pred"], ssub["wild_AF_pol"])
            slope_, intc = np.polyfit(ssub["qstar_pred"],
                                         ssub["wild_AF_pol"], 1)
            rows_results.append({"stringency": label, "subset": set_label,
                                  "n": len(ssub),
                                  "spearman": float(rs), "spearman_p": float(ps),
                                  "pearson": float(rp), "pearson_p": float(pp),
                                  "slope": float(slope_)})

    R = pd.DataFrame(rows_results)
    R.to_csv(f"{OUT_BASE}_sweep.tsv", sep="\t", index=False)


    fig = plt.figure(figsize=(12.5, 5.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.2], wspace=0.30)


    ax = fig.add_subplot(gs[0, 0])
    set_styles = {
        "all genome-wide": ("o", "#1a1a1a", "All GW"),
        "chr_439 LD block": ("s", "#C84A45", "chr_439 LD"),
        "rest of genome": ("^", "#3d6cb5", "rest of genome"),
    }
    for set_key, (marker, color, label) in set_styles.items():
        sub = R[R["subset"].str.contains(set_key)]
        ax.plot(range(len(sub_thresholds)), sub["spearman"].values,
                  marker=marker, color=color, linewidth=1.4, markersize=8,
                  markeredgecolor="black", markeredgewidth=0.5,
                  label=label)
        for i, (_, r) in enumerate(sub.iterrows()):
            if r["spearman_p"] < 0.001:
                ax.text(i, r["spearman"] + 0.01, "***",
                          ha="center", fontsize=9, color=color,
                          fontweight="bold")
            elif r["spearman_p"] < 0.01:
                ax.text(i, r["spearman"] + 0.01, "**",
                          ha="center", fontsize=9, color=color,
                          fontweight="bold")
            elif r["spearman_p"] < 0.05:
                ax.text(i, r["spearman"] + 0.01, "*",
                          ha="center", fontsize=9, color=color,
                          fontweight="bold")
    ax.axhline(0, color="#888", linestyle=":", linewidth=0.5)
    ax.set_xticks(range(len(sub_thresholds)))
    ax.set_xticklabels([t[0] for t in sub_thresholds],
                          fontsize=8, rotation=15, ha="right")
    ax.set_ylabel("Spearman ρ\n(q*_predicted vs wild_AF_polarized)",
                   fontsize=9)
    ax.set_title("Per-SNP correlation vs LRT stringency", fontsize=10)
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    ax.tick_params(labelsize=8)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)

    ax = fig.add_subplot(gs[0, 1])
    thr = df_antag["lrt"].quantile((0.99 - 0.95) / (1 - 0.95))
    sub = df_antag[df_antag["lrt"] >= thr]
    rng = np.random.default_rng(42)

    if len(sub) > 5000:
        idx_show = rng.choice(len(sub), size=5000, replace=False)
        sub_show = sub.iloc[idx_show]
    else:
        sub_show = sub

    rest = sub_show[~sub_show["in_chr439_LD"]]
    chr439 = sub_show[sub_show["in_chr439_LD"]]
    ax.scatter(rest["qstar_pred"], rest["wild_AF_pol"],
                s=4, color="#bbbbbb", alpha=0.4, edgecolor="none",
                label=f"Rest of genome (n={len(rest)})")
    ax.scatter(chr439["qstar_pred"], chr439["wild_AF_pol"],
                s=18, color="#C84A45", alpha=0.85,
                edgecolor="black", linewidth=0.4,
                label=f"chr_439 LD (n={len(chr439)})")

    ax.plot([0, 1], [0, 1], color="#999", linestyle=":", linewidth=0.7,
              label="y = x")

    s_all, i_all = np.polyfit(sub["qstar_pred"], sub["wild_AF_pol"], 1)
    rs_all, ps_all = spearmanr(sub["qstar_pred"], sub["wild_AF_pol"])
    x_fit = np.linspace(0, 1, 100)
    ax.plot(x_fit, s_all * x_fit + i_all,
              color="black", linewidth=1.4,
              label=f"All top-1% slope={s_all:+.2f}")
    ax.text(0.03, 0.97,
              f"Top 1% LRT genome-wide (n={len(sub):,})\n"
              f"Spearman ρ = {rs_all:+.3f}\n"
              f"p = {ps_all:.2e}",
              transform=ax.transAxes, ha="left", va="top",
              fontsize=8.5, family="monospace",
              bbox=dict(facecolor="white", alpha=0.9, edgecolor="#bbb",
                          linewidth=0.5, boxstyle="round,pad=0.3"))
    ax.axhline(0.5, color="#888", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.axvline(0.5, color="#888", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("q*_predicted from lab Dempster\n(per-SNP, additive h=0.5)",
                   fontsize=9)
    ax.set_ylabel("Wild allele frequency (polarized to lab B-favored)",
                   fontsize=9)
    ax.set_title("Top 1% LRT genome-wide: q*_pred vs wild AF",
                  fontsize=10)
    ax.legend(loc="lower right", fontsize=7.5, frameon=False)
    ax.tick_params(labelsize=8)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)

    fig.suptitle("Genome-wide direct test: per-SNP lab q*_predicted vs wild AF",
                  fontsize=11.5)
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.13, top=0.86,
                          wspace=0.25)
    fig.savefig(f"{OUT_BASE}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.svg", bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
