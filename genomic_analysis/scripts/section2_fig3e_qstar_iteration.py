#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(".")
AD_TSV = ROOT / "variance_analysis/merged_ad.tsv"
SAMPLES = ROOT / "variance_analysis/sample_list.txt"
GLM_CSV = ROOT / "glm_lrt_gw_final/glmV1full.csv"
OUT_DIR = ROOT / "final_plots/wild"

CHROM = "chr_ScDA7r2_439_HRSCAF_779"
LD_S, LD_E = 2_640_000, 3_610_000

W_AA_B, W_Aa_B, W_aa_B = 1.000, 0.913, 0.612
W_AA_T, W_Aa_T, W_aa_T = 0.611, 0.806, 0.902
sB_P = (W_AA_B - W_aa_B) / W_AA_B
hB_P = (W_Aa_B - W_aa_B) / (W_AA_B - W_aa_B)
sT_P = (W_aa_T - W_AA_T) / W_aa_T
hT_P = (W_Aa_T - W_AA_T) / (W_aa_T - W_AA_T)
qstar_pheno = (0.5 * sB_P * hB_P) / (0.5 * sB_P * hB_P + 0.5 * sT_P * hT_P)

LRT_TOP_PCT = 0.95
MIN_DEPTH = 30
MIN_AF, MAX_AF = 0.05, 0.95
MIN_GENS = 2
G_MAX_VALUES = [2, 3, 4, 5, 6, 7, 8, 9, 10]


def parse_samples():
    names = open(SAMPLES).read().splitlines()
    cols = {n: i + 4 for i, n in enumerate(names)}
    return cols


def af_at(f, idx_list):
    vals = []
    for ci in idx_list:
        ad = f[ci]
        if "," not in ad: continue
        try: r, a = (int(x) for x in ad.split(","))
        except: continue
        if r + a < MIN_DEPTH: continue
        vals.append(a / (r + a))
    return np.mean(vals) if len(vals) >= 2 else np.nan


def logit(p, eps=1e-3):
    p = np.clip(p, eps, 1-eps)
    return np.log(p / (1 - p))


def slope_from_traj(af_per_gen, gens, min_gens=MIN_GENS):
    if len(af_per_gen) < min_gens: return np.nan
    x = np.array(gens, dtype=float)
    y = logit(np.array(af_per_gen))
    if not np.all(np.isfinite(y)): return np.nan
    if np.std(x) == 0: return np.nan
    slope, *_ = stats.linregress(x, y)
    return slope


def collect_per_snp_traj(cols):
    glm = pd.read_csv(GLM_CSV)
    glm = glm[(glm["chrom"] == CHROM) & (glm["pos"] >= LD_S) & (glm["pos"] <= LD_E)]
    cutoff = glm["LRT_chisq"].quantile(LRT_TOP_PCT)
    top_set = set(glm.loc[glm["LRT_chisq"] >= cutoff, "pos"].astype(int))

    F_idx = [cols[f"F{r}G00"] for r in range(1, 5) if f"F{r}G00" in cols]
    B_per_gen = {g: [cols[f"B{r}G{g:02d}"] for r in range(1, 5)
                      if f"B{r}G{g:02d}" in cols] for g in range(1, 11)}
    T_per_gen = {g: [cols[f"T{r}G{g:02d}"] for r in range(1, 5)
                      if f"T{r}G{g:02d}" in cols] for g in range(1, 11)}
    M_per_gen = {g: [cols[f"M{r}G{g:02d}"] for r in range(1, 5)
                      if f"M{r}G{g:02d}" in cols] for g in range(1, 11)}

    snps = {}
    n_seen = 0
    with open(AD_TSV) as fh:
        for line in fh:
            f = line.rstrip().split("\t")
            chrom, pos = f[0], int(f[1])
            if chrom != CHROM: continue
            if not (LD_S <= pos <= LD_E): continue

            af_F = af_at(f, F_idx)
            if np.isnan(af_F) or not (MIN_AF <= af_F <= MAX_AF): continue

            B_traj, T_traj, M_traj = {}, {}, {}
            for g in range(1, 11):
                v = af_at(f, B_per_gen[g])
                if not np.isnan(v): B_traj[g] = v
                v = af_at(f, T_per_gen[g])
                if not np.isnan(v): T_traj[g] = v
                v = af_at(f, M_per_gen[g])
                if not np.isnan(v): M_traj[g] = v
            if len(B_traj) < MIN_GENS or len(T_traj) < MIN_GENS:
                continue

            snps[pos] = {
                "founder_af": af_F,
                "B_traj": B_traj, "T_traj": T_traj, "M_traj": M_traj,
                "in_top": pos in top_set,
                "lrt": float(glm.loc[glm["pos"] == pos, "LRT_chisq"].iloc[0])
                       if pos in glm["pos"].values else np.nan,
            }
            n_seen += 1
    return snps


def slopes_at_cutoff(snps, g_max):
    rows = []
    for pos, d in snps.items():
        B_g = sorted(g for g in d["B_traj"] if g <= g_max)
        T_g = sorted(g for g in d["T_traj"] if g <= g_max)
        M_g = sorted(g for g in d.get("M_traj", {}) if g <= g_max)
        B_af = [d["B_traj"][g] for g in B_g]
        T_af = [d["T_traj"][g] for g in T_g]
        M_af = [d["M_traj"][g] for g in M_g] if M_g else []
        sB = slope_from_traj(B_af, B_g)
        sT = slope_from_traj(T_af, T_g)
        sM = (slope_from_traj(M_af, M_g) if len(M_g) >= MIN_GENS else np.nan)
        if not (np.isfinite(sB) and np.isfinite(sT)):
            continue
        # Polarize so sB>0 (apply same flip to sT, sM, and AF values)
        polarized = False
        if sB < 0:
            sB, sT = -sB, -sT
            if np.isfinite(sM): sM = -sM
            polarized = True
        founder_af = d["founder_af"] if not polarized else (1 - d["founder_af"])
        # M_avg = mean polarized AF over the trajectory
        if M_af:
            M_af_pol = [(1 - x) if polarized else x for x in M_af]
            M_avg = float(np.mean(M_af_pol))
        else:
            M_avg = np.nan
        rows.append({"pos": pos, "founder_af": founder_af,
                     "sB": sB, "sT": sT, "sM": sM, "M_avg": M_avg,
                     "in_top": d["in_top"], "lrt": d["lrt"]})
    return pd.DataFrame(rows)


def qstar(sB, sT, h=0.5, c=0.5):
    if sB <= 0 or sT >= 0:
        return np.nan
    return (c * sB * h) / (c * sB * h + (1 - c) * abs(sT) * h)


def bootstrap_qstar(top_df, n_boot=1000, seed=0):
    rng = np.random.default_rng(seed)
    n = len(top_df)
    sB_v = top_df["sB"].values
    sT_v = top_df["sT"].values
    out = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        out[i] = qstar(np.median(sB_v[idx]), np.median(sT_v[idx]))
    return out


def main():
    plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42,
                          "font.family": "sans-serif",
                          "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
                          "axes.linewidth": 0.6})

    cols = parse_samples()
    snps = collect_per_snp_traj(cols)
    n_top_def = sum(1 for d in snps.values() if d["in_top"])
    res = []
    for g_max in G_MAX_VALUES:
        df = slopes_at_cutoff(snps, g_max)
        top = df[df["in_top"]]
        if len(top) < 20:
            continue
        med_sB = float(top["sB"].median())
        med_sT = float(top["sT"].median())
        med_sM = float(top["sM"].dropna().median()) if top["sM"].notna().any() else np.nan
        n_M = int(top["sM"].notna().sum())
        sM_pred_add = 0.5 * (med_sB + med_sT)
        dom_ratio = (med_sM / sM_pred_add
                       if (np.isfinite(med_sM) and abs(sM_pred_add) > 1e-6)
                       else np.nan)

        q_pt = qstar(med_sB, med_sT)
        boot = bootstrap_qstar(top, n_boot=1500, seed=g_max)
        valid = ~np.isnan(boot)
        q_lo = float(np.nanquantile(boot, 0.025))
        q_hi = float(np.nanquantile(boot, 0.975))
        res.append({"g_max": g_max, "n_top": len(top), "n_M_snp": n_M,
                     "med_sB": med_sB, "med_sT": med_sT, "med_sM": med_sM,
                     "sM_pred_add": sM_pred_add, "dom_ratio": dom_ratio,
                     "qstar": q_pt, "q_lo": q_lo, "q_hi": q_hi,
                     "valid_frac": valid.mean()})

    R = pd.DataFrame(res)
    R.to_csv(OUT_DIR / "section2_fig3e_qstar_iteration.tsv", sep="\t", index=False)

    g10_df = slopes_at_cutoff(snps, 10)
    g10_df.to_csv(OUT_DIR / "section2_fig3e_per_snp_g10_full.tsv",
                    sep="\t", index=False)

    fig = plt.figure(figsize=(9.5, 4.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.3], wspace=0.30)

    ax = fig.add_subplot(gs[0, 0])
    ax.axhline(qstar_pheno, color="#a04848", linestyle="-", linewidth=1.2,
                alpha=0.85, label=f"q*_pheno = {qstar_pheno:.3f}")
    ax.axhline(0.5, color="#888", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.errorbar(R["g_max"], R["qstar"],
                  yerr=[R["qstar"] - R["q_lo"], R["q_hi"] - R["qstar"]],
                  fmt="o-", color="#1a5e1a", markersize=7, capsize=3,
                  linewidth=1.4, label="q*_geno (top-LRT, h=0.5)")
    for _, r in R.iterrows():
        ax.text(r["g_max"], r["qstar"] + 0.05, f"{r['qstar']:.3f}",
                 ha="center", fontsize=7, color="#1a5e1a")
    ax.set_xlabel("included generations (G1 → G_max)", fontsize=8)
    ax.set_ylabel("q*", fontsize=8)
    ax.set_xticks(G_MAX_VALUES)
    ax.set_xticklabels([f"G1–G{g}" for g in G_MAX_VALUES], fontsize=7, rotation=20)
    ax.set_ylim(0, 1)
    ax.set_title("q*_geno convergence as more generations are added",
                  fontsize=10)
    ax.tick_params(labelsize=7)
    ax.legend(loc="lower right", fontsize=7, frameon=False)

    ax = fig.add_subplot(gs[0, 1])
    ax.axhline(0, color="#888", linewidth=0.5)
    ax.plot(R["g_max"], R["med_sB"], "o-", color="#a04848",
             markersize=6, linewidth=1.4, label="median s_B (B host)")
    ax.plot(R["g_max"], R["med_sT"], "s-", color="#3d6cb5",
             markersize=6, linewidth=1.4, label="median s_T (T host)")

    ax.set_xlabel("included generations (G1 → G_max)", fontsize=8)
    ax.set_ylabel("median per-SNP s (top-LRT)", fontsize=8)
    ax.set_xticks(G_MAX_VALUES)
    ax.set_xticklabels([f"G1–G{g}" for g in G_MAX_VALUES], fontsize=7, rotation=20)
    ax.set_title("Aggregate per-SNP selection coefficients",
                  fontsize=10)
    ax.tick_params(labelsize=7)
    ax.legend(loc="upper right", fontsize=7, frameon=False)

    for _, r in R.iterrows():
        ratio = r["med_sB"] / (r["med_sB"] + abs(r["med_sT"]))
        ax.text(r["g_max"], r["med_sB"] + 0.005,
                  f"q*={ratio:.2f}", ha="center", fontsize=6, color="#1a5e1a")

    fig.suptitle("Sensitivity of q*_geno to the time window of allele-frequency dynamics",
                  fontsize=11)
    fig.tight_layout()
    save(fig, "section2_fig3e_qstar_iteration")


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
