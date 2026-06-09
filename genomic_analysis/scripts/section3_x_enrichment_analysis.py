#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, fisher_exact

ROOT = Path(".")
QSTAR_TSV = ROOT / "final_plots/wild/section2_wild_qstar_correlation_gw_top50_perrep.tsv"
AF_MATRIX = ROOT / "final_plots/wild/af_matrix_22pools.csv"
DIVERSITY_TSV = ROOT / "grenfst/diversity/trajectory_pi_200000_nonoverlap.csv"
SCAFF_MULLER = ROOT / "final_plots/wild/sfla_v2_scaffold_muller.tsv"
SNP_MULLER = ROOT / "final_plots/wild/sfla_v2_snp_muller.tsv"

OUT_PER_MULLER = ROOT / "final_plots/wild/section3_x_enrichment_per_muller.tsv"
OUT_CONTRAST = ROOT / "final_plots/wild/section3_x_enrichment_contrast.tsv"

WILD_B = ["AVB", "PSB", "RMB"]
WILD_T = ["AVT", "PST", "RMT"]
LAB_Z = 2.0
WITHIN_HOST_SD_MAX = 0.07
M_VAR_QUANTILE = 0.75

N_BOOT = 2000
N_PERM = 5000
SEED = 42


def beta_MOM(x: np.ndarray, y: np.ndarray, var_eps_x: float) -> float:
    x, y = np.asarray(x), np.asarray(y)
    if len(x) < 3:
        return np.nan
    cov_xy = np.cov(x, y, ddof=1)[0, 1]
    var_x_obs = np.var(x, ddof=1)
    var_true = var_x_obs - var_eps_x
    if var_true <= 0:
        return np.nan
    return float(cov_xy / var_true)


def boot_rho_ci(x: np.ndarray, y: np.ndarray, n_boot: int,
                rng: np.random.Generator) -> tuple[float, float]:
    n = len(x)
    if n < 5:
        return (np.nan, np.nan)
    rhos = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        rhos[i] = spearmanr(x[idx], y[idx]).correlation
    return float(np.nanpercentile(rhos, 2.5)), float(np.nanpercentile(rhos, 97.5))


def load_q_pred_with_filters() -> pd.DataFrame:
    df = pd.read_csv(QSTAR_TSV, sep="\t")
    df = df[df["antag"]
            & df["qstar_pred"].notna()
            & df["wild_AF_pol"].notna()
            & df["sB_SE"].notna()
            & df["sT_SE"].notna()].copy()

    af = pd.read_csv(AF_MATRIX, usecols=["chrom_pos"] + WILD_B + WILD_T)
    af[["chrom", "pos"]] = af["chrom_pos"].str.split(":", n=1, expand=True)
    af["pos"] = af["pos"].astype(int)
    df = df.merge(af, on=["chrom", "pos"], how="left")
    Bp = df[WILD_B].to_numpy(dtype=float)
    Tp = df[WILD_T].to_numpy(dtype=float)
    df["within_host_SD"] = 0.5 * (np.nanstd(Bp, axis=1, ddof=1)
                                   + np.nanstd(Tp, axis=1, ddof=1))

    div = pd.read_csv(DIVERSITY_TSV)
    div["start"] = div["start"].astype(int)
    div["end"] = div["end"].astype(int)
    for arm in ["B", "T", "M"]:
        rep_cols = [f"{arm}{r}G10.theta_watterson" for r in range(1, 5)
                    if f"{arm}{r}G10.theta_watterson" in div.columns]
        div[f"{arm}_SD_G10"] = div[rep_cols].std(axis=1, ddof=1)
    div = div.dropna(subset=["B_SD_G10", "T_SD_G10", "M_SD_G10"]).copy()
    div["M_var_excess"] = (div["M_SD_G10"]
                            - 0.5 * (div["B_SD_G10"] + div["T_SD_G10"]))
    snp_to_window: dict[tuple[str, int], int] = {}
    for chrom, sub_div in div.groupby("chrom"):
        sub = sub_div.sort_values("start").reset_index(drop=True)
        starts = sub["start"].values
        ends = sub["end"].values
        snps_chrom = df[df["chrom"] == chrom]["pos"].values
        for pos in snps_chrom:
            idx = np.searchsorted(starts, pos, side="right") - 1
            if 0 <= idx < len(starts) and starts[idx] <= pos <= ends[idx]:
                snp_to_window[(chrom, int(pos))] = int(starts[idx])
    div_lookup = {(r["chrom"], int(r["start"])): float(r["M_var_excess"])
                  for _, r in div.iterrows()}
    df["div_window_start"] = [snp_to_window.get((c, int(p)), -1)
                               for c, p in zip(df["chrom"], df["pos"])]
    df["M_var_excess"] = [div_lookup.get((c, int(b)), np.nan)
                          for c, b in zip(df["chrom"], df["div_window_start"])]
    thrM = div["M_var_excess"].quantile(M_VAR_QUANTILE)

    df["zB"] = df["sB"].abs() / df["sB_SE"].clip(lower=1e-9)
    df["zT"] = df["sT"].abs() / df["sT_SE"].clip(lower=1e-9)
    INPUT_TOP_PCT = 0.50
    q_local_5 = max(0.0, (0.95 - INPUT_TOP_PCT) / (1 - INPUT_TOP_PCT))
    thr5 = df["lrt"].quantile(q_local_5)
    df["bg_pass"] = (df["lrt"] >= thr5) & (df["zB"] >= LAB_Z) & (df["zT"] >= LAB_Z)
    df["fg_pass"] = (df["bg_pass"]
                     & (df["within_host_SD"] <= WITHIN_HOST_SD_MAX)
                     & (df["M_var_excess"] >= thrM))


    return df, thrM


def per_muller_table(df: pd.DataFrame, set_col: str,
                     rng: np.random.Generator) -> pd.DataFrame:
    sub = df[df[set_col] & df["muller_assigned"].notna()].copy()

    rows = []
    for muller, grp in sub.groupby("muller_assigned"):
        x = grp["qstar_pred"].to_numpy(dtype=float)
        y = grp["wild_AF_pol"].to_numpy(dtype=float)
        if len(x) < 5:
            rows.append({"muller": muller, "n": len(x),
                         "rho": np.nan, "rho_p": np.nan,
                         "ci_lo": np.nan, "ci_hi": np.nan,
                         "beta_eiv": np.nan})
            continue
        r = spearmanr(x, y)
        ci_lo, ci_hi = boot_rho_ci(x, y, N_BOOT, rng)
        var_eps_x = np.nanmean(grp["qstar_pred_SE"].to_numpy(dtype=float) ** 2)
        beta_eiv = beta_MOM(x, y, var_eps_x)
        rows.append({"muller": muller, "n": len(x),
                     "rho": float(r.correlation), "rho_p": float(r.pvalue),
                     "ci_lo": ci_lo, "ci_hi": ci_hi,
                     "beta_eiv": beta_eiv})
    return pd.DataFrame(rows).sort_values("muller").reset_index(drop=True)


def x_vs_autosomes_perm(df: pd.DataFrame, set_col: str, n_perm: int,
                         rng: np.random.Generator) -> dict:
    AUTOSOMES = {"B", "C", "D", "E", "F"}
    sub = df[df[set_col] & df["muller_assigned"].isin(
        {"A", "B", "C", "D", "E", "F"})].copy()
    if len(sub) < 20:
        return {"rho_A": np.nan, "rho_auto": np.nan,
                "delta_obs": np.nan, "p_perm": np.nan}

    def stat(s: pd.DataFrame) -> tuple[float, float, float]:
        a = s[s["muller_assigned"] == "A"]
        au = s[s["muller_assigned"].isin(AUTOSOMES)]
        if len(a) < 5 or len(au) < 5:
            return (np.nan, np.nan, np.nan)
        rA = spearmanr(a["qstar_pred"], a["wild_AF_pol"]).correlation
        rAu = spearmanr(au["qstar_pred"], au["wild_AF_pol"]).correlation
        return (float(rA), float(rAu), float(rA - rAu))

    rho_A, rho_auto, delta_obs = stat(sub)
    n_A_obs = int((sub["muller_assigned"] == "A").sum())

    qstar = sub["qstar_pred"].to_numpy()
    waf = sub["wild_AF_pol"].to_numpy()
    n = len(sub)
    deltas = np.empty(n_perm)
    for i in range(n_perm):
        idx_A = rng.choice(n, size=n_A_obs, replace=False)
        is_A = np.zeros(n, dtype=bool)
        is_A[idx_A] = True
        rA = spearmanr(qstar[is_A], waf[is_A]).correlation
        rAu = spearmanr(qstar[~is_A], waf[~is_A]).correlation
        deltas[i] = rA - rAu
    p_perm = float(((deltas >= delta_obs) & np.isfinite(deltas)).sum()
                    / np.isfinite(deltas).sum())
    return {"rho_A": rho_A, "rho_auto": rho_auto,
            "delta_obs": delta_obs, "p_perm": p_perm,
            "n_A": n_A_obs,
            "n_auto": int(sub["muller_assigned"].isin(AUTOSOMES).sum())}


def compositional_test(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    bg = df[df["bg_pass"] & df["muller_assigned"].notna()]
    fg = df[df["fg_pass"] & df["muller_assigned"].notna()]
    for muller in ["A", "B", "C", "D", "E", "F"]:
        a = int((fg["muller_assigned"] == muller).sum())
        b = int(len(fg) - a)
        c = int((bg["muller_assigned"] == muller).sum())
        d = int(len(bg) - c)
        if a + b == 0 or c + d == 0:
            continue
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({"muller": muller,
                     "fg_in": a, "fg_out": b, "fg_pct": a / (a + b),
                     "bg_in": c, "bg_out": d, "bg_pct": c / (c + d),
                     "odds_ratio": float(odds), "fisher_p": float(p)})
    return pd.DataFrame(rows)


def main():
    rng = np.random.default_rng(SEED)
    if not SNP_MULLER.exists():
        print(f"ERROR: {SNP_MULLER} not built yet. Run "
              f"`build_snp_muller_assignment.py` after the diamond "
              f"blastp completes on the cluster.", file=sys.stderr)
        sys.exit(2)

    df, thrM = load_q_pred_with_filters()
    snp_muller = pd.read_csv(SNP_MULLER, sep="\t")
    snp_muller = snp_muller.rename(columns={"muller": "muller_assigned"})
    df = df.merge(snp_muller[["chrom", "pos", "muller_assigned",
                              "dist_to_gene"]],
                  on=["chrom", "pos"], how="left")
    n_assigned = df["muller_assigned"].notna().sum()
    bg_table = per_muller_table(df, "bg_pass", rng)
    bg_table["set"] = "background"
    fg_table = per_muller_table(df, "fg_pass", rng)
    fg_table["set"] = "orange"
    pm = pd.concat([bg_table, fg_table], ignore_index=True)
    pm = pm[["set", "muller", "n", "rho", "rho_p", "ci_lo", "ci_hi", "beta_eiv"]]
    pm.to_csv(OUT_PER_MULLER, sep="\t", index=False)
    contrast_rows = []
    for set_col, label in [("bg_pass", "background"), ("fg_pass", "orange")]:
        res = x_vs_autosomes_perm(df, set_col, N_PERM, rng)
        res["set"] = label
        contrast_rows.append(res)
        print(f"  {label}: rho_A={res['rho_A']:.3f} (n={res.get('n_A','?')}) "
              f"rho_autosomes={res['rho_auto']:.3f} (n={res.get('n_auto','?')}) "
              f"delta={res['delta_obs']:+.3f}  p_perm={res['p_perm']:.4f}")
    pd.DataFrame(contrast_rows).to_csv(OUT_CONTRAST, sep="\t", index=False)

    comp = compositional_test(df)
    comp.to_csv(OUT_CONTRAST.with_name(
        "section3_x_enrichment_compositional.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
