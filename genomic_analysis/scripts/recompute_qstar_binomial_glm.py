#!/usr/bin/env python3
from __future__ import annotations
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

ROOT = Path(".")
GLM_CSV = ROOT / "glm_lrt_gw_final/glmV1full.csv"
AD_TSV = ROOT / "variance_analysis/merged_ad.tsv"
SAMPLES = ROOT / "variance_analysis/sample_list.txt"
OUT = ROOT / "final_plots/wild/qstar_binomial_glm.tsv"

GENS = list(range(1, 11))
MIN_DEPTH = 30
TOP_PCT = 0.50

def parse_samples():
    return {n: i + 4
            for i, n in enumerate(open(SAMPLES).read().splitlines())}


def parse_ad(s):
    if "," not in s:
        return None, None
    try:
        r, a = (int(x) for x in s.split(","))
    except (ValueError, AttributeError):
        return None, None
    return r, a


def fit_binomial_glm(alt, total, gen_arr, treat_arr, rep_arr):
    if len(alt) < 8:
        return None

    n = len(alt)
    X_cols = [np.ones(n), gen_arr, treat_arr, gen_arr * treat_arr]

    unique_reps = np.unique(rep_arr)
    if len(unique_reps) > 1:
        for r in unique_reps[1:]:
            X_cols.append((rep_arr == r).astype(float))

    X = np.column_stack(X_cols)

    ref = total - alt
    y = np.column_stack([alt, ref])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            model = sm.GLM(y, X, family=sm.families.Binomial())
            res = model.fit(maxiter=50, atol=1e-6, disp=0)
        except Exception:
            return None

    if not res.converged:
        return None

    beta = res.params
    cov = res.cov_params()

    sB = beta[1]
    sT = beta[1] + beta[3]
    sB_var = cov[1, 1]
    sT_var = cov[1, 1] + cov[3, 3] + 2 * cov[1, 3]
    if sB_var <= 0 or sT_var <= 0:
        return None
    sB_SE = np.sqrt(sB_var)
    sT_SE = np.sqrt(sT_var)
    return float(sB), float(sT), float(sB_SE), float(sT_SE), True, n


def main():
    g = pd.read_csv(GLM_CSV, usecols=["chrom", "pos", "LRT_chisq",
                                         "converged", "error"])
    g = g[(g["converged"] == True) & (g["error"] == "OK")
          & g["LRT_chisq"].notna()]
    thr = g["LRT_chisq"].quantile(TOP_PCT)
    g_top = g[g["LRT_chisq"] >= thr][["chrom", "pos", "LRT_chisq"]].copy()
    g_top["pos"] = g_top["pos"].astype(int)
    target_set = set(zip(g_top["chrom"], g_top["pos"]))
    chrom_pos_set = {}
    for chrom, pos in target_set:
        chrom_pos_set.setdefault(chrom, set()).add(pos)
    lrt_lookup = dict(zip(zip(g_top["chrom"], g_top["pos"]),
                            g_top["LRT_chisq"]))

    cols = parse_samples()
    B_per_rep_gen = {(r, gn): cols[f"B{r}G{gn:02d}"]
                       for r in range(1, 5) for gn in GENS
                       if f"B{r}G{gn:02d}" in cols}
    T_per_rep_gen = {(r, gn): cols[f"T{r}G{gn:02d}"]
                       for r in range(1, 5) for gn in GENS
                       if f"T{r}G{gn:02d}" in cols}

    t0 = time.time()
    rows_done = 0
    fits = 0
    out_rows = []
    with open(AD_TSV) as fh:
        header = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            chrom = f[0]
            try:
                pos = int(f[1])
            except ValueError:
                continue
            rows_done += 1
            if rows_done % 50000 == 0:
                rate = fits / max(time.time() - t0, 1e-6)
                pct = fits / max(len(target_set), 1) * 100
            if pos not in chrom_pos_set.get(chrom, ()):
                continue

            alt_list, total_list, gen_list, treat_list, rep_list = ([] for _ in range(5))
            for (r, gn), ci in B_per_rep_gen.items():
                if ci >= len(f):
                    continue
                ref, alt = parse_ad(f[ci])
                if ref is None or ref + alt < MIN_DEPTH:
                    continue
                alt_list.append(alt)
                total_list.append(ref + alt)
                gen_list.append(gn)
                treat_list.append(0)
                rep_list.append(r)
            for (r, gn), ci in T_per_rep_gen.items():
                if ci >= len(f):
                    continue
                ref, alt = parse_ad(f[ci])
                if ref is None or ref + alt < MIN_DEPTH:
                    continue
                alt_list.append(alt)
                total_list.append(ref + alt)
                gen_list.append(gn)
                treat_list.append(1)
                rep_list.append(r + 4)  # offset to keep B/T reps distinct

            if len(alt_list) < 16:
                continue
            result = fit_binomial_glm(
                np.array(alt_list, dtype=float),
                np.array(total_list, dtype=float),
                np.array(gen_list, dtype=float),
                np.array(treat_list, dtype=float),
                np.array(rep_list))
            if result is None:
                continue
            sB, sT, sB_SE, sT_SE, conv, n_obs = result
            fits += 1

            if sB < 0:
                sB, sT = -sB, -sT
            antag = (sB > 0) and (sT < 0)
            if antag:
                qstar = sB / (sB + abs(sT))
                denom4 = (sB - sT) ** 4
                var_q = (sT ** 2 * sB_SE ** 2 + sB ** 2 * sT_SE ** 2) / denom4
                qstar_SE = float(np.sqrt(max(var_q, 0)))
            else:
                qstar = np.nan
                qstar_SE = np.nan

            out_rows.append({
                "chrom": chrom, "pos": pos,
                "lrt": float(lrt_lookup[(chrom, pos)]),
                "sB": sB, "sT": sT,
                "sB_SE": sB_SE, "sT_SE": sT_SE,
                "n_obs": n_obs, "antag": antag,
                "qstar_pred": qstar, "qstar_pred_SE": qstar_SE,
            })

    df = pd.DataFrame(out_rows)
    df.to_csv(OUT, sep="\t", index=False)

if __name__ == "__main__":
    main()
