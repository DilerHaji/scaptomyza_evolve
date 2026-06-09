#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(".")
DIV_CSV = ROOT / "grenfst/diversity_attrition/attrition_pi_390000diversity.csv"
OUT_S1  = ROOT / "final_plots/wild/s3_table_S3_1_diversity_summary.tsv"
OUT_S2  = ROOT / "final_plots/wild/s3_table_S3_2_wild_paired_wilcoxon.tsv"

METRICS = [("theta_pi", "θ_π"), ("theta_watterson", "θ_W"), ("tajimas_d", "Tajima's D")]

GROUPS = {
    "Wild":      ["AVB", "AVT", "PSB", "PST", "RMB", "RMT"],
    "Founder":   ["F1G00", "F2G00", "F3G00", "F4G00"],
    "G10 B":     ["B1G10", "B2G10", "B3G10", "B4G10"],
    "G10 T":     ["T1G10", "T2G10", "T3G10", "T4G10"],
    "G10 B+T":   ["M1G10", "M2G10", "M3G10", "M4G10"],
}

WILD_SITES = {"AV": ("AVB", "AVT"), "PS": ("PSB", "PST"), "RM": ("RMB", "RMT")}


def pool_series(raw: pd.DataFrame, pool: str, metric: str) -> pd.Series:
    for col in (f"{pool}.{metric}", f"{pool}.1.{metric}"):
        if col in raw.columns:
            return pd.to_numeric(raw[col], errors="coerce")
    raise KeyError(f"{pool}/{metric} not found in columns")


def summarise(raw: pd.DataFrame, group: str, pools: list[str], key: str,
              label: str) -> dict:
    per_pool_medians = [float(np.nanmedian(pool_series(raw, p, key).values))
                         for p in pools]
    pooled = np.concatenate([pool_series(raw, p, key).dropna().values for p in pools])
    return {
        "group": group,
        "metric": label,
        "n_pools": len(pools),
        "n_windows_per_pool": int(np.mean([
            pool_series(raw, p, key).notna().sum() for p in pools])),
        "per_pool_median_mean": round(float(np.mean(per_pool_medians)), 4),
        "per_pool_median_sd":   round(float(np.std(per_pool_medians, ddof=1)), 4),
        "per_pool_median_min":  round(float(np.min(per_pool_medians)), 4),
        "per_pool_median_max":  round(float(np.max(per_pool_medians)), 4),
        "pooled_q25":           round(float(np.nanpercentile(pooled, 25)), 4),
        "pooled_median":        round(float(np.nanmedian(pooled)), 4),
        "pooled_q75":           round(float(np.nanpercentile(pooled, 75)), 4),
    }


def main() -> None:
    raw = pd.read_csv(DIV_CSV)

    rows = []
    for group, pools in GROUPS.items():
        for key, label in METRICS:
            rows.append(summarise(raw, group, pools, key, label))
    s1 = pd.DataFrame(rows)
    OUT_S1.parent.mkdir(parents=True, exist_ok=True)
    s1.to_csv(OUT_S1, sep="\t", index=False)

    rows = []
    for site, (pool_B, pool_T) in WILD_SITES.items():
        for key, label in METRICS:
            v_B = pool_series(raw, pool_B, key)
            v_T = pool_series(raw, pool_T, key)
            mask = v_B.notna() & v_T.notna()
            b = v_B[mask].values
            t = v_T[mask].values
            med_B = float(np.median(b))
            med_T = float(np.median(t))

            median_diff = float(np.median(b - t))
            try:
                stat = stats.wilcoxon(b, t, alternative="two-sided",
                                       zero_method="wilcox")
                W, p = float(stat.statistic), float(stat.pvalue)
            except ValueError:
                W, p = float("nan"), float("nan")
            rows.append({
                "site": site,
                "metric": label,
                "n_paired_windows": int(mask.sum()),
                f"median_B": round(med_B, 4),
                f"median_T": round(med_T, 4),
                "median(B − T)": round(median_diff, 4),
                "wilcoxon_W": round(W, 2),
                "p_value": f"{p:.2e}",
            })
    s2 = pd.DataFrame(rows)
    s2.to_csv(OUT_S2, sep="\t", index=False)


if __name__ == "__main__":
    main()
