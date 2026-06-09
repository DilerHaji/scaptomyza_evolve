#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as stats

ROOT = Path(".")
MASTER_TSV = ROOT / "final_plots/wild/section2_candidate_master_v2.tsv"
HV_BLOCKS_TSV = ROOT / "final_plots/wild/section2_hv_blocks_filtered.tsv"
GLMM_PER_WINDOW_TSV = ROOT / "final_plots/wild/section2_glmm_per_window.tsv"
OUT_P_TSV = ROOT / "final_plots/wild/section2_per_test_perm_p.tsv"
OUT_THR_TSV = ROOT / "final_plots/wild/section2_per_test_perm_thresholds.tsv"

WIN = 200_000
SCAFF_MIN_WINDOWS = 50
N_PERM = 5000
ALPHA_LEVELS = (0.05, 0.01)

SEGMENT_TRACKS = [
    ("hv_B", "B"),
    ("hv_T", "T"),
    ("hv_M", "M"),
]
STRIPE_TRACKS = [
    ("cvtk_antag_n",  "high"),
    ("cov_BT",        "low"),
    ("slope_div_z",   "high"),
    ("permFST_emp_p", "p"),
    ("glmm_lrt",      "high"),
    ("wild_C2_max",   "high"),
]


def hv_score_per_window(df_p, blocks, treatment_letter):
    blk = blocks[blocks["treatment"] == treatment_letter]
    s = np.zeros(len(df_p), dtype=float)
    if blk.empty:
        return s
    by_chr = {c: g for c, g in blk.groupby("chr")}
    for i, (_, r) in enumerate(df_p.iterrows()):
        g = by_chr.get(r["chrom"])
        if g is None:
            continue
        m = (g["start"] <= r["end"]) & (g["end"] >= r["start"])
        if m.any():
            s[i] = float(g.loc[m, "n_snps"].max())
    return s


def signed_score(values: np.ndarray, direction: str) -> np.ndarray:
    s = np.asarray(values, dtype=float)
    if direction == "low":
        return -s
    if direction == "p":
        return stats.norm.ppf(1 - np.clip(s, 1e-10, 1 - 1e-10))
    return s


def per_test_perm_p(scores: np.ndarray, df_p: pd.DataFrame, n_perm: int,
                      seed: int = 42):
    rng = np.random.default_rng(seed)
    df_p = df_p.reset_index(drop=True)
    s = np.nan_to_num(scores, nan=0.0)
    n = len(s)

    chrom_groups = []
    for chrom, g in df_p.groupby("chrom", sort=False):
        chrom_groups.append((chrom, g.index.to_numpy()))

    null_max = np.empty(n_perm)

    perwin_hits = np.zeros(n, dtype=np.int64)

    for i in range(n_perm):
        s_perm = np.empty_like(s)
        for chrom, idx in chrom_groups:
            offset = int(rng.integers(0, len(idx)))
            s_perm[idx] = np.roll(s[idx], offset)
        null_max[i] = s_perm.max()
        perwin_hits += (s_perm >= s).astype(np.int64)

    pwise = (perwin_hits + 1) / (n_perm + 1)
    thresholds = {a: float(np.quantile(null_max, 1 - a)) for a in ALPHA_LEVELS}
    return pwise, thresholds


def main() -> None:
    df = pd.read_csv(MASTER_TSV, sep="\t")
    blocks = pd.read_csv(HV_BLOCKS_TSV, sep="\t")
    glmm = pd.read_csv(GLMM_PER_WINDOW_TSV, sep="\t")
    df = df.merge(glmm[["chrom", "start", "glmm_lrt"]],
                   on=["chrom", "start"], how="left")

    chrom_sizes = df.groupby("chrom").size().sort_values(ascending=False)
    keep = chrom_sizes[chrom_sizes >= SCAFF_MIN_WINDOWS].index.tolist()
    df_p = df[df["chrom"].isin(keep)].copy().reset_index(drop=True)

    track_scores = {}
    for key, treat in SEGMENT_TRACKS:
        track_scores[key] = signed_score(
            hv_score_per_window(df_p, blocks, treat), "high"
        )
    for col, dirn in STRIPE_TRACKS:
        track_scores[col] = signed_score(df_p[col].values, dirn)

    out_p = pd.DataFrame({"chrom": df_p["chrom"], "start": df_p["start"]})
    threshold_rows = []
    for name, s in track_scores.items():
        pwise, thr = per_test_perm_p(s, df_p, N_PERM)
        out_p[f"p_{name}"] = pwise
        threshold_rows.append({
            "track": name,
            "threshold_p05_score": thr[0.05],
            "threshold_p01_score": thr[0.01],
        })

    out_p.to_csv(OUT_P_TSV, sep="\t", index=False)
    pd.DataFrame(threshold_rows).to_csv(OUT_THR_TSV, sep="\t", index=False)


if __name__ == "__main__":
    main()
