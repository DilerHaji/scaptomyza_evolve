#!/usr/bin/env python3
from __future__ import annotations
from itertools import combinations
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(".")
OUT_DIR = ROOT / "final_plots/wild"
OUT_TBL_MASTER = OUT_DIR / "section2_candidate_master_v2.tsv"
OUT_TBL_BYVOTES = OUT_DIR / "section2_candidate_byvotes_v2.tsv"
OUT_TBL_REGIONS = OUT_DIR / "section2_candidate_regions_v2.tsv"
OUT_BLOCKS_FILT = OUT_DIR / "section2_hv_blocks_filtered.tsv"

WIN = 200_000
TOP_FRAC = 0.05
MERGE_GAP = 200_000
MIN_WIN_NSNPS = 100    # Pre-filter: drop 200kb windows with < 100 SNPs (noisy)

MERGED_TSV = ROOT / "variance_analysis/ascertainment_validation_windows_merged.tsv"
CVTK_ANTAG = ROOT / "variance_analysis/cvtkpy_final/antagonistic_windows.tsv"
SLOPE_DIV = ROOT / ("grenfst/divergence_mixture_lmm_old/"
                     "e10Ffe9wGREN_btwTB_200000_1000/divergence_stats.csv")
PERM_FST = ROOT / "final_plots/wild/section2_lab_BvT_permFST_windows.tsv"
AF_MATRIX = ROOT / "final_plots/wild/af_matrix_22pools.csv"
TRAJ_DIV = ROOT / "grenfst/diversity/trajectory_pi_200000_nonoverlap.csv"
HV_ROOT = ROOT / "hv_results_cluster"

WILD_B_POOLS = ["AVB", "PSB", "RMB"]
WILD_T_POOLS = ["AVT", "PST", "RMT"]
B_G10 = ["B1G10", "B2G10", "B3G10", "B4G10"]
T_G10 = ["T1G10", "T2G10", "T3G10", "T4G10"]
REPS = [1, 2, 3, 4]


def compute_diversity_slope(traj: pd.DataFrame, trt: str, stat: str) -> pd.DataFrame:
    cols = [f"{trt}{r}G{g:02d}.{stat}"
            for r in REPS for g in range(1, 11)
            if f"{trt}{r}G{g:02d}.{stat}" in traj.columns]
    keep = ["chrom", "start", "end"] + cols
    sub = traj[keep].copy()
    long = sub.melt(id_vars=["chrom", "start", "end"],
                     var_name="pool_stat", value_name="theta_w")
    long = long.dropna(subset=["theta_w"])
    long["rep"] = long["pool_stat"].str.extract(r"(\d)G\d{2}")
    long["gen"] = long["pool_stat"].str.extract(r"G(\d{2})").astype(int)
    per_gen = (long.groupby(["chrom", "start", "end", "gen"])["theta_w"]
                    .mean().reset_index())
    def slope(df):
        if len(df) < 4:
            return np.nan
        x = df["gen"].values.astype(float)
        y = df["theta_w"].values.astype(float)
        return np.polyfit(x, y, 1)[0]
    short = {"theta_pi": "thetaPi", "theta_watterson": "thetaW",
             "tajimas_d": "tajD"}[stat]
    colname = f"{trt}_{short}_slope"
    slopes = (per_gen.groupby(["chrom", "start", "end"])
                     .apply(slope).reset_index(name=colname))
    slopes["bin_start"] = (slopes["start"] // WIN) * WIN
    binned = (slopes.groupby(["chrom", "bin_start"])[colname]
                    .mean().reset_index().rename(columns={"bin_start": "start"}))
    return binned


def compute_rep_convergence_FST(af: pd.DataFrame, trt_cols: list[str]) -> pd.DataFrame:
    n_eff = 58
    P = af[trt_cols].to_numpy(dtype=float)
    pairs = list(combinations(range(P.shape[1]), 2))
    total_num = np.zeros(P.shape[0])
    total_den = np.zeros(P.shape[0])
    for (i, j) in pairs:
        p1 = P[:, i]; p2 = P[:, j]
        num = (p1 - p2) ** 2 - (p1 * (1 - p1)) / (n_eff - 1) - (p2 * (1 - p2)) / (n_eff - 1)
        den = p1 * (1 - p2) + p2 * (1 - p1)
        total_num += np.nan_to_num(num)
        total_den += np.nan_to_num(den)
    af_out = pd.DataFrame({
        "chrom": af["chrom"].values,
        "start": (af["pos"] // WIN) * WIN,
        "pair_num_sum": total_num / len(pairs),
        "pair_den_sum": total_den / len(pairs),
    })
    grp = af_out.groupby(["chrom", "start"], as_index=False).agg(
        pair_num_sum=("pair_num_sum", "sum"),
        pair_den_sum=("pair_den_sum", "sum"),
    )
    grp["mean_within_FST"] = grp["pair_num_sum"] / np.maximum(grp["pair_den_sum"], 1e-9)
    return grp[["chrom", "start", "mean_within_FST"]]


HV_MIN_NSNPS = 1000   # high-information filter (top blocks by SNP support)


def load_hv_blocks_effective(trt: str) -> pd.DataFrame:
    fp = HV_ROOT / trt / "B_effective" / "dominant_blocks.tsv"
    df = pd.read_csv(fp, sep="\t")
    df = df[df["n_snps"] >= HV_MIN_NSNPS].copy()
    df["treatment"] = trt
    return df


def hv_blocks_center_in_window(blocks: pd.DataFrame, win_df: pd.DataFrame) -> np.ndarray:
    flag = np.zeros(len(win_df), dtype=int)
    blk = blocks.copy()
    blk["center"] = (blk["start"] + blk["end"]) / 2
    win_keyed = win_df.reset_index(drop=True).copy()
    for _, b in blk.iterrows():
        hit = win_keyed[(win_keyed["chrom"] == b["chr"]) &
                         (win_keyed["start"] <= b["center"]) &
                         (win_keyed["end"] > b["center"])]
        if not hit.empty:
            for idx in hit.index:
                flag[idx] = 1
    return flag


def load_and_join() -> pd.DataFrame:
    merged = pd.read_csv(MERGED_TSV, sep="\t")
    keep = ["chrom", "start", "end", "n_snps",
            "C2_max", "C2_mean",
            "hv_any_block", "hv_any_antag", "hv_mean_deltaAF_product",
            "cov_BT", "cov_BM", "cov_TM", "across_rep_var_M",
            "ratio_M_over_F"]
    keep = [c for c in keep if c in merged.columns]
    df = merged[keep].copy().rename(columns={"C2_max": "wild_C2_max"})

    antag = pd.read_csv(CVTK_ANTAG, sep="\t")
    antag["bin_start"] = (antag["start"] // WIN) * WIN
    antag_w = antag.groupby(["chrom", "bin_start"]).size().reset_index(
        name="cvtk_antag_n").rename(columns={"bin_start": "start"})
    df = df.merge(antag_w, on=["chrom", "start"], how="left").fillna({"cvtk_antag_n": 0})
    df["cvtk_antag_n"] = df["cvtk_antag_n"].astype(int)

    sd = pd.read_csv(SLOPE_DIV, usecols=["chrom", "start", "end",
                                           "slope_divergence", "z_divergence"])
    sd_thin = sd[sd["start"] % WIN == 1].copy()
    sd_thin["start"] = sd_thin["start"] - 1
    sd_thin = sd_thin.rename(columns={"z_divergence": "slope_div_z",
                                        "slope_divergence": "slope_div"})
    df = df.merge(sd_thin[["chrom", "start", "slope_div", "slope_div_z"]],
                   on=["chrom", "start"], how="left")


    pf = pd.read_csv(PERM_FST, sep="\t").rename(columns={"bin_start": "start"})
    df = df.merge(pf[["chrom", "start", "fst_BvT", "emp_p", "fst_ratio"]],
                   on=["chrom", "start"], how="left").rename(
                   columns={"fst_BvT": "permFST_BvT",
                             "emp_p": "permFST_emp_p",
                             "fst_ratio": "permFST_ratio"})


    af = pd.read_csv(AF_MATRIX, usecols=["chrom_pos"] + WILD_B_POOLS + WILD_T_POOLS
                                         + B_G10 + T_G10)
    af[["chrom", "pos"]] = af["chrom_pos"].str.split(":", expand=True)
    af["pos"] = af["pos"].astype(int)
    af["bin_start"] = (af["pos"] // WIN) * WIN

    sd_AV = (af["AVB"] - af["AVT"]).abs()
    sd_PS = (af["PSB"] - af["PST"]).abs()
    sd_RM = (af["RMB"] - af["RMT"]).abs()
    af["wild_per_snp_mean_abs"] = (sd_AV + sd_PS + sd_RM) / 3
    wild_agg = af.groupby(["chrom", "bin_start"]).agg(
        wild_AFdiff_mean=("wild_per_snp_mean_abs", "mean"),
    ).reset_index().rename(columns={"bin_start": "start"})
    df = df.merge(wild_agg, on=["chrom", "start"], how="left")

    B_conv = compute_rep_convergence_FST(af, B_G10).rename(
        columns={"mean_within_FST": "B_rep_FST"})
    T_conv = compute_rep_convergence_FST(af, T_G10).rename(
        columns={"mean_within_FST": "T_rep_FST"})
    df = df.merge(B_conv, on=["chrom", "start"], how="left")
    df = df.merge(T_conv, on=["chrom", "start"], how="left")

    traj = pd.read_csv(TRAJ_DIV)
    for trt in ("B", "T"):
        for stat in ("theta_pi", "theta_watterson", "tajimas_d"):
            slp = compute_diversity_slope(traj, trt, stat)
            df = df.merge(slp, on=["chrom", "start"], how="left")

    hv_B = load_hv_blocks_effective("B")
    hv_T = load_hv_blocks_effective("T")
    hv_M = load_hv_blocks_effective("M")
    df["IND_HV_blocks_B"] = hv_blocks_center_in_window(hv_B, df.reset_index(drop=True))
    df["IND_HV_blocks_T"] = hv_blocks_center_in_window(hv_T, df.reset_index(drop=True))
    df["IND_HV_blocks_M"] = hv_blocks_center_in_window(hv_M, df.reset_index(drop=True))

    pd.concat([hv_B, hv_T, hv_M], ignore_index=True).to_csv(
        OUT_BLOCKS_FILT, sep="\t", index=False)

    n_pre = len(df)
    df = df[df["n_snps"] >= MIN_WIN_NSNPS].reset_index(drop=True)


    return df


def add_indicators(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    out = df.copy()

    out["IND_cvtk_antag"] = (out["cvtk_antag_n"] >= 1).astype(int)
    thr = out["slope_div_z"].quantile(1 - TOP_FRAC)
    out["IND_slope_div"] = (out["slope_div_z"] >= thr).fillna(False).astype(int)
    out["IND_permFST"] = (out["permFST_emp_p"] <= TOP_FRAC).fillna(False).astype(int)

    thr = out["cov_BT"].quantile(TOP_FRAC)
    out["IND_cov_BT_neg"] = (out["cov_BT"] <= thr).fillna(False).astype(int)

    thr = out["wild_C2_max"].quantile(1 - TOP_FRAC)
    out["IND_wild_C2"] = (out["wild_C2_max"] >= thr).fillna(False).astype(int)

    ind_cols = [c for c in out.columns if c.startswith("IND_")]
    out["votes"] = out[ind_cols].sum(axis=1)

    return out, ind_cols


def merge_adjacent(df: pd.DataFrame, gap=MERGE_GAP) -> pd.DataFrame:
    df = df.sort_values(["chrom", "start"]).reset_index(drop=True)
    if df.empty:
        return df
    regions = []
    cur = {"chrom": df.iloc[0]["chrom"], "start": int(df.iloc[0]["start"]),
           "end": int(df.iloc[0]["end"]),
           "max_votes": int(df.iloc[0]["votes"]),
           "votes_sum": int(df.iloc[0]["votes"]), "n_windows": 1}
    for _, r in df.iloc[1:].iterrows():
        if (r["chrom"] == cur["chrom"] and r["start"] - cur["end"] <= gap):
            cur["end"] = int(r["end"])
            cur["max_votes"] = max(cur["max_votes"], int(r["votes"]))
            cur["votes_sum"] += int(r["votes"])
            cur["n_windows"] += 1
        else:
            regions.append(cur)
            cur = {"chrom": r["chrom"], "start": int(r["start"]),
                   "end": int(r["end"]),
                   "max_votes": int(r["votes"]),
                   "votes_sum": int(r["votes"]), "n_windows": 1}
    regions.append(cur)
    out = pd.DataFrame(regions)
    out["span_Mb"] = (out["end"] - out["start"]) / 1e6
    return out


def main() -> None:
    df = load_and_join()
    df, ind_cols = add_indicators(df)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TBL_MASTER, sep="\t", index=False)

    rows = []
    for v in range(0, len(ind_cols) + 1):
        sub = df[df["votes"] >= v]
        if sub.empty:
            continue
        merged = merge_adjacent(sub)
        rows.append(dict(min_votes=v, n_windows=len(sub),
                          n_regions_after_merge=len(merged),
                          median_span_Mb=merged["span_Mb"].median(),
                          total_span_Mb=merged["span_Mb"].sum()))
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_TBL_BYVOTES, sep="\t", index=False)

    cands = df[df["votes"] >= 2].copy()
    regs = merge_adjacent(cands).sort_values("max_votes", ascending=False)
    regs.to_csv(OUT_TBL_REGIONS, sep="\t", index=False)
    for _, r in regs.head(10).iterrows():
        sc = r["chrom"].split("_")[2] if "_" in r["chrom"] else r["chrom"]

if __name__ == "__main__":
    main()
