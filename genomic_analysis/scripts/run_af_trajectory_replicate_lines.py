#!/usr/bin/env python3
from __future__ import annotations
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from section2_fig3c_diversity_strengthening import (
    load_master_with_votes,
    setup_rc, save,
    GENS, B_REPS, T_REPS, M_REPS, F_REPS,
    PEAK_CHROM, PEAK_START, WIN,
    COL_B, COL_T, COL_M, PREFIX, ROOT,
)


PB_THRESHOLD = 0.01
DEPTH_MIN = 5
BG_VOTE = 1
MAX_BG_SNPS = 30_000
N_RESAMPLES = 1000
MAF_TOL = 0.05
DEPTH_TOL_REL = 0.20


def load_af_data(df_p):
    AD_TSV = ROOT / "variance_analysis/merged_ad.tsv"
    SAMPLE_LIST = ROOT / "variance_analysis/sample_list.txt"
    GLMM = ROOT / "glm_lrt_gw_final/glmV1full.csv"
    samples = [s.strip() for s in open(SAMPLE_LIST)]

    bg_wins = df_p[df_p["votes_v3"] == BG_VOTE][["chrom", "start"]].copy()
    bg_chroms = set(bg_wins["chrom"])
    keep_chroms = bg_chroms | {PEAK_CHROM}
    bg_keys = set(zip(bg_wins["chrom"], bg_wins["start"]))

    glmm = pd.read_csv(GLMM)
    qc_base = ((glmm["converged"] == True) & (glmm["error"] == "OK") &
               (glmm["singular"] == False) & glmm["PB_p_val"].notna())
    g_qc = glmm[qc_base & glmm["chrom"].isin(keep_chroms)].copy()
    g_qc["bin"] = (g_qc["pos"] // WIN) * WIN

    g_peak = g_qc[(g_qc["chrom"] == PEAK_CHROM) &
                   (g_qc["bin"] == PEAK_START) &
                   (g_qc["PB_p_val"] < PB_THRESHOLD)]
    peak_pos = set(g_peak["pos"])

    g_bg_w = g_qc[g_qc.apply(
        lambda r: (r["chrom"], r["bin"]) in bg_keys, axis=1)]
    g_bg = g_bg_w[g_bg_w["PB_p_val"] >= PB_THRESHOLD]
    if len(g_bg) > MAX_BG_SNPS:
        g_bg = g_bg.sample(n=MAX_BG_SNPS, random_state=0)
    bg_pos_by_chrom = g_bg.groupby("chrom")["pos"].apply(set).to_dict()

    cols = ["chrom", "pos", "ref", "alt"] + samples
    chunks = []
    for chunk in pd.read_csv(AD_TSV, sep="\t", header=None, names=cols,
                              chunksize=200_000, dtype=str, low_memory=False):
        chunk = chunk[chunk["chrom"].isin(keep_chroms)].reset_index(drop=True)
        if chunk.empty:
            continue
        chunk["pos"] = chunk["pos"].astype(int)
        chrom_arr = chunk["chrom"].values
        pos_arr = chunk["pos"].values
        keep_mask = (chrom_arr == PEAK_CHROM) & np.isin(pos_arr, list(peak_pos))
        for chrom, pset in bg_pos_by_chrom.items():
            keep_mask |= (chrom_arr == chrom) & np.isin(pos_arr, list(pset))
        chunk = chunk[keep_mask]
        if not chunk.empty:
            chunks.append(chunk)
    df_ad = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()

    af = {}
    depth_arrs = []
    for s in samples:
        col = df_ad[s].astype(str).str.split(",", expand=True)
        ref = pd.to_numeric(col[0], errors="coerce").values
        alt = pd.to_numeric(col[1], errors="coerce").values
        depth = ref + alt
        af[s] = np.where(depth >= DEPTH_MIN, alt / np.maximum(depth, 1),
                          np.nan)
        depth_arrs.append(depth)
    af_df = pd.DataFrame(af)
    af_df["chrom"] = df_ad["chrom"].values
    af_df["pos"] = df_ad["pos"].values
    af_df["mean_depth"] = np.nanmean(np.stack(depth_arrs), axis=0)

    F_cols = [f"F{r}G00" for r in F_REPS if f"F{r}G00" in samples]
    f_alts, f_tots = [], []
    for s in F_cols:
        col = df_ad[s].astype(str).str.split(",", expand=True)
        ref = pd.to_numeric(col[0], errors="coerce").values
        alt = pd.to_numeric(col[1], errors="coerce").values
        f_alts.append(alt); f_tots.append(ref + alt)
    fa = np.nansum(np.stack(f_alts), axis=0)
    ft = np.nansum(np.stack(f_tots), axis=0)
    founder_af = np.where(ft > 0, fa / np.maximum(ft, 1), np.nan)
    af_df["founder_maf"] = np.minimum(founder_af, 1 - founder_af)
    af_df["bin"] = (af_df["pos"] // WIN) * WIN

    is_peak = (af_df["chrom"] == PEAK_CHROM) & af_df["pos"].isin(peak_pos)
    is_bg = af_df.apply(
        lambda r: r["chrom"] in bg_pos_by_chrom and
                  r["pos"] in bg_pos_by_chrom[r["chrom"]],
        axis=1) & ~is_peak

    af_peak = af_df[is_peak].reset_index(drop=True)
    af_bg = af_df[is_bg].reset_index(drop=True)
    n_target = len(af_peak)

    focal_med_maf = float(np.nanmedian(af_peak["founder_maf"]))
    focal_med_depth = float(np.nanmedian(af_peak["mean_depth"]))

    bg_win_stats = af_bg.groupby(["chrom", "bin"]).agg(
        n_snps=("pos", "size"),
        med_maf=("founder_maf", "median"),
        med_depth=("mean_depth", "median"),
    ).reset_index()
    bg_win_stats["maf_ok"] = (bg_win_stats["med_maf"] - focal_med_maf).abs() <= MAF_TOL
    bg_win_stats["depth_ok"] = (
        (bg_win_stats["med_depth"] - focal_med_depth).abs() / focal_med_depth
        <= DEPTH_TOL_REL)
    bg_win_stats["snps_ok"] = bg_win_stats["n_snps"] >= n_target
    cand = bg_win_stats[bg_win_stats["maf_ok"] & bg_win_stats["depth_ok"] &
                          bg_win_stats["snps_ok"]].reset_index(drop=True)

    bg_by_window = {}
    for (chrom, bin_), g in af_bg.groupby(["chrom", "bin"]):
        bg_by_window[(chrom, bin_)] = g.index.to_numpy()

    sample_cols = list(samples)
    sample_idx = {s: i for i, s in enumerate(sample_cols)}
    bg_arr = af_bg[sample_cols].values.astype(float)
    peak_arr = af_peak[sample_cols].values.astype(float)

    F_idx = [sample_idx[f"F{r}G00"] for r in F_REPS if f"F{r}G00" in sample_idx]
    B10_all_idx = [sample_idx[f"B{r}G10"] for r in B_REPS
                    if f"B{r}G10" in sample_idx]
    B10_loo_idx = {r: [sample_idx[f"B{rr}G10"] for rr in B_REPS
                         if rr != r and f"B{rr}G10" in sample_idx]
                    for r in B_REPS}

    def _build_flips(af_arr):
        founder = np.nanmean(af_arr[:, F_idx], axis=1)
        b10_all = np.nanmean(af_arr[:, B10_all_idx], axis=1)
        flip_all = b10_all < founder
        flip_loo = {}
        for r in B_REPS:
            b10_oth = np.nanmean(af_arr[:, B10_loo_idx[r]], axis=1)
            flip_loo[r] = b10_oth < founder
        return flip_all, flip_loo

    def _aggregate_per_trt(af_arr, flip_all, flip_loo):
        out = {}
        for trt, reps in [("B", B_REPS), ("T", T_REPS), ("M", M_REPS)]:
            rep_traj = np.full((len(GENS), len(reps)), np.nan)
            for ri, r in enumerate(reps):
                flip = flip_loo[r] if trt == "B" else flip_all
                for gi, g in enumerate(GENS):
                    key = f"{trt}{r}G{g:02d}"
                    if key in sample_idx:
                        vals = af_arr[:, sample_idx[key]].copy()
                        vals[flip] = 1.0 - vals[flip]
                        rep_traj[gi, ri] = np.nanmean(vals)
            out[trt] = rep_traj
        return out

    flip_all_pk, flip_loo_pk = _build_flips(peak_arr)
    rep_trajs_focal = _aggregate_per_trt(peak_arr, flip_all_pk, flip_loo_pk)
    focal_per_rep = {}
    for trt, _ in [("B", B_REPS), ("T", T_REPS), ("M", M_REPS)]:
        rep_traj = rep_trajs_focal[trt]
        g01 = rep_traj[0:1, :]
        focal_per_rep[trt] = rep_traj - g01

    rng = np.random.default_rng(0)
    bg_cloud = {trt: np.full((N_RESAMPLES, len(GENS)), np.nan)
                 for trt in ("B", "T", "M")}
    for i in range(N_RESAMPLES):
        wi = int(rng.integers(0, len(cand)))
        win_key = (cand.iloc[wi]["chrom"], cand.iloc[wi]["bin"])
        snp_idx_in_win = bg_by_window[win_key]
        replace = len(snp_idx_in_win) < n_target
        sel = rng.choice(snp_idx_in_win, size=n_target, replace=replace)
        af_sub = bg_arr[sel]
        flip_all_s, flip_loo_s = _build_flips(af_sub)
        rep_trajs = _aggregate_per_trt(af_sub, flip_all_s, flip_loo_s)
        for trt, _ in [("B", B_REPS), ("T", T_REPS), ("M", M_REPS)]:
            rep_traj = rep_trajs[trt]
            g01 = rep_traj[0:1, :]
            bg_cloud[trt][i, :] = np.nanmean(rep_traj - g01, axis=1)

    return focal_per_rep, bg_cloud, n_target, len(cand)


def main():
    setup_rc()
    df_p, _ = load_master_with_votes()
    focal_per_rep, bg_cloud, n_target, n_cand = load_af_data(df_p)

    fig, ax = plt.subplots(1, 1, figsize=(3.0, 3.0))
    trt_offsets = {"B": -0.25, "T": 0.0, "M": 0.25}
    cmap = {"B": COL_B, "T": COL_T, "M": COL_M}
    rng = np.random.default_rng(0)

    for trt in ("B", "T", "M"):
        col = cmap[trt]
        off = trt_offsets[trt]
        x_all, y_all = [], []
        for gi, g in enumerate(GENS):
            vals = bg_cloud[trt][:, gi]
            x_jit = g + off + rng.uniform(-0.08, 0.08, len(vals))
            x_all.append(x_jit)
            y_all.append(vals)
        ax.scatter(np.concatenate(x_all), np.concatenate(y_all),
                   s=1.0, color=col, alpha=0.22, edgecolor="none", zorder=2)

    for trt in ("B", "T", "M"):
        col = cmap[trt]
        off = trt_offsets[trt]
        per_rep = focal_per_rep[trt]   # (n_gens, n_reps)
        x_focal = np.array(GENS) + off
        for ri in range(per_rep.shape[1]):
            y = per_rep[:, ri]
            mask = ~np.isnan(y)
            ax.plot(x_focal[mask], y[mask],
                    color=col, lw=0.7, alpha=0.85, zorder=4,
                    marker="o", markersize=2.2,
                    markeredgewidth=0)

    ax.axhline(0, color="#888888", linestyle=":", linewidth=0.5, alpha=0.7,
               zorder=1)
    ax.set_xlabel("Experimental Generation\n(t, Wright-Fisher Process)",
                  fontsize=6.5)
    ax.set_ylabel(r"B-polarized $\Delta$AF relative to $t_1$", fontsize=6.5)
    ax.set_xticks(GENS)
    ax.set_xlim(0.4, 10.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="both", labelsize=5.5, length=1.5)

    fig.tight_layout()
    save(fig, f"{PREFIX}_af_trajectory_replicate_lines")


if __name__ == "__main__":
    main()
