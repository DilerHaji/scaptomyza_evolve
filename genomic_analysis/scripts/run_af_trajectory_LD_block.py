#!/usr/bin/env python3
from __future__ import annotations
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sp_stats

sys.path.insert(0, str(Path(__file__).parent))
from section2_fig3c_diversity_strengthening import (
    load_master_with_votes,
    setup_rc, save,
    GENS, B_REPS, T_REPS, M_REPS, F_REPS,
    PEAK_CHROM, PEAK_START, WIN,
    COL_B, COL_T, COL_M, PREFIX, ROOT,
)

LD_BLOCK_START = 2_640_000
LD_BLOCK_END   = 3_610_000

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
                   (g_qc["pos"] >= LD_BLOCK_START) &
                   (g_qc["pos"] <  LD_BLOCK_END) &
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

    cand = bg_win_stats[bg_win_stats["maf_ok"] & bg_win_stats["depth_ok"]].reset_index(drop=True)

    if cand.empty:
        bg_win_stats["depth_ok"] = (
            (bg_win_stats["med_depth"] - focal_med_depth).abs() /
            focal_med_depth <= 0.50)
        cand = bg_win_stats[bg_win_stats["maf_ok"] & bg_win_stats["depth_ok"]].reset_index(drop=True)
    if cand.empty:
        return None

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

    pooled_idx = np.concatenate(
        [bg_by_window[(cand.iloc[wi]["chrom"], cand.iloc[wi]["bin"])]
         for wi in range(len(cand))])

    for i in range(N_RESAMPLES):
        sel = rng.choice(pooled_idx, size=n_target,
                         replace=len(pooled_idx) < n_target)
        af_sub = bg_arr[sel]
        flip_all_s, flip_loo_s = _build_flips(af_sub)
        rep_trajs = _aggregate_per_trt(af_sub, flip_all_s, flip_loo_s)
        for trt, _ in [("B", B_REPS), ("T", T_REPS), ("M", M_REPS)]:
            rep_traj = rep_trajs[trt]
            g01 = rep_traj[0:1, :]
            bg_cloud[trt][i, :] = np.nanmean(rep_traj - g01, axis=1)


    return focal_per_rep, bg_cloud, n_target, len(cand)


def per_gen_pairwise_welch(focal_per_rep):
    rows = []
    for gi, g in enumerate(GENS):
        a = focal_per_rep["B"][gi, :]
        b = focal_per_rep["T"][gi, :]
        c = focal_per_rep["M"][gi, :]
        for label, x, y in [("B-T", a, b), ("B-M", a, c), ("T-M", b, c)]:
            x = x[~np.isnan(x)]
            y = y[~np.isnan(y)]
            if len(x) < 2 or len(y) < 2:
                rows.append({"gen": g, "pair": label, "t": np.nan, "p": np.nan,
                             "n_x": len(x), "n_y": len(y),
                             "mean_x": np.nanmean(x) if len(x) else np.nan,
                             "mean_y": np.nanmean(y) if len(y) else np.nan})
                continue
            t, p = sp_stats.ttest_ind(x, y, equal_var=False)
            rows.append({"gen": g, "pair": label, "t": t, "p": p,
                         "n_x": len(x), "n_y": len(y),
                         "mean_x": float(np.mean(x)),
                         "mean_y": float(np.mean(y))})
    return pd.DataFrame(rows)


def whole_series_slope_test(focal_per_rep):
    g_arr = np.array(GENS, dtype=float)
    slope_dict = {}
    for trt, mat in focal_per_rep.items():
        slopes = np.full(mat.shape[1], np.nan)
        for ri in range(mat.shape[1]):
            y = mat[:, ri]
            m = ~np.isnan(y)
            if m.sum() >= 3:
                slopes[ri] = np.polyfit(g_arr[m], y[m], 1)[0]
        slope_dict[trt] = slopes

    rows = []
    for label, a, b in [("B-T", "B", "T"), ("B-M", "B", "M"), ("T-M", "T", "M")]:
        x = slope_dict[a]; y = slope_dict[b]
        x = x[~np.isnan(x)]; y = y[~np.isnan(y)]
        if len(x) < 2 or len(y) < 2:
            rows.append({"pair": label, "t": np.nan, "p": np.nan,
                         "mean_x": np.nan, "mean_y": np.nan})
            continue
        t, p = sp_stats.ttest_ind(x, y, equal_var=False)
        rows.append({"pair": label, "t": t, "p": p,
                     "mean_x": float(np.mean(x)),
                     "mean_y": float(np.mean(y))})
    return slope_dict, pd.DataFrame(rows)


def main():
    setup_rc()
    df_p, _ = load_master_with_votes()
    res = load_af_data(df_p)
    if res is None:
        return
    focal_per_rep, bg_cloud, n_target, n_cand = res

    pgen = per_gen_pairwise_welch(focal_per_rep)

    slope_dict, pwhole = whole_series_slope_test(focal_per_rep)

    fig, ax = plt.subplots(1, 1, figsize=(3.6, 3.6))
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
                   s=1.0, color=col, alpha=0.18, edgecolor="none", zorder=2)

    for trt in ("B", "T", "M"):
        col = cmap[trt]
        off = trt_offsets[trt]
        per_rep = focal_per_rep[trt]
        n = np.sum(~np.isnan(per_rep), axis=1)
        mean_g = np.nanmean(per_rep, axis=1)
        sd_g = np.nanstd(per_rep, axis=1, ddof=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            se_g = sd_g / np.sqrt(np.maximum(n, 1))
        ci = 1.96 * se_g
        x_focal = np.array(GENS) + off
        ax.errorbar(x_focal, mean_g, yerr=ci, fmt="-o", color=col,
                    ecolor=col, elinewidth=0.7, capsize=1.5, capthick=0.4,
                    markersize=3, linewidth=1.0,
                    markeredgecolor="black", markeredgewidth=0.3,
                    zorder=5, label="B+T" if trt == "M" else trt)

    for _, row in pgen[pgen["pair"] == "B-T"].iterrows():
        if pd.isna(row["p"]):
            continue
        if row["p"] < 0.001:
            star = "***"
        elif row["p"] < 0.01:
            star = "**"
        elif row["p"] < 0.05:
            star = "*"
        else:
            continue
        ax.text(row["gen"], 0.058, star, fontsize=7, ha="center", va="center",
                color="black")

    ax.axhline(0, color="#888888", linestyle=":", linewidth=0.5, alpha=0.7,
               zorder=1)
    ax.set_xlabel("Experimental Generation\n(t, Wright-Fisher Process)",
                  fontsize=6.5)
    ax.set_ylabel(r"B-polarized $\Delta$AF relative to $t_1$", fontsize=6.5)
    ax.set_xticks(GENS)
    ax.set_xlim(0.4, 10.6)
    ax.legend(frameon=False, fontsize=6, loc="upper left",
              handletextpad=0.3)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="both", labelsize=5.5, length=1.5)

    bt_p = float(pwhole.loc[pwhole["pair"] == "B-T", "p"].iloc[0])
    bm_p = float(pwhole.loc[pwhole["pair"] == "B-M", "p"].iloc[0])
    tm_p = float(pwhole.loc[pwhole["pair"] == "T-M", "p"].iloc[0])
    txt = ("whole-series slope Welch (4v4):\n"
           f"  B vs T: p={bt_p:.3f}\n"
           f"  B vs B+T: p={bm_p:.3f}\n"
           f"  T vs B+T: p={tm_p:.3f}\n"
           f"focal: {n_target} SNPs in LD block\n"
           f"bg: {n_cand} matched windows")
    ax.text(0.98, 0.02, txt, transform=ax.transAxes,
            fontsize=5, va="bottom", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", fc="white",
                      ec="#cccccc", lw=0.4))

    fig.tight_layout()
    save(fig, f"{PREFIX}_af_trajectory_LD_block")

    pgen.to_csv(ROOT / "final_plots/wild/section2_fig3c_diversity_af_trajectory_LD_block_pgen.tsv",
                sep="\t", index=False)
    pwhole.to_csv(ROOT / "final_plots/wild/section2_fig3c_diversity_af_trajectory_LD_block_whole.tsv",
                  sep="\t", index=False)


if __name__ == "__main__":
    main()
