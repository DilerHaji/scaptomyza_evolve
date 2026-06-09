#!/usr/bin/env python3
from __future__ import annotations
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from section2_fig3c_diversity_strengthening import (
    load_master_with_votes, load_trajectory,
    per_rep_gen_full,
    setup_rc, save,
    GENS, B_REPS, T_REPS, M_REPS, F_REPS,
    PEAK_CHROM, PEAK_START, WIN,
    COL_B, COL_T, COL_M, PREFIX,
)
from run_diversity_trajectory_robust_null import (
    far_from_candidates,
    DIST_MB, N_SNPS_TOL, THETA_TOL, CANDIDATE_VOTE,
    STAT_LABELS,
)


LD_BLOCK_START = 2_640_000
LD_BLOCK_END   = 3_610_000
N_BOOT = 1000


def _per_rep_slopes(mat_full, gens_arr):
    n_w, _, n_r = mat_full.shape
    out = np.full((n_w, n_r), np.nan)
    for ri in range(n_r):
        y = mat_full[:, :, ri]
        for wi in range(n_w):
            yi = y[wi]
            m = ~np.isnan(yi)
            if m.sum() >= 3:
                out[wi, ri] = np.polyfit(gens_arr[m], yi[m], 1)[0]
    return out


def compute_panel(df_p, traj, stat):
    F_cols_stat = [f"F{r}G00.{stat}" for r in F_REPS
                   if f"F{r}G00.{stat}" in traj.columns]
    traj = traj.copy()
    traj["founder_stat"] = traj[F_cols_stat].mean(axis=1, skipna=True).values

    keys = list(zip(traj["chrom"], traj["bin"]))
    key_to_idx = {k: i for i, k in enumerate(keys)}
    df_p = df_p.copy()
    df_p["traj_idx"] = [key_to_idx.get((c, s), -1)
                        for c, s in zip(df_p["chrom"], df_p["start"])]
    valid = (df_p["traj_idx"] >= 0).values
    df_p["founder_stat"] = np.nan
    df_p.loc[valid, "founder_stat"] = (
        traj.iloc[df_p.loc[valid, "traj_idx"].values]["founder_stat"].values
    )

    far_mask = far_from_candidates(df_p, DIST_MB * 1e6, CANDIDATE_VOTE)

    focal = df_p[(df_p["chrom"] == PEAK_CHROM) & (df_p["start"] == PEAK_START)]
    focal_n = float(focal["n_snps"].iloc[0])
    focal_t = float(focal["founder_stat"].iloc[0])
    n_lo, n_hi = focal_n * (1 - N_SNPS_TOL), focal_n * (1 + N_SNPS_TOL)
    t_lo, t_hi = focal_t * (1 - THETA_TOL), focal_t * (1 + THETA_TOL)
    n_ok = (df_p["n_snps"] >= n_lo) & (df_p["n_snps"] <= n_hi)
    t_ok = (df_p["founder_stat"] >= t_lo) & (df_p["founder_stat"] <= t_hi)

    focal_mask = ((df_p["chrom"] == PEAK_CHROM) &
                  (df_p["start"] >= LD_BLOCK_START) &
                  (df_p["start"] < LD_BLOCK_END))
    focal_traj_idx = df_p.loc[focal_mask & valid, "traj_idx"].values.astype(int)

    is_focal_in_scale = focal_mask.values
    null_mask = (far_mask & (n_ok & t_ok).values & valid &
                 ~is_focal_in_scale)
    null_traj_idx = df_p.loc[null_mask, "traj_idx"].values.astype(int)
    K = len(focal_traj_idx)

    pi_B_full = per_rep_gen_full(traj, "B", GENS, B_REPS, stat)
    pi_T_full = per_rep_gen_full(traj, "T", GENS, T_REPS, stat)
    pi_M_full = per_rep_gen_full(traj, "M", GENS, M_REPS, stat)

    gens_arr = np.array(GENS, dtype=float)
    slope_by_lab = {
        "B": _per_rep_slopes(pi_B_full, gens_arr),
        "T": _per_rep_slopes(pi_T_full, gens_arr),
        "M": _per_rep_slopes(pi_M_full, gens_arr),
    }

    rng = np.random.default_rng(0)
    focal_d, null_d, p_d = {}, {}, {}

    for lab, slope_mat in slope_by_lab.items():
        focal_per_rep = np.nanmean(slope_mat[focal_traj_idx, :], axis=0)
        null_dist = np.full(N_BOOT, np.nan)
        for b in range(N_BOOT):
            bs = rng.choice(null_traj_idx, size=K, replace=True)
            bs_per_rep = np.nanmean(slope_mat[bs, :], axis=0)
            null_dist[b] = np.nanmean(bs_per_rep)
        null_dist = null_dist[~np.isnan(null_dist)]
        null_median = float(np.median(null_dist))

        focal_d[lab] = focal_per_rep - null_median
        null_d[lab] = null_dist - null_median

        f_mean = float(np.nanmean(focal_per_rep))
        p = ((null_dist <= f_mean).sum() + 1) / (len(null_dist) + 1)
        p_d[lab] = p

    return focal_d, null_d, p_d, K


def _strip(ax, focal_d, null_d, p_d, ylab, title, rng):
    rep_ids = {"B": B_REPS, "T": T_REPS, "M": M_REPS}
    cmap = {"B": COL_B, "T": COL_T, "M": COL_M}
    x_pos = {"B": 1, "T": 2, "M": 3}

    for lab, x in x_pos.items():
        col = cmap[lab]
        vals = null_d[lab]
        x_jit = x + rng.uniform(-0.22, 0.22, len(vals))
        ax.scatter(x_jit, vals, s=3, color=col, alpha=0.35,
                   edgecolor="none", zorder=1)

        f_r = focal_d[lab]
        keep = ~np.isnan(f_r)
        reps_kept = [r for r, ok in zip(rep_ids[lab], keep) if ok]
        f_r = f_r[keep]
        x_focal = x + rng.uniform(-0.10, 0.10, len(f_r))
        ax.scatter(x_focal, f_r, s=20, color=col, edgecolor="black",
                   linewidth=0.5, zorder=4)
        for xi, yi, ri in zip(x_focal, f_r, reps_kept):
            ax.text(xi + 0.07, yi, str(ri), fontsize=5.5,
                    color="black", ha="left", va="center", zorder=5)

        p = p_d[lab]
        p_str = f"p={p:.3f}" if p >= 1e-3 else "p<.001"
        f_mean = float(np.mean(f_r))
        weight = "bold" if p < 0.05 else "normal"
        ax.text(x, f_mean, p_str, fontsize=6, ha="center", va="bottom",
                fontweight=weight)

    ax.axhline(0, color="#888888", linestyle=":", linewidth=0.5, alpha=0.7)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["B", "T", "B+T"], fontsize=6.5)
    ax.set_ylabel(ylab, fontsize=6.5)
    ax.set_title(title, fontsize=7.5, loc="left", pad=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="both", labelsize=5.5, length=2)
    ax.set_box_aspect(1)


def main():
    setup_rc()
    df_p, keep = load_master_with_votes()
    traj = load_trajectory(keep)

    fpi, npi, ppi, K_pi = compute_panel(df_p, traj, "theta_pi")

    fw, nw, pw, K_w = compute_panel(df_p, traj, "theta_watterson")

    fig, axes = plt.subplots(1, 2, figsize=(4.6, 2.6),
                             gridspec_kw={"wspace": 0.55})
    rng = np.random.default_rng(0)

    pi_tex, _ = STAT_LABELS["theta_pi"]
    w_tex, _ = STAT_LABELS["theta_watterson"]

    _strip(axes[0], fpi, npi, ppi,
           ylab=(rf"per-rep ${pi_tex}$ slope (G1-G10)"
                 "\n− matched-bg bootstrap median"),
           title=rf"${pi_tex}$",
           rng=rng)
    _strip(axes[1], fw, nw, pw,
           ylab=(rf"per-rep ${w_tex}$ slope (G1-G10)"
                 "\n− matched-bg bootstrap median"),
           title=rf"${w_tex}$",
           rng=rng)

    fig.tight_layout()
    save(fig, f"{PREFIX}_panel_combined")


if __name__ == "__main__":
    main()
