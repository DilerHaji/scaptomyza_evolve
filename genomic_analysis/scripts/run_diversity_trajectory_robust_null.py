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
    load_master_with_votes, load_trajectory, per_rep_gen_mean,
    setup_rc, save,
    GENS, B_REPS, T_REPS, M_REPS, F_REPS,
    PEAK_CHROM, PEAK_START, WIN,
    COL_B, COL_T, COL_M, PREFIX,
)


DIST_MB = 1.0
N_SNPS_TOL = 0.20
THETA_TOL = 0.20
CANDIDATE_VOTE = 3

STAT_LABELS = {
    "theta_pi": (r"\theta_\pi", "θπ"),
    "theta_watterson": (r"\theta_W", "θW"),
}


def far_from_candidates(df, dist_bp, vote_min=CANDIDATE_VOTE):
    cand = df[df["votes_v3"] >= vote_min]
    cand_centers = {c: (g["start"].values + WIN / 2)
                    for c, g in cand.groupby("chrom")}
    centers = df["start"].values + WIN / 2
    chrom_arr = df["chrom"].values
    far = np.ones(len(df), dtype=bool)
    for i in range(len(df)):
        c = chrom_arr[i]
        if c in cand_centers:
            if np.min(np.abs(cand_centers[c] - centers[i])) < dist_bp:
                far[i] = False
    return far


def fig_robust(df_p, traj, stat: str):
    tex_lab, ascii_lab = STAT_LABELS[stat]
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

    n_cand = (df_p["votes_v3"] >= CANDIDATE_VOTE).sum()
    far_mask = far_from_candidates(df_p, DIST_MB * 1e6, CANDIDATE_VOTE)

    focal = df_p[(df_p["chrom"] == PEAK_CHROM) & (df_p["start"] == PEAK_START)]
    if focal.empty:
        return
    focal_n = float(focal["n_snps"].iloc[0])
    focal_t = float(focal["founder_stat"].iloc[0])


    n_lo, n_hi = focal_n * (1 - N_SNPS_TOL), focal_n * (1 + N_SNPS_TOL)
    t_lo, t_hi = focal_t * (1 - THETA_TOL), focal_t * (1 + THETA_TOL)
    n_ok = (df_p["n_snps"] >= n_lo) & (df_p["n_snps"] <= n_hi)
    t_ok = (df_p["founder_stat"] >= t_lo) & (df_p["founder_stat"] <= t_hi)
    match_mask = (n_ok & t_ok).values

    null_mask = far_mask & match_mask & valid

    is_focal = ((df_p["chrom"] == PEAK_CHROM) &
                (df_p["start"] == PEAK_START)).values
    null_mask = null_mask & ~is_focal

    null_df_idx = np.where(null_mask)[0]
    null_traj_idx = df_p.iloc[null_df_idx]["traj_idx"].values.astype(int)



    pi_B = per_rep_gen_mean(traj, "B", GENS, B_REPS, stat)
    pi_T = per_rep_gen_mean(traj, "T", GENS, T_REPS, stat)
    pi_M = per_rep_gen_mean(traj, "M", GENS, M_REPS, stat)
    peak_idx = key_to_idx[(PEAK_CHROM, PEAK_START)]


    gw_g01 = {}    # mean θ at G01
    gw_g10 = {}    # mean θ at G10
    gw_fold = {}   # mean (G10 / G01) — i.e., genome-wide fold-change
    for lab, mat in [("B", pi_B), ("T", pi_T), ("M", pi_M)]:
        gw_g01[lab] = float(np.nanmean(mat[:, 0]))
        gw_g10[lab] = float(np.nanmean(mat[:, -1]))
        gw_fold[lab] = gw_g10[lab] / gw_g01[lab]


    fig, axes = plt.subplots(1, 4, figsize=(11.0, 3.0),
                             gridspec_kw={"width_ratios": [1.3, 1.3, 1, 1]})
    rng = np.random.default_rng(0)
    x_pos = {"B": 1, "T": 2, "M": 3}
    cmap = {"B": COL_B, "T": COL_T, "M": COL_M}

    def _strip(ax, focal_d, null_d, ylab, title):
        for lab, x in x_pos.items():
            col = cmap[lab]
            vals = null_d[lab]
            x_jit = x + rng.uniform(-0.22, 0.22, len(vals))
            ax.scatter(x_jit, vals, s=5, color=col, alpha=0.45,
                       edgecolor="none", zorder=1)
            f = focal_d[lab]
            ax.scatter([x], [f], s=30, color=col, edgecolor="black",
                       linewidth=0.5, zorder=4)
            p = ((vals <= f).sum() + 1) / (len(vals) + 1)
            p_str = f"p={p:.3f}" if p >= 1e-3 else "p<.001"
            ax.text(x, f, p_str, fontsize=6, ha="center", va="bottom")

        ax.axhline(1, color="#888888", linestyle=":", linewidth=0.5, alpha=0.7)
        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(["B", "T", "M"], fontsize=7)
        ax.set_ylabel(ylab, fontsize=7.5)
        ax.set_title(title, fontsize=8, loc="left", pad=2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(axis="both", labelsize=6.5, length=2)

    focal_d, null_d = {}, {}
    for lab, mat in [("B", pi_B), ("T", pi_T), ("M", pi_M)]:
        focal_g01 = mat[peak_idx, 0]
        focal_g10 = mat[peak_idx, -1]
        focal_fold = focal_g10 / focal_g01 if focal_g01 > 0 else np.nan
        focal_d[lab] = focal_fold / gw_fold[lab]

        null_g01 = mat[null_traj_idx, 0]
        null_g10 = mat[null_traj_idx, -1]
        with np.errstate(invalid="ignore", divide="ignore"):
            null_fold = np.where(null_g01 > 0, null_g10 / null_g01, np.nan)
            nv = null_fold / gw_fold[lab]
        null_d[lab] = nv[~np.isnan(nv)]
    _strip(axes[0], focal_d, null_d,
           rf"$({{{tex_lab}}}_{{G10}} / {{{tex_lab}}}_{{G01}})$ /  genome-wide fold-change",
           "ratio of G10/G01 fold-changes (low-tail vs matched bg)")

    focal_e, null_e = {}, {}
    for lab, mat in [("B", pi_B), ("T", pi_T), ("M", pi_M)]:
        focal_e[lab] = mat[peak_idx, -1] / gw_g10[lab]
        with np.errstate(invalid="ignore", divide="ignore"):
            nv = mat[null_traj_idx, -1] / gw_g10[lab]
        null_e[lab] = nv[~np.isnan(nv)]
    _strip(axes[1], focal_e, null_e,
           rf"${tex_lab}$ at G10  /  genome-wide mean",
           "G10 endpoint ratio (low-tail vs matched bg)")

    null_n = df_p.iloc[null_df_idx]["n_snps"].values
    all_n = df_p.loc[far_mask & valid, "n_snps"].values
    axes[2].hist(all_n, bins=30, color="#dddddd", alpha=0.7,
                 edgecolor="white", linewidth=0.3, label="all far-bg")
    axes[2].hist(null_n, bins=30, color="#666666", alpha=0.85,
                 edgecolor="white", linewidth=0.3, label="matched")
    axes[2].axvline(focal_n, color="red", linewidth=1.2, label="focal")
    axes[2].set_xlabel("n_snps per window", fontsize=7)
    axes[2].set_title(f"n_snps match  ±{int(N_SNPS_TOL*100)}%",
                      fontsize=7.5, loc="left", pad=2)
    axes[2].legend(frameon=False, fontsize=5.5, loc="upper right")


    null_t = df_p.iloc[null_df_idx]["founder_stat"].values
    all_t = df_p.loc[far_mask & valid, "founder_stat"].values
    axes[3].hist(all_t, bins=30, color="#dddddd", alpha=0.7,
                 edgecolor="white", linewidth=0.3, label="all far-bg")
    axes[3].hist(null_t, bins=30, color="#666666", alpha=0.85,
                 edgecolor="white", linewidth=0.3, label="matched")
    axes[3].axvline(focal_t, color="red", linewidth=1.2, label="focal")
    axes[3].set_xlabel(rf"founder ${tex_lab}$", fontsize=7)
    axes[3].set_title(f"founder {ascii_lab} match  +/-{int(THETA_TOL*100)}%",
                      fontsize=7.5, loc="left", pad=2)
    axes[3].legend(frameon=False, fontsize=5.5, loc="upper right")

    for ax in axes[2:]:
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(axis="both", labelsize=6.5, length=2)
        ax.set_ylabel("count", fontsize=7)

    fig.suptitle(f"chr_439:{PEAK_START/1e6:.1f}-{(PEAK_START+WIN)/1e6:.1f} Mb  "
                 f"({stat}; n_null={len(null_traj_idx)}; >={DIST_MB} Mb from vote>={CANDIDATE_VOTE} "
                 f"+ covariate-matched on n_snps & founder {ascii_lab})",
                 fontsize=8, y=1.01)
    fig.tight_layout()
    save(fig, f"{PREFIX}_trajectory_robust_{stat}")


def main():
    stat = sys.argv[1] if len(sys.argv) > 1 else "theta_pi"
    if stat not in STAT_LABELS:
        sys.exit(f"unknown STAT '{stat}'; choose one of {list(STAT_LABELS)}")
    setup_rc()
    df_p, keep = load_master_with_votes()
    traj = load_trajectory(keep)
    fig_robust(df_p, traj, stat)


if __name__ == "__main__":
    main()
