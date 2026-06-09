#!/usr/bin/env python3
from __future__ import annotations
import re
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(".")
MASTER_TSV = ROOT / "final_plots/wild/section2_candidate_master_v2.tsv"
GLMM_PER_WINDOW_TSV = ROOT / "final_plots/wild/section2_glmm_per_window.tsv"
TRAJ_CSV = ROOT / "grenfst/diversity/trajectory_pi_200000_nonoverlap.csv"
OUT_DIR = ROOT / "final_plots/wild"
PREFIX = OUT_DIR / "section2_fig3c_diversity"

WIN = 200_000
SCAFF_MIN_WINDOWS = 50
TOP_FRAC = 0.05
PEAK_CHROM = "chr_ScDA7r2_439_HRSCAF_779"
PEAK_START = 2_800_000
N_PERM = 5000

GREY = "#bdbdbd"
RED = "#C84A45"
BLUE = "#3a6fb0"
ORANGE = "#d99440"

GENS = list(range(1, 11))
B_REPS = [1, 2, 3, 4]
T_REPS = [1, 2, 3, 4]
M_REPS = [1, 2, 3, 4]
F_REPS = [1, 2, 3, 4]

# treatment palette (peak lines)
COL_B = "#C84A45"   # red
COL_T = "#d99440"   # orange
COL_M = "#3a6fb0"   # blue


def emp_p(value, null_vals, direction):
    null_vals = np.asarray(null_vals, dtype=float)
    null_vals = null_vals[~np.isnan(null_vals)]
    if len(null_vals) == 0 or np.isnan(value):
        return np.nan, 0
    if direction == "high":
        n_extreme = int((null_vals >= value).sum())
    else:
        n_extreme = int((null_vals <= value).sum())
    return (n_extreme + 1) / (len(null_vals) + 1), len(null_vals)


def fmt_p(p):
    if np.isnan(p):
        return "p=NA"
    return f"p={p:.3f}" if p >= 1e-3 else "p<0.001"


def setup_rc():
    plt.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.6,
    })


def save(fig, base):
    fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{base}.svg", bbox_inches="tight")
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    svg = Path(f"{base}.svg")
    txt = svg.read_text()
    txt = re.sub(r"<clipPath[^>]*>.*?</clipPath>", "", txt, flags=re.DOTALL)
    txt = re.sub(r'\s*clip-path="url\([^)]+\)"', "", txt)
    svg.write_text(txt)
    plt.close(fig)


def load_master_with_votes():
    df = pd.read_csv(MASTER_TSV, sep="\t")
    glmm = pd.read_csv(GLMM_PER_WINDOW_TSV, sep="\t")
    df = df.merge(glmm[["chrom", "start", "glmm_lrt"]],
                  on=["chrom", "start"], how="left")
    chrom_sizes = df.groupby("chrom").size().sort_values(ascending=False)
    keep = chrom_sizes[chrom_sizes >= SCAFF_MIN_WINDOWS].index.tolist()
    df_p = df[df["chrom"].isin(keep)].copy().reset_index(drop=True)
    glmm_thr = df_p["glmm_lrt"].quantile(1 - TOP_FRAC)
    df_p["IND_glmm_lrt"] = (df_p["glmm_lrt"] >= glmm_thr).fillna(False).astype(int)
    IND_COLS = [
        "IND_HV_blocks_B", "IND_HV_blocks_T", "IND_HV_blocks_M",
        "IND_cov_BT_neg", "IND_slope_div", "IND_permFST",
        "IND_glmm_lrt", "IND_wild_C2",
    ]
    df_p["votes_v3"] = df_p[IND_COLS].fillna(0).astype(int).sum(axis=1)
    return df_p, keep


def load_trajectory(keep):
    raw = pd.read_csv(TRAJ_CSV)
    raw = raw[raw["chrom"].isin(set(keep))].copy()
    raw["bin"] = (raw["start"] // WIN) * WIN
    raw = raw.rename(columns={"chrom": "chrom"})
    return raw


def per_rep_gen_mean(df, prefix_pattern, gens, reps, suffix):
    out = np.full((len(df), len(gens)), np.nan)
    for gi, g in enumerate(gens):
        cols = []
        for r in reps:
            c = f"{prefix_pattern}{r}G{g:02d}.{suffix}"
            if c in df.columns:
                cols.append(c)
        if cols:
            out[:, gi] = df[cols].mean(axis=1, skipna=True).values
    return out


def per_rep_gen_full(df, prefix_pattern, gens, reps, suffix):
    out = np.full((len(df), len(gens), len(reps)), np.nan)
    for gi, g in enumerate(gens):
        for ri, r in enumerate(reps):
            c = f"{prefix_pattern}{r}G{g:02d}.{suffix}"
            if c in df.columns:
                out[:, gi, ri] = df[c].values
    return out


def founder_mean(df, suffix="theta_pi"):
    cols = [f"F{r}G00.{suffix}" for r in F_REPS if f"F{r}G00.{suffix}" in df.columns]
    return df[cols].mean(axis=1, skipna=True).values


def fig_endpoint_ratio(df_p, traj):
    pi_F = founder_mean(traj, "theta_pi")
    pi_B10 = per_rep_gen_mean(traj, "B", [10], B_REPS, "theta_pi")[:, 0]
    pi_T10 = per_rep_gen_mean(traj, "T", [10], T_REPS, "theta_pi")[:, 0]
    ratio_B = np.where(pi_F > 0, pi_B10 / pi_F, np.nan)
    ratio_T = np.where(pi_F > 0, pi_T10 / pi_F, np.nan)
    ratio_BT = np.nanmean(np.stack([ratio_B, ratio_T]), axis=0)
    keys = list(zip(traj["chrom"], traj["bin"]))
    map_B = dict(zip(keys, ratio_B))
    map_T = dict(zip(keys, ratio_T))
    map_BT = dict(zip(keys, ratio_BT))
    df_p = df_p.copy()
    df_p["pi_ratio_B"] = [map_B.get((c, s), np.nan) for c, s in
                          zip(df_p["chrom"], df_p["start"])]
    df_p["pi_ratio_T"] = [map_T.get((c, s), np.nan) for c, s in
                          zip(df_p["chrom"], df_p["start"])]
    df_p["pi_ratio_BT"] = [map_BT.get((c, s), np.nan) for c, s in
                           zip(df_p["chrom"], df_p["start"])]
    _pointcloud_three(df_p, ["pi_ratio_B", "pi_ratio_T", "pi_ratio_BT"],
                      ["B", "T", "B+T"],
                      r"$\theta_\pi^{G10} / \theta_\pi^{F}$",
                      "low",
                      f"{PREFIX}_endpoint_ratio")



def fig_trajectory(df_p, traj):
    pi_B = per_rep_gen_mean(traj, "B", GENS, B_REPS, "theta_pi")
    pi_T = per_rep_gen_mean(traj, "T", GENS, T_REPS, "theta_pi")
    pi_M = per_rep_gen_mean(traj, "M", GENS, M_REPS, "theta_pi")
    keys = list(zip(traj["chrom"], traj["bin"]))
    key_to_idx = {k: i for i, k in enumerate(keys)}

    peak_key = (PEAK_CHROM, PEAK_START)
    peak_idx = key_to_idx.get(peak_key, -1)
    if peak_idx < 0:
        return

    df_p = df_p.copy()
    df_p["traj_idx"] = [key_to_idx.get((c, s), -1)
                         for c, s in zip(df_p["chrom"], df_p["start"])]
    null_idx = df_p.loc[(df_p["votes_v3"] == 0) & (df_p["traj_idx"] >= 0),
                          "traj_idx"].values

    fig, axes = plt.subplots(1, 2, figsize=(5.4, 3.0))
    rng = np.random.default_rng(0)
    x_pos = {"B": 1, "T": 2, "M": 3}
    cmap = {"B": COL_B, "T": COL_T, "M": COL_M}

    def _strip(ax, focal_dict, null_dict, ylab, title, direction):
        for lab, x in x_pos.items():
            col = cmap[lab]
            vals = null_dict[lab]
            x_jit = x + rng.uniform(-0.22, 0.22, len(vals))
            ax.scatter(x_jit, vals, s=3, color=col, alpha=0.30,
                       edgecolor="none", zorder=1)
            f = focal_dict[lab]
            ax.scatter([x], [f], s=28, color=col, edgecolor="black",
                       linewidth=0.5, zorder=4)
            if direction == "low":
                n_ext = int((vals <= f).sum())
            else:
                n_ext = int((vals >= f).sum())
            p = (n_ext + 1) / (len(vals) + 1)
            p_str = f"p={p:.3f}" if p >= 1e-3 else "p<.001"
            ax.text(x, f, p_str, fontsize=6, ha="center", va="bottom",
                    color="black")
        ax.axhline(0 if "change" in title.lower() else
                    float(np.nanmedian(np.concatenate(list(null_dict.values())))),
                    color="#888888", linestyle=":", linewidth=0.5,
                    alpha=0.7)
        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(["B", "T", "M"], fontsize=7)
        ax.set_ylabel(ylab, fontsize=7.5)
        ax.set_title(title, fontsize=8, loc="left", pad=2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(axis="both", labelsize=6.5, length=2)

    focal_d = {}
    null_d = {}
    for lab, mat in [("B", pi_B), ("T", pi_T), ("M", pi_M)]:
        focal_d[lab] = mat[peak_idx, -1] - mat[peak_idx, 0]
        nv = mat[null_idx, -1] - mat[null_idx, 0]
        null_d[lab] = nv[~np.isnan(nv)]
    _strip(axes[0], focal_d, null_d,
           r"$\theta_\pi(G10) - \theta_\pi(G01)$",
           "change G01 to G10 (low-tail vs vote=0)",
           "low")

    focal_e = {}
    null_e = {}
    for lab, mat in [("B", pi_B), ("T", pi_T), ("M", pi_M)]:
        focal_e[lab] = mat[peak_idx, -1]
        nv = mat[null_idx, -1]
        null_e[lab] = nv[~np.isnan(nv)]
    _strip(axes[1], focal_e, null_e,
           r"$\theta_\pi$ at G10",
           "G10 endpoint (low-tail vs vote=0)",
           "low")

    fig.suptitle(f"chr_439:{PEAK_START/1e6:.1f}-{(PEAK_START+WIN)/1e6:.1f} Mb "
                  f"(n_null={len(null_idx)})",
                  fontsize=8, y=0.99)
    fig.tight_layout()
    save(fig, f"{PREFIX}_trajectory")


def fig_trajectory_lines(df_p, traj):
    pi_B = per_rep_gen_full(traj, "B", GENS, B_REPS, "theta_pi")
    pi_T = per_rep_gen_full(traj, "T", GENS, T_REPS, "theta_pi")
    pi_M = per_rep_gen_full(traj, "M", GENS, M_REPS, "theta_pi")
    keys = list(zip(traj["chrom"], traj["bin"]))
    key_to_idx = {k: i for i, k in enumerate(keys)}
    peak_idx = key_to_idx.get((PEAK_CHROM, PEAK_START), -1)
    if peak_idx < 0:
        return

    fig, ax = plt.subplots(1, 1, figsize=(4.6, 3.0))
    for arr, col, lab in [(pi_B, COL_B, "B"),
                           (pi_T, COL_T, "T"),
                           (pi_M, COL_M, "M")]:
        rep_vals = arr[peak_idx, :, :]              # (n_gens, n_reps)
        g01 = rep_vals[0:1, :]
        delta_rep = rep_vals - g01
        n_per_gen = np.sum(~np.isnan(delta_rep), axis=1)
        mean_g = np.nanmean(delta_rep, axis=1)
        sd_g = np.nanstd(delta_rep, axis=1, ddof=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            se_g = sd_g / np.sqrt(np.maximum(n_per_gen, 1))
        ax.fill_between(GENS, mean_g - se_g, mean_g + se_g,
                        color=col, alpha=0.18, edgecolor="none", zorder=2)
        ax.plot(GENS, mean_g, color=col, linewidth=1.5, marker="o",
                markersize=4, zorder=3, label=lab)

    ax.axhline(0, color="#888888", linestyle=":", linewidth=0.5, alpha=0.7,
               zorder=1)
    ax.axvspan(2.6, 5.4, color="#fafafa", zorder=0)
    ax.text(4.0, 0.0, "2-rep gens", fontsize=5.5, ha="center", va="bottom",
            color="#888888", transform=ax.get_xaxis_transform())

    ax.set_xlabel("generation", fontsize=8)
    ax.set_ylabel(r"$\theta_\pi$ change vs G01 (per-rep)", fontsize=8)
    ax.set_xticks(GENS)
    ax.set_xlim(0.5, 10.5)
    ax.set_title(f"chr_439:{PEAK_START/1e6:.1f}-{(PEAK_START+WIN)/1e6:.1f} Mb "
                  f"(band = ±1 SE across reps per gen)",
                  fontsize=7.5, loc="left", pad=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="both", labelsize=7, length=2)
    ax.legend(loc="best", frameon=False, fontsize=7)
    fig.tight_layout()
    save(fig, f"{PREFIX}_trajectory_lines")


def _polarize_and_aggregate(af_df, samples):
    F_cols = [f"F{r}G00" for r in F_REPS if f"F{r}G00" in samples]
    B10_cols = [f"B{r}G10" for r in B_REPS if f"B{r}G10" in samples]
    founder_af = af_df[F_cols].mean(axis=1, skipna=True).values
    b10_af = af_df[B10_cols].mean(axis=1, skipna=True).values
    flip = b10_af < founder_af
    n_kept = (~flip & ~np.isnan(b10_af) & ~np.isnan(founder_af)).sum()
    n_flip = (flip & ~np.isnan(b10_af) & ~np.isnan(founder_af)).sum()
    for s in af_df.columns:
        v = af_df[s].values.copy()
        v[flip] = 1 - v[flip]
        af_df[s] = v

    out = {}
    for trt, reps in [("B", B_REPS), ("T", T_REPS), ("M", M_REPS)]:
        rep_traj = np.full((len(GENS), len(reps)), np.nan)
        for gi, g in enumerate(GENS):
            for ri, r in enumerate(reps):
                key = f"{trt}{r}G{g:02d}"
                if key in af_df.columns:
                    rep_traj[gi, ri] = np.nanmean(af_df[key].values)
        out[trt] = rep_traj
    return out


def _plot_traj_panel(ax, rep_trajs, title, n_snps):
    for trt, col in [("B", COL_B), ("T", COL_T), ("M", COL_M)]:
        rep_traj = rep_trajs[trt]
        g01 = rep_traj[0:1, :]
        delta = rep_traj - g01
        n = np.sum(~np.isnan(delta), axis=1)
        mean_g = np.nanmean(delta, axis=1)
        sd_g = np.nanstd(delta, axis=1, ddof=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            se_g = sd_g / np.sqrt(np.maximum(n, 1))
        ax.fill_between(GENS, mean_g - se_g, mean_g + se_g,
                        color=col, alpha=0.18, edgecolor="none", zorder=2)
        ax.plot(GENS, mean_g, color=col, linewidth=1.5, marker="o",
                markersize=4, zorder=3, label=trt)
    ax.axhline(0, color="#888888", linestyle=":", linewidth=0.5, alpha=0.7,
               zorder=1)
    ax.axvspan(2.6, 5.4, color="#fafafa", zorder=0)
    ax.set_xlabel("generation", fontsize=8)
    ax.set_xticks(GENS)
    ax.set_xlim(0.5, 10.5)
    ax.set_title(f"{title}  (n SNPs = {n_snps})",
                 fontsize=7.5, loc="left", pad=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="both", labelsize=7, length=2)


def fig_af_bt_trajectory(df_p):
    AD_TSV = ROOT / "variance_analysis/merged_ad.tsv"
    SAMPLE_LIST = ROOT / "variance_analysis/sample_list.txt"
    GLMM = ROOT / "glm_lrt_gw_final/glmV1full.csv"
    LRT_THRESHOLD = 8.06
    DEPTH_MIN = 5

    samples = [s.strip() for s in open(SAMPLE_LIST)]

    vote0 = df_p[df_p["votes_v3"] == 0][["chrom", "start"]].copy()
    bg_chroms = set(vote0["chrom"])
    keep_chroms = bg_chroms | {PEAK_CHROM}
    vote0_keys = set(zip(vote0["chrom"], vote0["start"]))

    glmm = pd.read_csv(GLMM)
    qc = ((glmm["converged"] == True) & (glmm["error"] == "OK") &
          (glmm["singular"] == False) & (glmm["LRT_chisq"] >= LRT_THRESHOLD))
    g_keep = glmm[qc & glmm["chrom"].isin(keep_chroms)].copy()
    g_keep["bin"] = (g_keep["pos"] // WIN) * WIN

    g_peak = g_keep[(g_keep["chrom"] == PEAK_CHROM) &
                     (g_keep["bin"] == PEAK_START)]
    g_bg = g_keep[g_keep.apply(
        lambda r: (r["chrom"], r["bin"]) in vote0_keys, axis=1)]
    peak_pos = set(g_peak["pos"])
    bg_pos_by_chrom = g_bg.groupby("chrom")["pos"].apply(set).to_dict()

    cols = ["chrom", "pos", "ref", "alt"] + samples
    chunks = []
    for chunk in pd.read_csv(AD_TSV, sep="\t", header=None, names=cols,
                              chunksize=200_000, dtype=str, low_memory=False):
        chunk = chunk[chunk["chrom"].isin(keep_chroms)].reset_index(drop=True)
        if chunk.empty:
            continue
        chunk["pos"] = chunk["pos"].astype(int)
        keep_mask = np.zeros(len(chunk), dtype=bool)
        chrom_arr = chunk["chrom"].values
        pos_arr = chunk["pos"].values
        keep_mask |= (chrom_arr == PEAK_CHROM) & np.isin(pos_arr, list(peak_pos))
        for chrom, pset in bg_pos_by_chrom.items():
            keep_mask |= (chrom_arr == chrom) & np.isin(pos_arr, list(pset))
        chunk = chunk[keep_mask]
        if not chunk.empty:
            chunks.append(chunk)
    df_ad = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
    if df_ad.empty:
        return

    af = {}
    for s in samples:
        col = df_ad[s].astype(str).str.split(",", expand=True)
        ref = pd.to_numeric(col[0], errors="coerce").values
        alt = pd.to_numeric(col[1], errors="coerce").values
        depth = ref + alt
        af[s] = np.where(depth >= DEPTH_MIN, alt / np.maximum(depth, 1), np.nan)
    af_df = pd.DataFrame(af)
    af_df["chrom"] = df_ad["chrom"].values
    af_df["pos"] = df_ad["pos"].values

    is_peak = (af_df["chrom"] == PEAK_CHROM) & af_df["pos"].isin(peak_pos)
    is_bg = af_df.apply(
        lambda r: r["chrom"] in bg_pos_by_chrom and
                  r["pos"] in bg_pos_by_chrom[r["chrom"]],
        axis=1) & ~is_peak

    af_peak = af_df[is_peak].drop(columns=["chrom", "pos"]).reset_index(drop=True)
    af_bg = af_df[is_bg].drop(columns=["chrom", "pos"]).reset_index(drop=True)

    def _bt_trajectory(af_sub):
        """Returns (mean_g, se_g, n_g) — n_g = number of (B,T) combos used."""
        mean_g = np.full(len(GENS), np.nan)
        se_g = np.full(len(GENS), np.nan)
        n_g = np.zeros(len(GENS), dtype=int)
        for gi, g in enumerate(GENS):
            combos = []
            for br in B_REPS:
                for tr in T_REPS:
                    bk = f"B{br}G{g:02d}"
                    tk = f"T{tr}G{g:02d}"
                    if bk in af_sub.columns and tk in af_sub.columns:
                        d = np.abs(af_sub[bk].values - af_sub[tk].values)
                        d = d[~np.isnan(d)]
                        if len(d):
                            combos.append(d.mean())
            if combos:
                arr = np.array(combos)
                mean_g[gi] = arr.mean()
                se_g[gi] = arr.std(ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 else np.nan
                n_g[gi] = len(arr)
        return mean_g, se_g, n_g

    pk_m, pk_se, pk_n = _bt_trajectory(af_peak)
    bg_m, bg_se, bg_n = _bt_trajectory(af_bg)

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.0), sharey=True)
    for ax, m, se, n_combos, title, n_snps in [
        (axes[0], pk_m, pk_se, pk_n,
         f"chr_439:{PEAK_START/1e6:.1f}-{(PEAK_START+WIN)/1e6:.1f} Mb peak",
         len(af_peak)),
        (axes[1], bg_m, bg_se, bg_n,
         "vote=0 background", len(af_bg)),
    ]:
        ax.fill_between(GENS, m - se, m + se, color="#444444", alpha=0.18,
                        edgecolor="none", zorder=2)
        ax.plot(GENS, m, color="#222222", linewidth=1.5, marker="o",
                markersize=4, zorder=3)
        ax.axvspan(2.6, 5.4, color="#fafafa", zorder=0)
        ax.text(4.0, 0.0, "2-rep gens", fontsize=5.5, ha="center",
                va="bottom", color="#888888",
                transform=ax.get_xaxis_transform())
        ax.set_xlabel("generation", fontsize=8)
        ax.set_xticks(GENS)
        ax.set_xlim(0.5, 10.5)
        ax.set_title(f"{title}  (n SNPs = {n_snps})",
                     fontsize=7.5, loc="left", pad=2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(axis="both", labelsize=7, length=2)
    axes[0].set_ylabel(r"mean $|$AF$_B$ − AF$_T|$" + "\n(across SNPs and B×T combos)",
                       fontsize=8)
    fig.suptitle(f"unpolarized B vs T differentiation at "
                  f"GLMM-LRT>={LRT_THRESHOLD} SNPs",
                  fontsize=8, y=1.0)
    fig.tight_layout()
    save(fig, f"{PREFIX}_af_bt_trajectory")


def fig_af_trajectory(df_p, maf_tol=0.05, depth_tol_rel=0.20,
                       sampling_mode="within_window", suffix=""):
    AD_TSV = ROOT / "variance_analysis/merged_ad.tsv"
    SAMPLE_LIST = ROOT / "variance_analysis/sample_list.txt"
    GLMM = ROOT / "glm_lrt_gw_final/glmV1full.csv"
    PB_THRESHOLD = 0.01         # parametric-bootstrap p cutoff (calibrated)
    DEPTH_MIN = 5
    BG_VOTE = 1
    MAX_BG_SNPS = 30_000   # subsample if huge

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
        keep_mask = np.zeros(len(chunk), dtype=bool)
        chrom_arr = chunk["chrom"].values
        pos_arr = chunk["pos"].values
        peak_mask = (chrom_arr == PEAK_CHROM) & np.isin(pos_arr, list(peak_pos))
        keep_mask |= peak_mask
        for chrom, pset in bg_pos_by_chrom.items():
            mm = (chrom_arr == chrom) & np.isin(pos_arr, list(pset))
            keep_mask |= mm
        chunk = chunk[keep_mask]
        if not chunk.empty:
            chunks.append(chunk)
    df_ad = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()

    if df_ad.empty:
        return


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


    N_RESAMPLES = 1000
    n_target = len(af_peak)
    MAF_TOL = maf_tol
    DEPTH_TOL_REL = depth_tol_rel


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

    if cand.empty:
        bg_win_stats["depth_ok"] = (
            (bg_win_stats["med_depth"] - focal_med_depth).abs() /
            focal_med_depth <= 0.50)
        cand = bg_win_stats[bg_win_stats["maf_ok"] & bg_win_stats["depth_ok"] &
                              bg_win_stats["snps_ok"]].reset_index(drop=True)
    if cand.empty:
        return

    bg_by_window = {}
    for (chrom, bin_), g in af_bg.groupby(["chrom", "bin"]):
        bg_by_window[(chrom, bin_)] = g.index.to_numpy()

    sample_cols = [s for s in samples]
    sample_idx = {s: i for i, s in enumerate(sample_cols)}
    bg_arr = af_bg[sample_cols].values.astype(float)        # (n_bg, n_samples)
    peak_arr = af_peak[sample_cols].values.astype(float)    # (n_peak, n_samples)

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

    TRT_REPS = [("B", B_REPS), ("T", T_REPS), ("M", M_REPS)]

    flip_all_pk, flip_loo_pk = _build_flips(peak_arr)
    rep_trajs_focal = _aggregate_per_trt(peak_arr, flip_all_pk, flip_loo_pk)
    focal_per_rep = {}                      # trt -> (n_gens, n_reps) Δ vs G01
    for trt, _ in TRT_REPS:
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
        if sampling_mode == "across_windows":
            sel = rng.choice(pooled_idx, size=n_target,
                              replace=len(pooled_idx) < n_target)
        else:  # within_window
            wi = int(rng.integers(0, len(cand)))
            win_key = (cand.iloc[wi]["chrom"], cand.iloc[wi]["bin"])
            snp_idx_in_win = bg_by_window[win_key]
            if len(snp_idx_in_win) >= n_target:
                sel = rng.choice(snp_idx_in_win, size=n_target, replace=False)
            else:
                sel = rng.choice(snp_idx_in_win, size=n_target, replace=True)
        af_sub = bg_arr[sel]
        flip_all_s, flip_loo_s = _build_flips(af_sub)
        rep_trajs = _aggregate_per_trt(af_sub, flip_all_s, flip_loo_s)
        for trt, _ in TRT_REPS:
            rep_traj = rep_trajs[trt]
            g01 = rep_traj[0:1, :]
            bg_cloud[trt][i, :] = np.nanmean(rep_traj - g01, axis=1)

    fig, ax = plt.subplots(1, 1, figsize=(3.0, 3.0))
    trt_offsets = {"B": -0.25, "T": 0.0, "M": 0.25}
    cmap = {"B": COL_B, "T": COL_T, "M": COL_M}
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
                   s=1.2, color=col, alpha=0.50, edgecolor="none", zorder=2)
        per_rep = focal_per_rep[trt]
        n = np.sum(~np.isnan(per_rep), axis=1)
        mean_g = np.nanmean(per_rep, axis=1)
        sd_g = np.nanstd(per_rep, axis=1, ddof=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            se_g = sd_g / np.sqrt(np.maximum(n, 1))
        ci = 1.96 * se_g
        x_focal = np.array(GENS) + off
        ax.errorbar(x_focal, mean_g, yerr=ci, fmt="o", color=col,
                    ecolor=col, elinewidth=0.7, capsize=1.5, capthick=0.4,
                    markersize=3, markeredgecolor="black",
                    markeredgewidth=0.3, zorder=5, label=trt)

    ax.axhline(0, color="#888888", linestyle=":", linewidth=0.5, alpha=0.7,
               zorder=1)
    ax.set_xlabel("generation", fontsize=6.5)
    ax.set_ylabel("polarized ΔAF vs G01", fontsize=6.5)
    ax.set_xticks(GENS)
    ax.set_xlim(0.4, 10.6)
    ax.set_title(f"chr_439:{PEAK_START/1e6:.1f}-{(PEAK_START+WIN)/1e6:.1f} Mb "
                  f"(n={n_target}, PB<{PB_THRESHOLD})\n"
                  f"{N_RESAMPLES} bg ({sampling_mode}, "
                  f"MAF±{MAF_TOL}, depth±{int(DEPTH_TOL_REL*100)}%, "
                  f"{len(cand)} wins)",
                  fontsize=5.5, loc="left", pad=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="both", labelsize=5.5, length=1.5)
    ax.legend(loc="best", frameon=False, fontsize=6,
               handletextpad=0.3, borderaxespad=0.2)
    fig.tight_layout()
    save(fig, f"{PREFIX}_af_trajectory{suffix}")






def fig_regional_zoom(df_p, traj):
    pi_F = founder_mean(traj, "theta_pi")
    pi_B10 = per_rep_gen_mean(traj, "B", [10], B_REPS, "theta_pi")[:, 0]
    pi_T10 = per_rep_gen_mean(traj, "T", [10], T_REPS, "theta_pi")[:, 0]

    rB = np.where(pi_F > 0, pi_B10 / pi_F, np.nan)
    rT = np.where(pi_F > 0, pi_T10 / pi_F, np.nan)

    sub = traj[traj["chrom"] == PEAK_CHROM].copy()
    sub_idx = sub.index.values
    pos_mb = (sub["bin"].values + WIN / 2) / 1e6
    rB_chr = rB[sub_idx]
    rT_chr = rT[sub_idx]

    keys = list(zip(traj["chrom"], traj["bin"]))
    df_idx_map = {k: i for i, k in enumerate(keys)}
    df_p = df_p.copy()
    df_p["traj_idx"] = [df_idx_map.get((c, s), -1)
                         for c, s in zip(df_p["chrom"], df_p["start"])]
    null_idx = df_p.loc[(df_p["votes_v3"] == 0) & (df_p["traj_idx"] >= 0),
                          "traj_idx"].values
    null_rB = rB[null_idx]; null_rT = rT[null_idx]

    fig, axes = plt.subplots(2, 1, figsize=(7.0, 4.2), sharex=True)
    for ax, ratio, null_r, title, color in [
        (axes[0], rB_chr, null_rB, "B at G10", BLUE),
        (axes[1], rT_chr, null_rT, "T at G10", ORANGE),
    ]:
        med = np.nanmedian(null_r)
        lo  = np.nanpercentile(null_r, 10)
        hi  = np.nanpercentile(null_r, 90)
        ax.axhspan(lo, hi, color=GREY, alpha=0.30,
                    label=f"genome-wide vote=0 (10-90%)")
        ax.axhline(med, color="#666666", linestyle="--", linewidth=0.6,
                    label="vote=0 median")
        ax.axhline(1.0, color="black", linestyle=":", linewidth=0.5,
                    alpha=0.6)
        ax.plot(pos_mb, ratio, color=color, linewidth=1.0, marker="o",
                 markersize=3, zorder=4)

        ax.axvspan((PEAK_START)/1e6, (PEAK_START + WIN)/1e6,
                    color=RED, alpha=0.15, zorder=0)

        ax.set_title(title, fontsize=8, loc="left")
        ax.set_ylabel(r"$\theta_\pi^{G10} / \theta_\pi^{F}$", fontsize=7)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(axis="both", labelsize=6.5, length=2)
        ax.legend(loc="lower right", frameon=False, fontsize=5.5)
    axes[1].set_xlabel(f"position on {PEAK_CHROM} (Mb)", fontsize=7)
    fig.tight_layout()
    save(fig, f"{PREFIX}_regional_zoom")



def fig_tajD_endpoint(df_p, traj):
    D_B10 = per_rep_gen_mean(traj, "B", [10], B_REPS, "tajimas_d")[:, 0]
    D_T10 = per_rep_gen_mean(traj, "T", [10], T_REPS, "tajimas_d")[:, 0]
    D_BT = np.nanmean(np.stack([D_B10, D_T10]), axis=0)
    keys = list(zip(traj["chrom"], traj["bin"]))
    map_B  = dict(zip(keys, D_B10))
    map_T  = dict(zip(keys, D_T10))
    map_BT = dict(zip(keys, D_BT))
    df_p = df_p.copy()
    df_p["D_B_G10"]  = [map_B.get((c, s), np.nan)  for c, s in zip(df_p["chrom"], df_p["start"])]
    df_p["D_T_G10"]  = [map_T.get((c, s), np.nan)  for c, s in zip(df_p["chrom"], df_p["start"])]
    df_p["D_BT_G10"] = [map_BT.get((c, s), np.nan) for c, s in zip(df_p["chrom"], df_p["start"])]
    _pointcloud_three(df_p, ["D_B_G10", "D_T_G10", "D_BT_G10"],
                      ["B", "T", "B+T"],
                      r"Tajima's $D$ at G10",
                      "low",
                      f"{PREFIX}_tajD_endpoint")




def fig_perm_fwer(df_p, traj):
    pi_F = founder_mean(traj, "theta_pi")
    pi_B10 = per_rep_gen_mean(traj, "B", [10], B_REPS, "theta_pi")[:, 0]
    pi_T10 = per_rep_gen_mean(traj, "T", [10], T_REPS, "theta_pi")[:, 0]
    ratio_B = np.where(pi_F > 0, pi_B10 / pi_F, np.nan)
    ratio_T = np.where(pi_F > 0, pi_T10 / pi_F, np.nan)
    keys = list(zip(traj["chrom"], traj["bin"]))
    map_B  = dict(zip(keys, ratio_B))
    map_T  = dict(zip(keys, ratio_T))
    df_p = df_p.copy()
    df_p["pi_ratio_B"] = [map_B.get((c, s), np.nan) for c, s in zip(df_p["chrom"], df_p["start"])]
    df_p["pi_ratio_T"] = [map_T.get((c, s), np.nan) for c, s in zip(df_p["chrom"], df_p["start"])]


    tracks = [
        ("pi_ratio_B",      "B θπ G10/F"),
        ("pi_ratio_T",      "T θπ G10/F"),
        ("B_thetaPi_slope", "B θπ slope"),
        ("T_thetaPi_slope", "T θπ slope"),
    ]
    rng = np.random.default_rng(42)
    rows = []
    chrom_groups = []
    for chrom, g in df_p.groupby("chrom", sort=False):
        chrom_groups.append((chrom, g.index.to_numpy()))

    peak = df_p[(df_p["chrom"] == PEAK_CHROM) & (df_p["start"] == PEAK_START)]
    for col, label in tracks:
        s = df_p[col].fillna(np.nanmedian(df_p[col])).values

        peak_v = float(peak[col].iloc[0]) if len(peak) and peak[col].notna().any() else np.nan
        if np.isnan(peak_v):
            rows.append((label, np.nan, np.nan, np.nan))
            continue
        null_min = np.empty(N_PERM)
        for i in range(N_PERM):
            sp = np.empty_like(s)
            for chrom, idx in chrom_groups:
                off = int(rng.integers(0, len(idx)))
                sp[idx] = np.roll(s[idx], off)
            null_min[i] = sp.min()

        p_fwer = ((null_min <= peak_v).sum() + 1) / (N_PERM + 1)
        thr05 = float(np.quantile(null_min, 0.05))
        thr01 = float(np.quantile(null_min, 0.01))
        rows.append((label, peak_v, p_fwer, thr05))

    fig, ax = plt.subplots(figsize=(5.0, 0.5 + 0.35*len(rows)))
    ax.axis("off")
    cells = [["track", "peak value", "FWER p (low-tail, 5000 perm)", "5% null threshold"]]
    for label, v, p, thr in rows:
        cells.append([label, f"{v:.4f}",
                      "NA" if np.isnan(p) else (f"{p:.4f}" if p >= 1e-3 else "<0.001"),
                      "NA" if np.isnan(thr) else f"{thr:.4f}"])
    tbl = ax.table(cellText=cells, loc="center", cellLoc="left")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.scale(1.0, 1.3)
    for j in range(len(cells[0])):
        tbl[(0, j)].set_text_props(weight="bold")
        tbl[(0, j)].set_facecolor("#eeeeee")
    fig.tight_layout()
    save(fig, f"{PREFIX}_perm_fwer")





def fig_joint_stouffer(df_p):
    cols = ["B_thetaPi_slope", "B_thetaW_slope", "B_tajD_slope"]
    sub = df_p[df_p[cols].notna().all(axis=1)].copy()
    z = np.zeros((len(sub), len(cols)))
    for j, c in enumerate(cols):
        v = sub[c].values
        z[:, j] = -(v - v.mean()) / v.std(ddof=1)
    z_combined = z.sum(axis=1) / np.sqrt(len(cols))
    sub["z_diversity_combined"] = z_combined


    null_vals = sub.loc[sub["votes_v3"] == 0, "z_diversity_combined"].values
    peak = sub[(sub["chrom"] == PEAK_CHROM) & (sub["start"] == PEAK_START)]
    peak_v = float(peak["z_diversity_combined"].iloc[0]) if len(peak) else np.nan
    p, n = emp_p(peak_v, null_vals, "high")  # high z = candidate
    base = np.nanmean(null_vals)

    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(1, 2, figsize=(3.6, 2.6))
    bg = sub[sub["votes_v3"] < 3]
    x_jit = 1 + rng.uniform(-0.18, 0.18, len(bg))
    axes[0].scatter(x_jit, bg["z_diversity_combined"], s=4, color=GREY,
                    alpha=0.55, edgecolor="none", zorder=1)
    if not np.isnan(peak_v):
        axes[0].scatter([1], [peak_v], s=24, color=RED, edgecolor="black",
                        linewidth=0.4, zorder=5)
    axes[0].axhline(base, color="#555555", linestyle="--", linewidth=0.5,
                    alpha=0.7)
    axes[0].set_xlim(0.4, 1.6); axes[0].set_xticks([])
    axes[0].set_ylabel("composite z\n(B θπ + θW + D slopes)", fontsize=7)
    axes[0].set_title(f"raw\n{fmt_p(p)} (n_null={n}, high-tail)",
                      fontsize=6.0, loc="left", pad=2)

    dev = bg["z_diversity_combined"] - base
    axes[1].scatter(x_jit, dev, s=4, color=GREY, alpha=0.55,
                    edgecolor="none", zorder=1)
    if not np.isnan(peak_v):
        axes[1].scatter([1], [peak_v - base], s=24, color=RED,
                        edgecolor="black", linewidth=0.4, zorder=5)
    axes[1].axhline(0, color="#555555", linestyle="--", linewidth=0.5,
                    alpha=0.7)
    axes[1].set_xlim(0.4, 1.6); axes[1].set_xticks([])
    axes[1].set_title("Δ vs vote=0", fontsize=6.5, loc="left", pad=2)

    for ax in axes:
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(axis="both", labelsize=5.5, length=2)
    fig.tight_layout()
    save(fig, f"{PREFIX}_joint_stouffer")




def _pointcloud_three(df_p, cols, sublabels, ylab_top, direction, base_path):
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(1, len(cols) * 2,
                             figsize=(len(cols) * 2 * 1.6, 2.6))
    peak = df_p[(df_p["chrom"] == PEAK_CHROM) & (df_p["start"] == PEAK_START)]
    for s_i, (col, sublabel) in enumerate(zip(cols, sublabels)):
        ax_raw = axes[s_i * 2]
        ax_dev = axes[s_i * 2 + 1]
        sub = df_p[df_p[col].notna()].copy()
        null_vals = sub.loc[sub["votes_v3"] == 0, col].dropna().values
        base = null_vals.mean() if len(null_vals) else np.nan
        peak_v = float(peak[col].iloc[0]) if len(peak) and peak[col].notna().any() else np.nan
        p, n = emp_p(peak_v, null_vals, direction)
        stat_str = f"{fmt_p(p)} (n_null={n})"

        x_jit = 1 + rng.uniform(-0.18, 0.18, len(sub))
        ax_raw.scatter(x_jit, sub[col], s=3, color=GREY, alpha=0.55,
                       edgecolor="none", zorder=1)
        if not np.isnan(peak_v):
            ax_raw.scatter([1], [peak_v], s=22, color=RED,
                           edgecolor="black", linewidth=0.4, zorder=5)
        if not np.isnan(base):
            ax_raw.axhline(base, color="#555555", linestyle="--",
                           linewidth=0.5, alpha=0.7)

        dev = sub[col] - base
        ax_dev.scatter(x_jit, dev, s=3, color=GREY, alpha=0.55,
                       edgecolor="none", zorder=1)
        if not np.isnan(peak_v) and not np.isnan(base):
            ax_dev.scatter([1], [peak_v - base], s=22, color=RED,
                           edgecolor="black", linewidth=0.4, zorder=5)
        ax_dev.axhline(0, color="#555555", linestyle="--",
                       linewidth=0.5, alpha=0.7)

        for ax in (ax_raw, ax_dev):
            ax.set_xlim(0.4, 1.6); ax.set_xticks([])
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            ax.tick_params(axis="both", labelsize=5.5, length=2)

        if s_i == 0:
            ax_raw.set_ylabel(ylab_top, fontsize=7)
        ax_raw.set_title(f"{sublabel} · raw\n{stat_str}",
                         fontsize=6.0, loc="left", pad=2)
        ax_dev.set_title(f"{sublabel} · Δ vs vote=0",
                         fontsize=6.0, loc="left", pad=2)
    fig.tight_layout()
    save(fig, base_path)



def main() -> None:
    setup_rc()
    df_p, keep = load_master_with_votes()
    traj = load_trajectory(keep)

    fig_endpoint_ratio(df_p, traj)
    fig_trajectory(df_p, traj)
    fig_regional_zoom(df_p, traj)
    fig_tajD_endpoint(df_p, traj)
    fig_perm_fwer(df_p, traj)
    fig_joint_stouffer(df_p)


if __name__ == "__main__":
    main()
