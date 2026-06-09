#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
from matplotlib.gridspec import GridSpec

ROOT = Path(".")
MASTER_TSV = ROOT / "final_plots/wild/section2_candidate_master_v2.tsv"
HV_BLOCKS_TSV = ROOT / "final_plots/wild/section2_hv_blocks_filtered.tsv"
GLMM_CSV = ROOT / "glm_lrt_gw_final/glmV1full.csv"
GLMM_PER_WINDOW_TSV = ROOT / "final_plots/wild/section2_glmm_per_window.tsv"
OUT_BASE = ROOT / "final_plots/wild/section2_candidate_manhattan_v3_combined"

WIN = 200_000
GAP = 4_000_000
ROLL_RADIUS = 1            # ±1 windows (3-window rolling sum, 600 kb total)
Z_CAP = 5.0                # y-axis cap for per-test z-score Manhattans
Z_SIG = 1.645              # one-sided z corresponding to p = 0.05
SCAFF_MIN_WINDOWS = 50
TOP_FRAC = 0.05            # threshold for the new IND_glmm_lrt indicator

COL_B  = "#499FFF"
COL_T  = "#EDB72D"
COL_M  = "#B01754"
COL_BLACK = "#222222"
HIGHLIGHT_COLOR = "#C84A45"
COMPARE_COLOR  = "#6E8E3D"   # for Fig 3E Test 3 comparison regions

COMPARE_REGIONS = [
    ("chr_ScDA7r2_597_HRSCAF_953", 3_600_000, 4_000_000, "chr_597:3.7"),
    ("chr_ScDA7r2_110_HRSCAF_295", 5_000_000, 5_200_000, "chr_110:5.1"),
    ("chr_ScDA7r2_439_HRSCAF_779", 7_800_000, 8_000_000, "chr_439:7.9"),
    ("chr_ScDA7r2_439_HRSCAF_779", 16_400_000, 16_800_000, "chr_439:16.6"),
]

SEGMENT_TRACKS = [
    ("hv_B", "HV blocks (Barbarea-evolved)",     COL_B),
    ("hv_T", "HV blocks (Turritis-evolved)",     COL_T),
    ("hv_M", "HV blocks (B+T mixture-evolved)",  COL_M),
]
STRIPE_TRACKS = [
    ("cov_BT",        "Lab cov(B,T) — antagonistic",    "low"),
    ("slope_div_z",   "Lab BvT temporal FST slope",     "high"),
    ("permFST_BvT",   "Lab BvT G10 FST",                "high"),
    ("glmm_lrt",      "Lab GLMM gen×trt LRT",           "high"),
    ("wild_C2_max",   "Wild BayPass C2",                "high"),
]


def load_glmm_per_window(force: bool = False) -> pd.DataFrame:
    if GLMM_PER_WINDOW_TSV.exists() and not force:
        return pd.read_csv(GLMM_PER_WINDOW_TSV, sep="\t")
    g = pd.read_csv(GLMM_CSV,
                     usecols=["chrom", "pos", "LRT_chisq", "PB_p_val",
                                "converged", "error"])
    qc = (g["converged"] == True) & (g["error"] == "OK") & \
         g["PB_p_val"].notna() & g["LRT_chisq"].notna()
    g = g[qc].copy()
    g["start"] = (g["pos"] // WIN) * WIN
    agg = (g.groupby(["chrom", "start"])
              .agg(glmm_lrt=("LRT_chisq", "mean"),
                    glmm_n_snps=("LRT_chisq", "size"))
              .reset_index())
    agg.to_csv(GLMM_PER_WINDOW_TSV, sep="\t", index=False)
    return agg


def emp_neglogp(scores, direction: str) -> np.ndarray:
    s = np.asarray(scores, dtype=float)
    s = np.nan_to_num(s, nan=0.0)
    N = len(s)
    floor = 1.0 / (2.0 * N)
    if direction == "p":
        return -np.log10(np.clip(s, floor, 1.0))
    sorted_s = np.sort(s)
    if direction == "high":
        counts = N - np.searchsorted(sorted_s, s, side="left")
    elif direction == "low":
        counts = np.searchsorted(sorted_s, s, side="right")
    else:
        raise ValueError(direction)
    p = np.clip(counts / N, floor, 1.0)
    return -np.log10(p)


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



def standardize_test(scores: np.ndarray, direction: str) -> np.ndarray:
    s = np.asarray(scores, dtype=float)
    if direction == "low":
        s = -s
    elif direction == "p":
        s = stats.norm.ppf(1 - np.clip(s, 1e-10, 1 - 1e-10))
    elif direction != "high":
        raise ValueError(direction)
    valid = np.isfinite(s)
    if valid.sum() == 0:
        return np.full_like(s, np.nan, dtype=float)
    mu = np.nanmean(s[valid])
    sd = np.nanstd(s[valid], ddof=1)
    z = np.full_like(s, np.nan, dtype=float)
    z[valid] = (s[valid] - mu) / sd
    return z


def rolling_mean_per_scaffold(z: np.ndarray, df_p: pd.DataFrame,
                                radius: int) -> np.ndarray:
    out = np.copy(z).astype(float)
    df_p = df_p.reset_index(drop=True)
    for chrom, g in df_p.groupby("chrom", sort=False):
        idx = g.index.to_numpy()
        sub = pd.Series(z[idx])
        out[idx] = (sub.rolling(window=2 * radius + 1, center=True,
                                  min_periods=1)
                       .mean().to_numpy())
    return out


def stouffer_combine(z_smoothed: np.ndarray) -> np.ndarray:
    valid = ~np.isnan(z_smoothed)
    n_valid = valid.sum(axis=1).astype(float)
    z_sum = np.nansum(z_smoothed, axis=1)
    out = np.full(z_smoothed.shape[0], np.nan)
    nz = n_valid > 0
    out[nz] = z_sum[nz] / np.sqrt(n_valid[nz])
    return out


def permutation_null_max(z_unsmoothed: np.ndarray, df_p: pd.DataFrame,
                            radius: int, n_perm: int, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n_windows, n_tests = z_unsmoothed.shape
    null_max = np.empty(n_perm)
    z = np.nan_to_num(z_unsmoothed, nan=0.0)
    for i in range(n_perm):
        z_perm = np.empty_like(z)
        for t in range(n_tests):
            z_perm[:, t] = z[rng.permutation(n_windows), t]
        z_smoothed = np.empty_like(z_perm)
        for t in range(n_tests):
            z_smoothed[:, t] = rolling_mean_per_scaffold(z_perm[:, t], df_p, radius)
        comb = stouffer_combine(z_smoothed)
        null_max[i] = np.nanmax(comb)
    return null_max


def main() -> None:
    plt.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.7,
    })

    df = pd.read_csv(MASTER_TSV, sep="\t")
    blocks = pd.read_csv(HV_BLOCKS_TSV, sep="\t")

    glmm = load_glmm_per_window()
    df = df.merge(glmm[["chrom", "start", "glmm_lrt"]],
                   on=["chrom", "start"], how="left")

    chrom_sizes = df.groupby("chrom").size().sort_values(ascending=False)
    keep = chrom_sizes[chrom_sizes >= SCAFF_MIN_WINDOWS].index.tolist()

    offsets = {}
    cum = 0
    for s in keep:
        offsets[s] = cum
        slen = df.loc[df["chrom"] == s, "end"].max()
        cum += slen + GAP
    df_p = df[df["chrom"].isin(keep)].copy().reset_index(drop=True)
    df_p["x"] = df_p["chrom"].map(offsets) + df_p["start"]
    n_windows = len(df_p)

    test_specs = []
    for key, _, _ in SEGMENT_TRACKS:
        treat = {"hv_B": "B", "hv_T": "T", "hv_M": "M"}[key]
        s = hv_score_per_window(df_p, blocks, treat)
        df_p[f"raw_{key}"] = s
        test_specs.append((key, f"raw_{key}", "high"))
    for col, _, direction in STRIPE_TRACKS:
        df_p[f"raw_{col}"] = df_p[col].astype(float).values
        test_specs.append((col, f"raw_{col}", direction))

    for name, col, dirn in test_specs:
        df_p[f"nlp_{name}"] = emp_neglogp(df_p[col].values, dirn)
        df_p[f"z_{name}"] = standardize_test(df_p[col].values, dirn)


    glmm_thr = df_p["raw_glmm_lrt"].quantile(1 - TOP_FRAC)
    df_p["IND_glmm_lrt"] = (df_p["raw_glmm_lrt"] >= glmm_thr).fillna(False).astype(int)

    IND_COLS = [
        "IND_cov_BT_neg", "IND_slope_div",
        "IND_permFST", "IND_glmm_lrt", "IND_wild_C2",
    ]
    n_tests = len(IND_COLS)
    df_p["votes_pw"] = df_p[IND_COLS].fillna(0).astype(int).sum(axis=1)

    df_p["votes_total"] = df_p["votes_pw"].astype(int)
    max_possible = n_tests

    n_seg = len(SEGMENT_TRACKS)
    n_str = len(STRIPE_TRACKS)
    n_panels = 1 + n_seg + n_str

    fig = plt.figure(figsize=(12, 3))
    gs = GridSpec(n_panels, 1, figure=fig,
                  height_ratios=([1.10] +
                                  [0.18] * n_seg +
                                  [0.55] * n_str),
                  hspace=0.45)

    xmin = df_p["x"].min() - WIN
    xmax = df_p["x"].max() + WIN

    ax_top = fig.add_subplot(gs[0])
    headline_idx = int(df_p["votes_total"].idxmax())
    headline_row = df_p.iloc[headline_idx]
    HEADLINE_X = headline_row["x"] + WIN / 2
    HEADLINE_HALFWIDTH = 5 * WIN
    headline_total = int(headline_row["votes_total"])

    ax_top.axvspan(HEADLINE_X - HEADLINE_HALFWIDTH,
                    HEADLINE_X + HEADLINE_HALFWIDTH,
                    color=HIGHLIGHT_COLOR, alpha=0.07, zorder=0)

    nonzero = df_p[df_p["votes_total"] >= 1].copy()
    for _, r in nonzero.iterrows():
        x = r["x"] + WIN / 2
        v = int(r["votes_total"])
        is_peak = (r["chrom"] == headline_row["chrom"]
                    and r["start"] == headline_row["start"])

        if is_peak:
            col, a, lw = HIGHLIGHT_COLOR, 1.0, 1.1
        elif v >= 3:
            col, a, lw = COL_BLACK, 1.0, 0.9
        else:
            col, a, lw = "#888888", 0.6, 0.5
        ax_top.vlines(x, 0, v, color=col, alpha=a, linewidth=lw, zorder=3)
        ax_top.scatter(x, v, s=10 + 6*v, color=col, alpha=a,
                       edgecolor="none", zorder=4)

    sc = headline_row["chrom"].split("_")[2] if "_" in headline_row["chrom"] else headline_row["chrom"]
    ax_top.annotate(f"chr_{sc}:{headline_row['start']/1e6:.1f} Mb "
                    f"({headline_total}/{n_tests})",
                    xy=(HEADLINE_X, headline_total),
                    xytext=(0, 5), textcoords="offset points",
                    fontsize=9, ha="center", va="bottom",
                    color=HIGHLIGHT_COLOR, fontweight="bold")

    for chrom_c, start_c, end_c, label_c in COMPARE_REGIONS:
        if chrom_c not in offsets:
            continue
        x_center = offsets[chrom_c] + (start_c + end_c) / 2
        in_reg = (df_p["chrom"] == chrom_c) & \
                  (df_p["start"] >= start_c) & (df_p["start"] < end_c)
        votes_at_reg = int(df_p.loc[in_reg, "votes_total"].max()) \
                          if in_reg.any() else 0
        y_top = headline_total + 0.6
        ax_top.annotate(label_c,
                          xy=(x_center, votes_at_reg + 0.05),
                          xytext=(x_center, y_top),
                          fontsize=6.5, ha="center", va="bottom",
                          color=COMPARE_COLOR,
                          arrowprops=dict(arrowstyle="-|>",
                                            color=COMPARE_COLOR,
                                            shrinkA=0, shrinkB=2,
                                            lw=0.7, alpha=0.9))

    ax_top.set_ylabel(f"Votes\n(of {n_tests} tests)",
                       fontsize=8, labelpad=8)
    ax_top.set_xlim(xmin, xmax)
    ax_top.set_ylim(0, headline_total + 2.4)
    yt_max = headline_total + 1
    yt_step = max(1, yt_max // 4) if yt_max > 4 else 1
    ax_top.set_yticks(range(0, yt_max + 1, yt_step))
    ax_top.set_xticks([])
    ax_top.tick_params(axis="x", which="both", bottom=False, top=False,
                        labelbottom=False)
    plt.setp(ax_top.get_xticklabels(), visible=False)
    for sp in ("top", "right", "bottom"):
        ax_top.spines[sp].set_visible(False)
    ax_top.tick_params(axis="y", labelsize=7)

    seg_axes = []
    for ti, (key, label, colour) in enumerate(SEGMENT_TRACKS):
        ax = fig.add_subplot(gs[1 + ti], sharex=ax_top)
        seg_axes.append(ax)
        for i, s in enumerate(keep):
            xs = df_p.loc[df_p["chrom"] == s, "x"]
            if i % 2 == 1 and len(xs):
                ax.axvspan(xs.min() - WIN/2, xs.max() + WIN/2,
                           color="#F2F2F2", zorder=0)
        treat = {"hv_B": "B", "hv_T": "T", "hv_M": "M"}[key]
        blk_t = blocks[blocks["treatment"] == treat]
        if not blk_t.empty:
            log_n = np.log10(blk_t["n_snps"].astype(float))
            spread = max(log_n.max() - log_n.min(), 0.1)
            norm = ((log_n - log_n.min()) / spread) ** 2
            blk_t = blk_t.assign(_alpha=0.20 + 0.80 * norm)
        for _, blk in blk_t.iterrows():
            if blk["chr"] not in offsets:
                continue
            x0 = offsets[blk["chr"]] + blk["start"]
            x1 = offsets[blk["chr"]] + blk["end"]
            ax.fill_betweenx([-0.4, 0.4], x0, x1, color=colour,
                             alpha=float(blk["_alpha"]),
                             edgecolor="none", zorder=3)
        ax.text(1.005, 0.5, label, transform=ax.transAxes,
                ha="left", va="center", fontsize=7, color=colour)
        ax.set_yticks([]); ax.set_ylim(-0.5, 0.5)
        ax.tick_params(axis="x", which="both", bottom=False, top=False,
                        labelbottom=False)
        plt.setp(ax.get_xticklabels(), visible=False)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_color("#DDDDDD")


    str_axes = []
    for ti, (col, label, direction) in enumerate(STRIPE_TRACKS):
        ax = fig.add_subplot(gs[1 + n_seg + ti], sharex=ax_top)
        str_axes.append(ax)
        for i, s in enumerate(keep):
            xs = df_p.loc[df_p["chrom"] == s, "x"]
            if i % 2 == 1 and len(xs):
                ax.axvspan(xs.min() - WIN/2, xs.max() + WIN/2,
                           color="#F2F2F2", zorder=0)

        z = df_p[f"z_{col}"].values.astype(float)
        z = np.nan_to_num(z, nan=0.0)
        x = df_p["x"].values + WIN / 2

        z_show = np.maximum(z, 0)

        high_vote = df_p["votes_pw"].values >= 3
        if (~high_vote).any():
            ax.vlines(x[~high_vote], 0, np.minimum(z_show[~high_vote], Z_CAP),
                      color="#888888", alpha=0.55, linewidth=0.5, zorder=2)
        if high_vote.any():
            ax.vlines(x[high_vote], 0, np.minimum(z_show[high_vote], Z_CAP),
                      color=COL_BLACK, alpha=1.0, linewidth=0.9, zorder=3)

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(0, Z_CAP * 1.05)
        ax.set_yticks([0, Z_CAP])
        ax.set_yticklabels(["0", f"{Z_CAP:.0f}"], fontsize=6)
        ax.tick_params(axis="y", labelsize=6, length=2, pad=1)
        ax.text(1.005, 0.5, label, transform=ax.transAxes,
                ha="left", va="center", fontsize=7, color=COL_BLACK)
        ax.tick_params(axis="x", which="both", bottom=False, top=False,
                        labelbottom=False)
        plt.setp(ax.get_xticklabels(), visible=False)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.spines["left"].set_color("#BBBBBB")
        ax.spines["bottom"].set_color("#BBBBBB")

    fig.text(0.05, 0.30,
             "Per-test z-score (standardized within track; positive tail only)",
             rotation=90, ha="center", va="center", fontsize=8)

    last_ax = str_axes[-1]
    tick_x = []
    for s in keep:
        xs = df_p.loc[df_p["chrom"] == s, "x"]
        tick_x.append(xs.mean())
    short_labels = [s.split("_")[2] if "_" in s else s for s in keep]
    last_ax.set_xticks(tick_x)
    last_ax.set_xticklabels(short_labels, fontsize=8)
    last_ax.tick_params(axis="x", which="both", bottom=True,
                        labelbottom=True, length=2)
    last_ax.set_xlabel("Scaffold", fontsize=8.5)

    fig.savefig(f"{OUT_BASE}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.svg", bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.pdf", bbox_inches="tight")
    svg = Path(f"{OUT_BASE}.svg")
    txt = svg.read_text()
    txt = re.sub(r"<clipPath[^>]*>.*?</clipPath>", "", txt, flags=re.DOTALL)
    txt = re.sub(r'\s*clip-path="url\([^)]+\)"', "", txt)
    svg.write_text(txt)


if __name__ == "__main__":
    main()
