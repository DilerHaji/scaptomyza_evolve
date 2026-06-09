#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

ROOT = Path(".")
MASTER_TSV = ROOT / "final_plots/wild/section2_candidate_master_v2.tsv"
GLMM_PER_WINDOW_TSV = ROOT / "final_plots/wild/section2_glmm_per_window.tsv"
OUT_BASE = ROOT / "final_plots/wild/section2_vote_sweep"

WIN = 200_000
GAP = 4_000_000
SCAFF_MIN_WINDOWS = 50
TOP_FRAC = 0.05

WINDOW_SIZES = [1, 3, 5, 11, 25]      # in 200 kb bins → 0.2/0.6/1.0/2.2/5.0 Mb
STRIDE_MODES = [
    ("stride = 1 (full overlap)", lambda w: 1),
    ("stride = win/2 (half overlap)", lambda w: max(1, w // 2)),
    ("stride = win (non-overlapping)", lambda w: max(1, w)),
]

COL_BLACK = "#222222"
HIGHLIGHT_COLOR = "#C84A45"
ZOOM_CHROM = "chr_ScDA7r2_439_HRSCAF_779"


def rolling_total_votes(df_p: pd.DataFrame, ind_cols, win_bins: int,
                          stride: int) -> pd.DataFrame:
    out_rows = []
    df_p = df_p.reset_index(drop=True)
    df_p = df_p.copy()
    df_p["_pw_votes"] = df_p[ind_cols].fillna(0).astype(int).sum(axis=1)
    for chrom, g in df_p.groupby("chrom", sort=False):
        g = g.sort_values("start").reset_index(drop=True)
        n = len(g)
        if n == 0:
            continue
        s = (g["_pw_votes"]
                .rolling(window=win_bins, center=True, min_periods=1).sum()
                .to_numpy().astype(int))
        idx = np.arange(0, n, stride)
        for i in idx:
            x = g.iloc[i]["x"] + WIN / 2
            out_rows.append({"chrom": chrom, "x": x, "total_votes": int(s[i])})
    return pd.DataFrame(out_rows)


def main() -> None:
    plt.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.6,
    })

    df = pd.read_csv(MASTER_TSV, sep="\t")
    glmm = pd.read_csv(GLMM_PER_WINDOW_TSV, sep="\t")
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

    glmm_thr = df_p["glmm_lrt"].quantile(1 - TOP_FRAC)
    df_p["IND_glmm_lrt"] = (df_p["glmm_lrt"] >= glmm_thr).fillna(False).astype(int)

    IND_COLS = [
        "IND_HV_blocks_B", "IND_HV_blocks_T", "IND_HV_blocks_M",
        "IND_cvtk_antag", "IND_cov_BT_neg", "IND_slope_div",
        "IND_permFST", "IND_glmm_lrt", "IND_wild_C2",
    ]
    n_tests = len(IND_COLS)

    xmin = df_p["x"].min() - WIN
    xmax = df_p["x"].max() + WIN

    headline_chrom = ZOOM_CHROM
    HEADLINE_HALFWIDTH = 5 * WIN

    chr439 = df_p[df_p["chrom"] == headline_chrom]
    peak_pw = chr439.loc[chr439[IND_COLS].sum(axis=1).idxmax()]
    HEADLINE_X = peak_pw["x"] + WIN / 2

    n_rows = len(WINDOW_SIZES)
    n_cols = len(STRIDE_MODES)
    fig = plt.figure(figsize=(12, 9))
    gs = GridSpec(n_rows, n_cols, figure=fig, hspace=0.55, wspace=0.18,
                   left=0.06, right=0.985, top=0.94, bottom=0.05)

    all_traces = {}
    for ri, wbins in enumerate(WINDOW_SIZES):
        for ci, (stride_label, stride_fn) in enumerate(STRIDE_MODES):
            stride = stride_fn(wbins)
            tr = rolling_total_votes(df_p, IND_COLS, wbins, stride)
            all_traces[(ri, ci)] = (tr, stride)
            max_v = tr["total_votes"].max() if len(tr) else 0
    global_ymax = max(t["total_votes"].max() for (t, _) in all_traces.values())

    for ri, wbins in enumerate(WINDOW_SIZES):
        for ci, (stride_label, stride_fn) in enumerate(STRIDE_MODES):
            ax = fig.add_subplot(gs[ri, ci])
            tr, stride = all_traces[(ri, ci)]
            max_possible = wbins * n_tests

            for i, s in enumerate(keep):
                xs = df_p.loc[df_p["chrom"] == s, "x"]
                if i % 2 == 1 and len(xs):
                    ax.axvspan(xs.min() - WIN/2, xs.max() + WIN/2,
                                color="#F5F5F5", zorder=0)

            ax.axvspan(HEADLINE_X - HEADLINE_HALFWIDTH,
                        HEADLINE_X + HEADLINE_HALFWIDTH,
                        color=HIGHLIGHT_COLOR, alpha=0.10, zorder=0.5)

            for chrom in keep:
                t_chr = tr[tr["chrom"] == chrom].sort_values("x")
                if t_chr.empty:
                    continue
                x = t_chr["x"].values
                y = t_chr["total_votes"].values
                ax.fill_between(x, 0, y, color=COL_BLACK, alpha=0.20,
                                  step="mid" if stride > 1 else None,
                                  linewidth=0, zorder=2)
                ax.plot(x, y, color=COL_BLACK, linewidth=0.6,
                         drawstyle="steps-mid" if stride > 1 else "default",
                         zorder=3)

                if chrom == headline_chrom:
                    in_band = (x >= HEADLINE_X - HEADLINE_HALFWIDTH) & \
                              (x <= HEADLINE_X + HEADLINE_HALFWIDTH)
                    if in_band.any():
                        ax.fill_between(x[in_band], 0, y[in_band],
                                          color=HIGHLIGHT_COLOR, alpha=0.55,
                                          step="mid" if stride > 1 else None,
                                          linewidth=0, zorder=4)
                        ax.plot(x[in_band], y[in_band], color=HIGHLIGHT_COLOR,
                                 linewidth=1.0,
                                 drawstyle="steps-mid" if stride > 1 else "default",
                                 zorder=5)

            ax.set_xlim(xmin, xmax)
            ax.set_ylim(0, global_ymax * 1.05)
            ax.set_yticks([0, max_possible])
            ax.set_yticklabels(["0", f"{max_possible}"], fontsize=6)

            chr439_max = (tr[tr["chrom"] == headline_chrom]["total_votes"].max()
                            if (tr["chrom"] == headline_chrom).any() else 0)
            global_max = tr["total_votes"].max()
            ax.set_title(f"win={wbins} bins ({wbins*WIN/1e6:.1f} Mb), stride={stride}\n"
                          f"chr_439 peak={int(chr439_max)} / "
                          f"genome-wide max={int(global_max)} / "
                          f"max possible={max_possible}",
                          fontsize=7, loc="left", pad=3)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            ax.spines["left"].set_color("#BBBBBB")
            ax.spines["bottom"].set_color("#BBBBBB")

            if ri == n_rows - 1:
                tick_x = [df_p.loc[df_p["chrom"] == s, "x"].mean() for s in keep]
                short_labels = [s.split("_")[2] if "_" in s else s for s in keep]
                ax.set_xticks(tick_x)
                ax.set_xticklabels(short_labels, fontsize=6)
                ax.tick_params(axis="x", length=2)
            else:
                ax.set_xticks([])

            if ri == 0:
                ax.text(0.5, 1.30, stride_label, transform=ax.transAxes,
                         ha="center", va="bottom", fontsize=9,
                         fontweight="bold", color="#222222")

    fig.text(0.005, 0.5, "Total votes (rolling sum)",
              rotation=90, ha="left", va="center", fontsize=9)
    fig.text(0.5, 0.005, "Scaffold", ha="center", va="bottom", fontsize=9)

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
