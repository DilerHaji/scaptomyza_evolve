#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_s3_fig1_wild_diversity import (  # type: ignore
    DIV_FILE,
    POOL_ORDER,
    POOL_LABELS,
    GROUP_OF,
    SITE_OF,
    METRICS,
    C_B,
    C_T,
    C_FOUNDER,
    C_BT,
    C_GREY_MID,
    C_GREY_LIGHT,
    _adj_positions,
    _color_of,
    load_wild_long,
)

ROOT = Path(".")
OUT_DIR = ROOT / "final_plots/wild"
G_DIR = ROOT / "variance_analysis/cvtkpy_results"


def load_g_k2() -> dict[str, float]:
    out: dict[str, float] = {}
    for treat in ("B", "T", "M"):
        path = G_DIR / f"{treat}_G.tsv"
        df = pd.read_csv(path, sep="\t")
        df = df[df["replicate"].astype(str).str.startswith("rep")]
        for _, row in df.iterrows():
            rep = int(str(row["replicate"]).replace("rep", ""))
            out[f"{treat}{rep}G10"] = float(row["G_k2_validation"])
    return out

CLOUD_COLOR = "#D9D9D9"   # light grey
CLOUD_ALPHA = 0.30
CLOUD_SIZE = 2.0           # cloud points (small but visible)
JITTER_HALF_WIDTH = 0.22
POINT_SIZE = 9             # median marker — small to read against the cloud
ERR_LW = 1.0
RNG = np.random.default_rng(42)
SAVE_DPI = 600

Y_CONFIG: dict[str, dict] = {
    "theta_pi": {
        "broken": True,
        "bot_ylim": (0.0, 0.35),
        "top_ylim": (0.40, 0.70),
        "ratio": (1, 5),
    },
    "theta_watterson": {
        "broken": True,
        "bot_ylim": (0.18, 0.35),
        "top_ylim": (0.40, 0.96),
        "ratio": (1, 4),
    },
    "tajimas_d": {
        "broken": False,
        "ylim": (-5.0, 4.0),
    },
}


def median_iqr(values: np.ndarray) -> tuple[float, float, float]:
    if values.size == 0:
        return (np.nan, np.nan, np.nan)
    med = float(np.median(values))
    q25, q75 = np.percentile(values, [25, 75])
    return med, float(q25), float(q75)


def wild_baseline(df: pd.DataFrame, metric: str) -> float:
    wild_pools = [p for p in POOL_ORDER if GROUP_OF[p] == "wild"]
    medians = []
    for p in wild_pools:
        v = df.loc[(df["pool"] == p) & (df["metric"] == metric), "value"].values
        if v.size:
            medians.append(np.median(v))
    return float(np.median(medians)) if medians else float("nan")


def _draw_separators(ax: plt.Axes, adj_positions: np.ndarray) -> None:
    super_of = lambda q: ("wild" if GROUP_OF[q] == "wild"
                          else "founder" if GROUP_OF[q] == "founder"
                          else GROUP_OF[q])
    for i, p in enumerate(POOL_ORDER[:-1]):
        nxt = POOL_ORDER[i + 1]
        if (GROUP_OF[p] == "wild" and GROUP_OF[nxt] == "wild"
                and SITE_OF[p] == SITE_OF[nxt]):
            continue
        xmid = (adj_positions[i] + adj_positions[i + 1]) / 2
        if super_of(p) == super_of(nxt):
            ax.axvline(xmid, color=C_GREY_LIGHT, linestyle=":", linewidth=0.6, zorder=0)
        else:
            ax.axvline(xmid, color=C_GREY_MID, linestyle="-", linewidth=0.5,
                       zorder=0, alpha=0.45)


def _scatter_cloud_and_medians(ax: plt.Axes, df: pd.DataFrame, metric: str,
                                adj_positions: np.ndarray) -> None:
    for x, p in zip(adj_positions, POOL_ORDER):
        vals = df.loc[(df["pool"] == p) & (df["metric"] == metric), "value"].values
        if vals.size == 0:
            continue
        jitter = RNG.uniform(-JITTER_HALF_WIDTH, JITTER_HALF_WIDTH, size=vals.size)
        ax.scatter(x + jitter, vals, s=CLOUD_SIZE, c=CLOUD_COLOR,
                   alpha=CLOUD_ALPHA, edgecolors="none", linewidths=0,
                   marker=".", rasterized=True, zorder=1)
    for x, p in zip(adj_positions, POOL_ORDER):
        vals = df.loc[(df["pool"] == p) & (df["metric"] == metric), "value"].values
        med, lo, hi = median_iqr(vals)
        if np.isnan(med):
            continue
        c = _color_of(p)
        ax.errorbar(x, med, yerr=[[med - lo], [hi - med]],
                    fmt="none", ecolor=c, elinewidth=ERR_LW,
                    capsize=1.6, capthick=ERR_LW, zorder=3)
        ax.scatter(x, med, s=POINT_SIZE, c=c, edgecolors="none",
                   linewidths=0, zorder=4)


def _format_x(ax: plt.Axes, adj_positions: np.ndarray) -> None:
    ax.set_xticks(adj_positions)
    ax.set_xticklabels([POOL_LABELS[p] for p in POOL_ORDER], fontsize=6)
    for tick, p in zip(ax.get_xticklabels(), POOL_ORDER):
        tick.set_color(_color_of(p))
        tick.set_fontweight("bold")
    ax.tick_params(axis="x", length=0, pad=2)


def _add_super_labels(ax: plt.Axes, adj_positions: np.ndarray) -> None:
    g10_label = {"g10_B": "B", "g10_T": "T", "g10_M": "B+T"}
    group_centers: dict[str, list[float]] = {}
    super_centers: dict[str, list[float]] = {"Wild": [], "Founder": [], "G10": []}
    for i, p in enumerate(POOL_ORDER):
        if GROUP_OF[p] == "wild":
            key = SITE_OF[p]
            super_centers["Wild"].append(adj_positions[i])
        elif GROUP_OF[p] == "founder":
            key = "Founder"
            super_centers["Founder"].append(adj_positions[i])
        else:
            key = g10_label[GROUP_OF[p]]
            super_centers["G10"].append(adj_positions[i])
        group_centers.setdefault(key, []).append(adj_positions[i])

    y_max = ax.get_ylim()[1]
    y_min = ax.get_ylim()[0]
    rng = y_max - y_min
    y_grp = y_max + 0.04 * rng
    y_super = y_max + 0.45 * rng
    for grp, positions in group_centers.items():
        x_center = float(np.mean(positions))
        ax.text(x_center, y_grp, grp, ha="center", va="bottom",
                fontsize=7, color=C_GREY_MID, fontweight="bold", clip_on=False)
    for sgrp, positions in super_centers.items():
        x_lo, x_hi = min(positions), max(positions)
        x_c = (x_lo + x_hi) / 2
        ax.plot([x_lo - 0.25, x_hi + 0.25], [y_super, y_super],
                color=C_GREY_MID, linewidth=0.55, clip_on=False)
        ax.text(x_c, y_super + 0.05 * rng, sgrp, ha="center", va="bottom",
                fontsize=8.5, color="black", fontweight="bold", clip_on=False)


def _draw_break_marks(ax_top: plt.Axes, ax_bot: plt.Axes) -> None:
    d = 0.012  # length of the diagonal lines (in axes coords)
    kw = dict(color="black", clip_on=False, linewidth=0.7)
    # Bottom of top axis: small diagonals
    ax_top.plot([-d, +d], [-d * 5, +d * 5], transform=ax_top.transAxes, **kw)
    ax_top.plot([1 - d, 1 + d], [-d * 5, +d * 5], transform=ax_top.transAxes, **kw)
    # Top of bottom axis: small diagonals
    ax_bot.plot([-d, +d], [1 - d, 1 + d], transform=ax_bot.transAxes, **kw)
    ax_bot.plot([1 - d, 1 + d], [1 - d, 1 + d], transform=ax_bot.transAxes, **kw)


def _plot_g_strip(ax: plt.Axes, adj_positions: np.ndarray,
                   g_values: dict[str, float]) -> None:
    bar_w = 0.55
    for x, p in zip(adj_positions, POOL_ORDER):
        if p not in g_values:
            continue
        g = g_values[p]
        c = _color_of(p)
        ax.bar(x, g, width=bar_w, color=c, alpha=0.85, edgecolor="none",
               linewidth=0, zorder=2)
    ax.axhline(0, color="black", linewidth=0.6, zorder=1)
    ax.set_ylabel("G$_{k\\geq2}$", fontsize=8, labelpad=4)
    ax.tick_params(axis="y", labelsize=6.5, length=2)
    ax.locator_params(axis="y", nbins=3)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    # match x-axis range to data axis
    ax.set_xlim(adj_positions.min() - 0.6, adj_positions.max() + 0.6)


def plot_metric_figure(df: pd.DataFrame, metric: str, label: str,
                       show_zero: bool, out_stem: Path,
                       g_values: dict[str, float]) -> None:
    cfg = Y_CONFIG[metric]
    adj_positions = _adj_positions()
    wild_y = wild_baseline(df, metric)
    fig = plt.figure(figsize=(6.4, 3.1))

    if cfg["broken"]:
        rt, rb = cfg["ratio"]
        gs = gridspec.GridSpec(3, 1, figure=fig,
                                height_ratios=[rt, rb, max(2, (rt + rb) / 4)],
                                hspace=0.08,
                                left=0.10, right=0.99, top=0.80, bottom=0.18)
        ax_top = fig.add_subplot(gs[0])
        ax_bot = fig.add_subplot(gs[1], sharex=ax_top)
        ax_g = fig.add_subplot(gs[2], sharex=ax_top)

        for ax in (ax_top, ax_bot):
            _scatter_cloud_and_medians(ax, df, metric, adj_positions)
            _draw_separators(ax, adj_positions)
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)

        ax_top.set_ylim(*cfg["top_ylim"])
        ax_bot.set_ylim(*cfg["bot_ylim"])
        ax_top.spines["bottom"].set_visible(False)
        ax_bot.spines["top"].set_visible(False)
        ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False, top=False)
        ax_bot.tick_params(axis="x", which="both", bottom=False, labelbottom=False, top=False)
        ax_top.tick_params(axis="y", labelsize=7, length=2.5)
        ax_bot.tick_params(axis="y", labelsize=7, length=2.5)
        ax_top.locator_params(axis="y", nbins=3)
        ax_bot.locator_params(axis="y", nbins=5)

        if show_zero:
            ax_bot.axhline(0, color=C_GREY_MID, linestyle="--", linewidth=0.7, zorder=0)
        bot_lo, bot_hi = cfg["bot_ylim"]
        if bot_lo <= wild_y <= bot_hi:
            ax_bot.axhline(wild_y, color="black", linestyle="--", linewidth=0.7,
                            alpha=0.55, zorder=0)

        _draw_separators(ax_g, adj_positions)
        _plot_g_strip(ax_g, adj_positions, g_values)
        _format_x(ax_g, adj_positions)

        fig.text(0.022, 0.55, label, ha="center", va="center",
                 rotation="vertical", fontsize=10)

        _add_super_labels(ax_top, adj_positions)
        _draw_break_marks(ax_top, ax_bot)
    else:
        gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[4, 1.3], hspace=0.08,
                                left=0.10, right=0.99, top=0.80, bottom=0.18)
        ax = fig.add_subplot(gs[0])
        ax_g = fig.add_subplot(gs[1], sharex=ax)
        _scatter_cloud_and_medians(ax, df, metric, adj_positions)
        _draw_separators(ax, adj_positions)
        ax.set_ylim(*cfg["ylim"])
        if show_zero:
            ax.axhline(0, color=C_GREY_MID, linestyle="--", linewidth=0.7, zorder=0)
        ax.axhline(wild_y, color="black", linestyle="--", linewidth=0.7,
                   alpha=0.55, zorder=0)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        ax.tick_params(axis="y", labelsize=7, length=2.5)
        ax.locator_params(axis="y", nbins=5)
        ax.set_ylabel(label, fontsize=10)
        _add_super_labels(ax, adj_positions)

        _draw_separators(ax_g, adj_positions)
        _plot_g_strip(ax_g, adj_positions, g_values)
        _format_x(ax_g, adj_positions)

    out_png = out_stem.with_suffix(".png")
    out_svg = out_stem.with_suffix(".svg")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=SAVE_DPI, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    plt.rcParams.update({
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.7,
    })
    df = load_wild_long()
    g_values = load_g_k2()
    for metric, label in METRICS:
        out_stem = OUT_DIR / f"s3_main_diversity_{metric}"
        show_zero = (metric == "tajimas_d")
        plot_metric_figure(df, metric, label, show_zero=show_zero,
                            out_stem=out_stem, g_values=g_values)


if __name__ == "__main__":
    main()
