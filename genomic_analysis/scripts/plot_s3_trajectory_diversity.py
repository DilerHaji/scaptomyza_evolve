#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(".")
DEFAULT_INPUT = ROOT / "grenfst/diversity/trajectory_pi_200000_nonoverlap.csv"
OUT_INTERMEDIATE_PNG = ROOT / "final_plots/wild/s3_trajectory_intermediate.png"
OUT_INTERMEDIATE_SVG = ROOT / "final_plots/wild/s3_trajectory_intermediate.svg"
OUT_ENDPOINTS_PNG    = ROOT / "final_plots/wild/s3_trajectory_endpoints.png"
OUT_ENDPOINTS_SVG    = ROOT / "final_plots/wild/s3_trajectory_endpoints.svg"

WILD_POOLS    = ["AVB", "AVT", "PSB", "PST", "RMB", "RMT"]
FOUNDER_POOLS = [f"F{i}G00" for i in range(1, 5)]

TREATMENTS = ["B", "T", "M"]
TREAT_LABEL = {"B": "B", "T": "T", "M": "B+T"}
TREAT_COLOR = {"B": "#499FFF", "T": "#EDB72D", "M": "#9BAB96"}
REP_MARKER  = {1: "o", 2: "s", 3: "^", 4: "D"}
REP_MFC_FRAC = {1: 1.0, 2: 0.75, 3: 0.55, 4: 0.35}   # dim within-treatment shade per rep

METRICS = [
    ("theta_pi",        r"$\theta_\pi$"),
    ("theta_watterson", r"$\theta_W$"),
    ("tajimas_d",       r"Tajima's $D$"),
]

C_WILD = "#555555"


def load_trajectory(csv_path: Path) -> tuple[pd.DataFrame, dict]:
    raw = pd.read_csv(csv_path)
    metric_suffixes = ("theta_pi", "theta_watterson", "tajimas_d")
    rows = []
    per_window: dict[tuple[str, str], np.ndarray] = {}
    for col in raw.columns:
        metric = next((m for m in metric_suffixes if col.endswith("." + m)), None)
        if metric is None:
            continue
        prefix = col[: -(len(metric) + 1)]
        pool = prefix.rstrip(".1") if prefix.endswith(".1") else prefix
        vals = pd.to_numeric(raw[col], errors="coerce").dropna().values
        per_window[(pool, metric)] = vals
        rows.append({"pool": pool, "metric": metric, "median": float(np.median(vals))})
    df = pd.DataFrame(rows)

    def decode(p):
        if p in WILD_POOLS:
            return pd.Series({"treatment": "wild", "rep": np.nan, "generation": np.nan})
        if p in FOUNDER_POOLS:
            return pd.Series({"treatment": "founder", "rep": np.nan, "generation": 0})
        treat = p[0]
        rep = int(p[1])
        gen = int(p[3:])
        return pd.Series({"treatment": treat, "rep": rep, "generation": gen})
    decoded = df["pool"].apply(decode)
    df = pd.concat([df, decoded], axis=1)
    return df, per_window


def plot_intermediate(df: pd.DataFrame, out_png: Path, out_svg: Path) -> None:
    import matplotlib.colors as mcolors
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    fig, axes = plt.subplots(3, 3, figsize=(12.5, 9.5), sharex="col", sharey="row",
                              gridspec_kw={"wspace": 0.16, "hspace": 0.18})

    for ci, treat in enumerate(TREATMENTS):
        for ri, (metric, ylabel) in enumerate(METRICS):
            ax = axes[ri, ci]

            color_line   = mcolors.to_rgba(TREAT_COLOR[treat], alpha=0.55)
            color_marker = mcolors.to_rgba(TREAT_COLOR[treat], alpha=0.85)
            color_reg    = mcolors.to_rgba(TREAT_COLOR[treat], alpha=0.95)

            pooled_mid_x: list[float] = []
            pooled_mid_y: list[float] = []
            per_rep_slopes: list[tuple[int, float]] = []
            for rep in sorted(df.loc[df["treatment"] == treat, "rep"].dropna().unique()):
                sub = df[(df["treatment"] == treat) & (df["rep"] == rep)
                         & (df["metric"] == metric)
                         & (df["generation"] >= 1) & (df["generation"] <= 9)
                         ].sort_values("generation")
                xs = sub["generation"].astype(int).tolist()
                ys = sub["median"].astype(float).tolist()
                if not xs:
                    continue

                ax.plot(xs, ys, "-", color=color_line, linewidth=0.7, zorder=2)
                ax.scatter(xs, ys, s=11, c=[color_marker], marker="o",
                            edgecolors="none", linewidths=0, zorder=3)

                for g, v in zip(xs, ys):
                    pooled_mid_x.append(float(g))
                    pooled_mid_y.append(float(v))

                if len(xs) >= 3:
                    s, _ = np.polyfit(xs, ys, 1)
                    per_rep_slopes.append((int(rep), float(s)))

                label_x = xs[-1] + 0.15
                ax.text(label_x, ys[-1], str(int(rep)),
                        color=color_reg, fontsize=9, fontweight="bold",
                        va="center", ha="left", zorder=5)

            if len(pooled_mid_x) >= 3:
                slope, intercept = np.polyfit(pooled_mid_x, pooled_mid_y, 1)
                x_fit = np.array([1.0, 9.0])
                y_fit = slope * x_fit + intercept
                ax.plot(x_fit, y_fit, "--", color=color_reg, linewidth=1.8, zorder=4)

            if metric == "tajimas_d":
                ax.axhline(0, color=C_WILD, linestyle=":", linewidth=0.6,
                           alpha=0.5, zorder=0)

            ax.set_xlim(0.5, 9.8)
            ax.set_xticks(range(1, 10))
            if ri == 2:
                ax.set_xlabel("Generation", fontsize=10)
            if ci == 0:
                ax.set_ylabel(ylabel, fontsize=11)
            if ri == 0:
                ax.set_title(TREAT_LABEL[treat], fontsize=12, fontweight="bold")
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
            ax.tick_params(axis="both", labelsize=8, length=2.5)

            if per_rep_slopes:
                iax = inset_axes(ax, width="22%", height="35%",
                                  loc="upper right", borderpad=0.7)
                reps = [r for r, _ in per_rep_slopes]
                slopes = [s for _, s in per_rep_slopes]
                bar_colors = [color_reg] * len(reps)
                iax.bar(range(len(reps)), slopes, color=bar_colors,
                         edgecolor="none", width=0.75)
                iax.axhline(0, color="black", linewidth=0.5, zorder=0)
                iax.set_xticks(range(len(reps)))
                iax.set_xticklabels([str(r) for r in reps], fontsize=7)
                iax.tick_params(axis="x", length=0, pad=1)
                iax.tick_params(axis="y", labelsize=6, length=1.5, pad=1)
                iax.locator_params(axis="y", nbins=3)
                for spine in ("top", "right"):
                    iax.spines[spine].set_visible(False)
                iax.spines["left"].set_linewidth(0.5)
                iax.spines["bottom"].set_linewidth(0.5)
                iax.set_title("slope", fontsize=7, pad=1)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    fig.savefig(out_svg.with_suffix(".pdf"), bbox_inches="tight")
    _strip_svg_clips(out_svg)
    plt.close(fig)


def plot_endpoints(df: pd.DataFrame, per_window: dict, out_png: Path,
                    out_svg: Path) -> None:
    import matplotlib.colors as mcolors
    fig, axes = plt.subplots(3, 3, figsize=(9.0, 9.5), sharey="row",
                              gridspec_kw={"wspace": 0.06, "hspace": 0.18})

    rng = np.random.default_rng(42)
    CLOUD_COLOR = "#CFCFCF"
    CLOUD_ALPHA = 0.25
    CLOUD_SIZE  = 1.5
    CLOUD_JITTER = 0.14
    MEDIAN_JITTER = 0.10

    founder_pools = [p for p in FOUNDER_POOLS]

    for ci, treat in enumerate(TREATMENTS):
        for ri, (metric, ylabel) in enumerate(METRICS):
            ax = axes[ri, ci]

            color_founder = "black"
            color_g10     = mcolors.to_rgba(TREAT_COLOR[treat], alpha=0.90)

            for p in founder_pools:
                vals = per_window.get((p, metric))
                if vals is None or len(vals) == 0:
                    continue
                jitter = rng.uniform(-CLOUD_JITTER, CLOUD_JITTER, size=vals.size)
                ax.scatter(0 + jitter, vals, s=CLOUD_SIZE, c=CLOUD_COLOR,
                           alpha=CLOUD_ALPHA, edgecolors="none", linewidths=0,
                           marker=".", rasterized=True, zorder=1)
            f_vals = df[(df["treatment"] == "founder") & (df["metric"] == metric)
                         ]["median"].astype(float).tolist()
            f_x = rng.uniform(-MEDIAN_JITTER, MEDIAN_JITTER, size=len(f_vals))
            ax.scatter(f_x, f_vals, s=26, c=color_founder, marker="o",
                       edgecolors="none", linewidths=0, zorder=3)

            g10 = df[(df["treatment"] == treat) & (df["metric"] == metric)
                     & (df["generation"] == 10)
                     ].sort_values("rep")
            reps = g10["rep"].astype(int).tolist()
            ys = g10["median"].astype(float).tolist()
            for r in reps:
                pool = f"{treat}{r}G10"
                vals = per_window.get((pool, metric))
                if vals is None or len(vals) == 0:
                    continue
                jitter = rng.uniform(-CLOUD_JITTER, CLOUD_JITTER, size=vals.size)
                ax.scatter(1 + jitter, vals, s=CLOUD_SIZE, c=CLOUD_COLOR,
                           alpha=CLOUD_ALPHA, edgecolors="none", linewidths=0,
                           marker=".", rasterized=True, zorder=1)
            g_x = 1 + rng.uniform(-MEDIAN_JITTER, MEDIAN_JITTER, size=len(reps))
            ax.scatter(g_x, ys, s=26, c=[color_g10], marker="o",
                       edgecolors="none", linewidths=0, zorder=3)
            for x, y, r in zip(g_x, ys, reps):
                ax.text(x + 0.03, y, str(r), color=color_g10, fontsize=9,
                        fontweight="bold", va="center", ha="left", zorder=5)

            if metric == "tajimas_d":
                ax.axhline(0, color=C_WILD, linestyle=":", linewidth=0.6,
                           alpha=0.5, zorder=0)

            ax.set_xlim(-0.4, 1.5)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["F", "G10"] if ri == 2 else ["", ""], fontsize=9)
            if ci == 0:
                ax.set_ylabel(ylabel, fontsize=11)
            if ri == 0:
                ax.set_title(TREAT_LABEL[treat], fontsize=12, fontweight="bold")
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
            ax.tick_params(axis="y", labelsize=8, length=2.5)
            ax.tick_params(axis="x", labelsize=9, length=0, pad=2)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    fig.savefig(out_svg.with_suffix(".pdf"), bbox_inches="tight")
    _strip_svg_clips(out_svg)
    plt.close(fig)
    

def _strip_svg_clips(svg_path: Path) -> None:
    import re
    txt = svg_path.read_text()
    txt = re.sub(r"<clipPath[^>]*>.*?</clipPath>", "", txt, flags=re.DOTALL)
    txt = re.sub(r'\s*clip-path="url\([^)]+\)"', "", txt)
    svg_path.write_text(txt)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(DEFAULT_INPUT),
                    help="path to grenedalf trajectory diversity CSV")
    args = ap.parse_args()

    plt.rcParams.update({
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.7,
    })
    df, per_window = load_trajectory(Path(args.input))
    plot_intermediate(df, OUT_INTERMEDIATE_PNG, OUT_INTERMEDIATE_SVG)
    plot_endpoints(df, per_window, OUT_ENDPOINTS_PNG, OUT_ENDPOINTS_SVG)


if __name__ == "__main__":
    main()
