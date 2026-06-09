#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(".")
PER_WIN_CSV = ROOT / "grenfst/fst_wild/wild_fst_390000fst.csv"
SUMMARY_TSV = ROOT / "final_plots/wild/wild_fst_summary.tsv"
OUT_PNG = ROOT / "final_plots/wild/wild_fst_twopanel.png"
OUT_SVG = ROOT / "final_plots/wild/wild_fst_twopanel.svg"

POOL_ORDER = ["AVB", "AVT", "PSB", "PST", "RMB", "RMT"]
HOST_OF = {"AVB": "B", "AVT": "T", "PSB": "B", "PST": "T", "RMB": "B", "RMT": "T"}

C_B = "#499FFF"   # blue   — B-host
C_T = "#EDB72D"   # gold   — T-host
C_WITHIN = "#CC3333"     # red dashed outline for within-site B-vs-T cells
C_GREY = "#7F7F7F"
C_LIGHT = "#CCCCCC"
CATEGORY_COLORS = {
    "within_site_BvsT":         "#CC3333",   # red — host axis, within-site
    "between_site_same_host":   "#7F7F7F",   # grey — geography, same host
    "between_site_diff_host":   "#BBBBBB",   # light grey — geography, different host
}
CATEGORY_LABEL = {
    "within_site_BvsT":       "Within-site\nB vs T",
    "between_site_same_host": "Between-site\nsame host",
    "between_site_diff_host": "Between-site\ndifferent host",
}
CATEGORY_ORDER = [
    "within_site_BvsT",
    "between_site_same_host",
    "between_site_diff_host",
]

def load_summary(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["poolA"] = df["pair"].str.split(" vs ").str[0]
    df["poolB"] = df["pair"].str.split(" vs ").str[1]
    return df


def load_per_window(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path)
    fst_cols = [c for c in raw.columns if c.endswith(".fst")]
    long_rows = []
    for col in fst_cols:
        stem = col[:-4]                    # drop ".fst"
        a, b = stem.split(":")             # "AVB.1", "AVT.1"
        pool_a = a.split(".")[0]
        pool_b = b.split(".")[0]
        if pool_a not in POOL_ORDER or pool_b not in POOL_ORDER:
            continue
        sub = pd.DataFrame({
            "pair": f"{pool_a} vs {pool_b}",
            "fst": pd.to_numeric(raw[col], errors="coerce"),
        })
        long_rows.append(sub)
    long = pd.concat(long_rows, ignore_index=True)
    long = long[long["fst"].notna()].copy()
    return long


def categorise_pair(pool_a: str, pool_b: str) -> str:
    site_a, site_b = pool_a[:2], pool_b[:2]
    host_a, host_b = HOST_OF[pool_a], HOST_OF[pool_b]
    if site_a == site_b:
        return "within_site_BvsT"
    if host_a == host_b:
        return "between_site_same_host"
    return "between_site_diff_host"


def normalise_pair_key(pair: str) -> tuple[str, str] | None:
    parts = [p.strip() for p in pair.split("vs")]
    if len(parts) != 2:
        return None
    a, b = parts
    if a not in POOL_ORDER or b not in POOL_ORDER:
        return None
    return (a, b) if POOL_ORDER.index(a) <= POOL_ORDER.index(b) else (b, a)

def build_fst_matrix(summary: pd.DataFrame) -> np.ndarray:
    n = len(POOL_ORDER)
    mat = np.full((n, n), np.nan)
    idx = {p: i for i, p in enumerate(POOL_ORDER)}
    for _, row in summary.iterrows():
        i, j = idx[row["poolA"]], idx[row["poolB"]]
        # use median_fst
        mat[i, j] = row["median_fst"]
        mat[j, i] = row["median_fst"]
    np.fill_diagonal(mat, 0.0)
    return mat


def plot_heatmap(ax: plt.Axes, mat: np.ndarray) -> None:
    show = np.tril(np.ones_like(mat, dtype=bool))
    display = np.where(show, mat, np.nan)

    vmax = float(np.nanmax(np.abs(display)))
    vmax = max(vmax, 0.011)
    cmap = plt.cm.Greys
    im = ax.imshow(display, cmap=cmap, vmin=0.0, vmax=vmax, aspect="equal")

    for i in range(len(POOL_ORDER)):
        for j in range(len(POOL_ORDER)):
            if not show[i, j]:
                continue
            if i == j:
                ax.text(j, i, POOL_ORDER[i], ha="center", va="center",
                        fontsize=8, color=C_GREY)
                continue
            val = mat[i, j]
            disp_val = 0.0 if abs(val) < 5e-4 else val
            txt_color = "white" if val > 0.55 * vmax else "black"
            ax.text(j, i, f"{disp_val:.3f}", ha="center", va="center",
                    fontsize=7.5, color=txt_color)

    ax.set_xticks(range(len(POOL_ORDER)))
    ax.set_yticks(range(len(POOL_ORDER)))
    ax.set_xticklabels(POOL_ORDER, rotation=0)
    ax.set_yticklabels(POOL_ORDER)
    for tick, pool in zip(ax.get_xticklabels(), POOL_ORDER):
        tick.set_color(C_B if HOST_OF[pool] == "B" else C_T)
        tick.set_fontweight("bold")
    for tick, pool in zip(ax.get_yticklabels(), POOL_ORDER):
        tick.set_color(C_B if HOST_OF[pool] == "B" else C_T)
        tick.set_fontweight("bold")

    within_site_pairs = [("AVB", "AVT"), ("PSB", "PST"), ("RMB", "RMT")]
    idx = {p: i for i, p in enumerate(POOL_ORDER)}
    for a, b in within_site_pairs:
        i, j = idx[a], idx[b]
        i, j = max(i, j), min(i, j)
        rect = mpatches.Rectangle(
            (j - 0.5, i - 0.5), 1, 1,
            fill=False, linestyle="--", linewidth=1.8,
            edgecolor=C_WITHIN, zorder=5, clip_on=False,
        )
        ax.add_patch(rect)

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlim(-0.5, len(POOL_ORDER) - 0.5)
    ax.set_ylim(len(POOL_ORDER) - 0.5, -0.5)

    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.68, pad=0.02)
    cbar.set_label("Median pairwise $F_{ST}$", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax.set_title("A   Median pairwise $F_{ST}$ among wild pools",
                 loc="left", fontsize=11, fontweight="bold")


def plot_violins(ax: plt.Axes, per_window: pd.DataFrame, summary: pd.DataFrame) -> None:
    per_window = per_window.copy()
    pair_keys = per_window["pair"].map(normalise_pair_key)
    per_window["poolA"] = pair_keys.map(lambda t: t[0] if t else None)
    per_window["poolB"] = pair_keys.map(lambda t: t[1] if t else None)
    per_window = per_window.dropna(subset=["poolA", "poolB"])
    per_window["category"] = per_window.apply(
        lambda r: categorise_pair(r["poolA"], r["poolB"]), axis=1
    )

    lo, hi = -0.04, 0.10
    pw_clip = per_window.copy()
    pw_clip["fst"] = pw_clip["fst"].clip(lo, hi)

    data = [pw_clip.loc[pw_clip["category"] == c, "fst"].values
            for c in CATEGORY_ORDER]

    parts = ax.violinplot(data, positions=range(len(CATEGORY_ORDER)),
                          widths=0.78, showextrema=False, showmedians=False)
    for pc, cat in zip(parts["bodies"], CATEGORY_ORDER):
        pc.set_facecolor(CATEGORY_COLORS[cat])
        pc.set_alpha(0.35)
        pc.set_edgecolor("black")
        pc.set_linewidth(0.8)

    for i, c in enumerate(CATEGORY_ORDER):
        vals = per_window.loc[per_window["category"] == c, "fst"].values
        q25, q50, q75 = np.quantile(vals, [0.25, 0.5, 0.75])
        ax.plot([i - 0.18, i + 0.18], [q50, q50], color="black", lw=2.0, zorder=5)
        ax.plot([i, i], [q25, q75], color="black", lw=1.4, zorder=5)

    cat_x = {c: i for i, c in enumerate(CATEGORY_ORDER)}
    rng = np.random.default_rng(12345)
    for _, row in summary.iterrows():
        cat = row["type"]
        if cat not in cat_x:
            continue
        x = cat_x[cat] + rng.uniform(-0.22, 0.22)
        y = row["median_fst"]
        if cat == "within_site_BvsT":
            marker_color = C_WITHIN
            edge = "black"
        elif HOST_OF[row["poolA"]] == "B":
            marker_color = C_B if cat == "between_site_same_host" else "#AAAAAA"
            edge = "black"
        else:
            marker_color = C_T if cat == "between_site_same_host" else "#AAAAAA"
            edge = "black"
        ax.scatter(x, y, s=60, color=marker_color, edgecolor=edge,
                   linewidth=0.9, zorder=6)
        if cat == "within_site_BvsT":
            ax.annotate(row["poolA"][:2], (x, y),
                        textcoords="offset points", xytext=(8, -2),
                        fontsize=8, color="black")

    ax.axhline(0, color=C_GREY, linestyle=":", linewidth=1, zorder=1)
    ax.set_xticks(range(len(CATEGORY_ORDER)))
    ax.set_xticklabels([CATEGORY_LABEL[c] for c in CATEGORY_ORDER], fontsize=9)
    ax.set_ylabel("Per-window $F_{ST}$", fontsize=10)
    ax.set_ylim(lo, hi)
    ax.tick_params(axis="y", labelsize=8)

    legend_handles = [
        mpatches.Patch(facecolor=C_WITHIN, edgecolor="black", alpha=0.5,
                       label="Within-site B vs T"),
        mpatches.Patch(facecolor="#7F7F7F", edgecolor="black", alpha=0.5,
                       label="Between-site (geographic)"),
        plt.Line2D([0], [0], marker="o", linestyle="none",
                   markerfacecolor=C_B, markeredgecolor="black", markersize=8,
                   label="Median, B-host pair"),
        plt.Line2D([0], [0], marker="o", linestyle="none",
                   markerfacecolor=C_T, markeredgecolor="black", markersize=8,
                   label="Median, T-host pair"),
    ]
    ax.legend(handles=legend_handles, loc="upper right",
              fontsize=7.5, frameon=False, handletextpad=0.4)

    ax.set_title("B   Per-window $F_{ST}$ by pair category",
                 loc="left", fontsize=11, fontweight="bold")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

def main() -> None:
    plt.rcParams.update({
        "svg.fonttype": "none",           # editable text in Figma / Illustrator
        "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.8,
    })

    summary = load_summary(SUMMARY_TSV)
    per_window = load_per_window(PER_WIN_CSV)

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(11.5, 4.6),
        gridspec_kw={"width_ratios": [1.0, 1.15], "wspace": 0.28},
    )

    mat = build_fst_matrix(summary)
    plot_heatmap(ax1, mat)
    plot_violins(ax2, per_window, summary)

    fig.tight_layout()

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")

if __name__ == "__main__":
    main()
