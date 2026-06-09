#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(".")
DIV_FILE = ROOT / "grenfst/diversity_attrition/attrition_pi_390000diversity.csv"
OUT_PNG  = ROOT / "final_plots/wild/s3_fig1_wild_diversity.png"
OUT_SVG  = ROOT / "final_plots/wild/s3_fig1_wild_diversity.svg"

WILD_POOLS = ["AVB", "AVT", "PSB", "PST", "RMB", "RMT"]
FOUNDER_POOLS = ["F1G00", "F2G00", "F3G00", "F4G00"]
G10_B_POOLS = ["B1G10", "B2G10", "B3G10", "B4G10"]
G10_T_POOLS = ["T1G10", "T2G10", "T3G10", "T4G10"]
G10_M_POOLS = ["M1G10", "M2G10", "M3G10", "M4G10"]
G10_POOLS = G10_B_POOLS + G10_T_POOLS + G10_M_POOLS
POOL_ORDER = WILD_POOLS + FOUNDER_POOLS + G10_POOLS

POOL_LABELS = {p: p for p in WILD_POOLS}
POOL_LABELS.update({"F1G00": "F1", "F2G00": "F2", "F3G00": "F3", "F4G00": "F4"})
POOL_LABELS.update({p: p[:2] for p in G10_POOLS})  # B1, B2, T1, M1, ...

HOST_OF = {"AVB": "B", "AVT": "T", "PSB": "B", "PST": "T", "RMB": "B", "RMT": "T"}
SITE_OF = {"AVB": "AV", "AVT": "AV", "PSB": "PS", "PST": "PS", "RMB": "RM", "RMT": "RM"}
GROUP_OF = {**{p: "wild" for p in WILD_POOLS},
            **{p: "founder" for p in FOUNDER_POOLS},
            **{p: "g10_B" for p in G10_B_POOLS},
            **{p: "g10_T" for p in G10_T_POOLS},
            **{p: "g10_M" for p in G10_M_POOLS}}

C_B = "#499FFF"        # blue  — B host
C_T = "#EDB72D"        # gold  — T host
C_FOUNDER = "#444444"  # dark grey — founders
C_BT = "#9BAB96"       # sage green — B+T (Levene multi-host)
C_GREY_MID = "#555555"
C_GREY_LIGHT = "#BBBBBB"

METRICS = [
    ("theta_pi",        r"$\theta_\pi$"),
    ("theta_watterson", r"$\theta_W$"),
    ("tajimas_d",       r"Tajima's $D$"),
]


def load_wild_long() -> pd.DataFrame:
    raw = pd.read_csv(DIV_FILE)
    rows = []
    for pool in POOL_ORDER:
        for key, _label in METRICS:
            col = f"{pool}.1.{key}"
            if col not in raw.columns:
                continue
            sub = pd.DataFrame({
                "pool": pool,
                "group": GROUP_OF[pool],
                "site": SITE_OF.get(pool, "Founder"),
                "host": HOST_OF.get(pool, "F"),
                "metric": key,
                "value": pd.to_numeric(raw[col], errors="coerce"),
            })
            rows.append(sub)
    df = pd.concat(rows, ignore_index=True)
    df = df.dropna(subset=["value"])
    return df


def _adj_positions() -> np.ndarray:
    site_groups = {"AV": 0, "PS": 1, "RM": 2}
    block_gap = 1.4   # gap between super-groups (wild | founder | G10-B | G10-T | G10-B+T)
    within = 0.85     # spacing between adjacent pools within a super-group
    pos = []
    cursor = 0.0
    last_super = None
    for p in POOL_ORDER:
        grp = GROUP_OF[p]
        if grp == "wild":
            super_grp = "wild"
        elif grp == "founder":
            super_grp = "founder"
        else:
            super_grp = grp  # g10_B / g10_T / g10_M each their own super-group
        if last_super is None:
            cursor = 0.0
        elif super_grp != last_super:
            cursor += block_gap
        else:
            cursor += within
            if super_grp == "wild":
                prev_p = POOL_ORDER[POOL_ORDER.index(p) - 1]
                if SITE_OF[prev_p] != SITE_OF[p]:
                    cursor += 0.4
        pos.append(cursor)
        last_super = super_grp
    return np.array(pos)


def _color_of(pool: str) -> str:
    grp = GROUP_OF[pool]
    if grp == "founder":
        return C_FOUNDER
    if grp == "g10_M":
        return C_BT
    if grp in ("g10_B",):
        return C_B
    if grp in ("g10_T",):
        return C_T
    return C_B if HOST_OF[pool] == "B" else C_T


def plot_panel(ax: plt.Axes, df: pd.DataFrame, metric: str, label: str,
               show_zero: bool = False) -> None:
    adj_positions = _adj_positions()
    data = [df.loc[(df["pool"] == p) & (df["metric"] == metric), "value"].values
            for p in POOL_ORDER]
    colors = [_color_of(p) for p in POOL_ORDER]

    bp = ax.boxplot(
        data,
        positions=adj_positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.4),
        whiskerprops=dict(color="black", linewidth=0.7),
        capprops=dict(color="black", linewidth=0.7),
        boxprops=dict(edgecolor="black", linewidth=0.6),
    )
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(mcolors.to_rgba(c, alpha=0.85))

    if show_zero:
        ax.axhline(0, color=C_GREY_MID, linestyle="--", linewidth=0.8, zorder=1)

    super_of = lambda q: ("wild" if GROUP_OF[q] == "wild"
                          else "founder" if GROUP_OF[q] == "founder"
                          else GROUP_OF[q])
    for i, p in enumerate(POOL_ORDER[:-1]):
        nxt = POOL_ORDER[i + 1]
        if (GROUP_OF[p] == "wild" and GROUP_OF[nxt] == "wild"
                and SITE_OF[p] == SITE_OF[nxt]):
            continue
        if super_of(p) == super_of(nxt):
            xmid = (adj_positions[i] + adj_positions[i + 1]) / 2
            ax.axvline(xmid, color=C_GREY_LIGHT, linestyle=":", linewidth=0.8, zorder=0)
        else:
            xmid = (adj_positions[i] + adj_positions[i + 1]) / 2
            ax.axvline(xmid, color=mcolors.to_rgba(C_GREY_MID, alpha=0.45),
                       linestyle="-", linewidth=0.7, zorder=0)

    ax.set_xticks(adj_positions)
    ax.set_xticklabels([POOL_LABELS[p] for p in POOL_ORDER], fontsize=9)
    for tick, p in zip(ax.get_xticklabels(), POOL_ORDER):
        tick.set_color(_color_of(p))
        tick.set_fontweight("bold")
    ax.set_ylabel(label, fontsize=11)
    ax.tick_params(axis="y", labelsize=9)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def main() -> None:
    plt.rcParams.update({
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.8,
    })

    df = load_wild_long()

    fig, axes = plt.subplots(1, 3, figsize=(20, 4.3), gridspec_kw={"wspace": 0.22})
    for ax, (metric, label) in zip(axes, METRICS):
        show_zero = (metric == "tajimas_d")
        plot_panel(ax, df, metric, label, show_zero=show_zero)

    axes[0].set_title("A", loc="left", fontsize=16, fontweight="bold", pad=8)
    axes[1].set_title("B", loc="left", fontsize=16, fontweight="bold", pad=8)
    axes[2].set_title("C", loc="left", fontsize=16, fontweight="bold", pad=8)

    adj_positions = _adj_positions()
    group_centers: dict[str, list[float]] = {}
    g10_label = {"g10_B": "B", "g10_T": "T", "g10_M": "B+T"}
    for i, p in enumerate(POOL_ORDER):
        if GROUP_OF[p] == "wild":
            key = SITE_OF[p]
        elif GROUP_OF[p] == "founder":
            key = "Founder"
        else:
            key = g10_label[GROUP_OF[p]]
        group_centers.setdefault(key, []).append(adj_positions[i])
    super_centers = {"Wild": [], "Founder": [], "G10": []}
    for i, p in enumerate(POOL_ORDER):
        if GROUP_OF[p] == "wild":
            super_centers["Wild"].append(adj_positions[i])
        elif GROUP_OF[p] == "founder":
            super_centers["Founder"].append(adj_positions[i])
        else:
            super_centers["G10"].append(adj_positions[i])

    for ax in axes:
        y_max = ax.get_ylim()[1]
        y_min = ax.get_ylim()[0]
        rng = y_max - y_min
        y_grp = y_max + 0.02 * rng
        y_super = y_max + 0.13 * rng
        for grp, positions in group_centers.items():
            x_center = float(np.mean(positions))
            ax.text(x_center, y_grp, grp, ha="center", va="bottom",
                    fontsize=9, color=C_GREY_MID, fontweight="bold")
        for sgrp, positions in super_centers.items():
            x_lo, x_hi = min(positions), max(positions)
            x_c = (x_lo + x_hi) / 2
            ax.plot([x_lo - 0.25, x_hi + 0.25], [y_super, y_super],
                    color=C_GREY_MID, linewidth=0.7, clip_on=False)
            ax.text(x_c, y_super + 0.012 * rng, sgrp, ha="center", va="bottom",
                    fontsize=10.5, color="black", fontweight="bold")
        ax.set_ylim(y_min, y_super + 0.09 * rng)

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")
    fig.savefig(OUT_SVG.with_suffix(".pdf"), bbox_inches="tight")
    _strip_svg_clips(OUT_SVG)


def _strip_svg_clips(svg_path: Path) -> None:
    import re
    txt = svg_path.read_text()
    txt = re.sub(r"<clipPath[^>]*>.*?</clipPath>", "", txt, flags=re.DOTALL)
    txt = re.sub(r'\s*clip-path="url\([^)]+\)"', "", txt)
    svg_path.write_text(txt)

if __name__ == "__main__":
    main()
