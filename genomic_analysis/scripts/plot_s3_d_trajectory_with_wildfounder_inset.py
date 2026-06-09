#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(".")
DIV_FILE  = ROOT / "grenfst/diversity_attrition/attrition_pi_390000diversity.csv"
TRAJ_FILE = ROOT / "grenfst/diversity/trajectory_pi_200000_nonoverlap.csv"

WILD_POOLS    = ["AVB", "AVT", "PSB", "PST", "RMB", "RMT"]
FOUNDER_POOLS = ["F1G00", "F2G00", "F3G00", "F4G00"]
GENS = list(range(1, 10))
REPS = [1, 2, 3, 4]

METRIC = "tajimas_d"
METRIC_LABEL = r"Tajima's $D$"

C_B  = "#499FFF"
C_T  = "#EDB72D"
C_BT = "#9BAB96"
C_WILD    = "#6C7B8B"
C_FOUNDER = "#444444"
C_GREY    = "#555555"

TREATMENTS = [("B", C_B, "B"), ("T", C_T, "T"), ("M", C_BT, "B+T")]


def _trajectory_pooled() -> dict[tuple[str, int], np.ndarray]:
    raw = pd.read_csv(TRAJ_FILE)
    out: dict[tuple[str, int], np.ndarray] = {}
    for trt, _, _ in TREATMENTS:
        for gen in GENS:
            chunks: list[np.ndarray] = []
            for rep in REPS:
                col = f"{trt}{rep}G{gen:02d}.{METRIC}"
                if col in raw.columns:
                    v = pd.to_numeric(raw[col], errors="coerce").dropna().to_numpy()
                    if len(v) > 0:
                        chunks.append(v)
            if chunks:
                out[(trt, gen)] = np.concatenate(chunks)
    return out


def _trajectory_per_rep_medians() -> dict[tuple[str, int], dict[int, float]]:
    raw = pd.read_csv(TRAJ_FILE)
    out: dict[tuple[str, int], dict[int, float]] = {}
    for trt, _, _ in TREATMENTS:
        for rep in REPS:
            d: dict[int, float] = {}
            for gen in GENS:
                col = f"{trt}{rep}G{gen:02d}.{METRIC}"
                if col in raw.columns:
                    v = pd.to_numeric(raw[col], errors="coerce").dropna().to_numpy()
                    if len(v) > 0:
                        d[gen] = float(np.median(v))
            if d:
                out[(trt, rep)] = d
    return out


def _inset_distributions() -> dict[str, np.ndarray]:
    raw = pd.read_csv(DIV_FILE)
    wild_all, founder_all = [], []
    for p in WILD_POOLS:
        col = f"{p}.1.{METRIC}"
        if col in raw.columns:
            wild_all.append(pd.to_numeric(raw[col], errors="coerce").dropna().to_numpy())
    for p in FOUNDER_POOLS:
        col = f"{p}.1.{METRIC}"
        if col in raw.columns:
            founder_all.append(pd.to_numeric(raw[col], errors="coerce").dropna().to_numpy())
    return {
        "Wild":    np.concatenate(wild_all)   if wild_all    else np.array([]),
        "Founder": np.concatenate(founder_all) if founder_all else np.array([]),
    }


def _strip_svg_clips(svg_path: Path) -> None:
    import re
    txt = svg_path.read_text()
    txt = re.sub(r"<clipPath[^>]*>.*?</clipPath>", "", txt, flags=re.DOTALL)
    txt = re.sub(r'\s*clip-path="url\([^)]+\)"', "", txt)
    svg_path.write_text(txt)


def main() -> None:
    plt.rcParams.update({
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.8,
    })

    per_rep = _trajectory_per_rep_medians()
    inset = _inset_distributions()

    fig = plt.figure(figsize=(13, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[5.5, 1.0], wspace=0.20)
    ax = fig.add_subplot(gs[0, 0])
    ax_in = fig.add_subplot(gs[0, 1], sharey=ax)

    for trt, color, label in TREATMENTS:
        xs, means, sems = [], [], []
        for gen in GENS:
            rep_vals = [per_rep[(trt, r)][gen]
                        for r in REPS
                        if (trt, r) in per_rep and gen in per_rep[(trt, r)]]
            if len(rep_vals) < 2:
                continue
            xs.append(gen)
            means.append(float(np.mean(rep_vals)))
            sems.append(float(np.std(rep_vals, ddof=1) / np.sqrt(len(rep_vals))))
        if not xs:
            continue
        xs_a = np.array(xs)
        m_a  = np.array(means)
        s_a  = np.array(sems)
        ax.fill_between(xs_a, m_a - s_a, m_a + s_a,
                        color=color, alpha=0.20, linewidth=0, zorder=1)
        ax.plot(xs_a, m_a, color=color, linewidth=2.2,
                marker="o", markersize=5, markeredgewidth=0.8,
                markeredgecolor="white", label=label, zorder=3)

    ax.axhline(0, color=C_GREY, linestyle=":", linewidth=0.8, zorder=1)
    ax.set_xlabel("Generation", fontsize=11)
    ax.set_ylabel(METRIC_LABEL, fontsize=11)
    ax.set_xlim(0.5, 9.9)
    ax.set_ylim(-1.7, 1.9) 
    ax.tick_params(axis="both", labelsize=9)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    leg = ax.legend(
        title="Treatment", loc="lower left",
        frameon=False, fontsize=9, title_fontsize=9,
        handlelength=2.0, handletextpad=0.6,
    )
    leg._legend_box.align = "left"

    labels = ["Wild", "Founder"]
    data = [inset[k] for k in labels]
    colors = [C_WILD, C_FOUNDER]
    bp = ax_in.boxplot(
        data, positions=[0, 1], widths=0.55,
        patch_artist=True, showfliers=False,
        medianprops=dict(color="black", linewidth=1.2),
        whiskerprops=dict(color="black", linewidth=0.6),
        capprops=dict(color="black", linewidth=0.6),
        boxprops=dict(edgecolor="black", linewidth=0.6),
    )
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(mcolors.to_rgba(c, alpha=0.85))
    ax_in.axhline(0, color=C_GREY, linestyle=":", linewidth=0.8, zorder=0)
    ax_in.set_xticks([0, 1])
    ax_in.set_xticklabels(labels, fontsize=9)
    ax_in.set_xlim(-0.6, 1.6)
    ax_in.tick_params(axis="x", labelsize=9)
    ax_in.tick_params(axis="y", which="both", left=False, labelleft=False)
    ax_in.set_title("Wild vs Founder", fontsize=10, fontweight="bold", pad=4)
    for spine in ("top", "right", "left"):
        ax_in.spines[spine].set_visible(False)

    from matplotlib.ticker import FixedLocator, FixedFormatter
    ax.xaxis.set_major_locator(FixedLocator(GENS))
    ax.xaxis.set_major_formatter(FixedFormatter([str(g) for g in GENS]))

    out_base = ROOT / "final_plots/wild/s3_d_trajectory_with_wildfounder_inset"
    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_base}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_base}.svg", bbox_inches="tight")
    fig.savefig(f"{out_base}.pdf", bbox_inches="tight")
    _strip_svg_clips(Path(f"{out_base}.svg"))
    
if __name__ == "__main__":
    main()
