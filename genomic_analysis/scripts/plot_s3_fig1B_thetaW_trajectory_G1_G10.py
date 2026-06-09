#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(".")
TRAJ_FILE = ROOT / "grenfst/diversity/trajectory_pi_200000_nonoverlap.csv"

METRIC = "theta_watterson"
METRIC_LABEL = r"$\theta_W$"

GENS = list(range(1, 11))   # G01 through G10
REPS = [1, 2, 3, 4]

TREATMENTS = [
    ("M", "#9BAB96", "B+T"),
    ("B", "#499FFF", "B"),
    ("T", "#EDB72D", "T"),
]

C_GREY = "#555555"


def _per_rep_means(raw: pd.DataFrame) -> dict:
    out: dict[tuple[str, int], dict[int, float]] = {}
    for trt, _, _ in TREATMENTS:
        for rep in REPS:
            d: dict[int, float] = {}
            for gen in GENS:
                col = f"{trt}{rep}G{gen:02d}.{METRIC}"
                if col in raw.columns:
                    v = pd.to_numeric(raw[col], errors="coerce").dropna().to_numpy()
                    if len(v) > 0:
                        d[gen] = float(np.mean(v))
            if d:
                out[(trt, rep)] = d
    return out


def _strip_svg_clips(svg_path: Path) -> None:
    import re
    txt = svg_path.read_text()
    txt = re.sub(r"<clipPath[^>]*>.*?</clipPath>", "", txt, flags=re.DOTALL)
    txt = re.sub(r'\s*clip-path="url\([^)]+\)"', "", txt)
    svg_path.write_text(txt)


def main() -> None:
    plt.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.8,
    })

    raw = pd.read_csv(TRAJ_FILE)
    per_rep = _per_rep_means(raw)

    fig, ax = plt.subplots(figsize=(5.8, 4.3))

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
        ribbon = ax.fill_between(xs_a, m_a - s_a, m_a + s_a,
                                  color=color, alpha=0.22, linewidth=0, zorder=1)
        ribbon.set_gid(f"ribbon_{trt}")
        ln, = ax.plot(xs_a, m_a, color=color, linewidth=2.2,
                       marker="o", markersize=5, markeredgewidth=0.8,
                       markeredgecolor="white", label=label, zorder=3)
        ln.set_gid(f"mean_{trt}")

    ax.set_xlabel("Generation", fontsize=11)
    ax.set_ylabel(METRIC_LABEL, fontsize=11)
    ax.set_xlim(0.5, 10.5)
    ax.set_xticks(GENS)
    ax.tick_params(axis="both", labelsize=9)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    leg = ax.legend(title="Treatment", loc="upper right",
                    frameon=False, fontsize=9, title_fontsize=9,
                    handlelength=2.0, handletextpad=0.6)
    leg._legend_box.align = "left"

    out_base = ROOT / "final_plots/wild/s3_fig1B_thetaW_trajectory_G1_G10"
    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_base}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_base}.svg", bbox_inches="tight")
    fig.savefig(f"{out_base}.pdf", bbox_inches="tight")
    _strip_svg_clips(Path(f"{out_base}.svg"))

if __name__ == "__main__":
    main()
