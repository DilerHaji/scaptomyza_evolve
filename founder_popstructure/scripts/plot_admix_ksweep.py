#!/usr/bin/env python3

import argparse
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none" 
import matplotlib.pyplot as plt


def save_png_and_svg(fig, out_png, dpi):
    base, _ = os.path.splitext(out_png)
    fig.savefig(base + ".png", dpi=dpi, bbox_inches="tight")
    for ax in fig.get_axes():
        for artist in ax.get_children():
            try:
                artist.set_clip_on(False)
            except (AttributeError, TypeError):
                pass
    fig.savefig(base + ".svg", bbox_inches=None)


HOST_COLORS = {
    "B":       "#D55E00",   # vermillion
    "T":       "#0072B2",   # deep blue
    "both":    "#CC79A7",   # reddish purple
    "unknown": "#BBBBBB",   # light grey
}


COMP_PALETTE = [
    "#0072B2",  # deep blue
    "#D55E00",  # vermillion
    "#F0E442",  # bright yellow
    "#CC79A7",  # reddish purple
    "#56B4E9",  # sky blue
    "#E69F00",  # amber orange
    "#9E8C4D",  # olive (warm tan)
    "#BBBBBB",  # light grey
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--admix-files", nargs="+", required=True,
                   help="K=<N>=<path> entries, one per K")
    p.add_argument("--metadata", required=True)
    p.add_argument("--output",   required=True)
    p.add_argument("--dpi",      type=int, default=200)
    return p.parse_args()


def load_admix(entries):
    """Return dict {K: np.array(n, K)} sorted by K."""
    out = {}
    for e in entries:
        m = re.match(r"^(\d+)=(.+)$", e)
        if not m:
            raise SystemExit(f"Bad --admix-files entry: {e}")
        k = int(m.group(1))
        q = np.loadtxt(m.group(2))
        if q.ndim == 1:
            q = q[:, None]
        out[k] = q
    return dict(sorted(out.items()))


def main():
    args = parse_args()
    ks = load_admix(args.admix_files)
    meta = pd.read_csv(args.metadata, sep="\t")

    n = len(meta)
    order = np.arange(n)

    n_k = len(ks)
    fig, axes = plt.subplots(
        n_k + 1, 1,
        figsize=(max(10, 0.08 * n), 1.4 * n_k + 1.2),
        sharex=True,
        gridspec_kw={"height_ratios": [0.4] + [1] * n_k, "hspace": 0.12},
    )
    axes = np.atleast_1d(axes)
    ax_top   = axes[0]
    ax_admix = axes[1:]

    host = meta["host_plant"].values[order]
    for i, h in enumerate(host):
        ax_top.add_patch(plt.Rectangle(
            (i - 0.5, 0), 1, 1,
            color=HOST_COLORS.get(h, "#CCCCCC"),
            linewidth=0,
        ))
    ax_top.set_xlim(-0.5, n - 0.5)
    ax_top.set_ylim(0, 1)
    ax_top.set_yticks([])
    ax_top.set_ylabel("host", rotation=0, ha="right", va="center", fontsize=8)

    for ax, (k, q) in zip(ax_admix, ks.items()):
        q_o = q[order]
        bottom = np.zeros(n)
        x = np.arange(n)
        for c in range(k):
            ax.bar(x, q_o[:, c], bottom=bottom,
                   color=COMP_PALETTE[c % len(COMP_PALETTE)],
                   width=1.0, linewidth=0)
            bottom += q_o[:, c]
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.5, 1])
        ax.set_ylabel(f"K={k}", rotation=0, ha="right", va="center", fontsize=10)

    sorted_tubes = meta["tube"].values[order]
    bounds = [0]
    for i in range(1, n):
        if sorted_tubes[i] != sorted_tubes[i - 1]:
            bounds.append(i)
    bounds.append(n)
    tick_pos = [(bounds[i] + bounds[i + 1] - 1) / 2 for i in range(len(bounds) - 1)]
    tick_labels = [sorted_tubes[bounds[i]] for i in range(len(bounds) - 1)]
    ax_admix[-1].set_xticks(tick_pos)
    ax_admix[-1].set_xticklabels(tick_labels, rotation=90, fontsize=6)
    ax_admix[-1].set_xlabel(f"Individuals (n={n}), in bamlist order")

    fig.suptitle("PCAngsd admixture — K-sweep", fontsize=12)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    save_png_and_svg(fig, args.output, args.dpi)


if __name__ == "__main__":
    main()
