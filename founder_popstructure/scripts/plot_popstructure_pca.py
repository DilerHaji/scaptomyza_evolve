#!/usr/bin/env python3


import argparse
import os
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

EDGE_COLOR = "#333333"       # dark grey 
POOLTYPE_MARKERS = {
    "single_host":  "o",
    "mixed":        "s",
    "all_combined": "^",
    "exp_evolution":"D",
    "unknown":      "X",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cov",      required=True, help="PCAngsd covariance matrix (.cov)")
    p.add_argument("--metadata", required=True, help="Metadata TSV (row-aligned to cov)")
    p.add_argument("--color-by", default="host_plant",
                   choices=["host_plant", "pool_type", "tube"])
    p.add_argument("--output",   required=True, help="Output PNG")
    p.add_argument("--dpi",      type=int, default=200)
    p.add_argument("--label-tubes", action="store_true",
                   help="Draw tube code next to each point (small font)")
    return p.parse_args()


def pca_from_cov(cov):
    eigvals, eigvecs = np.linalg.eigh(cov)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    pcs     = eigvecs * np.sqrt(np.abs(eigvals))
    pct     = 100 * eigvals / np.abs(eigvals).sum()
    return pcs, pct


def main():
    args = parse_args()
    cov  = np.loadtxt(args.cov)
    meta = pd.read_csv(args.metadata, sep="\t")
    assert len(meta) == cov.shape[0], (
        f"Metadata rows ({len(meta)}) != cov dim ({cov.shape[0]})"
    )

    pcs, pct = pca_from_cov(cov)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2),
                             gridspec_kw={"width_ratios": [1, 1, 0.5]})
    ax12, ax13, axv = axes

    groups = meta[args.color_by].fillna("unknown").astype(str)
    unique = sorted(groups.unique())
    if args.color_by == "host_plant":
        cmap = {g: HOST_COLORS.get(g, "#888888") for g in unique}
    else:
        base = plt.get_cmap("tab10").colors
        cmap = {g: base[i % len(base)] for i, g in enumerate(unique)}

    markers = meta["pool_type"].fillna("unknown").map(
        lambda x: POOLTYPE_MARKERS.get(x, "o")).values

    for g in unique:
        sel = (groups == g).values
        for m in np.unique(markers[sel]):
            idx = sel & (markers == m)
            ax12.scatter(pcs[idx, 0], pcs[idx, 1], c=[cmap[g]], marker=m,
                         s=55, alpha=0.85, edgecolors=EDGE_COLOR, linewidths=0.4,
                         label=f"{g} / {m}" if m != "o" else g)
            ax13.scatter(pcs[idx, 0], pcs[idx, 2], c=[cmap[g]], marker=m,
                         s=55, alpha=0.85, edgecolors=EDGE_COLOR, linewidths=0.4)

    if args.label_tubes:
        for i in range(len(meta)):
            ax12.annotate(meta["tube"].iloc[i], (pcs[i, 0], pcs[i, 1]),
                          fontsize=3.5, alpha=0.6,
                          xytext=(2, 2), textcoords="offset points")
            ax13.annotate(meta["tube"].iloc[i], (pcs[i, 0], pcs[i, 2]),
                          fontsize=3.5, alpha=0.6,
                          xytext=(2, 2), textcoords="offset points")

    for ax in (ax12, ax13):
        ax.axhline(0, lw=0.5, color="grey", ls="--")
        ax.axvline(0, lw=0.5, color="grey", ls="--")
    ax12.set_xlabel(f"PC1 ({pct[0]:.2f}%)")
    ax12.set_ylabel(f"PC2 ({pct[1]:.2f}%)")
    ax13.set_xlabel(f"PC1 ({pct[0]:.2f}%)")
    ax13.set_ylabel(f"PC3 ({pct[2]:.2f}%)")
    ax12.set_title("PC1 vs PC2")
    ax13.set_title("PC1 vs PC3")


    axv.axis("off")
    handles = []
    for g in unique:
        handles.append(plt.Line2D([0], [0], marker="o", linestyle="",
                                  markerfacecolor=cmap[g], markeredgecolor=EDGE_COLOR,
                                  markersize=8, label=str(g)))
    for ptype, mk in POOLTYPE_MARKERS.items():
        if ptype in meta["pool_type"].astype(str).unique():
            handles.append(plt.Line2D([0], [0], marker=mk, linestyle="",
                                      markerfacecolor="white", markeredgecolor=EDGE_COLOR,
                                      markersize=8, label=ptype))
    axv.legend(handles=handles, loc="center left", frameon=False,
               title=args.color_by, fontsize=9)


    axin = fig.add_axes([0.72, 0.58, 0.22, 0.28])
    n_show = min(10, len(pct))
    axin.bar(np.arange(1, n_show + 1), pct[:n_show], color="#555555")
    axin.set_xlabel("PC", fontsize=8)
    axin.set_ylabel("% var", fontsize=8)
    axin.tick_params(labelsize=7)
    axin.set_title("Variance explained", fontsize=9)

    n = cov.shape[0]
    fig.suptitle(f"PCAngsd — n={n}  colored by {args.color_by}", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    save_png_and_svg(fig, args.output, args.dpi)


if __name__ == "__main__":
    main()
