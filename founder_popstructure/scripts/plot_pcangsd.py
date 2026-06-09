#!/usr/bin/env python3


import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cov",     required=True, help="PCAngsd covariance matrix (.cov)")
    p.add_argument("--admix",   required=True, help="PCAngsd admixture matrix (*.Q)")
    p.add_argument("--bamlist", required=True, help="BAM list used for ANGSD (one per line)")
    p.add_argument("--n-pcs",   type=int, default=10, help="PCs computed (default: 10)")
    p.add_argument("--output",  required=True, help="Output PNG path")
    p.add_argument("--dpi",     type=int, default=200)
    return p.parse_args()


def read_sample_names(bamlist_path):
    names = []
    with open(bamlist_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            basename = os.path.basename(line)
            name = basename.replace(".bam", "").replace(".BAM", "")
            names.append(name)
    return names


def pca_from_cov(cov):
    eigvals, eigvecs = np.linalg.eigh(cov)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    pcs = eigvecs * np.sqrt(np.abs(eigvals))
    pct_var = 100 * eigvals / eigvals.sum()
    return pcs, pct_var


def main():
    args = parse_args()
    cov   = np.loadtxt(args.cov)
    admix = np.loadtxt(args.admix)
    n_samples, K = admix.shape

    sample_names = read_sample_names(args.bamlist)
    if len(sample_names) != n_samples:
        while len(sample_names) < n_samples:
            sample_names.append(f"sample_{len(sample_names)}")
        sample_names = sample_names[:n_samples]

    pcs, pct_var = pca_from_cov(cov)

    sort_idx = np.argsort(admix[:, 0])

    COLORS = ["#E07B54", "#4C86A8", "#6AAB6E", "#C17BC4", "#E8CC5B"]

    fig = plt.figure(figsize=(16, 12))
    gs  = gridspec.GridSpec(2, 2, figure=fig,
                            hspace=0.38, wspace=0.32,
                            left=0.07, right=0.97, top=0.93, bottom=0.07)

    ax_pc12   = fig.add_subplot(gs[0, 0])
    ax_pc13   = fig.add_subplot(gs[0, 1])
    ax_admix  = fig.add_subplot(gs[1, :]) 

    sc = ax_pc12.scatter(
        pcs[:, 0], pcs[:, 1],
        c=admix[:, 0], cmap="RdBu_r",
        vmin=0, vmax=1,
        s=30, alpha=0.8, linewidths=0.5, edgecolors="black",
        zorder=3,
    )

    for i, name in enumerate(sample_names):
        ax_pc12.annotate(
            name, (pcs[i, 0], pcs[i, 1]),
            fontsize=3.5, alpha=0.7,
            xytext=(2, 2), textcoords="offset points",
        )
    ax_pc12.axhline(0, lw=0.5, color="grey", ls="--")
    ax_pc12.axvline(0, lw=0.5, color="grey", ls="--")
    ax_pc12.set_xlabel(f"PC1 ({pct_var[0]:.1f}% var)")
    ax_pc12.set_ylabel(f"PC2 ({pct_var[1]:.1f}% var)")
    ax_pc12.set_title("PC1 vs PC2")
    plt.colorbar(sc, ax=ax_pc12, label="Ancestry component 1", shrink=0.85)



    sc2 = ax_pc13.scatter(
        pcs[:, 0], pcs[:, 2],
        c=admix[:, 0], cmap="RdBu_r",
        vmin=0, vmax=1,
        s=30, alpha=0.8, linewidths=0.5, edgecolors="black",
        zorder=3,
    )
    for i, name in enumerate(sample_names):
        ax_pc13.annotate(
            name, (pcs[i, 0], pcs[i, 2]),
            fontsize=3.5, alpha=0.7,
            xytext=(2, 2), textcoords="offset points",
        )
    ax_pc13.axhline(0, lw=0.5, color="grey", ls="--")
    ax_pc13.axvline(0, lw=0.5, color="grey", ls="--")
    ax_pc13.set_xlabel(f"PC1 ({pct_var[0]:.1f}% var)")
    ax_pc13.set_ylabel(f"PC3 ({pct_var[2]:.1f}% var)")
    ax_pc13.set_title("PC1 vs PC3")
    plt.colorbar(sc2, ax=ax_pc13, label="Ancestry component 1", shrink=0.85)




    sorted_admix = admix[sort_idx]
    sorted_names = [sample_names[i] for i in sort_idx]
    bottom = np.zeros(n_samples)
    x = np.arange(n_samples)
    for k in range(K):
        ax_admix.bar(
            x, sorted_admix[:, k], bottom=bottom,
            color=COLORS[k % len(COLORS)],
            width=1.0, linewidth=0,
            label=f"Component {k+1}",
        )
        bottom += sorted_admix[:, k]

    ax_admix.set_xlim(-0.5, n_samples - 0.5)
    ax_admix.set_ylim(0, 1)
    ax_admix.set_xlabel(f"Individuals (n={n_samples}), sorted by component 1 proportion")
    ax_admix.set_ylabel("Admixture proportion")
    ax_admix.set_title(f"K={K} admixture proportions")
    ax_admix.legend(loc="upper right", framealpha=0.9, ncol=K)

    ax_admix.set_xticks(x)
    ax_admix.set_xticklabels(sorted_names, rotation=90, fontsize=4, ha="center")


    mean_q1 = admix[:, 0].mean()
    std_q1  = admix[:, 0].std()
    n_mostly_comp1 = (admix[:, 0] > 0.7).sum()
    n_mostly_comp2 = (admix[:, 0] < 0.3).sum()
    n_admixed      = n_samples - n_mostly_comp1 - n_mostly_comp2

    summary = (
        f"n={n_samples}  |  mean Q1={mean_q1:.3f} +/- {std_q1:.3f}  |  "
        f">70% comp1: {n_mostly_comp1}  |  "
        f">70% comp2: {n_mostly_comp2}  |  "
        f"admixed (30-70%): {n_admixed}"
    )
    fig.suptitle(f"PCAngsd — Founder individual genotype structure\n{summary}",
                 fontsize=11)

    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")


if __name__ == "__main__":
    main()
