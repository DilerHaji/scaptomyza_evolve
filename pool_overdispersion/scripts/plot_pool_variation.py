"""
Generate diagnostic figures for inter-pool variation analysis.

Figures:
  1. AF scatter: p1 vs p2 for each pool group (+ R^2, RMSE)
  2. Overdispersion (phi) distribution across groups
  3. Phi vs MAF bin to check for structure contamination
  4. Summary bar chart of n_eff across all groups
"""

import polars as pl
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import os
import itertools


def main():
    summary_path = snakemake.input.summary
    per_site_path = snakemake.input.per_site
    af_table_path = snakemake.input.af_table
    pool_groups = dict(snakemake.params.pool_groups)
    pair_groups = dict(snakemake.params.pair_groups)
    outdir = snakemake.params.outdir

    os.makedirs(outdir, exist_ok=True)

    summary = pl.read_csv(summary_path, separator="\t")
    per_site = pl.read_csv(per_site_path, separator="\t")
    af_table = pl.read_csv(af_table_path, separator="\t")

    # ---- Figure 1: AF scatter plots (p1 vs p2) for each pair group ----
    n_pairs = len(pair_groups) + 1  # +1 for founders
    ncols = min(4, n_pairs)
    nrows = (n_pairs + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows), squeeze=False)
    axes_flat = axes.flatten()

    idx = 0
    for group, samples in sorted(pair_groups.items()):
        s1, s2 = sorted(samples)
        af1_col = f"{s1}_af"
        af2_col = f"{s2}_af"
        if af1_col in af_table.columns and af2_col in af_table.columns:
            sub = af_table.select([af1_col, af2_col]).drop_nulls()
            x = sub[af1_col].to_numpy()
            y = sub[af2_col].to_numpy()
            mask = np.isfinite(x) & np.isfinite(y)
            x, y = x[mask], y[mask]

            ax = axes_flat[idx]
            ax.hexbin(x, y, gridsize=80, cmap="inferno", mincnt=1)
            ax.plot([0, 1], [0, 1], "r--", alpha=0.5)
            r2 = np.corrcoef(x, y)[0, 1] ** 2
            rmse = np.sqrt(np.mean((x - y) ** 2))
            ax.set_title(f"{group}\nR²={r2:.3f}, RMSE={rmse:.4f}", fontsize=10)
            ax.set_xlabel(s1)
            ax.set_ylabel(s2)
            idx += 1

    # Founders: pick one representative pair
    founder_samples = sorted(pool_groups.get("founders", []))
    if len(founder_samples) >= 2:
        s1, s2 = founder_samples[0], founder_samples[1]
        af1_col = f"{s1}_af"
        af2_col = f"{s2}_af"
        if af1_col in af_table.columns and af2_col in af_table.columns:
            sub = af_table.select([af1_col, af2_col]).drop_nulls()
            x = sub[af1_col].to_numpy()
            y = sub[af2_col].to_numpy()
            mask = np.isfinite(x) & np.isfinite(y)
            x, y = x[mask], y[mask]

            ax = axes_flat[idx]
            ax.hexbin(x, y, gridsize=80, cmap="inferno", mincnt=1)
            ax.plot([0, 1], [0, 1], "r--", alpha=0.5)
            r2 = np.corrcoef(x, y)[0, 1] ** 2
            rmse = np.sqrt(np.mean((x - y) ** 2))
            ax.set_title(f"Founders ({s1} vs {s2})\nR²={r2:.3f}, RMSE={rmse:.4f}", fontsize=10)
            ax.set_xlabel(s1)
            ax.set_ylabel(s2)
            idx += 1

    # Hide unused axes
    for i in range(idx, len(axes_flat)):
        axes_flat[i].set_visible(False)

    plt.tight_layout()
    fig.savefig(os.path.join(outdir, "af_scatter_p1_vs_p2.png"), dpi=150)
    plt.close(fig)

    # ---- Figure 2: Phi distribution by group ----
    fig, ax = plt.subplots(figsize=(10, 5))
    groups = per_site["group"].unique().to_list()
    phi_data = []
    labels = []
    for g in sorted(groups):
        vals = per_site.filter(pl.col("group") == g)["phi_twostage"].to_numpy()
        vals = vals[np.isfinite(vals)]
        # Clip extreme values for visualization
        vals = vals[(vals > 0) & (vals < np.percentile(vals, 99))]
        phi_data.append(vals)
        labels.append(g)

    ax.boxplot(phi_data, labels=labels, showfliers=False)
    ax.axhline(y=1.0, color="red", linestyle="--", alpha=0.5, label="No overdispersion")
    ax.set_ylabel("Phi (overdispersion)")
    ax.set_title("Inter-pool overdispersion by group")
    ax.legend()
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, "phi_distribution_by_group.png"), dpi=150)
    plt.close(fig)

    # ---- Figure 3: Phi vs MAF bin ----
    fig, ax = plt.subplots(figsize=(8, 5))
    per_site_np = per_site.to_pandas()
    per_site_np["maf_bin"] = np.minimum(per_site_np["p_mean"], 1 - per_site_np["p_mean"])
    per_site_np["maf_bin"] = (per_site_np["maf_bin"] * 20).astype(int) / 20  # 0.05 bins

    maf_phi = per_site_np.groupby("maf_bin")["phi_twostage"].median().reset_index()
    ax.plot(maf_phi["maf_bin"], maf_phi["phi_twostage"], "o-", markersize=4)
    ax.axhline(y=1.0, color="red", linestyle="--", alpha=0.5)
    ax.set_xlabel("Minor allele frequency bin")
    ax.set_ylabel("Median phi (two-stage)")
    ax.set_title("Overdispersion vs MAF\n(flat = pure technical; rising at extremes = structure)")
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, "phi_vs_maf.png"), dpi=150)
    plt.close(fig)

    # ---- Figure 4: Summary n_eff bar chart ----
    fig, ax = plt.subplots(figsize=(10, 5))
    comparisons = summary["comparison"].to_list()
    n_eff = summary["n_eff_diploid"].to_numpy()
    colors = ["#2196F3" if "founders" not in c else "#FF9800" for c in comparisons]

    ax.bar(range(len(comparisons)), n_eff, color=colors)
    ax.axhline(y=80, color="red", linestyle="--", alpha=0.5, label="Nominal (80 diploid)")
    ax.set_xticks(range(len(comparisons)))
    ax.set_xticklabels(comparisons, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Effective pool size (diploid)")
    ax.set_title("Effective pool size across comparisons")
    ax.legend()
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, "n_eff_summary.png"), dpi=150)
    plt.close(fig)


main()
