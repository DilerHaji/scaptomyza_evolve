#!/usr/bin/env python3
import argparse
import gzip
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pool-freq", required=True,
                   help="Pool-seq freq table TSV (chrom, pos, ref, alt, F1G00..F4G00)")
    p.add_argument("--pool-depth", required=True,
                   help="Pool-seq depth table TSV (chrom, pos, F1G00_DP..F4G00_DP)")
    p.add_argument("--pool-size", type=int, default=80,
                   help="Number of diploid individuals per pool (default: 80)")
    p.add_argument("--admix-q", default=None,
                   help="PCAngsd admixture Q matrix (to filter structured sites)")
    p.add_argument("--ind-mafs", default=None,
                   help="ANGSD ind mafs.gz (for AF comparison residual filter)")
    p.add_argument("--max-residual", type=float, default=0.15,
                   help="Max |pool - ind| residual to keep site (default: 0.15)")
    p.add_argument("--min-depth", type=int, default=30,
                   help="Min depth per pool to include site (default: 30)")
    p.add_argument("--max-depth-quantile", type=float, default=0.99,
                   help="Exclude sites above this depth quantile (default: 0.99)")
    p.add_argument("--min-maf", type=float, default=0.05,
                   help="Min MAF to include site (default: 0.05)")
    p.add_argument("--output-stats", required=True, help="Output stats TSV")
    p.add_argument("--output-fig", required=True, help="Output diagnostic figure")
    p.add_argument("--dpi", type=int, default=200)
    return p.parse_args()


def two_stage_expected_var(p, n_haploid, coverage):
    # Per-pool variance: C_j = 1/n + 1/R_j - 1/(n*R_j)
    C = 1.0 / n_haploid + 1.0 / coverage - 1.0 / (n_haploid * coverage)
    # Expected between-pool variance = p(1-p) * mean(C)
    return p * (1 - p) * np.mean(C)


def main():
    args = parse_args()
    n_haploid = 2 * args.pool_size  # pool haploid count

    freq_df = pd.read_csv(args.pool_freq, sep="\t")
    pool_cols = [c for c in freq_df.columns if c.startswith("F") and c.endswith("G00")]
    for col in pool_cols:
        freq_df[col] = pd.to_numeric(freq_df[col], errors="coerce")

    depth_df = pd.read_csv(args.pool_depth, sep="\t")
    depth_cols = [c for c in depth_df.columns if "DP" in c or c in pool_cols]
    # Merge on chrom + pos
    df = freq_df.merge(depth_df, on=["chrom", "pos"], how="inner",
                       suffixes=("_freq", "_depth"))

    dp_cols = [c for c in depth_df.columns if c not in ("chrom", "pos")]
    if len(dp_cols) != len(pool_cols):
        print(f"{len(dp_cols)} depth cols vs {len(pool_cols)} freq cols",
              file=sys.stderr)

    freq_matrix = df[pool_cols].values  # (n_sites, n_pools)
    depth_matrix = df[dp_cols].values.astype(float)

    mean_af = np.nanmean(freq_matrix, axis=1)
    maf = np.minimum(mean_af, 1 - mean_af)
    obs_var = np.nanvar(freq_matrix, axis=1, ddof=1)

    exp_var = np.array([
        two_stage_expected_var(mean_af[i], n_haploid, depth_matrix[i])
        for i in range(len(df))
    ])

    naive_var = mean_af * (1 - mean_af) / n_haploid

    df["mean_af"] = mean_af
    df["maf"] = maf
    df["obs_var"] = obs_var
    df["exp_var_two_stage"] = exp_var
    df["exp_var_naive"] = naive_var
    df["mean_depth"] = np.nanmean(depth_matrix, axis=1)
    df["min_depth"] = np.nanmin(depth_matrix, axis=1)
    
    n_total = len(df)
    mask = np.ones(n_total, dtype=bool)

    mask &= maf >= args.min_maf
    mask &= df["min_depth"].values >= args.min_depth
    depth_cutoff = np.quantile(df["mean_depth"].values[mask], args.max_depth_quantile)
    mask &= df["mean_depth"].values <= depth_cutoff
    mask &= np.isfinite(obs_var) & np.isfinite(exp_var) & (exp_var > 0)
    filtered = df[mask].copy()
    phi_mean = filtered["obs_var"].mean() / filtered["exp_var_two_stage"].mean()
    per_site_phi = filtered["obs_var"] / filtered["exp_var_two_stage"]
    phi_median = per_site_phi.median()

    from numpy.linalg import lstsq
    X = filtered["exp_var_two_stage"].values.reshape(-1, 1)
    y = filtered["obs_var"].values
    phi_regression = float(lstsq(X, y, rcond=None)[0][0])

    phi_naive = filtered["obs_var"].mean() / filtered["exp_var_naive"].mean()

    n_eff_mean = n_haploid / phi_mean
    n_eff_median = n_haploid / phi_median
    n_eff_regression = n_haploid / phi_regression

    print(f"\n--- Overdispersion estimates ---", file=sys.stderr)
    print(f"  Naive (pooling only):        phi = {phi_naive:.3f}  "
          f"-> n_eff = {n_haploid/phi_naive:.0f} haploid", file=sys.stderr)
    print(f"  Two-stage mean ratio:        phi = {phi_mean:.3f}  "
          f"-> n_eff = {n_eff_mean:.0f} haploid", file=sys.stderr)
    print(f"  Two-stage median per-site:   phi = {phi_median:.3f}  "
          f"-> n_eff = {n_eff_median:.0f} haploid", file=sys.stderr)
    print(f"  Two-stage regression:        phi = {phi_regression:.3f}  "
          f"-> n_eff = {n_eff_regression:.0f} haploid", file=sys.stderr)

    maf_bins = pd.cut(filtered["maf"], bins=[0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
    cov_bins = pd.cut(filtered["mean_depth"],
                      bins=np.quantile(filtered["mean_depth"], [0, 0.25, 0.5, 0.75, 1.0]))

    maf_phi = filtered.groupby(maf_bins, observed=True).apply(
        lambda g: g["obs_var"].mean() / g["exp_var_two_stage"].mean()
        if g["exp_var_two_stage"].mean() > 0 else np.nan
    )

    cov_phi = filtered.groupby(cov_bins, observed=True).apply(
        lambda g: g["obs_var"].mean() / g["exp_var_two_stage"].mean()
        if g["exp_var_two_stage"].mean() > 0 else np.nan
    )

    print(f"\n--- Overdispersion by MAF bin ---", file=sys.stderr)
    for b, v in maf_phi.items():
        print(f"  {b}: phi = {v:.3f}", file=sys.stderr)

    print(f"\n--- Overdispersion by coverage quartile ---", file=sys.stderr)
    for b, v in cov_phi.items():
        print(f"  {b}: phi = {v:.3f}", file=sys.stderr)

    stats = [
        {"metric": "n_sites_used", "value": len(filtered)},
        {"metric": "n_sites_total", "value": n_total},
        {"metric": "nominal_haploid_pool_size", "value": n_haploid},
        {"metric": "phi_naive_pooling_only", "value": phi_naive},
        {"metric": "phi_two_stage_mean_ratio", "value": phi_mean},
        {"metric": "phi_two_stage_median_per_site", "value": phi_median},
        {"metric": "phi_two_stage_regression", "value": phi_regression},
        {"metric": "n_eff_haploid_mean", "value": n_eff_mean},
        {"metric": "n_eff_haploid_median", "value": n_eff_median},
        {"metric": "n_eff_haploid_regression", "value": n_eff_regression},
        {"metric": "n_eff_naive", "value": n_haploid / phi_naive},
        {"metric": "depth_cutoff_upper", "value": depth_cutoff},
    ]
    for b, v in maf_phi.items():
        stats.append({"metric": f"phi_maf_{b}", "value": v})
    for b, v in cov_phi.items():
        stats.append({"metric": f"phi_cov_{b}", "value": v})

    pd.DataFrame(stats).to_csv(args.output_stats, sep="\t", index=False)
    print(f"\nStats saved to {args.output_stats}", file=sys.stderr)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle(
        f"Pool-seq overdispersion estimation\n"
        f"phi(two-stage, regression) = {phi_regression:.3f}  |  "
        f"n_eff = {n_eff_regression:.0f} haploid  |  "
        f"n_sites = {len(filtered):,}",
        fontsize=11,
    )

    ax = axes[0, 0]
    ax.scatter(filtered["exp_var_two_stage"], filtered["obs_var"],
               s=1, alpha=0.1, rasterized=True)
    lim = max(filtered["exp_var_two_stage"].quantile(0.99),
              filtered["obs_var"].quantile(0.99))
    ax.plot([0, lim], [0, lim], "r--", lw=1, label="1:1 (no overdispersion)")
    ax.plot([0, lim], [0, lim * phi_regression], "orange", lw=1.5,
            label=f"phi = {phi_regression:.2f}")
    ax.set_xlabel("Expected two-stage variance")
    ax.set_ylabel("Observed between-pool variance")
    ax.set_title("Observed vs expected variance")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim * phi_regression * 1.5)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    clipped = per_site_phi.clip(0, 20)
    ax.hist(clipped, bins=100, color="#4C86A8", edgecolor="none", alpha=0.8)
    ax.axvline(phi_median, color="red", lw=1.5, ls="--",
               label=f"Median = {phi_median:.2f}")
    ax.axvline(phi_mean, color="orange", lw=1.5, ls="--",
               label=f"Mean ratio = {phi_mean:.2f}")
    ax.axvline(1.0, color="grey", lw=1, ls=":")
    ax.set_xlabel("Per-site overdispersion (obs/exp)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of per-site phi")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    maf_labels = [str(b) for b in maf_phi.index]
    maf_vals = maf_phi.values
    ax.bar(range(len(maf_labels)), maf_vals, color="#E07B54")
    ax.axhline(phi_regression, color="orange", lw=1.5, ls="--",
               label=f"Genome-wide = {phi_regression:.2f}")
    ax.set_xticks(range(len(maf_labels)))
    ax.set_xticklabels(maf_labels, fontsize=8, rotation=15)
    ax.set_xlabel("MAF bin")
    ax.set_ylabel("Overdispersion (phi)")
    ax.set_title("Overdispersion by MAF")
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    cov_labels = [str(b) for b in cov_phi.index]
    cov_vals = cov_phi.values
    ax.bar(range(len(cov_labels)), cov_vals, color="#6AAB6E")
    ax.axhline(phi_regression, color="orange", lw=1.5, ls="--",
               label=f"Genome-wide = {phi_regression:.2f}")
    ax.set_xticks(range(len(cov_labels)))
    ax.set_xticklabels(cov_labels, fontsize=8, rotation=15)
    ax.set_xlabel("Mean coverage quartile")
    ax.set_ylabel("Overdispersion (phi)")
    ax.set_title("Overdispersion by coverage")
    ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(args.output_fig, dpi=args.dpi, bbox_inches="tight")


if __name__ == "__main__":
    main()
