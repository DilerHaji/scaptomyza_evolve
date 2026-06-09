#!/usr/bin/env python3


import argparse
import gzip
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats


# -- Argument parsing ----------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pool-freq",    required=True,
                   help="Pool-seq freq table TSV (from pool_freq_table rule)")
    p.add_argument("--ind-mafs",     required=True,
                   help="ANGSD .mafs.gz from ind_pseudopool_freq rule")
    p.add_argument("--output-fig",   required=True, help="Output PNG figure")
    p.add_argument("--output-stats", required=True, help="Output stats TSV")
    p.add_argument("--min-depth",    type=int, default=20,
                   help="Min nInd at a site in ANGSD output (default: 20)")
    p.add_argument("--dpi",          type=int, default=200)
    return p.parse_args()


# -- Data loading --------------------------------------------------------------

def load_ind_mafs(path, min_depth=20):
    rows = []
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom      = parts[0]
            pos        = int(parts[1])
            major      = parts[2]
            minor      = parts[3]
            ref_allele = parts[4]
            freq_minor = float(parts[5])
            n_ind      = int(parts[7])
            if n_ind < min_depth:
                continue
            rows.append((chrom, pos, major, minor, ref_allele, freq_minor, n_ind))
    df = pd.DataFrame(rows, columns=["chrom", "pos", "major_ind", "minor_ind",
                                     "ref_ind", "freq_minor", "n_ind"])
    df["site_id"] = df["chrom"] + ":" + df["pos"].astype(str)
    return df.set_index("site_id")


def load_pool_freq(path):
    df = pd.read_csv(path, sep="\t")
    pool_cols = [c for c in df.columns if c.startswith("F") and c.endswith("G00")]

    # Convert freq strings to float; "." -> NaN
    for col in pool_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["site_id"] = df["chrom"] + ":" + df["pos"].astype(str)
    df = df.set_index("site_id")
    df["freq_pool_mean"] = df[pool_cols].mean(axis=1)
    df["freq_pool_var"]  = df[pool_cols].var(axis=1, ddof=1)
    return df, pool_cols


def binomial_variance(p, n):
    return p * (1 - p) / (2 * n)


def regression_stats(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 10:
        return dict(r=np.nan, r2=np.nan, slope=np.nan, intercept=np.nan, rmse=np.nan, n=len(x))
    slope, intercept, r, p, se = stats.linregress(x, y)
    resid = y - (slope * x + intercept)
    rmse  = np.sqrt(np.mean(resid ** 2))
    return dict(r=r, r2=r**2, slope=slope, intercept=intercept, rmse=rmse, n=len(x))



def main():
    args = parse_args()
    ind_df = load_ind_mafs(args.ind_mafs, min_depth=args.min_depth)
    pool_df, pool_cols = load_pool_freq(args.pool_freq)
    join_cols = ["chrom", "pos", "ref", "alt", "freq_pool_mean",
                 "freq_pool_var"] + pool_cols
    merged = ind_df.join(pool_df[join_cols], how="inner", rsuffix="_pool")

    if len(merged) == 0:
        print("error",
              file=sys.stderr)
        sys.exit(1)

    minor_is_alt = merged["minor_ind"] == merged["alt"]
    major_is_alt = merged["major_ind"] == merged["alt"]

    merged["freq_ind"] = np.where(
        minor_is_alt, merged["freq_minor"],
        np.where(major_is_alt, 1.0 - merged["freq_minor"], np.nan)
    )

    n_before = len(merged)
    n_matched_minor = minor_is_alt.sum()
    n_matched_major = major_is_alt.sum()
    merged = merged.dropna(subset=["freq_ind"])
    n_dropped = n_before - len(merged)


    POOL_N = 80

    merged["binom_var_expected"] = merged["freq_pool_mean"].apply(
        lambda p: binomial_variance(p, POOL_N) if np.isfinite(p) else np.nan
    )
    merged["residual"]      = merged["freq_pool_mean"] - merged["freq_ind"]
    merged["abs_residual"]  = merged["residual"].abs()

    mean_freq = (merged["freq_ind"] + merged["freq_pool_mean"]) / 2
    merged["maf"] = np.minimum(mean_freq, 1 - mean_freq)
    merged["maf_bin"] = pd.cut(
        merged["maf"],
        bins=[0, 0.1, 0.2, 0.3, 0.4, 0.5],
        labels=["0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5"],
    )

    reg = regression_stats(
        merged["freq_ind"].values,
        merged["freq_pool_mean"].values,
    )

    per_pool_stats = {}
    for col in pool_cols:
        r = regression_stats(merged["freq_ind"].values,
                             merged[col].values)
        r["pool"] = col
        per_pool_stats[col] = r

    stats_rows = []
    stats_rows.append({"metric": "n_shared_sites",     "value": len(merged)})
    stats_rows.append({"metric": "n_minor_is_alt",     "value": int(n_matched_minor)})
    stats_rows.append({"metric": "n_major_is_alt",     "value": int(n_matched_major)})
    stats_rows.append({"metric": "n_allele_mismatch",  "value": n_dropped})
    stats_rows.append({"metric": "overall_r2",         "value": reg["r2"]})
    stats_rows.append({"metric": "overall_slope",      "value": reg["slope"]})
    stats_rows.append({"metric": "overall_intercept",  "value": reg["intercept"]})
    stats_rows.append({"metric": "overall_rmse",       "value": reg["rmse"]})
    stats_rows.append({"metric": "mean_abs_residual",  "value": merged["abs_residual"].mean()})
    stats_rows.append({"metric": "median_abs_residual","value": merged["abs_residual"].median()})
    mean_obs_var  = merged["freq_pool_var"].mean(skipna=True)
    mean_exp_var  = merged["binom_var_expected"].mean(skipna=True)
    stats_rows.append({"metric": "mean_observed_pool_variance",   "value": mean_obs_var})
    stats_rows.append({"metric": "mean_expected_binomial_variance","value": mean_exp_var})
    stats_rows.append({"metric": "variance_inflation_factor",
                       "value": mean_obs_var / mean_exp_var if mean_exp_var > 0 else np.nan})
    for col, r in per_pool_stats.items():
        stats_rows.append({"metric": f"r2_{col}",    "value": r["r2"]})
        stats_rows.append({"metric": f"rmse_{col}",  "value": r["rmse"]})
        stats_rows.append({"metric": f"slope_{col}", "value": r["slope"]})

    pd.DataFrame(stats_rows).to_csv(args.output_stats, sep="\t", index=False)

    if mean_exp_var > 0:
        print(f"variance inflation {mean_obs_var/mean_exp_var:.2f}x")

    fig = plt.figure(figsize=(14, 12))
    gs  = gridspec.GridSpec(2, 2, figure=fig,
                            hspace=0.42, wspace=0.35,
                            left=0.08, right=0.97, top=0.92, bottom=0.07)

    ax1 = fig.add_subplot(gs[0, 0])   # pool mean vs ind
    ax2 = fig.add_subplot(gs[0, 1])   # residual distribution
    ax3 = fig.add_subplot(gs[1, 0])   # observed vs expected variance
    ax4 = fig.add_subplot(gs[1, 1])   # per-pool R2 bar chart

    plot_n = min(50_000, len(merged))
    plot_df = merged.sample(n=plot_n, random_state=42) if len(merged) > plot_n else merged

    sc = ax1.scatter(
        plot_df["freq_ind"],
        plot_df["freq_pool_mean"],
        c=plot_df["n_ind"],
        cmap="viridis_r",
        s=2, alpha=0.3, linewidths=0,
        rasterized=True,
    )
    lims = [0, 1]
    ax1.plot(lims, lims, "r--", lw=1, label="y = x")
    x_line = np.array(lims)
    ax1.plot(x_line, reg["slope"] * x_line + reg["intercept"],
             "orange", lw=1.5, ls="-",
             label=f"OLS slope={reg['slope']:.3f}, R2={reg['r2']:.3f}")
    ax1.set_xlim(0, 1); ax1.set_ylim(0, 1)
    ax1.set_xlabel("Individual-seq ALT allele frequency")
    ax1.set_ylabel("Pool-seq mean ALT allele frequency (F1-F4G00)")
    ax1.set_title(f"Individual vs pool-seq AF\n(n={len(merged):,} sites, subsample shown)")
    ax1.legend(fontsize=8)
    plt.colorbar(sc, ax=ax1, label="nInd (ANGSD)", shrink=0.85)

   
   
   
    resid = merged["residual"].dropna().values
    ax2.hist(resid, bins=100, color="#4C86A8", edgecolor="none", alpha=0.8)
    ax2.axvline(0, color="red", lw=1, ls="--")
    ax2.axvline(np.percentile(resid, 2.5),  color="grey", lw=1, ls=":")
    ax2.axvline(np.percentile(resid, 97.5), color="grey", lw=1, ls=":",
                label=f"95% CI: [{np.percentile(resid,2.5):.3f}, {np.percentile(resid,97.5):.3f}]")
    ax2.set_xlabel("Pool-seq - individual-seq frequency")
    ax2.set_ylabel("Count")
    ax2.set_title(f"Residual distribution\nRMSE={reg['rmse']:.4f}, "
                  f"median|d|={np.median(np.abs(resid)):.4f}")
    ax2.legend(fontsize=8)

   
   
   
    var_data = merged[["maf_bin", "freq_pool_var", "binom_var_expected"]].dropna()
    maf_bins = var_data["maf_bin"].cat.categories.tolist()
    obs_means = [var_data.loc[var_data["maf_bin"] == b, "freq_pool_var"].mean()
                 for b in maf_bins]
    exp_means = [var_data.loc[var_data["maf_bin"] == b, "binom_var_expected"].mean()
                 for b in maf_bins]
    ns        = [int((var_data["maf_bin"] == b).sum()) for b in maf_bins]

    x_pos = np.arange(len(maf_bins))
    w = 0.35
    ax3.bar(x_pos - w/2, obs_means, width=w, color="#E07B54",
            label="Observed (between-pool var)")
    ax3.bar(x_pos + w/2, exp_means, width=w, color="#4C86A8",
            label=f"Expected binomial (N=80)")
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels([f"{b}\n(n={n:,})" for b, n in zip(maf_bins, ns)],
                        fontsize=8)
    ax3.set_xlabel("MAF bin (minor allele frequency)")
    ax3.set_ylabel("Allele frequency variance")
    ax3.set_title("Observed vs expected pool-seq variance\n"
                  f"Inflation factor: {mean_obs_var/mean_exp_var:.2f}x  "
                  f"(obs/exp, genome-wide mean)")
    ax3.legend(fontsize=8)

    
    
    
    
    for i, (o, e) in enumerate(zip(obs_means, exp_means)):
        if e and e > 0 and np.isfinite(o) and np.isfinite(e):
            ax3.text(i, max(o, e) * 1.05, f"{o/e:.1f}x",
                     ha="center", va="bottom", fontsize=7, color="black")




    pool_names = pool_cols
    r2_vals  = [per_pool_stats[p]["r2"]   for p in pool_names]
    rmse_vals = [per_pool_stats[p]["rmse"] for p in pool_names]

    colors = ["#E07B54", "#4C86A8", "#6AAB6E", "#C17BC4"]
    bars = ax4.bar(pool_names, r2_vals, color=colors[:len(pool_names)], alpha=0.85)
    ax4_r = ax4.twinx()
    ax4_r.plot(pool_names, rmse_vals, "k^--", ms=8, lw=1.5, label="RMSE")
    ax4_r.set_ylabel("RMSE (allele frequency)")

    for bar, r2 in zip(bars, r2_vals):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                 f"R2={r2:.3f}", ha="center", va="bottom", fontsize=9)

    ax4.set_ylim(0, 1.1)
    ax4.set_xlabel("Founder pool")
    ax4.set_ylabel("R2 vs individual-seq pseudo-pool")
    ax4.set_title("Per-pool correlation with\nindividual-seq frequencies")
    ax4_r.legend(loc="upper right", fontsize=8)

    fig.suptitle(
        "Pool-seq vs individual-seq allele frequency comparison — Founder population\n"
        f"N_ind=192 (~0.5x each)  |  N_pool=80/pool  |  "
        f"Shared sites: {len(merged):,}  |  Overall R2={reg['r2']:.3f}",
        fontsize=11,
    )

    fig.savefig(args.output_fig, dpi=args.dpi, bbox_inches="tight")


if __name__ == "__main__":
    main()
