"""
Compute inter-pool overdispersion from allele frequency comparisons.

For each pool group (e.g., B2G10: p1 vs p2), computes:
  - Between-pool AF variance at each site
  - Expected binomial variance: p(1-p) / (2*N) where N = pool_size_diploid
  - Overdispersion ratio: observed / expected
  - Corrected phi accounting for finite read depth (two-stage sampling)

For founders (F1-F4G00): 4-way comparison, all pairwise.

Output:
  - summary: one row per group with mean phi, median phi, n_eff
  - per_site: per-site overdispersion for diagnostic plots
"""

import polars as pl
import numpy as np
import json
import itertools


def compute_two_stage_expected_var(p, n_haploid, depth):
    """
    Expected AF variance under two-stage sampling (pooling + sequencing).
    E[Var] = p(1-p) * [1/n_haploid + 1/depth - 1/(n_haploid * depth)]
    """
    return p * (1 - p) * (1.0 / n_haploid + 1.0 / depth - 1.0 / (n_haploid * depth))


def compute_binomial_expected_var(p, n_haploid):
    """Expected AF variance from pooling alone: p(1-p) / n_haploid."""
    return p * (1 - p) / n_haploid


def main():
    af_table_path = snakemake.input.af_table
    pool_groups = dict(snakemake.params.pool_groups)
    pair_groups = dict(snakemake.params.pair_groups)
    pool_size_diploid = snakemake.params.pool_size
    n_haploid = 2 * pool_size_diploid  # 160

    summary_path = snakemake.output.summary
    per_site_path = snakemake.output.per_site

    df = pl.read_csv(af_table_path, separator="\t", has_header=True)

    summary_rows = []
    per_site_dfs = []

    # Process pairwise pool groups (G10 samples with p1 vs p2)
    for group, samples in pair_groups.items():
        s1, s2 = sorted(samples)
        af1_col = f"{s1}_af"
        af2_col = f"{s2}_af"
        d1_col = f"{s1}_depth"
        d2_col = f"{s2}_depth"

        if af1_col not in df.columns or af2_col not in df.columns:
            continue

        sub = df.select(["chrom", "pos", af1_col, af2_col, d1_col, d2_col]).drop_nulls()

        af1 = sub[af1_col].to_numpy()
        af2 = sub[af2_col].to_numpy()
        d1 = sub[d1_col].to_numpy().astype(float)
        d2 = sub[d2_col].to_numpy().astype(float)

        # Mean AF across the two pools
        p_mean = (af1 + af2) / 2.0

        # Filter: require both pools have depth >= 10 and 0.05 < p_mean < 0.95
        mask = (d1 >= 10) & (d2 >= 10) & (p_mean > 0.05) & (p_mean < 0.95)
        af1 = af1[mask]
        af2 = af2[mask]
        d1 = d1[mask]
        d2 = d2[mask]
        p_mean = p_mean[mask]
        chrom_arr = sub.filter(pl.lit(mask))["chrom"].to_numpy()
        pos_arr = sub.filter(pl.lit(mask))["pos"].to_numpy()

        if len(af1) == 0:
            continue

        # Observed variance between two pools = (af1 - af2)^2 / 2
        obs_var = (af1 - af2) ** 2 / 2.0

        # Expected: two-stage (accounts for sequencing noise in both pools)
        mean_depth = (d1 + d2) / 2.0
        exp_var_twostage = compute_two_stage_expected_var(p_mean, n_haploid, mean_depth)
        exp_var_binomial = compute_binomial_expected_var(p_mean, n_haploid)

        # Per-site overdispersion
        phi_twostage = np.where(exp_var_twostage > 0, obs_var / exp_var_twostage, np.nan)
        phi_binomial = np.where(exp_var_binomial > 0, obs_var / exp_var_binomial, np.nan)

        # Summary statistics
        valid = np.isfinite(phi_twostage)
        mean_phi = np.nanmean(phi_twostage[valid])
        median_phi = np.nanmedian(phi_twostage[valid])
        n_eff_haploid = n_haploid / mean_phi if mean_phi > 0 else np.nan
        mean_phi_naive = np.nanmean(phi_binomial[valid])

        summary_rows.append({
            "group": group,
            "comparison": f"{s1}_vs_{s2}",
            "n_sites": int(valid.sum()),
            "phi_twostage_mean": round(mean_phi, 4),
            "phi_twostage_median": round(median_phi, 4),
            "phi_naive_mean": round(mean_phi_naive, 4),
            "n_eff_haploid": round(n_eff_haploid, 1),
            "n_eff_diploid": round(n_eff_haploid / 2, 1),
        })

        # Per-site output
        site_df = pl.DataFrame({
            "chrom": chrom_arr,
            "pos": pos_arr,
            "group": [group] * len(af1),
            "af1": af1,
            "af2": af2,
            "p_mean": p_mean,
            "obs_var": obs_var,
            "exp_var_twostage": exp_var_twostage,
            "phi_twostage": phi_twostage,
        })
        per_site_dfs.append(site_df)

    # Process founders: all pairwise comparisons among F1-F4G00
    founder_samples = sorted(pool_groups.get("founders", []))
    if len(founder_samples) >= 2:
        for s1, s2 in itertools.combinations(founder_samples, 2):
            af1_col = f"{s1}_af"
            af2_col = f"{s2}_af"
            d1_col = f"{s1}_depth"
            d2_col = f"{s2}_depth"

            if af1_col not in df.columns or af2_col not in df.columns:
                continue

            sub = df.select(["chrom", "pos", af1_col, af2_col, d1_col, d2_col]).drop_nulls()
            af1 = sub[af1_col].to_numpy()
            af2 = sub[af2_col].to_numpy()
            d1 = sub[d1_col].to_numpy().astype(float)
            d2 = sub[d2_col].to_numpy().astype(float)
            p_mean = (af1 + af2) / 2.0

            mask = (d1 >= 10) & (d2 >= 10) & (p_mean > 0.05) & (p_mean < 0.95)
            af1 = af1[mask]
            af2 = af2[mask]
            d1 = d1[mask]
            d2 = d2[mask]
            p_mean = p_mean[mask]
            chrom_arr = sub.filter(pl.lit(mask))["chrom"].to_numpy()
            pos_arr = sub.filter(pl.lit(mask))["pos"].to_numpy()

            if len(af1) == 0:
                continue

            obs_var = (af1 - af2) ** 2 / 2.0
            mean_depth = (d1 + d2) / 2.0
            exp_var_twostage = compute_two_stage_expected_var(p_mean, n_haploid, mean_depth)
            exp_var_binomial = compute_binomial_expected_var(p_mean, n_haploid)

            phi_twostage = np.where(exp_var_twostage > 0, obs_var / exp_var_twostage, np.nan)
            phi_binomial = np.where(exp_var_binomial > 0, obs_var / exp_var_binomial, np.nan)

            valid = np.isfinite(phi_twostage)
            mean_phi = np.nanmean(phi_twostage[valid])
            median_phi = np.nanmedian(phi_twostage[valid])
            n_eff_haploid = n_haploid / mean_phi if mean_phi > 0 else np.nan
            mean_phi_naive = np.nanmean(phi_binomial[valid])

            summary_rows.append({
                "group": "founders",
                "comparison": f"{s1}_vs_{s2}",
                "n_sites": int(valid.sum()),
                "phi_twostage_mean": round(mean_phi, 4),
                "phi_twostage_median": round(median_phi, 4),
                "phi_naive_mean": round(mean_phi_naive, 4),
                "n_eff_haploid": round(n_eff_haploid, 1),
                "n_eff_diploid": round(n_eff_haploid / 2, 1),
            })

            site_df = pl.DataFrame({
                "chrom": chrom_arr,
                "pos": pos_arr,
                "group": [f"founders_{s1}_vs_{s2}"] * len(af1),
                "af1": af1,
                "af2": af2,
                "p_mean": p_mean,
                "obs_var": obs_var,
                "exp_var_twostage": exp_var_twostage,
                "phi_twostage": phi_twostage,
            })
            per_site_dfs.append(site_df)

    # Write outputs
    summary_df = pl.DataFrame(summary_rows)
    summary_df.write_csv(summary_path, separator="\t")

    if per_site_dfs:
        per_site_all = pl.concat(per_site_dfs)
        per_site_all.write_csv(per_site_path, separator="\t")
    else:
        with open(per_site_path, "w") as f:
            f.write("chrom\tpos\tgroup\taf1\taf2\tp_mean\tobs_var\texp_var_twostage\tphi_twostage\n")


main()
