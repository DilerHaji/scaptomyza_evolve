#!/usr/bin/env python3
import argparse
import sys
import os
import warnings
import numpy as np
import pandas as pd
import polars as pl
import statsmodels.api as sm
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.tools.sm_exceptions import ConvergenceWarning
import concurrent.futures

warnings.simplefilter('ignore', ConvergenceWarning)
warnings.simplefilter('ignore', UserWarning)

def fit_metric_slope(pdf, metric_col):
    if metric_col not in pdf.columns: return None
    if pdf[metric_col].dtype == 'object':
        pdf[metric_col] = pd.to_numeric(pdf[metric_col], errors='coerce')

    sub_df = pdf.dropna(subset=[metric_col, "gen", "lineage"])
    

    sub_df = sub_df[~np.isinf(sub_df[metric_col])]
    
    if len(sub_df) < 5 or sub_df["gen"].nunique() < 2: return None

    gen_vals = sub_df["gen"].values
    lineage_vals = sub_df["lineage"].values
    y_vals = sub_df[metric_col].values

    gen_mean = np.mean(gen_vals); gen_std = np.std(gen_vals)
    if gen_std == 0: return None
    gen_scaled = (gen_vals - gen_mean) / gen_std
    
    if np.std(y_vals) == 0: return None
    
    try:
        exog = sm.add_constant(gen_scaled)
        model = MixedLM(endog=y_vals, exog=exog, groups=lineage_vals)
        res = model.fit(reml=False)
        if not res.converged: return None
        return (res.params[1] / gen_std), res.tvalues[1], res.pvalues[1]
    except:
        return None

def process_site_group(group_data):
    (chrom, start, end), pdf = group_data
    stats = {
        "chrom": chrom, "start": start, "end": end,
        "slope_divergence": np.nan, "z_divergence": np.nan, "p_divergence": np.nan,
        "slope_pbs": np.nan, "z_pbs": np.nan, "p_pbs": np.nan,
        "slope_pbs_ref": np.nan, "z_pbs_ref": np.nan, "p_pbs_ref": np.nan,
        "slope_pbsn1": np.nan, "z_pbsn1": np.nan, "p_pbsn1": np.nan,
        "slope_pbsn1_ref": np.nan, "z_pbsn1_ref": np.nan, "p_pbsn1_ref": np.nan,
        "slope_pbe": np.nan, "z_pbe": np.nan, "p_pbe": np.nan,
        "slope_pbe_ref": np.nan, "z_pbe_ref": np.nan, "p_pbe_ref": np.nan
    }
    
    has_data = False
    metrics = [
        ("total_divergence", "divergence"),
        ("pbs_target", "pbs"),
        ("pbs_ref", "pbs_ref"),
        ("pbsn1_target", "pbsn1"),
        ("pbsn1_ref", "pbsn1_ref"),
        ("pbe_target", "pbe"),
        ("pbe_ref", "pbe_ref")
    ]

    for col, suff in metrics:
        res = fit_metric_slope(pdf, col)
        if res:
            stats[f"slope_{suff}"], stats[f"z_{suff}"], stats[f"p_{suff}"] = res
            has_data = True

    return stats if has_data else None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--chroms", required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    header = "chrom,start,end," \
             "slope_divergence,z_divergence,p_divergence," \
             "slope_pbs,z_pbs,p_pbs," \
             "slope_pbs_ref,z_pbs_ref,p_pbs_ref," \
             "slope_pbsn1,z_pbsn1,p_pbsn1," \
             "slope_pbsn1_ref,z_pbsn1_ref,p_pbsn1_ref," \
             "slope_pbe,z_pbe,p_pbe," \
             "slope_pbe_ref,z_pbe_ref,p_pbe_ref\n"
    
    try:
        lf = pl.scan_csv(args.input)
    except Exception:
        with open(args.output, "w") as f: f.write(header)
        return

    with open(args.output, "w") as f: f.write(header)
    
    chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]
    
    for chrom in chroms:
        try:
            df_chrom = lf.filter(pl.col("chrom") == chrom).collect()
        except: continue
        
        if df_chrom.height == 0: continue

        pdf = df_chrom.to_pandas()
        grouped = list(pdf.groupby(["chrom", "start", "end"]))
        
        results = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
            for stats in executor.map(process_site_group, grouped, chunksize=200):
                if stats: results.append(stats)
        
        if results:
            res_df = pl.DataFrame(results).select([
                "chrom", "start", "end", 
                "slope_divergence", "z_divergence", "p_divergence", 
                "slope_pbs", "z_pbs", "p_pbs",
                "slope_pbs_ref", "z_pbs_ref", "p_pbs_ref",
                "slope_pbsn1", "z_pbsn1", "p_pbsn1",
                "slope_pbsn1_ref", "z_pbsn1_ref", "p_pbsn1_ref",
                "slope_pbe", "z_pbe", "p_pbe",
                "slope_pbe_ref", "z_pbe_ref", "p_pbe_ref"
            ])
            with open(args.output, "a") as f:
                res_df.write_csv(f, include_header=False)

if __name__ == "__main__":
    main()