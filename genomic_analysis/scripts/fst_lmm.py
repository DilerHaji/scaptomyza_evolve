import argparse
import sys
import warnings
import gc
import os
import numpy as np
import polars as pl
import statsmodels.api as sm
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.tools.sm_exceptions import ConvergenceWarning
import concurrent.futures

warnings.simplefilter('ignore', ConvergenceWarning)
warnings.simplefilter('ignore', UserWarning)

def fit_lmm_fast(group_data):
    (chrom, start, end), pdf = group_data
    
    pdf = pdf.dropna(subset=["fst", "gen", "rep_pair"])
    
    if len(pdf) < 5: return None
    if pdf["gen"].nunique() < 2: return None
    if pdf["rep_pair"].nunique() < 2: return None
    
    fst_vals = pdf["fst"].values
    gen_vals = pdf["gen"].values
    rep_vals = pdf["rep_pair"].values
    
    if np.std(fst_vals) == 0: return None

    gen_std = np.std(gen_vals)
    if gen_std == 0: return None
    gen_mean = np.mean(gen_vals)
    gen_scaled = (gen_vals - gen_mean) / gen_std

    try:
        exog = sm.add_constant(gen_scaled)
        
        model = MixedLM(endog=fst_vals, exog=exog, groups=rep_vals)
        result = model.fit(reml=False)

        if not result.converged:
            return None

        slope_raw = result.params[1] / gen_std
        se_raw = result.bse[1] / gen_std
        z = result.tvalues[1]
        p = result.pvalues[1]

        return {
            "chrom": chrom,
            "start": start,
            "end": end,
            "slope_lmm": slope_raw,
            "se_lmm": se_raw,
            "z_score": z,
            "p_value": p
        }

    except Exception:
        return None

def process_chromosome(chrom, lazy_frame, output_file, is_first_write, max_workers):
    try:
        df_chrom = lazy_frame.filter(pl.col("chrom") == chrom).collect()
    except Exception as e:
        return is_first_write

    if df_chrom.height == 0:
        return is_first_write

    pdf = df_chrom.to_pandas()
    
    grouped = list(pdf.groupby(["chrom", "start", "end"]))
    
    del df_chrom
    del pdf
    gc.collect()

    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        processed_stats = executor.map(fit_lmm_fast, grouped, chunksize=50)
        
        for res in processed_stats:
            if res:
                results.append(res)

    if results:
        res_df = pl.DataFrame(results)
        with open(output_file, "a") as f:
            if is_first_write:
                res_df.write_csv(f)
                is_first_write = False
            else:
                res_df.write_csv(f, include_header=False)
    
    del grouped
    del results
    gc.collect()
    
    return is_first_write

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--chroms", required=True)
    parser.add_argument("--threads", type=int, default=os.cpu_count(), help="Number of CPU cores to use")
    args = parser.parse_args()

    target_chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]

    if os.path.exists(args.output):
        os.remove(args.output)
    

    q = pl.scan_csv(args.input, null_values=["NA", "nan", ".", "inf", "-inf"])

    is_first_write = True

    for chrom in target_chroms:
        is_first_write = process_chromosome(
            chrom, q, args.output, is_first_write, args.threads
        )

    if is_first_write:
        with open(args.output, "w") as f:
            f.write("chrom,start,end,slope_lmm,se_lmm,z_score,p_value\n")

if __name__ == "__main__":
    main()
