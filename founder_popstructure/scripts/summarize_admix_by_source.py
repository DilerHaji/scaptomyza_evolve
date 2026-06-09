#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--admix",    required=True, help="PCAngsd .Q file for a single K")
    p.add_argument("--metadata", required=True)
    p.add_argument("--output-prefix", required=True,
                   help="Output files: {prefix}_by_host.tsv, {prefix}_by_tube.tsv, {prefix}_per_sample.tsv")
    return p.parse_args()


def main():
    args = parse_args()
    q    = np.loadtxt(args.admix)
    if q.ndim == 1:
        q = q[:, None]
    meta = pd.read_csv(args.metadata, sep="\t")
    assert len(meta) == q.shape[0], (
        f"Metadata rows ({len(meta)}) != Q rows ({q.shape[0]})"
    )

    K = q.shape[1]
    comp_cols = [f"Q{c + 1}" for c in range(K)]
    df = meta.copy()
    for c, col in enumerate(comp_cols):
        df[col] = q[:, c]

    out_dir = os.path.dirname(args.output_prefix) or "."
    os.makedirs(out_dir, exist_ok=True)

    # Per-sample table
    df.to_csv(f"{args.output_prefix}_per_sample.tsv", sep="\t", index=False)

    # By host_plant
    by_host = (df.groupby("host_plant")[comp_cols]
                 .agg(["mean", "std", "count"]))
    by_host.columns = [f"{c}_{stat}" for c, stat in by_host.columns]
    by_host.to_csv(f"{args.output_prefix}_by_host.tsv", sep="\t")

    # By pool_type
    by_pool = (df.groupby("pool_type")[comp_cols]
                 .agg(["mean", "std", "count"]))
    by_pool.columns = [f"{c}_{stat}" for c, stat in by_pool.columns]
    by_pool.to_csv(f"{args.output_prefix}_by_pool_type.tsv", sep="\t")

    # By tube (sample sizes vary, but useful for spot checks)
    by_tube = (df.groupby(["tube", "host_plant", "pool_type"])[comp_cols]
                 .agg(["mean", "std", "count"]))
    by_tube.columns = [f"{c}_{stat}" for c, stat in by_tube.columns]
    by_tube.to_csv(f"{args.output_prefix}_by_tube.tsv", sep="\t")

    # Summary to stdout
    print(f"n samples = {len(df)}, K = {K}")
    print("Mean Q by host_plant:")
    print(df.groupby("host_plant")[comp_cols].mean().round(3))
    print("\nMean Q by pool_type:")
    print(df.groupby("pool_type")[comp_cols].mean().round(3))


if __name__ == "__main__":
    main()
