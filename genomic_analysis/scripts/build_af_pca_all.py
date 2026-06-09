#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf",          required=True, help="Input bgzipped VCF")
    ap.add_argument("--out-prefix",   required=True, help="Output path prefix")
    ap.add_argument("--min-cov",      type=int, default=10,
                    help="Minimum per-sample depth (REF+ALT) to accept a site (default 10)")
    ap.add_argument("--max-cov",      type=int, default=500,
                    help="Maximum per-sample depth to accept a site (default 500)")
    ap.add_argument("--min-af",       type=float, default=0.01,
                    help="Exclude SNPs where mean AF across kept samples is < min-af or > 1-min-af (default 0.01)")
    ap.add_argument("--drop-samples", default="",
                    help="Comma-separated substrings of sample names to exclude")
    ap.add_argument("--max-sites",    type=int, default=0,
                    help="Cap the number of SNPs used for PCA (0 = no cap). Randomly subsamples after filtering.")
    return ap.parse_args()


def vcf_open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def main():
    args = parse_args()
    samples = None
    with vcf_open(args.vcf) as f:
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break
    if samples is None:
        sys.exit("ERROR: #CHROM header not found in VCF")

    drop_patterns = [p.strip() for p in args.drop_samples.split(",") if p.strip()]
    keep_idx = [i for i, s in enumerate(samples)
                if not any(p in s for p in drop_patterns)]
    keep_names = [samples[i] for i in keep_idx]
    n_samples = len(keep_names)
    af_rows = []    
    chrom_pos = []
    t0 = time.time()
    n_lines = 0
    n_kept = 0
    with vcf_open(args.vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            n_lines += 1
            fields = line.rstrip("\n").split("\t")
            chrom, pos = fields[0], fields[1]
            fmt = fields[8].split(":")
            if "AD" not in fmt:
                continue
            ad_idx = fmt.index("AD")
            row = np.full(n_samples, np.nan, dtype=np.float32)
            all_depths = []
            for j, si in enumerate(keep_idx):
                sample_field = fields[9 + si]
                if sample_field == "." or sample_field == "./.":
                    continue
                parts = sample_field.split(":")
                if ad_idx >= len(parts):
                    continue
                ad = parts[ad_idx]
                if ad == "." or "," not in ad:
                    continue
                ad_parts = ad.split(",")
                try:
                    ref = int(ad_parts[0])
                    alt = int(ad_parts[1])
                except Exception:
                    continue
                d = ref + alt
                all_depths.append(d)
                if args.min_cov <= d <= args.max_cov and d > 0:
                    row[j] = alt / d
            if np.isfinite(row).all():
                mean_af = float(np.nanmean(row))
                if args.min_af <= mean_af <= 1 - args.min_af:
                    af_rows.append(row)
                    chrom_pos.append(f"{chrom}:{pos}")
                    n_kept += 1

    if n_kept == 0:
        sys.exit("no SNPs passed filters")

    AF = np.vstack(af_rows)                                          # n_sites x n_samples
    if args.max_sites > 0 and AF.shape[0] > args.max_sites:
        rng = np.random.default_rng(42)
        idx = rng.choice(AF.shape[0], args.max_sites, replace=False)
        AF = AF[idx]
        chrom_pos = [chrom_pos[i] for i in idx]

    p = AF.mean(axis=1, keepdims=True)
    denom = np.sqrt(p * (1 - p))
    denom = np.where(denom == 0, np.nan, denom)
    Z = (AF - p) / denom                                              # n_sites x n_samples
    Z = np.nan_to_num(Z, nan=0.0, posinf=0.0, neginf=0.0)


    X = Z.T                                                           # n_samples x n_sites
    pca = PCA(n_components=min(10, X.shape[0]-1, X.shape[1]), svd_solver="full")
    scores = pca.fit_transform(X)
    pve = pca.explained_variance_ratio_ * 100


    out_prefix = args.out_prefix
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)

    scores_df = pd.DataFrame(scores, columns=[f"PC{i+1}" for i in range(scores.shape[1])])
    scores_df.insert(0, "ind", keep_names)
    scores_df.to_csv(f"{out_prefix}_scores.csv", index=False)

    pd.DataFrame({"PC": np.arange(1, len(pve)+1), "pve": pve}) \
      .to_csv(f"{out_prefix}_eigenvals.csv", index=False)


if __name__ == "__main__":
    main()
