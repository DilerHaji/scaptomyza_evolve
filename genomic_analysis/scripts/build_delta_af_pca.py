#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf",          required=True)
    ap.add_argument("--out-prefix",   required=True)
    ap.add_argument("--min-cov",      type=int,   default=10)
    ap.add_argument("--max-cov",      type=int,   default=500)
    ap.add_argument("--min-af",       type=float, default=0.05,
                    help="Filter SNPs whose MEAN AF across kept samples is outside [min_af, 1-min_af]")
    ap.add_argument("--drop-samples", default="")
    ap.add_argument("--founder-regex", default=r"^F[1-4](G00)?$",
                    help="Regex identifying founder samples to use as baseline (default matches F1, F1G00, etc.)")
    ap.add_argument("--max-sites",    type=int,   default=0,
                    help="Cap SNP count used for PCA (0 = no cap). Random sub-sample after filter.")
    return ap.parse_args()


def vcf_open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def main():
    args = parse_args()
    founder_re = re.compile(args.founder_regex)

    with vcf_open(args.vcf) as f:
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break

    drop_patterns = [p.strip() for p in args.drop_samples.split(",") if p.strip()]
    keep_idx = [i for i, s in enumerate(samples)
                if not any(p in s for p in drop_patterns)]
    keep_names = [samples[i] for i in keep_idx]
    n_samples = len(keep_names)
    founder_sample_idx = [j for j, s in enumerate(keep_names) if founder_re.match(s)]
    evo_sample_idx     = [j for j in range(n_samples) if j not in founder_sample_idx]
    if len(founder_sample_idx) == 0:
        sys.exit("No founders identified — adjust --founder-regex")

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
            for j, si in enumerate(keep_idx):
                sf = fields[9 + si]
                if sf in (".", "./."):
                    continue
                parts = sf.split(":")
                if ad_idx >= len(parts):
                    continue
                ad = parts[ad_idx]
                if ad == "." or "," not in ad:
                    continue
                try:
                    ref, alt = int(ad.split(",")[0]), int(ad.split(",")[1])
                except Exception:
                    continue
                d = ref + alt
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

    AF = np.vstack(af_rows)  # n_sites × n_samples

    if args.max_sites > 0 and AF.shape[0] > args.max_sites:
        rng = np.random.default_rng(42)
        idx = rng.choice(AF.shape[0], args.max_sites, replace=False)
        AF = AF[idx]
        chrom_pos = [chrom_pos[i] for i in idx]

    founder_mean = AF[:, founder_sample_idx].mean(axis=1, keepdims=True)

    DAF = AF - founder_mean  # n_sites × n_samples; founders near 0 by construction

    evo_names = [keep_names[j] for j in evo_sample_idx]
    DAF_evo   = DAF[:, evo_sample_idx]

    p = founder_mean[:, 0]
    denom = np.sqrt(p * (1 - p))
    denom = np.where(denom == 0, np.nan, denom)
    Z = DAF_evo / denom[:, None]
    Z = np.nan_to_num(Z, nan=0.0, posinf=0.0, neginf=0.0)

    X = Z.T  # n_evo_samples × n_sites
    pca = PCA(n_components=min(10, X.shape[0]-1, X.shape[1]), svd_solver="full")
    scores = pca.fit_transform(X)
    pve = pca.explained_variance_ratio_ * 100

    out_prefix = args.out_prefix
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)

    scores_df = pd.DataFrame(scores, columns=[f"PC{i+1}" for i in range(scores.shape[1])])
    scores_df.insert(0, "ind", evo_names)

    founder_rows = pd.DataFrame(
        np.zeros((len(founder_sample_idx), scores.shape[1]), dtype=np.float32),
        columns=[f"PC{i+1}" for i in range(scores.shape[1])]
    )
    founder_rows.insert(0, "ind", [keep_names[j] for j in founder_sample_idx])
    scores_df = pd.concat([scores_df, founder_rows], ignore_index=True)

    scores_df.to_csv(f"{out_prefix}_scores.csv", index=False)
    pd.DataFrame({"PC": np.arange(1, len(pve)+1), "pve": pve}) \
      .to_csv(f"{out_prefix}_eigenvals.csv", index=False)

if __name__ == "__main__":
    main()
