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


WILD_SAMPLES = ["AVB", "AVT", "PSB", "PST", "RMB", "RMT"]


_RE_START = re.compile(r"\^.")
_RE_INDEL = re.compile(r"[+-](\d+)")


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf",          required=True)
    ap.add_argument("--mpileup-dir",  required=True,
                    help="Directory with {sample}.mpileup files for wild pools")
    ap.add_argument("--out-prefix",   required=True)
    ap.add_argument("--min-cov",      type=int,   default=10)
    ap.add_argument("--max-cov",      type=int,   default=500)
    ap.add_argument("--min-af",       type=float, default=0.01)
    ap.add_argument("--wild-min-cov", type=int,   default=5,
                    help="Minimum per-sample depth for wild pools (default 5; wild coverage is ~38x)")
    ap.add_argument("--wild-samples", default=",".join(WILD_SAMPLES),
                    help="Comma-separated wild sample names. Defaults to the 6 core wild pools.")
    ap.add_argument("--drop-samples", default="")
    ap.add_argument("--max-sites",    type=int, default=0)
    ap.add_argument("--joint-fit",    action="store_true",
                    help="Fit PCA on experimental + wild JOINTLY (default: fit on "
                         "experimental only and project wild). Joint fit reveals "
                         "wild-vs-lab divergence on its own axis.")
    return ap.parse_args()


def vcf_open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def strip_pileup_string(bases: str) -> str:
    bases = _RE_START.sub("", bases)
    bases = bases.replace("$", "")
    while True:
        m = _RE_INDEL.search(bases)
        if not m:
            break
        n = int(m.group(1))
        bases = bases[: m.start()] + bases[m.end() + n:]
    return bases


def count_alleles(ref: str, alt: str, raw_bases: str) -> tuple[int, int]:
    bases = strip_pileup_string(raw_bases)
    ref_count = bases.count(".") + bases.count(",")
    alt_count = bases.count(alt.upper()) + bases.count(alt.lower())
    return ref_count, alt_count


def parse_mpileup_at_positions(mp_path: str, pos_to_refalt: dict) -> dict:
    out = {}
    opener = gzip.open if mp_path.endswith(".gz") else open
    with opener(mp_path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, pos = parts[0], parts[1]
            key = (chrom, pos)
            if key not in pos_to_refalt:
                continue
            depth = int(parts[3])
            if depth == 0:
                out[key] = (0, 0)
                continue
            ref, alt = pos_to_refalt[key]
            r, a = count_alleles(ref, alt, parts[4])
            out[key] = (a, r + a)
    return out


def main():
    args = parse_args()
    t0 = time.time()
    with vcf_open(args.vcf) as f:
        samples = None
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break
    if samples is None:
        sys.exit("ERROR: #CHROM header not found")

    drop_patterns = [p.strip() for p in args.drop_samples.split(",") if p.strip()]
    keep_idx = [i for i, s in enumerate(samples) if not any(p in s for p in drop_patterns)]
    exp_names = [samples[i] for i in keep_idx]
    n_exp = len(exp_names)

    wild_names = [s.strip() for s in args.wild_samples.split(",") if s.strip()]
    n_wild = len(wild_names)

    pos_to_refalt = {}       # (chrom, pos_str) -> (ref, alt)
    exp_af_rows   = []       # per-site n_exp-length float array
    chrom_pos_list = []      # per-site (chrom, pos_str)

    n_lines = 0
    n_kept = 0
    with vcf_open(args.vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            n_lines += 1

            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            if len(ref) != 1 or "," in alt or len(alt) != 1:
                continue  # only biallelic SNPs (needed for wild mpileup allele counting)

            fmt = fields[8].split(":")
            if "AD" not in fmt:
                continue
            ad_idx = fmt.index("AD")

            row = np.full(n_exp, np.nan, dtype=np.float32)
            ok = True
            for j, si in enumerate(keep_idx):
                sf = fields[9 + si]
                if sf in (".", "./."):
                    ok = False
                    break
                parts = sf.split(":")
                if ad_idx >= len(parts):
                    ok = False
                    break
                ad = parts[ad_idx]
                if ad == "." or "," not in ad:
                    ok = False
                    break
                try:
                    r_, a_ = ad.split(",")[:2]
                    r_, a_ = int(r_), int(a_)
                except Exception:
                    ok = False
                    break
                d = r_ + a_
                if not (args.min_cov <= d <= args.max_cov) or d == 0:
                    ok = False
                    break
                row[j] = a_ / d
            if not ok:
                continue
            mean_af = float(row.mean())
            if not (args.min_af <= mean_af <= 1 - args.min_af):
                continue
            pos_to_refalt[(chrom, pos)] = (ref, alt)
            exp_af_rows.append(row)
            chrom_pos_list.append((chrom, pos))
            n_kept += 1

    if n_kept == 0:
        sys.exit("no SNPs passed experimental filters")

    AF_exp = np.vstack(exp_af_rows)   # n_sites x n_exp

    wild_counts_per_sample = {}
    mpileup_dir = Path(args.mpileup_dir)
    for s in wild_names:
        mp = mpileup_dir / f"{s}.mpileup"
        if not mp.exists():
            sys.exit(f"ERROR: wild mpileup not found: {mp}")
        wild_counts_per_sample[s] = parse_mpileup_at_positions(str(mp), pos_to_refalt)

    wild_af_rows = []
    keep_mask = np.zeros(n_kept, dtype=bool)
    for i, key in enumerate(chrom_pos_list):
        row = np.full(n_wild, np.nan, dtype=np.float32)
        ok = True
        for j, s in enumerate(wild_names):
            ad = wild_counts_per_sample[s].get(key)
            if ad is None:
                ok = False
                break
            alt_c, tot = ad
            if tot < args.wild_min_cov:
                ok = False
                break
            row[j] = alt_c / tot
        if ok:
            keep_mask[i] = True
            wild_af_rows.append(row)

    n_joint = int(keep_mask.sum())


    AF_exp   = AF_exp[keep_mask]              # n_joint x n_exp
    AF_wild  = np.vstack(wild_af_rows)        # n_joint x n_wild
    chrom_pos_list = [cp for cp, k in zip(chrom_pos_list, keep_mask) if k]

    if args.max_sites > 0 and AF_exp.shape[0] > args.max_sites:
        rng = np.random.default_rng(42)
        idx = rng.choice(AF_exp.shape[0], args.max_sites, replace=False)
        AF_exp = AF_exp[idx]
        AF_wild = AF_wild[idx]
        chrom_pos_list = [chrom_pos_list[i] for i in idx]

    if args.joint_fit:
        AF_all = np.hstack([AF_exp, AF_wild])                 # n_sites x (n_exp + n_wild)
        p = AF_all.mean(axis=1, keepdims=True)
        mode_str = "joint (exp + wild)"
    else:
        p = AF_exp.mean(axis=1, keepdims=True)
        mode_str = "experimental-only (wild projected)"
    scale = np.sqrt(p * (1 - p))
    scale = np.where(scale == 0, np.nan, scale)

    Z_exp  = np.nan_to_num((AF_exp  - p) / scale, nan=0.0, posinf=0.0, neginf=0.0)
    Z_wild = np.nan_to_num((AF_wild - p) / scale, nan=0.0, posinf=0.0, neginf=0.0)

    X_exp  = Z_exp.T
    X_wild = Z_wild.T

    if args.joint_fit:
        X = np.vstack([X_exp, X_wild])                        # (n_exp + n_wild) x n_sites
        n_comp = min(10, X.shape[0] - 1, X.shape[1])
        pca = PCA(n_components=n_comp, svd_solver="full")
        scores_all = pca.fit_transform(X)
        scores_exp  = scores_all[:n_exp]
        scores_wild = scores_all[n_exp:]
    else:
        n_comp = min(10, X_exp.shape[0] - 1, X_exp.shape[1])
        pca = PCA(n_components=n_comp, svd_solver="full")
        scores_exp  = pca.fit_transform(X_exp)
        scores_wild = pca.transform(X_wild)

    pve = pca.explained_variance_ratio_ * 100

    cols = [f"PC{i+1}" for i in range(scores_exp.shape[1])]
    df_exp  = pd.DataFrame(scores_exp,  columns=cols)
    df_exp.insert(0, "ind", exp_names)
    df_wild = pd.DataFrame(scores_wild, columns=cols)
    df_wild.insert(0, "ind", wild_names)
    df_all = pd.concat([df_exp, df_wild], ignore_index=True)

    out_prefix = args.out_prefix
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)

    df_all.to_csv(f"{out_prefix}_scores.csv", index=False)
    pd.DataFrame({"PC": np.arange(1, len(pve) + 1), "pve": pve}) \
        .to_csv(f"{out_prefix}_eigenvals.csv", index=False)

if __name__ == "__main__":
    main()
