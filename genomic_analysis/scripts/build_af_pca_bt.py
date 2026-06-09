#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf",            required=True)
    ap.add_argument("--out-prefix",     required=True)
    ap.add_argument("--treatments",     default="B,T",
                    help="Comma-separated treatment prefixes to keep (default 'B,T'). "
                         "Valid: B, T, M, F. E.g. 'B,T,F' to include founders as anchors.")
    ap.add_argument("--min-cov",        type=int,   default=10)
    ap.add_argument("--max-cov",        type=int,   default=500)
    ap.add_argument("--min-af",         type=float, default=0.01)
    ap.add_argument("--window-size",    type=int,   default=50000,
                    help="Window size in bp for SNP-density filter (default 50000).")
    ap.add_argument("--top-density-pct", type=float, default=25.0,
                    help="Keep SNPs in the top N%% of windows by SNP count "
                         "(default 25). Set 100 to disable density filter.")
    ap.add_argument("--drop-samples",   default="")
    ap.add_argument("--max-sites",      type=int, default=0)
    ap.add_argument("--project-treatments", default="",
                    help="Additional VCF-resident treatments to PROJECT onto the "
                         "B-T axes without including them in the fit (e.g. 'M').")
    ap.add_argument("--project-wild-mpileup-dir", default="",
                    help="If set, parses wild mpileups at the kept SNP positions "
                         "and projects those 6 pools onto the B-T axes.")
    ap.add_argument("--wild-samples", default="AVB,AVT,PSB,PST,RMB,RMT")
    ap.add_argument("--wild-min-cov", type=int, default=5)
    return ap.parse_args()


def vcf_open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def classify_treatment(name: str):
    m = re.match(r"^([BTM])([1-4])G(\d{2})$", name)
    if m:
        return m.group(1)
    m = re.match(r"^F([1-4])(G00)?$", name)
    if m:
        return "F"
    return None


_RE_START = re.compile(r"\^.")
_RE_INDEL = re.compile(r"[+-](\d+)")


def _strip_pileup_string(bases: str) -> str:
    bases = _RE_START.sub("", bases)
    bases = bases.replace("$", "")
    while True:
        m = _RE_INDEL.search(bases)
        if not m:
            break
        k = int(m.group(1))
        bases = bases[: m.start()] + bases[m.end() + k:]
    return bases


def _count_alleles(ref: str, alt: str, raw_bases: str):
    bases = _strip_pileup_string(raw_bases)
    rc = bases.count(".") + bases.count(",")
    ac = bases.count(alt.upper()) + bases.count(alt.lower())
    return rc, ac


def parse_wild_mpileup(mp_path: str, pos_to_refalt: dict) -> dict:
    out = {}
    opener = gzip.open if mp_path.endswith(".gz") else open
    with opener(mp_path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            key = (parts[0], parts[1])
            if key not in pos_to_refalt:
                continue
            depth = int(parts[3])
            if depth == 0:
                out[key] = (0, 0)
                continue
            ref, alt = pos_to_refalt[key]
            rc, ac = _count_alleles(ref, alt, parts[4])
            out[key] = (ac, rc + ac)
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

    wanted = {t.strip() for t in args.treatments.split(",") if t.strip()}
    project_set = {t.strip() for t in args.project_treatments.split(",") if t.strip()}
    drop_patterns = [p.strip() for p in args.drop_samples.split(",") if p.strip()]

    keep_idx, keep_names = [], []
    proj_idx, proj_names = [], []
    for i, s in enumerate(samples):
        if any(p in s for p in drop_patterns):
            continue
        trt = classify_treatment(s)
        if trt is None:
            continue
        if trt in wanted:
            keep_idx.append(i)
            keep_names.append(s)
        elif trt in project_set:
            proj_idx.append(i)
            proj_names.append(s)

    n = len(keep_idx)
    n_proj = len(proj_idx)
    if n < 3:
        sys.exit(f"ERROR: only {n} samples matched treatments={wanted}; need ≥3 for PCA")

    wild_samples = ([s.strip() for s in args.wild_samples.split(",") if s.strip()]
                    if args.project_wild_mpileup_dir else [])

    do_density = args.top_density_pct < 100.0
    kept_windows = None
    if do_density:
        window_counts: dict = defaultdict(int)
        n_scan = 0
        with vcf_open(args.vcf) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                n_scan += 1
                fields = line.split("\t", 5)
                chrom, pos_str, _, ref, alt = fields[0], fields[1], None, fields[3], fields[4].split(",", 1)[0]
                if len(ref) != 1 or len(alt) != 1:
                    continue
                pos = int(pos_str)
                window_counts[(chrom, pos // args.window_size)] += 1
        counts = np.fromiter(window_counts.values(), dtype=np.int32)
        threshold = float(np.percentile(counts, 100 - args.top_density_pct))
        kept_windows = {w for w, c in window_counts.items() if c >= threshold}

    af_fit_rows  = []
    af_proj_rows = []         # only allocated if n_proj > 0
    pos_to_refalt = {}        # (chrom, pos_str) -> (ref, alt) — needed for wild parse
    chrom_pos = []
    n_lines = 0
    n_kept  = 0

    with vcf_open(args.vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            n_lines += 1
            fields = line.rstrip("\n").split("\t")
            chrom, pos_str, ref, alt = fields[0], fields[1], fields[3], fields[4]
            if len(ref) != 1 or "," in alt or len(alt) != 1:
                continue
            if kept_windows is not None:
                pos = int(pos_str)
                if (chrom, pos // args.window_size) not in kept_windows:
                    continue

            fmt = fields[8].split(":")
            if "AD" not in fmt:
                continue
            ad_idx = fmt.index("AD")

            fit_row = np.full(n, np.nan, dtype=np.float32)
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
                fit_row[j] = a_ / d
            if not ok:
                continue
            mean_af = float(fit_row.mean())
            if not (args.min_af <= mean_af <= 1 - args.min_af):
                continue

            if n_proj > 0:
                proj_row = np.full(n_proj, np.nan, dtype=np.float32)
                proj_ok = True
                for j, si in enumerate(proj_idx):
                    sf = fields[9 + si]
                    if sf in (".", "./."):
                        proj_ok = False
                        break
                    parts = sf.split(":")
                    if ad_idx >= len(parts):
                        proj_ok = False
                        break
                    ad = parts[ad_idx]
                    if ad == "." or "," not in ad:
                        proj_ok = False
                        break
                    try:
                        r_, a_ = ad.split(",")[:2]
                        r_, a_ = int(r_), int(a_)
                    except Exception:
                        proj_ok = False
                        break
                    d = r_ + a_
                    if not (args.min_cov <= d <= args.max_cov) or d == 0:
                        proj_ok = False
                        break
                    proj_row[j] = a_ / d
                if not proj_ok:
                    continue
                af_proj_rows.append(proj_row)

            af_fit_rows.append(fit_row)
            pos_to_refalt[(chrom, pos_str)] = (ref, alt)
            chrom_pos.append((chrom, pos_str))
            n_kept += 1

    if n_kept == 0:
        sys.exit("no SNPs passed filters")

    AF_fit  = np.vstack(af_fit_rows)                         # n_sites x n
    AF_proj = np.vstack(af_proj_rows) if n_proj else None    # n_sites x n_proj

    if args.max_sites > 0 and AF_fit.shape[0] > args.max_sites:
        rng = np.random.default_rng(42)
        idx = rng.choice(AF_fit.shape[0], args.max_sites, replace=False)
        AF_fit = AF_fit[idx]
        if AF_proj is not None:
            AF_proj = AF_proj[idx]
        chrom_pos = [chrom_pos[i] for i in idx]
        pos_to_refalt = {chrom_pos[k]: pos_to_refalt[chrom_pos[k]]
                         for k in range(len(chrom_pos))}

    AF_wild = None
    if wild_samples:
        wild_counts = {}
        for s in wild_samples:
            mp = Path(args.project_wild_mpileup_dir) / f"{s}.mpileup"
            if not mp.exists():
                sys.exit(f"ERROR: wild mpileup not found: {mp}")
            wild_counts[s] = parse_wild_mpileup(str(mp), pos_to_refalt)

        n_before = len(chrom_pos)
        AF_wild_full = np.full((n_before, len(wild_samples)),
                               np.nan, dtype=np.float32)
        for i, key in enumerate(chrom_pos):
            for j, s in enumerate(wild_samples):
                ad = wild_counts[s].get(key)
                if ad is None:
                    continue
                alt_c, tot = ad
                if tot >= args.wild_min_cov:
                    AF_wild_full[i, j] = alt_c / tot
        keep_mask = (~np.isnan(AF_wild_full)).all(axis=1)
        n_after = int(keep_mask.sum())
        AF_fit  = AF_fit[keep_mask]
        if AF_proj is not None:
            AF_proj = AF_proj[keep_mask]
        AF_wild = AF_wild_full[keep_mask]
        chrom_pos = [cp for cp, k in zip(chrom_pos, keep_mask) if k]

    p = AF_fit.mean(axis=1, keepdims=True)
    scale = np.sqrt(p * (1 - p))
    scale = np.where(scale == 0, np.nan, scale)
    Z_fit = np.nan_to_num((AF_fit - p) / scale, nan=0.0, posinf=0.0, neginf=0.0)

    X_fit = Z_fit.T
    n_comp = min(10, X_fit.shape[0] - 1, X_fit.shape[1])
    pca = PCA(n_components=n_comp, svd_solver="full")
    scores_fit = pca.fit_transform(X_fit)
    pve = pca.explained_variance_ratio_ * 100

    cols = [f"PC{i+1}" for i in range(scores_fit.shape[1])]
    pieces = [pd.DataFrame(scores_fit, columns=cols).assign(ind=keep_names)]


    if AF_proj is not None:
        Z_proj = np.nan_to_num((AF_proj - p) / scale, nan=0.0, posinf=0.0, neginf=0.0)
        scores_proj = pca.transform(Z_proj.T)
        pieces.append(pd.DataFrame(scores_proj, columns=cols).assign(ind=proj_names))

    if AF_wild is not None:
        Z_wild = np.nan_to_num((AF_wild - p) / scale, nan=0.0, posinf=0.0, neginf=0.0)
        scores_wild = pca.transform(Z_wild.T)
        pieces.append(pd.DataFrame(scores_wild, columns=cols).assign(ind=wild_samples))

    df_all = pd.concat(pieces, ignore_index=True)
    df_all = df_all[["ind"] + cols]


    out_prefix = args.out_prefix
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    df_all.to_csv(f"{out_prefix}_scores.csv", index=False)
    pd.DataFrame({"PC": np.arange(1, len(pve) + 1), "pve": pve}) \
        .to_csv(f"{out_prefix}_eigenvals.csv", index=False)



if __name__ == "__main__":
    main()
