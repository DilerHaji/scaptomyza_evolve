#!/usr/bin/env python3

import sys
import os
import numpy as np
import pandas as pd
from scipy import stats

OUTDIR = "variance_analysis/section1_rigorous"
N_EFF = 29                    # effective diploids per pool (Poisson/founder estimate)
MIN_DEPTH = 10
MIN_MAF = 0.05
BLOCK_SIZE = 100_000          # bp
N_BOOT = 300
SEED = 42

GENS_4REP = [1, 2, 6, 7, 8, 9]
TREATMENTS = ["B", "T", "M"]


def load_data():
    samples = [l.strip() for l in open("variance_analysis/sample_list.txt") if l.strip()]
    chroms, positions, ad_ref, ad_alt = [], [], [], []
    with open("variance_analysis/merged_ad.tsv") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chroms.append(parts[0])
            positions.append(int(parts[1]))
            refs, alts = [], []
            for ad in parts[4:]:
                if ad == "." or ad == ".,.":
                    refs.append(0); alts.append(0)
                else:
                    r, a = ad.split(",")[:2]
                    refs.append(int(r)); alts.append(int(a))
            ad_ref.append(refs); ad_alt.append(alts)
    ad_ref = np.array(ad_ref, dtype=np.int32)
    ad_alt = np.array(ad_alt, dtype=np.int32)
    total = ad_ref + ad_alt
    freq = np.where(total > 0, ad_alt / total, np.nan)
    return samples, np.array(chroms), np.array(positions), freq, total


def build_common_filter(samples, freq, total):
    founder_idx = [samples.index(f"F{i}G00") for i in range(1, 5)]
    founder_freq = np.nanmean(freq[:, founder_idx], axis=1)
    maf_ok = (np.minimum(founder_freq, 1 - founder_freq) >= MIN_MAF) & np.isfinite(founder_freq)

    depth_ok = np.ones(freq.shape[0], dtype=bool)
    for trt in TREATMENTS:
        for gen in GENS_4REP:
            for rep in range(1, 5):
                sname = f"{trt}{rep}G{gen:02d}"
                if sname in samples:
                    depth_ok &= (total[:, samples.index(sname)] >= MIN_DEPTH)

    mask = maf_ok & depth_ok

    return mask


def precompute_site_contributions(samples, freq, total, mask):
    contribs = {}  # (trt, gen) -> dict with arrays
    idx_keep = np.where(mask)[0]
    pos_keep = None  # filled in caller
    for trt in TREATMENTS:
        for gen in GENS_4REP:
            rep_names = [f"{trt}{r}G{gen:02d}" for r in range(1, 5)]
            rep_idx = [samples.index(s) for s in rep_names if s in samples]
            if len(rep_idx) != 4:
                continue
            f_r = freq[idx_keep][:, rep_idx]  # n_keep x 4
            d_r = total[idx_keep][:, rep_idx]
            pbar = f_r.mean(axis=1)
            het = pbar * (1 - pbar)
            num_raw = ((f_r - pbar[:, None])**2).sum(axis=1) / (4 - 1)
            inv_d = (1.0 / d_r).mean(axis=1)
            num_noise = het * (1.0/(2*N_EFF) + inv_d)
            contribs[(trt, gen)] = dict(
                num_raw=num_raw.astype(np.float64),
                num_noise=num_noise.astype(np.float64),
                den=het.astype(np.float64),
            )
    return contribs, idx_keep


def build_blocks(chroms, positions, idx_keep):
    chroms_keep = chroms[idx_keep]
    pos_keep = positions[idx_keep]
    blocks = []  # each block = array of indices into the kept-sites axis
    offset = 0
    for chrom in np.unique(chroms_keep):
        sel = np.where(chroms_keep == chrom)[0]
        if len(sel) == 0:
            continue
        p = pos_keep[sel]
        lo = (p // BLOCK_SIZE).astype(np.int64)
        uniq_bins = np.unique(lo)
        for b in uniq_bins:
            in_bin = sel[lo == b]
            if len(in_bin) >= 10:
                blocks.append(in_bin)
    return blocks


def compute_ne_from_subset(contribs, site_idx):
    out = {}
    for trt in TREATMENTS:
        gens, F_raw, F_bio = [], [], []
        for gen in GENS_4REP:
            c = contribs.get((trt, gen))
            if c is None:
                continue
            num_raw = c["num_raw"][site_idx].sum()
            num_noise = c["num_noise"][site_idx].sum()
            den = c["den"][site_idx].sum()
            if den <= 0:
                continue
            F_raw.append(num_raw / den)
            F_bio.append((num_raw - num_noise) / den)
            gens.append(gen)
        if len(gens) < 3:
            out[trt] = dict(Ne_raw=np.nan, Ne_bio=np.nan,
                            slope_raw=np.nan, slope_bio=np.nan,
                            int_raw=np.nan, int_bio=np.nan)
            continue
        sr, ir, *_ = stats.linregress(gens, F_raw)
        sb, ib, *_ = stats.linregress(gens, F_bio)
        out[trt] = dict(
            Ne_raw=(1/(2*sr) if sr > 0 else np.inf),
            Ne_bio=(1/(2*sb) if sb > 0 else np.inf),
            slope_raw=sr, slope_bio=sb, int_raw=ir, int_bio=ib,
            F_raw=F_raw, F_bio=F_bio, gens=gens,
        )
    return out


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    samples, chroms, positions, freq, total = load_data()
    mask = build_common_filter(samples, freq, total)
    contribs, idx_keep = precompute_site_contributions(samples, freq, total, mask)
    blocks = build_blocks(chroms, positions, idx_keep)
    n_blocks = len(blocks)

    all_sites = np.arange(len(idx_keep))
    point = compute_ne_from_subset(contribs, all_sites)

    rows = []
    for trt in TREATMENTS:
        p = point[trt]
        for g, fr, fb in zip(p["gens"], p["F_raw"], p["F_bio"]):
            rows.append(dict(treatment=trt, generation=g, F_raw=fr, F_bio=fb))
    pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "F_noise_corrected_full.tsv"),
                              sep="\t", index=False)

    for trt in TREATMENTS:
        p = point[trt]

    rng = np.random.RandomState(SEED)
    boots_raw = {trt: [] for trt in TREATMENTS}
    boots_bio = {trt: [] for trt in TREATMENTS}
    for b in range(N_BOOT):
        samp = rng.choice(n_blocks, n_blocks, replace=True)
        site_idx = np.concatenate([blocks[i] for i in samp])
        res = compute_ne_from_subset(contribs, site_idx)
        for trt in TREATMENTS:
            r = res[trt]
            if np.isfinite(r["Ne_raw"]) and r["Ne_raw"] < 1e7:
                boots_raw[trt].append(r["Ne_raw"])
            if np.isfinite(r["Ne_bio"]) and r["Ne_bio"] < 1e7:
                boots_bio[trt].append(r["Ne_bio"])

    summary = []
    for trt in TREATMENTS:
        p = point[trt]
        br = np.array(boots_raw[trt])
        bb = np.array(boots_bio[trt])
        summary.append(dict(
            treatment=trt,
            Ne_raw_point=p["Ne_raw"],
            Ne_raw_lo=np.percentile(br, 2.5) if len(br) else np.nan,
            Ne_raw_hi=np.percentile(br, 97.5) if len(br) else np.nan,
            Ne_bio_point=p["Ne_bio"],
            Ne_bio_lo=np.percentile(bb, 2.5) if len(bb) else np.nan,
            Ne_bio_hi=np.percentile(bb, 97.5) if len(bb) else np.nan,
            slope_raw=p["slope_raw"], int_raw=p["int_raw"],
            slope_bio=p["slope_bio"], int_bio=p["int_bio"],
            n_boot=len(bb),
        ))
    sdf = pd.DataFrame(summary)
    sdf.to_csv(os.path.join(OUTDIR, "F_noise_corrected_bootstrap.tsv"),
               sep="\t", index=False)

if __name__ == "__main__":
    main()
