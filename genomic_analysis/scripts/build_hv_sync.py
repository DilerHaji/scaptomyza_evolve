#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path


ALLELE_INDEX = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
GENS_AVAILABLE = [0, 1, 2, 6, 7, 8, 9]  # intersection across all replicates
REPLICATES = [1, 2, 3, 4]


def parse_sample_list(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def get_sample_indices(sample_names, treatment, gens=GENS_AVAILABLE, reps=REPLICATES):
    name_to_idx = {name: i for i, name in enumerate(sample_names)}

    indices = []
    labels = []

    if len(treatment) == 1:
        trt_list = [treatment]
        reps_per_trt = reps
    elif len(treatment) == 2:
        trt_list = [treatment[0], treatment[1]]
        reps_per_trt = reps
    else:
        raise ValueError(f"Treatment must be 1 or 2 letters, got '{treatment}'")

    for trt in trt_list:
        for r in reps_per_trt:
            name = f"F{r}G00"
            if name not in name_to_idx:
                raise ValueError(f"Founder sample {name} not found in sample list")
            indices.append(name_to_idx[name])
            labels.append(name)

    for g in gens:
        if g == 0:
            continue
        for trt in trt_list:
            for r in reps_per_trt:
                name = f"{trt}{r}G{g:02d}"
                if name not in name_to_idx:
                    raise ValueError(f"Sample {name} not found in sample list")
                indices.append(name_to_idx[name])
                labels.append(name)

    return indices, labels


def counts_to_sync(ref_allele, alt_allele, ref_count, alt_count):
    counts = [0, 0, 0, 0, 0, 0]  # A, T, C, G, N, del
    ref_upper = ref_allele.upper()
    alt_upper = alt_allele.upper()

    if ref_upper in ALLELE_INDEX:
        counts[ALLELE_INDEX[ref_upper]] = ref_count
    if alt_upper in ALLELE_INDEX:
        counts[ALLELE_INDEX[alt_upper]] = alt_count

    return ':'.join(str(c) for c in counts)


def main():
    parser = argparse.ArgumentParser(description="Build sync file for haplovalidate")
    parser.add_argument('--treatment', required=True,
                        help='Treatment to extract: B, T, M (single) or BT, BM, TM (pairwise)')
    parser.add_argument('--scaffold', required=True,
                        help='Scaffold to subset (exact match on col 1), or "ALL" for all sites')
    parser.add_argument('--ad', required=True,
                        help='Path to merged_ad.tsv')
    parser.add_argument('--sample_list', required=True,
                        help='Path to sample_list.txt')
    parser.add_argument('--out', required=True,
                        help='Output sync file path')
    parser.add_argument('--gens', default=','.join(str(g) for g in GENS_AVAILABLE),
                        help='Comma-separated generation numbers (default: 0,1,2,6,7,8,9)')
    args = parser.parse_args()

    gens = [int(g) for g in args.gens.split(',')]

    sample_names = parse_sample_list(args.sample_list)
    indices, labels = get_sample_indices(sample_names, args.treatment, gens=gens)


    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    scaffold_set = set(args.scaffold.split(',')) if args.scaffold != "ALL" else set()

    n_sites = 0
    n_skipped = 0

    with open(args.ad) as fin, open(args.out, 'w') as fout:
        for line in fin:
            fields = line.rstrip('\n').split('\t')
            chrom = fields[0]

            if args.scaffold != "ALL" and chrom not in scaffold_set:
                continue

            pos = fields[1]
            ref_allele = fields[2]
            alt_allele = fields[3]
            sample_data = fields[4:]  # "ref_count,alt_count" per sample

            sync_cols = []
            skip = False
            for idx in indices:
                parts = sample_data[idx].split(',')
                ref_count = int(parts[0])
                alt_count = int(parts[1])
                sync_cols.append(counts_to_sync(ref_allele, alt_allele,
                                                ref_count, alt_count))

            fout.write(f"{chrom}\t{pos}\t{ref_allele.upper()}\t")
            fout.write('\t'.join(sync_cols))
            fout.write('\n')
            n_sites += 1

if __name__ == '__main__':
    main()
