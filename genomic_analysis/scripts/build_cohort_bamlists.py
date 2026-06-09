#!/usr/bin/env python3

import argparse
import os
import sys
import csv


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--samples_tsv", required=True)
    p.add_argument("--bam_dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--bam_suffix", default=".bam",
                   help="File suffix appended to sample name (default: .bam)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    samples = []
    with open(args.samples_tsv) as f:
        first = f.readline()
        delim = "\t" if "\t" in first else ","
        f.seek(0)
        reader = csv.reader(f, delimiter=delim)
        header = next(reader)
        sample_idx = 0
        for i, h in enumerate(header):
            if h.lower() in ("sample", "sample_id", "id", "name"):
                sample_idx = i
                break
        for row in reader:
            if not row or not row[sample_idx].strip():
                continue
            samples.append(row[sample_idx].strip())

    founders = []
    t2g07 = []
    for s in samples:
        if s.startswith("T2G07_"):
            t2g07.append(s)
        else:
            founders.append(s)

    def write_bamlist(cohort_name, sample_list):
        out_path = os.path.join(args.outdir, f"{cohort_name}.bamlist")
        missing = []
        with open(out_path, "w") as out:
            for s in sample_list:
                bam = os.path.join(args.bam_dir, f"{s}{args.bam_suffix}")
                if not os.path.exists(bam):
                    missing.append(bam)
                out.write(bam + "\n")

    write_bamlist("founders", founders)
    write_bamlist("T2G07", t2g07)


if __name__ == "__main__":
    main()
