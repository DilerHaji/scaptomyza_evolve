#!/usr/bin/env python3
"""
build_founder_metadata.py

Parse sample IDs from a BAM list into a metadata TSV aligned to the
row order of the BAM list (= row order of downstream ANGSD/PCAngsd output).

Each sample ID is of the form {TUBE}[_idx]_{A|B}_{pos}, where TUBE is the
first underscore-separated token (e.g., B1B7E2, T1T2E3m, FANH). TUBE codes
are looked up in the config tube_metadata map to get host_plant and pool_type.

Output columns:
  sample_id, tube, host_plant, pool_type, is_founder, row_index
"""

import argparse
import os
import sys
import yaml


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bamlist", required=True)
    p.add_argument("--config",  required=True)
    p.add_argument("--output",  required=True)
    return p.parse_args()


def sample_from_bampath(path):
    return os.path.basename(path).replace(".bam", "").replace(".BAM", "")


def tube_from_sample(sample_id):
    return sample_id.split("_")[0]


def main():
    args = parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    tube_meta = cfg.get("tube_metadata", {})

    rows = []
    unknown_tubes = set()

    with open(args.bamlist) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if not line:
                continue
            sample = sample_from_bampath(line)
            tube   = tube_from_sample(sample)
            meta   = tube_meta.get(tube)
            if meta is None:
                unknown_tubes.add(tube)
                host_plant = "unknown"
                pool_type  = "unknown"
            else:
                host_plant = meta.get("host_plant", "unknown")
                pool_type  = meta.get("pool_type",  "unknown")
            is_founder = pool_type != "exp_evolution"
            rows.append((sample, tube, host_plant, pool_type,
                         "TRUE" if is_founder else "FALSE", i))


    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as fh:
        fh.write("sample_id\ttube\thost_plant\tpool_type\tis_founder\trow_index\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    n_founder   = sum(1 for r in rows if r[4] == "TRUE")
    n_nonfound  = len(rows) - n_founder


if __name__ == "__main__":
    main()
