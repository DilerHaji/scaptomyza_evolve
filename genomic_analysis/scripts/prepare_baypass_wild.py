from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from pathlib import Path

CORE_SAMPLES = ["AVB", "AVT", "PSB", "PST", "RMB", "RMT"]
UNSPEC_SAMPLES = ["S.flavaMA", "S.flavaAZ"]

POOL_SIZES = {
    "AVB": 50, "AVT": 50,
    "PSB": 50, "PST": 50,
    "RMB": 50, "RMT": 50,
    "S.flavaMA": 50, "S.flavaAZ": 50,
}

TREATMENT_COV = {
    "AVB":  1, "PSB":  1, "RMB":  1,
    "AVT": -1, "PST": -1, "RMT": -1,
    "S.flavaMA": 0, "S.flavaAZ": 0,
}

_RE_START  = re.compile(r'\^.')          # read-start marker + mapq char
_RE_INDEL  = re.compile(r'[+-](\d+)')   # indel length marker


def _strip_pileup_string(bases: str) -> str:
    bases = _RE_START.sub('', bases)
    bases = bases.replace('$', '')
    while True:
        m = _RE_INDEL.search(bases)
        if not m:
            break
        n = int(m.group(1))
        bases = bases[:m.start()] + bases[m.end() + n:]
    return bases


def count_alleles(ref: str, alt: str, raw_bases: str) -> tuple[int, int]:
    bases = _strip_pileup_string(raw_bases)
    ref_count = bases.count('.') + bases.count(',')
    alt_up    = alt.upper()
    alt_lo    = alt.lower()
    alt_count = bases.count(alt_up) + bases.count(alt_lo)
    return ref_count, alt_count

def read_vcf_positions(vcf_path: str) -> dict:
    positions = {}
    opener = gzip.open if vcf_path.endswith('.gz') else open
    with opener(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.split('\t', 5)
            chrom, pos, _, ref, alt = fields[0], int(fields[1]), None, fields[3], fields[4].strip().split(',')[0]
            if len(ref) != 1 or len(alt) != 1:
                continue  # skip indels / multi-char alleles
            if ',' in fields[4]:
                continue  # skip multi-allelic
            positions[(chrom, pos)] = (ref, alt)
    return positions

def parse_mpileup(mpileup_path: str, positions: dict) -> dict:
    counts = {}
    opener = gzip.open if mpileup_path.endswith('.gz') else open
    with opener(mpileup_path, 'rt') as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            chrom, pos = parts[0], int(parts[1])
            key = (chrom, pos)
            if key not in positions:
                continue
            depth = int(parts[3])
            if depth == 0:
                counts[key] = (0, 0, 0)
                continue
            ref, alt = positions[key]
            raw_bases = parts[4]
            ref_c, alt_c = count_alleles(ref, alt, raw_bases)
            counts[key] = (ref_c, alt_c, depth)
    return counts

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--vcf', required=True,
                   help='VCF (gz or plain) used to define variant positions and REF/ALT alleles')
    p.add_argument('--mpileup-dir', required=True,
                   help='Directory containing {sample}.mpileup files')
    p.add_argument('--output-dir', required=True)
    p.add_argument('--prefix', default='wild')
    p.add_argument('--include-unspec', action='store_true',
                   help='Include S.flavaMA and S.flavaAZ as unspecialized populations')
    p.add_argument('--min-cov', type=int, default=5,
                   help='Minimum total read depth per sample per site (default: 5)')
    p.add_argument('--min-pop-cov', type=float, default=1.0,
                   help='Fraction of populations that must pass --min-cov to keep a site (default: 1.0 = all)')
    p.add_argument('--thin-step', type=int, default=16,
                   help='Keep every Nth SNP for the Omega estimation subset')
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    samples = CORE_SAMPLES + (UNSPEC_SAMPLES if args.include_unspec else [])
    n_pop = len(samples)

    positions = read_vcf_positions(args.vcf)
    ordered_pos = sorted(positions.keys())   # (chrom, pos) in genome order
    n_snps = len(ordered_pos)

    sample_counts = {}   # sample -> { (chrom,pos): (ref, alt, depth) }
    mpileup_dir = Path(args.mpileup_dir)
    for s in samples:
        mp_path = mpileup_dir / f"{s}.mpileup"
        if not mp_path.exists():
            sys.exit(f"ERROR: mpileup not found: {mp_path}")
        sample_counts[s] = parse_mpileup(str(mp_path), positions)
        found = len(sample_counts[s])

    min_pops_needed = max(1, round(args.min_pop_cov * n_pop))

    kept_pos = []
    ref_mat  = []   # list of n_pop-length lists
    alt_mat  = []

    for key in ordered_pos:
        row_ref = []
        row_alt = []
        pops_passing = 0
        for s in samples:
            counts = sample_counts[s].get(key, (0, 0, 0))
            r, a, d = counts
            if d >= args.min_cov:
                pops_passing += 1
            row_ref.append(r)
            row_alt.append(a)
        if pops_passing < min_pops_needed:
            continue
        kept_pos.append(key)
        ref_mat.append(row_ref)
        alt_mat.append(row_alt)

    n_keep = len(kept_pos)

    thin_idx = list(range(0, n_keep, args.thin_step))

    def write_geno(path, indices):
        with open(path, 'w') as fh:
            for i in indices:
                pairs = [f"{ref_mat[i][j]} {alt_mat[i][j]}" for j in range(n_pop)]
                fh.write(' '.join(pairs) + '\n')

    full_geno  = os.path.join(args.output_dir, f"{args.prefix}_pooldata.geno")
    omega_geno = os.path.join(args.output_dir, f"{args.prefix}_omega_pooldata.geno")
    write_geno(full_geno,  range(n_keep))
    write_geno(omega_geno, thin_idx)

    poolsize_path = os.path.join(args.output_dir, f"{args.prefix}_poolsize.txt")
    with open(poolsize_path, 'w') as fh:
        fh.write(' '.join(str(POOL_SIZES[s]) for s in samples) + '\n')

    cov_path = os.path.join(args.output_dir, f"{args.prefix}_treatment.cov")
    cov_vals = [str(TREATMENT_COV[s]) for s in samples]
    with open(cov_path, 'w') as fh:
        fh.write(' '.join(cov_vals) + '\n')

    contrast_path = os.path.join(args.output_dir, f"{args.prefix}_contrasts.con")
    has_unspec = args.include_unspec

    def make_contrast(rules: dict) -> list[str]:
        """rules: sample -> value; unlisted samples get 0."""
        return [str(rules.get(s, 0)) for s in samples]

    c1 = make_contrast({
        "AVB": 1, "PSB": 1, "RMB": 1,
        "AVT": -1, "PST": -1, "RMT": -1,
    })

    contrasts = [c1]
    contrast_names = ["B vs T"]

    if has_unspec:
        c2 = make_contrast({
            "AVB": 1, "PSB": 1, "RMB": 1,
            "S.flavaMA": -1, "S.flavaAZ": -1,
        })
        c3 = make_contrast({
            "AVT": 1, "PST": 1, "RMT": 1,
            "S.flavaMA": -1, "S.flavaAZ": -1,
        })
        contrasts += [c2, c3]
        contrast_names += ["B vs unspec", "T vs unspec"]

    with open(contrast_path, 'w') as fh:
        for c in contrasts:
            fh.write(' '.join(c) + '\n')

    pos_path = os.path.join(args.output_dir, f"{args.prefix}_snp_positions.csv")
    with open(pos_path, 'w') as fh:
        fh.write("mrk,chrom,pos\n")
        for i, (chrom, pos) in enumerate(kept_pos):
            fh.write(f"{i+1},{chrom},{pos}\n")
            
if __name__ == '__main__':
    main()
