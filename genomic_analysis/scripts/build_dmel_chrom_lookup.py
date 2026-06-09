#!/usr/bin/env python3

from __future__ import annotations
import gzip
from pathlib import Path

ROOT = Path(".")
SRC = ROOT / "data/Drosophila_melanogaster.gene_info.gz"
OUT = ROOT / "final_plots/wild/dmel_gene_chrom.tsv"

CHROM_TO_MULLER = {
    "X": "A",
    "2L": "B",
    "2R": "C",
    "3L": "D",
    "3R": "E",
    "4": "F",
    "Y": "Y",
    "mitochondrion_genome": "M",
}


def normalize_chrom(c: str) -> str:
    c = c.strip()
    if c in CHROM_TO_MULLER:
        return c

    for tok in c.replace("|", ";").replace(",", ";").split(";"):
        tok = tok.strip()
        if tok in CHROM_TO_MULLER:
            return tok
    return ""


def main():
    rows: list[tuple[str, str, str, str]] = []  # (symbol, chrom, muller, source)
    seen_canonical: set[str] = set()
    with gzip.open(SRC, "rt") as fh:
        header = fh.readline()
        assert header.startswith("#tax_id"), header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            symbol, synonyms, chrom = f[2], f[4], f[6]
            chrom_n = normalize_chrom(chrom)
            if not chrom_n:
                continue
            muller = CHROM_TO_MULLER[chrom_n]
            if symbol and symbol != "-":
                rows.append((symbol, chrom_n, muller, "canonical"))
                seen_canonical.add(symbol)
                if "\\" in symbol:
                    bare = symbol.split("\\", 1)[1]
                    if bare and bare not in seen_canonical:
                        rows.append((bare, chrom_n, muller, "canonical_bare"))
            for syn in synonyms.split("|") if synonyms != "-" else []:
                syn = syn.strip()
                if not syn or syn == symbol:
                    continue
                if "\\" in syn:
                    syn = syn.split("\\", 1)[1]
                if syn:
                    rows.append((syn, chrom_n, muller, "synonym"))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        fh.write("symbol\tchrom\tmuller\tsource\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    n_can = sum(1 for r in rows if r[3] in ("canonical", "canonical_bare"))
    n_syn = sum(1 for r in rows if r[3] == "synonym")
    chrom_counts: dict[str, int] = {}
    for r in rows:
        if r[3] == "canonical":
            chrom_counts[r[1]] = chrom_counts.get(r[1], 0) + 1

if __name__ == "__main__":
    main()
