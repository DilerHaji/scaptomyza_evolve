#!/usr/bin/env python3

from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import sys

from Bio import SeqIO
from Bio.Seq import Seq

GFF = Path("../Dissertation/Experimental_Evolution_Sflava/"
           "Sarah Lai/Fall 2023/sfla_v2.gff3")
FA = Path("../synteny/sfla_v2.fa")
OUT_FAA = Path("./"
               "final_plots/wild/sfla_v2_proteome.faa")
OUT_INDEX = OUT_FAA.with_suffix(".tsv")  # gene_id -> chrom (used for join later)


def parse_gff_cds():
    out: dict[str, dict] = {}
    with open(GFF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            chrom, cs, ce, strand, attr_str = f[0], int(f[3]), int(f[4]), f[6], f[8]
            attr = dict(kv.split("=", 1) for kv in attr_str.split(";") if "=" in kv)
            gid = attr.get("gene_id", "").strip()
            if not gid:
                gid = attr.get("ID", "").rsplit(".", 1)[0]
            if not gid:
                continue
            rec = out.setdefault(gid, {"chrom": chrom, "strand": strand, "cds": []})
            if rec["chrom"] != chrom:
                continue
            rec["cds"].append((cs, ce))
    return out


def main():
    genes = parse_gff_cds()
    chrom_counts: dict[str, int] = defaultdict(int)
    for g in genes.values():
        chrom_counts[g["chrom"]] += 1

    chrom_seq: dict[str, Seq] = {}
    needed = set(g["chrom"] for g in genes.values())
    for rec in SeqIO.parse(FA, "fasta"):
        if rec.id in needed:
            chrom_seq[rec.id] = rec.seq
    missing = needed - set(chrom_seq)


    OUT_FAA.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    n_skip_short = 0
    n_skip_intern_stop = 0
    with open(OUT_FAA, "w") as fh, open(OUT_INDEX, "w") as ix:
        ix.write("gene_id\tchrom\tstart\tend\tstrand\taa_len\n")
        for gid, info in genes.items():
            if info["chrom"] not in chrom_seq:
                continue
            cds = sorted(info["cds"])
            seq = chrom_seq[info["chrom"]]
            mrna = Seq("".join(str(seq[s - 1:e]) for s, e in cds))
            if info["strand"] == "-":
                mrna = mrna.reverse_complement()
            prot = str(mrna.translate(to_stop=False)).rstrip("*")
            if len(prot) < 30:
                n_skip_short += 1
                continue
            internal_stops = prot.count("*")
            if internal_stops > 0 and internal_stops / max(len(prot), 1) > 0.05:
                n_skip_intern_stop += 1
                continue
            prot = prot.replace("*", "X")
            fh.write(f">{gid} {info['chrom']}:{cds[0][0]}-{cds[-1][1]} "
                     f"strand={info['strand']} aa={len(prot)}\n")
            for i in range(0, len(prot), 60):
                fh.write(prot[i:i + 60] + "\n")
            ix.write(f"{gid}\t{info['chrom']}\t{cds[0][0]}\t{cds[-1][1]}\t"
                     f"{info['strand']}\t{len(prot)}\n")
            n_written += 1

if __name__ == "__main__":
    main()
