#!/usr/bin/env python3

from __future__ import annotations
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import sys

GFF = Path("../Dissertation/Experimental_Evolution_Sflava/Sarah Lai/Fall 2023/sfla_v2.gff3")
FA  = Path("../synteny/sfla_v2.fa")
OUT_FAA = Path("./final_plots/wild/chr439_region_proteins.faa")

CHROM = "chr_ScDA7r2_439_HRSCAF_779"
START = 2_500_000
END   = 3_750_000


def parse_gff_cds(chrom, start, end):
    out = {}
    with open(GFF) as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip().split("\t")
            if len(f) < 9 or f[0] != chrom or f[2] != "CDS": continue
            cs, ce = int(f[3]), int(f[4])
            if cs < start or ce > end: continue
            attr = dict(kv.split("=",1) for kv in f[8].split(";") if "=" in kv)
            gid = attr.get("gene_id", "").strip() or attr.get("ID","").rsplit(".",1)[0]
            if gid not in out:
                out[gid] = {"strand": f[6], "cds": []}
            out[gid]["cds"].append((cs, ce))
    return out


def main():
    chrom_seq = None
    for rec in SeqIO.parse(FA, "fasta"):
        if rec.id == CHROM:
            chrom_seq = rec.seq
            break
    if chrom_seq is None:
        sys.exit(1)

    cds_per_gene = parse_gff_cds(CHROM, START, END)

    OUT_FAA.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_FAA, "w") as fh:
        for gid, info in cds_per_gene.items():
            cds = sorted(info["cds"])
            parts = [chrom_seq[s-1:e] for s, e in cds]
            mrna = Seq("".join(str(p) for p in parts))
            if info["strand"] == "-":
                mrna = mrna.reverse_complement()
            prot = str(mrna.translate(to_stop=False)).rstrip("*")
            fh.write(f">{gid} {CHROM}:{cds[0][0]}-{cds[-1][1]} "
                     f"strand={info['strand']} aa={len(prot)}\n")
            for i in range(0, len(prot), 60):
                fh.write(prot[i:i+60] + "\n")

if __name__ == "__main__":
    main()
