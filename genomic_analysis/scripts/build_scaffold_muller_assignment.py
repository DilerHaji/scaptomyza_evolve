#!/usr/bin/env python3

from __future__ import annotations
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(".")
GENE_INDEX = ROOT / "final_plots/wild/sfla_v2_proteome.tsv"
BLASTP = ROOT / "final_plots/wild/sfla_v2_dmel_blastp.tsv"
DMEL_CHROM = ROOT / "final_plots/wild/dmel_gene_chrom.tsv"
OUT = ROOT / "final_plots/wild/sfla_v2_scaffold_muller.tsv"

MIN_HITS_CONFIDENT = 10
PCT_CONFIDENT = 0.70
MIN_HITS_AMBIG = 5
PCT_AMBIG = 0.50

GN_RE = re.compile(r"\bGN=(\S+)")
ID_RE = re.compile(r"^(?:sp|tr)\|([^|]+)\|(\S+)")  # accession, entry_name


def parse_blastp():
    with open(BLASTP) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            qid, sid, stitle = f[0], f[1], f[2]
            cands: list[str] = []
            m = GN_RE.search(stitle)
            if m:
                cands.append(m.group(1))
            m = ID_RE.match(sid)
            if m:
                acc, entry_name = m.group(1), m.group(2)
                if "_DROME" in entry_name:
                    cands.append(entry_name.split("_DROME")[0])
                cands.append(acc)
            yield qid, cands


def main():
    gene_idx = pd.read_csv(GENE_INDEX, sep="\t")
    chrom_lookup = pd.read_csv(DMEL_CHROM, sep="\t")
    priority = {"canonical": 0, "canonical_bare": 1, "synonym": 2}
    sym_to_chrom: dict[str, tuple[str, str, int]] = {}
    for _, r in chrom_lookup.iterrows():
        s = r["symbol"]
        p = priority.get(r["source"], 3)
        if s not in sym_to_chrom or p < sym_to_chrom[s][2]:
            sym_to_chrom[s] = (r["chrom"], r["muller"], p)

    sfla_to_muller: dict[str, str] = {}
    sfla_to_chrom: dict[str, str] = {}
    n_query = 0
    n_hit_with_chrom = 0
    n_hit_no_chrom = 0
    for sfla_gid, cands in parse_blastp():
        n_query += 1
        chrom = muller = None
        for c in cands:
            if c in sym_to_chrom:
                chrom, muller, _ = sym_to_chrom[c]
                break
        if muller is None:
            n_hit_no_chrom += 1
            continue
        sfla_to_muller[sfla_gid] = muller
        sfla_to_chrom[sfla_gid] = chrom
        n_hit_with_chrom += 1

    scaff_to_genes: dict[str, list[str]] = defaultdict(list)
    for _, r in gene_idx.iterrows():
        scaff_to_genes[r["chrom"]].append(r["gene_id"])

    rows = []
    for scaff, genes in scaff_to_genes.items():
        muller_hits = [sfla_to_muller[g] for g in genes if g in sfla_to_muller]
        chrom_hits = [sfla_to_chrom[g] for g in genes if g in sfla_to_chrom]
        n_genes_total = len(genes)
        n_blastp_resolved = len(muller_hits)
        c = Counter(muller_hits)
        if c:
            top_m, top_n = c.most_common(1)[0]
            top_pct = top_n / n_blastp_resolved
            second_m, second_pct = (None, 0.0)
            if len(c) >= 2:
                second_m, second_n = c.most_common(2)[1]
                second_pct = second_n / n_blastp_resolved
            if (n_blastp_resolved >= MIN_HITS_CONFIDENT
                    and top_pct >= PCT_CONFIDENT):
                cls = "confident"
            elif (n_blastp_resolved >= MIN_HITS_AMBIG
                  and top_pct >= PCT_AMBIG):
                cls = "ambiguous"
            else:
                cls = "unclassified"
        else:
            top_m = top_n = second_m = None
            top_pct = second_pct = 0.0
            cls = "unclassified"

        rows.append({
            "scaffold": scaff,
            "n_genes_total": n_genes_total,
            "n_blastp_resolved": n_blastp_resolved,
            "top_muller": top_m,
            "top_muller_n": top_n,
            "top_muller_pct": round(top_pct, 3) if top_pct else 0.0,
            "second_muller": second_m,
            "second_muller_pct": round(second_pct, 3) if second_pct else 0.0,
            "classification": cls,
            "muller_assigned": top_m if cls in ("confident", "ambiguous") else None,
        })

    df = pd.DataFrame(rows).sort_values("n_genes_total", ascending=False)
    df.to_csv(OUT, sep="\t", index=False)

    for s in ["chr_ScDA7r2_439_HRSCAF_779",   # chr_439, expected Muller A (X)
              "chr_ScDA7r2_597_HRSCAF_953",   # chr_597
              "chr_ScDA7r2_126_HRSCAF_325",   # chr_126
              "chr_ScDA7r2_110_HRSCAF_295"]:  # chr_110
        sub = df[df["scaffold"] == s]
        if not sub.empty:
            r = sub.iloc[0]

if __name__ == "__main__":
    main()
