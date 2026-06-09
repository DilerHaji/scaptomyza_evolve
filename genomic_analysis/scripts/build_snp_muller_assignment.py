#!/usr/bin/env python3

from __future__ import annotations
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(".")
GENE_INDEX = ROOT / "final_plots/wild/sfla_v2_proteome.tsv"
BLASTP = ROOT / "final_plots/wild/sfla_v2_dmel_blastp.tsv"
DMEL_CHROM = ROOT / "final_plots/wild/dmel_gene_chrom.tsv"
QSTAR_TSV = ROOT / "final_plots/wild/section2_wild_qstar_correlation_gw_top50_perrep.tsv"

OUT = ROOT / "final_plots/wild/sfla_v2_snp_muller.tsv"

MAX_DIST = 250_000 

GN_RE = re.compile(r"\bGN=(\S+)")
ID_RE = re.compile(r"^(?:sp|tr)\|([^|]+)\|(\S+)")


def main():
    cl = pd.read_csv(DMEL_CHROM, sep="\t")
    priority = {"canonical": 0, "canonical_bare": 1, "synonym": 2}
    sym_to_chrom: dict[str, tuple[str, str, int]] = {}
    for _, r in cl.iterrows():
        s = r["symbol"]
        p = priority.get(r["source"], 3)
        if s not in sym_to_chrom or p < sym_to_chrom[s][2]:
            sym_to_chrom[s] = (r["chrom"], r["muller"], p)

    sfla_to_muller: dict[str, str] = {}
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
                acc, en = m.group(1), m.group(2)
                if "_DROME" in en:
                    cands.append(en.split("_DROME")[0])
                cands.append(acc)
            for c in cands:
                if c in sym_to_chrom:
                    sfla_to_muller[qid] = sym_to_chrom[c][1]
                    break

    gi = pd.read_csv(GENE_INDEX, sep="\t")
    gi["mid"] = (gi["start"] + gi["end"]) // 2
    gi["muller"] = gi["gene_id"].map(sfla_to_muller)
    gi = gi.dropna(subset=["muller"]).copy()
    chrom_genes: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for ch, sub in gi.groupby("chrom"):
        sub = sub.sort_values("mid").reset_index(drop=True)
        chrom_genes[ch] = (sub["mid"].values.astype(np.int64),
                            sub["muller"].values)

    snps = pd.read_csv(QSTAR_TSV, sep="\t",
                       usecols=["chrom", "pos"])
    snps = snps.drop_duplicates().reset_index(drop=True)

    out_muller = []
    out_dist = []
    for ch, sub in snps.groupby("chrom"):
        if ch not in chrom_genes:
            out_muller.extend([None] * len(sub))
            out_dist.extend([np.nan] * len(sub))
            continue
        mids, mullers = chrom_genes[ch]
        positions = sub["pos"].values.astype(np.int64)
        idx_right = np.searchsorted(mids, positions, side="left")
        idx_right = np.clip(idx_right, 0, len(mids) - 1)
        idx_left = np.clip(idx_right - 1, 0, len(mids) - 1)
        d_right = np.abs(mids[idx_right] - positions)
        d_left = np.abs(mids[idx_left] - positions)
        use_right = d_right <= d_left
        idx_nearest = np.where(use_right, idx_right, idx_left)
        d_nearest = np.where(use_right, d_right, d_left)
        m_nearest = mullers[idx_nearest]
        m_final = np.where(d_nearest <= MAX_DIST, m_nearest, None)
        for i, _ in enumerate(sub.index):
            out_muller.append(m_final[i] if m_final[i] is not None else None)
            out_dist.append(int(d_nearest[i]))

    snps = snps.assign(muller=None, dist_to_gene=np.nan)
    cursor = 0
    for ch, sub in snps.groupby("chrom"):
        n = len(sub)
        snps.loc[sub.index, "muller"] = out_muller[cursor:cursor + n]
        snps.loc[sub.index, "dist_to_gene"] = out_dist[cursor:cursor + n]
        cursor += n

    snps.to_csv(OUT, sep="\t", index=False)

if __name__ == "__main__":
    main()
