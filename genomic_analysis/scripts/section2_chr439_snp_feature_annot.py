#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(".")
GLM = ROOT / "glm_lrt_gw_final/glmV1full.csv"
GFF = Path("../Dissertation/Experimental_Evolution_Sflava/Sarah Lai/Fall 2023/sfla_v2.gff3")
FA  = Path("../synteny/sfla_v2.fa")
ANNOT_TSV = ROOT / "final_plots/wild/chr439_region_annot.tsv"
OUT_TSV   = ROOT / "final_plots/wild/chr439_peak_snp_features.tsv"
OUT_SUMMARY = ROOT / "final_plots/wild/chr439_peak_snp_feature_summary.tsv"

CHROM = "chr_ScDA7r2_439_HRSCAF_779"
PEAK_START = 2_800_000
PEAK_END   = 3_000_000
LD_BLOCK_START = 2_640_000
LD_BLOCK_END   = 3_610_000
VIEW_START = 2_500_000   # background range for null comparison
VIEW_END   = 3_750_000

AD_TSV = ROOT / "variance_analysis/merged_ad.tsv"


def parse_gff(chrom):
    gene_intervals = []
    feats_per_mrna = defaultdict(list)  # mrna_id -> [(type, s, e)]
    mrna_to_gene   = {}
    with open(GFF) as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip().split("\t")
            if len(f) < 9 or f[0] != chrom: continue
            t = f[2]
            s, e = int(f[3]), int(f[4])
            strand = f[6]
            attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            gid = attr.get("gene_id", "").strip()
            ID  = attr.get("ID", "").strip()
            par = attr.get("Parent", "").strip()
            if t == "gene":
                gene_intervals.append((s, e, strand, gid))
            elif t in ("mRNA", "transcript"):
                mrna_to_gene[ID] = gid
                feats_per_mrna[ID]  # initialise key
            elif t in ("CDS", "five_prime_UTR", "three_prime_UTR"):
                # Parent may be CSV-list of mRNAs
                for p in par.split(","):
                    if p:
                        feats_per_mrna[p].append((t, s, e, strand))
    return gene_intervals, feats_per_mrna, mrna_to_gene


def annotate_snp(pos, gene_intervals, feats_per_mrna, mrna_to_gene):
    overlapping_genes = [(gs, ge, gst, gid) for (gs, ge, gst, gid) in gene_intervals
                          if gs <= pos <= ge]
    if not overlapping_genes:
        return ("intergenic", "", "", "")

    best = "intron"  # if we are inside a gene but not in any annotated feature
    best_gid = overlapping_genes[0][3]
    best_mrna = ""
    best_strand = overlapping_genes[0][2]
    rank = {"intergenic": 0, "intron": 1, "three_prime_UTR": 2,
            "five_prime_UTR": 3, "CDS": 4}
    for (gs, ge, gst, gid) in overlapping_genes:
        for mrna_id, feats in feats_per_mrna.items():
            if mrna_to_gene.get(mrna_id) != gid:
                continue
            for (ftype, fs, fe, fstrand) in feats:
                if fs <= pos <= fe and rank[ftype] >= rank.get(best, 0):
                    best = ftype
                    best_gid = gid
                    best_mrna = mrna_id
                    best_strand = fstrand
    return (best, best_gid, best_mrna, best_strand)




def codon_at(pos, mrna_id, feats_per_mrna, chrom_seq):
    cds = [(s, e, strand) for (t, s, e, strand) in feats_per_mrna[mrna_id] if t == "CDS"]
    if not cds:
        return None
    strand = cds[0][2]
    cds = sorted(cds, key=lambda x: x[0])
    cum = 0
    cds_pos = None
    for s, e, _ in cds:
        if s <= pos <= e:
            cds_pos = cum + (pos - s)
            break
        cum += (e - s + 1)
    if cds_pos is None:
        return None
    if strand == "-":
        total = sum(e - s + 1 for s, e, _ in cds)
        cds_pos = total - 1 - cds_pos
    codon_start = (cds_pos // 3) * 3
    codon_pos = cds_pos % 3
    parts = [chrom_seq[s-1:e] for s, e, _ in cds]
    mrna = Seq("".join(str(p) for p in parts))
    if strand == "-":
        mrna = mrna.reverse_complement()
    if codon_start + 3 > len(mrna):
        return None
    return (str(mrna[codon_start:codon_start+3]), codon_pos, strand)


def get_alt_allele(snp_id_lookup, chrom, pos):
    if snp_id_lookup is None:
        return None
    key = f"{chrom}__{pos}"
    return snp_id_lookup.get(key)


def main():
    glm = pd.read_csv(GLM)
    glm = glm[(glm["chrom"] == CHROM) & (glm["pos"] >= VIEW_START) & (glm["pos"] <= VIEW_END)].copy()
    glm["in_peak"] = (glm["pos"] >= PEAK_START) & (glm["pos"] <= PEAK_END)
    gene_intervals, feats_per_mrna, mrna_to_gene = parse_gff(CHROM)
    feats = []
    for _, row in glm.iterrows():
        feats.append(annotate_snp(int(row["pos"]), gene_intervals, feats_per_mrna, mrna_to_gene))
    glm["feature"], glm["gene_id"], glm["mrna_id"], glm["strand"] = zip(*feats)

    annot = pd.read_csv(ANNOT_TSV, sep="\t")
    annot_lookup = annot.set_index("sfla_id")[["dmel_symbol", "category"]]
    glm = glm.join(annot_lookup, on="gene_id")
    glm["category"] = glm["category"].fillna("OTHER")
    glm["dmel_symbol"] = glm["dmel_symbol"].fillna("")

    glm["in_ld"] = (glm["pos"] >= LD_BLOCK_START) & (glm["pos"] <= LD_BLOCK_END)
    cds_in_peak = glm[(glm["in_ld"]) & (glm["feature"] == "CDS")].copy()
    chrom_seq = None
    if len(cds_in_peak):
        for rec in SeqIO.parse(FA, "fasta"):
            if rec.id == CHROM:
                chrom_seq = rec.seq
                break

    snp_alt = {}
    if AD_TSV.exists():
        with open(AD_TSV) as fh:
            header = fh.readline().rstrip().split("\t")
            for line in fh:
                f = line.rstrip().split("\t")
                if len(f) < 4: continue
                if f[0] != CHROM: continue
                p = int(f[1])
                if VIEW_START <= p <= VIEW_END:
                    snp_alt[p] = (f[2], f[3])  # ref, alt

    cds_records = []
    for _, row in cds_in_peak.iterrows():
        pos = int(row["pos"])
        codon_info = codon_at(pos, row["mrna_id"], feats_per_mrna, chrom_seq)
        if codon_info is None:
            cds_records.append({"effect": "unknown", "ref_aa": "", "alt_aa": "",
                                 "ref": "", "alt": ""})
            continue
        codon, cpos, strand = codon_info
        ref_alt = snp_alt.get(pos, ("", ""))
        ref, alt = ref_alt
        if not ref:
            cds_records.append({"effect": "unknown_alt", "ref_aa": str(Seq(codon).translate()),
                                 "alt_aa": "", "ref": "", "alt": ""})
            continue
        if strand == "-":
            ref_c = str(Seq(ref).reverse_complement())
            alt_c = str(Seq(alt).reverse_complement())
        else:
            ref_c, alt_c = ref, alt
        ref_codon = codon
        if ref_codon[cpos] != ref_c:
            cds_records.append({"effect": "ref_mismatch",
                                 "ref_aa": str(Seq(ref_codon).translate()),
                                 "alt_aa": "",
                                 "ref": ref, "alt": alt})
            continue
        alt_codon = ref_codon[:cpos] + alt_c + ref_codon[cpos+1:]
        ref_aa = str(Seq(ref_codon).translate())
        alt_aa = str(Seq(alt_codon).translate())
        if ref_aa == alt_aa:
            eff = "synonymous"
        elif alt_aa == "*":
            eff = "nonsense"
        elif ref_aa == "*":
            eff = "stop_lost"
        else:
            eff = "missense"
        cds_records.append({"effect": eff, "ref_aa": ref_aa, "alt_aa": alt_aa,
                             "ref": ref, "alt": alt})
    cds_eff = pd.DataFrame(cds_records, index=cds_in_peak.index)
    glm = glm.join(cds_eff, how="left")
    glm["effect"] = glm["effect"].fillna("")
    glm["ref_aa"] = glm["ref_aa"].fillna("")
    glm["alt_aa"] = glm["alt_aa"].fillna("")

    glm.to_csv(OUT_TSV, sep="\t", index=False)

    in_peak = glm[glm["in_peak"]]
    out_peak = glm[~glm["in_peak"]]

    def show(df, label):
        c = df["feature"].value_counts(dropna=False)
        n = len(df)

    top_q = 0.95
    cutoff = in_peak["LRT_chisq"].quantile(top_q)
    top = in_peak[in_peak["LRT_chisq"] >= cutoff]
    show(top, "    top-LRT")

    if (in_peak["feature"] == "CDS").any():
        cds = in_peak[in_peak["feature"] == "CDS"]
        top_cds = cds[cds["LRT_chisq"] >= cutoff]

    in_ld = glm[glm["in_ld"]]
    cds_ld = in_ld[in_ld["feature"] == "CDS"]
    if (cds_ld["effect"].isin(["missense", "nonsense", "stop_lost"])).any():
        nonsyn = cds_ld[cds_ld["effect"].isin(["missense", "nonsense", "stop_lost"])].copy()
        nonsyn = nonsyn.sort_values("LRT_chisq", ascending=False)
        cols = ["pos", "LRT_chisq", "PB_p_val", "feature", "gene_id",
                "dmel_symbol", "category", "effect", "ref_aa", "alt_aa", "ref", "alt"]
        ld_top_q = in_ld["LRT_chisq"].quantile(0.95)
        ld_top1q = in_ld["LRT_chisq"].quantile(0.99)
        n_top5 = (nonsyn["LRT_chisq"] >= ld_top_q).sum()
        n_top1 = (nonsyn["LRT_chisq"] >= ld_top1q).sum()
        nonsyn[cols].to_csv(ROOT / "final_plots/wild/chr439_ldblock_nonsyn_snps.tsv",
                              sep="\t", index=False)


    rows = []
    for label, df in [("peak", in_peak), ("background", out_peak), ("peak_top5pct", top)]:
        n = len(df)
        for feat in ("CDS", "five_prime_UTR", "three_prime_UTR", "intron", "intergenic"):
            cnt = (df["feature"] == feat).sum()
            rows.append({"set": label, "feature": feat, "n": cnt, "frac": cnt/n if n else 0})
    pd.DataFrame(rows).to_csv(OUT_SUMMARY, sep="\t", index=False)


if __name__ == "__main__":
    main()
