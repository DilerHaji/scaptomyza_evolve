#!/usr/bin/env python3

from __future__ import annotations
from pathlib import Path
import re
import pandas as pd

ROOT = Path(".")
BLAST = ROOT / "final_plots/wild/chr439_region_dmel_blastp.tsv"
FAA   = ROOT / "final_plots/wild/chr439_region_proteins.faa"
OUT   = ROOT / "final_plots/wild/chr439_region_annot.tsv"


def parse_dmel_title(title: str):
    gn = re.search(r"\bGN=(\S+)", title)
    sym = gn.group(1) if gn else ""
    desc_match = re.match(r"\S+\s+(.+?)\s+OS=", title)
    desc = desc_match.group(1) if desc_match else ""
    return sym, desc


def classify(sym: str, desc: str) -> str:
    s = sym.lower()
    d = desc.lower()
    text = f"{s} {d}"

    # CHEMORECEPTION
    if re.match(r"^(or|gr|ir)\d", s) and "receptor" in d:
        return "CHEMORECEPTION"
    if re.match(r"^(or|gr|ir)\d+[a-z]?$", s):
        return "CHEMORECEPTION"
    if re.match(r"^obp\d", s) or "odorant-binding" in d or "odorant binding" in d:
        return "CHEMORECEPTION"
    if "gustatory receptor" in d or "olfactory receptor" in d or "ionotropic receptor" in d:
        return "CHEMORECEPTION"
    if re.match(r"^cheb", s):
        return "CHEMORECEPTION"

    # DETOX
    if re.match(r"^cyp\d", s) or "cytochrome p450" in d:
        return "DETOX"
    if re.match(r"^gst[a-z]?\d?", s) or "glutathione s-transferase" in d or "glutathione transferase" in d:
        return "DETOX"
    if re.match(r"^ugt", s) or "udp-glucuronosyl" in d or "udp-glycosyl" in d:
        return "DETOX"
    if re.match(r"^cce", s) or "carboxylesterase" in d or "carboxyl/cholinesterase" in d:
        return "DETOX"
    if re.match(r"^sod\d?$", s) or "superoxide dismutase" in d:
        return "DETOX"

    # CHANNEL / TRANSPORT
    if re.match(r"^trp", s) or "transient receptor potential" in d:
        return "CHANNEL"
    if re.match(r"^(para|nav|napilot)", s) and "channel" in d:
        return "CHANNEL"
    if "voltage-gated" in d and "channel" in d:
        return "CHANNEL"
    if re.match(r"^(slc|sglt|glut)\d?", s):
        return "CHANNEL"
    if re.match(r"^abc[a-g]?\d?", s):
        return "CHANNEL"
    if re.match(r"^ppk", s) or "pickpocket" in d:
        return "CHANNEL"
    if re.match(r"^aqp", s) or "aquaporin" in d:
        return "CHANNEL"
    if "ion channel" in d or "potassium channel" in d or "calcium channel" in d or "sodium channel" in d or "chloride channel" in d:
        return "CHANNEL"
    if re.match(r"^(eag|sk|shal|shab|shaw|shaker)$", s):
        return "CHANNEL"

    # IMMUNE
    if re.match(r"^(drs|dpt|atta|cec|def|mtk|listericin)", s) or "antimicrobial" in d:
        return "IMMUNE"
    if re.match(r"^toll", s) or re.match(r"^imd", s):
        return "IMMUNE"
    if "innate immune" in d or "antibacterial" in d or "antifungal peptide" in d:
        return "IMMUNE"

    # SIGNALING
    if re.match(r"^rh\d", s) or "rhodopsin" in d:
        return "SIGNALING"
    if re.match(r"^(dop|oct|tyr)r?\d?", s) and "receptor" in d:
        return "SIGNALING"
    if re.match(r"^nachr", s) or "acetylcholine receptor" in d:
        return "SIGNALING"
    if re.match(r"^mglur", s) or "metabotropic glutamate" in d:
        return "SIGNALING"
    if re.search(r"\bg[- ]protein[- ]coupled receptor\b", d) or "gpcr" in d:
        return "SIGNALING"
    if re.search(r"\bkinase\b", d) and "domain" not in d:
        return "SIGNALING"
    if re.search(r"\bphosphatase\b", d):
        return "SIGNALING"
    if re.match(r"^(ras|rho|rab|cdc42|arf)", s):
        return "SIGNALING"

    # REGULATION (transcription / translation / chromatin)
    if re.search(r"transcription factor", d) or re.search(r"\btf\b", d):
        return "REG"
    if "zinc finger" in d or "homeobox" in d or "helix-loop-helix" in d or "bhlh" in d:
        return "REG"
    if "chromatin" in d or "histone" in d or "nucleosome" in d:
        return "REG"
    if "translation initiation" in d or "ribosomal protein" in d or "ribosome" in d:
        return "REG"
    if "rna polymerase" in d or "splicing" in d or "spliceosome" in d:
        return "REG"

    # STRUCT (cytoskeleton, ECM, adhesion)
    if "actin" in d or "myosin" in d or "tubulin" in d or "kinesin" in d or "dynein" in d:
        return "STRUCT"
    if "cadherin" in d or "integrin" in d or "collagen" in d or "laminin" in d:
        return "STRUCT"
    if "extracellular matrix" in d:
        return "STRUCT"

    # METAB (broad metabolic enzyme)
    if re.search(r"\b(synthase|synthetase|reductase|dehydrogenase|hydrolase|hydratase|transferase|lyase|ligase|oxidase|peroxidase|isomerase|dehydratase|carboxylase|decarboxylase|aminotransferase|deaminase|aldolase|kinase|esterase|phosphorylase|mutase|epimerase|racemase)\b", d):
        return "METAB"
    if "metabolic" in d or "metabolism" in d:
        return "METAB"

    return "OTHER"


def main():
    rows = []
    with open(BLAST) as fh:
        for line in fh:
            f = line.rstrip().split("\t")
            if len(f) < 7: continue
            qid, sid, stitle, pid, length, evalue, bitscore = f[0], f[1], f[2], f[3], f[4], f[5], f[6]
            sym, desc = parse_dmel_title(stitle)
            cat = classify(sym, desc)
            rows.append({
                "sfla_id": qid,
                "dmel_acc": sid.split("|")[1] if "|" in sid else sid,
                "dmel_symbol": sym,
                "dmel_desc": desc,
                "category": cat,
                "pident": float(pid),
                "evalue": float(evalue),
                "bitscore": float(bitscore),
            })
    annot = pd.DataFrame(rows)

    all_qids = []
    with open(FAA) as fh:
        for line in fh:
            if line.startswith(">"):
                all_qids.append(line[1:].split()[0])
    missing = [q for q in all_qids if q not in set(annot["sfla_id"])]
    for q in missing:
        annot.loc[len(annot)] = {
            "sfla_id": q, "dmel_acc": "", "dmel_symbol": "", "dmel_desc": "",
            "category": "OTHER", "pident": 0.0, "evalue": float("inf"), "bitscore": 0.0,
        }

    annot = annot.sort_values("sfla_id").reset_index(drop=True)
    annot.to_csv(OUT, sep="\t", index=False)
    show = annot[annot["category"].isin(["CHEMORECEPTION", "DETOX", "CHANNEL", "IMMUNE"])]

if __name__ == "__main__":
    main()
