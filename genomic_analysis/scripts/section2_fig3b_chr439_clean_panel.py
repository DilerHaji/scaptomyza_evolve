#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle, FancyBboxPatch
import numpy as np
import pandas as pd

ROOT = Path(".")
GFF_PATH = Path("../Dissertation/Experimental_Evolution_Sflava/Sarah Lai/Fall 2023/sfla_v2.gff3")
ANNOT_PATH = ROOT / "final_plots/wild/chr439_region_annot.tsv"
PI_PATH = ROOT / "grenfst/diversity/trajectory_pi_200000_nonoverlap.csv"
HV_PATH = ROOT / "final_plots/wild/section2_hv_blocks_filtered.tsv"
LYNCH_TSV = ROOT / "final_plots/wild/section2_fig3e_lynch_chr439_wellbehaved.tsv"
OUT_BASE = ROOT / "final_plots/wild/section2_fig3b_chr439_clean_panel"

CHROM = "chr_ScDA7r2_439_HRSCAF_779"
VIEW_START = 2_500_000
VIEW_END   = 3_750_000
LD_BLOCK_START = 2_640_000
LD_BLOCK_END   = 3_610_000
PEAK_START = 2_800_000
PEAK_END   = 3_000_000

C_PEAK = "#7d2c2c"
C_LDBLOCK = "#888888"
C_PI = "#2c7a4d"        # green for π (nucleotide diversity)
C_TD = "#5a3a8a"        # purple for Tajima's D


HV_COLOR = {"B": "#a04848", "T": "#3d6cb5", "M": "#7e6d8a"}

CAT_COLOR = {
    "CHEMORECEPTION": ("#D55E00", "Chemoreception (Or/Gr/Ir/Obp)"),
    "DETOX":          ("#009E73", "Detoxification (Cyp/Gst/Ugt)"),
    "CHANNEL":        ("#0072B2", "Ion channel / transport"),
    "IMMUNE":         ("#CC79A7", "Immunity"),
    "SIGNALING":      ("#E69F00", "Signaling (GPCR/kinase)"),
    "REG":            ("#a8a8a8", "Transcription / translation"),
    "STRUCT":         ("#a8a8a8", "Cytoskeleton / ECM"),
    "METAB":          ("#a8a8a8", "Metabolism"),
    "OTHER":          ("#d8d8d8", "Other / unannotated"),
}
CAT_BOLD = {"CHEMORECEPTION", "DETOX", "CHANNEL", "IMMUNE", "SIGNALING"}


def load_genes():
    genes = []
    with open(GFF_PATH) as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip().split("\t")
            if len(f) < 9 or f[0] != CHROM or f[2] != "gene": continue
            s, e = int(f[3]), int(f[4])
            if e < VIEW_START or s > VIEW_END: continue
            attr = dict(kv.split("=",1) for kv in f[8].split(";") if "=" in kv)
            genes.append({"start": s, "end": e, "strand": f[6],
                          "gene_id": attr.get("gene_id","").strip()})
    return pd.DataFrame(genes)


def load_annot():
    return pd.read_csv(ANNOT_PATH, sep="\t")


def load_hv_blocks():
    df = pd.read_csv(HV_PATH, sep="\t")
    m = (df["chr"] == CHROM) & (df["end"] > VIEW_START) & (df["start"] < VIEW_END)
    return df[m].sort_values(["treatment", "start"])


def load_lynch_wellbehaved():
    if not LYNCH_TSV.exists():
        return pd.DataFrame(columns=["POS", "antagonistic"])
    df = pd.read_csv(LYNCH_TSV, sep="\t")
    df = df[(df["CHROM"] == CHROM) &
             (df["POS"] >= VIEW_START) & (df["POS"] <= VIEW_END)]
    return df


def load_diversity():
    df = pd.read_csv(PI_PATH)
    pi_cols = [f"F{i}G00.theta_pi"  for i in range(1, 5)]
    td_cols = [f"F{i}G00.tajimas_d" for i in range(1, 5)]
    df["pi_F"] = df[pi_cols].mean(axis=1)
    df["td_F"] = df[td_cols].mean(axis=1)
    df["mid"]  = (df["start"] + df["end"]) // 2
    gw_pi_med = df["pi_F"].median()
    gw_td_med = df["td_F"].median()
    sub = df[(df["chrom"] == CHROM) & (df["end"] > VIEW_START) & (df["start"] < VIEW_END)].copy()
    return sub, gw_pi_med, gw_td_med


def main():
    plt.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.6,
    })

    genes = load_genes()
    annot = load_annot()
    genes = genes.merge(annot[["sfla_id", "dmel_symbol", "dmel_desc", "category"]],
                         left_on="gene_id", right_on="sfla_id", how="left")
    genes["category"] = genes["category"].fillna("OTHER")
    div, gw_pi, gw_td = load_diversity()
    hv = load_hv_blocks()
    lynch_snps = load_lynch_wellbehaved()
    n_antag = int(lynch_snps["antagonistic"].sum()) if len(lynch_snps) else 0
    fig = plt.figure(figsize=(7.0, 4.8))
    gs = GridSpec(7, 1, figure=fig,
                   height_ratios=[0.20, 0.20, 0.30, 0.22, 0.55, 0.55, 1.55],
                   hspace=0.18)

    ax_p = fig.add_subplot(gs[0])
    ax_p.set_xlim(VIEW_START, VIEW_END)
    ax_p.set_ylim(0, 1)
    ax_p.axhspan(0.40, 0.70, xmin=(PEAK_START-VIEW_START)/(VIEW_END-VIEW_START),
                  xmax=(PEAK_END-VIEW_START)/(VIEW_END-VIEW_START),
                  color=C_PEAK, alpha=0.85)
    ax_p.text((PEAK_START+PEAK_END)/2, 0.55,
               f"GLMM peak ({(PEAK_END-PEAK_START)/1e3:.0f} kb)",
               ha="center", va="center", fontsize=7,
               color="white", fontweight="bold")
    ax_p.text(VIEW_START - (VIEW_END-VIEW_START)*0.005, 0.55,
               "Significant hit", ha="right", va="center", fontsize=7.5,
               color=C_PEAK, fontweight="bold")
    ax_p.set_xticks([])
    ax_p.set_yticks([])
    plt.setp(ax_p.get_xticklabels(), visible=False)
    plt.setp(ax_p.get_yticklabels(), visible=False)
    ax_p.tick_params(axis="both", which="both", length=0,
                      labelbottom=False, labelleft=False, labelright=False, labeltop=False)
    for sp in ("top","right","left","bottom"):
        ax_p.spines[sp].set_visible(False)

    ax_l = fig.add_subplot(gs[1], sharex=ax_p)
    ax_l.set_xlim(VIEW_START, VIEW_END)
    ax_l.set_ylim(0, 1)
    ax_l.axhspan(0.40, 0.70, xmin=(LD_BLOCK_START-VIEW_START)/(VIEW_END-VIEW_START),
                  xmax=(LD_BLOCK_END-VIEW_START)/(VIEW_END-VIEW_START),
                  color=C_LDBLOCK, alpha=0.7)
    ax_l.text((LD_BLOCK_START+LD_BLOCK_END)/2, 0.55,
               f"haplotype block ({(LD_BLOCK_END-LD_BLOCK_START)/1e3:.0f} kb)",
               ha="center", va="center", fontsize=7,
               color="white", fontweight="bold")
    ax_l.text(VIEW_START - (VIEW_END-VIEW_START)*0.005, 0.55,
               "LD block", ha="right", va="center", fontsize=7.5,
               color=C_LDBLOCK, fontweight="bold")
    ax_l.set_xticks([])
    ax_l.set_yticks([])
    plt.setp(ax_l.get_xticklabels(), visible=False)
    plt.setp(ax_l.get_yticklabels(), visible=False)
    ax_l.tick_params(axis="both", which="both", length=0,
                      labelbottom=False, labelleft=False, labelright=False, labeltop=False)
    for sp in ("top","right","left","bottom"):
        ax_l.spines[sp].set_visible(False)

    ax_hv = fig.add_subplot(gs[2], sharex=ax_p)
    ax_hv.set_xlim(VIEW_START, VIEW_END)
    treatments_present = list(dict.fromkeys(hv["treatment"].tolist()))
    n_t = max(len(treatments_present), 1)
    ax_hv.set_ylim(0, n_t)
    for i, trt in enumerate(treatments_present):
        rows = hv[hv["treatment"] == trt]
        for _, b in rows.iterrows():
            xmin = max(b["start"], VIEW_START)
            xmax = min(b["end"], VIEW_END)
            y0 = (n_t - 1 - i) + 0.20
            y1 = (n_t - 1 - i) + 0.80
            ax_hv.add_patch(Rectangle((xmin, y0), xmax - xmin, y1 - y0,
                                        facecolor=HV_COLOR.get(trt, "#888"),
                                        edgecolor="none", alpha=0.8, zorder=2))

            if b["start"] < VIEW_START:
                ax_hv.annotate("", xy=(VIEW_START + (VIEW_END-VIEW_START)*0.005, (y0+y1)/2),
                                xytext=(VIEW_START + (VIEW_END-VIEW_START)*0.025, (y0+y1)/2),
                                arrowprops=dict(arrowstyle="->", color="white",
                                                  lw=1.0))
            if b["end"] > VIEW_END:
                ax_hv.annotate("", xy=(VIEW_END - (VIEW_END-VIEW_START)*0.005, (y0+y1)/2),
                                xytext=(VIEW_END - (VIEW_END-VIEW_START)*0.025, (y0+y1)/2),
                                arrowprops=dict(arrowstyle="->", color="white",
                                                  lw=1.0))

            label = f"{trt}-host  ({b['span_Mb']:.1f} Mb, {b['n_snps']:,} SNPs)"
            ax_hv.text((xmin + xmax) / 2, (y0 + y1) / 2, label,
                        ha="center", va="center", fontsize=6.5,
                        color="white", fontweight="bold", zorder=3)
    ax_hv.text(VIEW_START - (VIEW_END-VIEW_START)*0.005, n_t / 2,
                "HV blocks", ha="right", va="center", fontsize=7.5,
                color="#444444", fontweight="bold")
    ax_hv.set_xticks([])
    ax_hv.set_yticks([])
    plt.setp(ax_hv.get_xticklabels(), visible=False)
    plt.setp(ax_hv.get_yticklabels(), visible=False)
    ax_hv.tick_params(axis="both", which="both", length=0,
                       labelbottom=False, labelleft=False, labelright=False, labeltop=False)
    for sp in ("top","right","left","bottom"):
        ax_hv.spines[sp].set_visible(False)

    ax_ly = fig.add_subplot(gs[3], sharex=ax_p)
    ax_ly.set_xlim(VIEW_START, VIEW_END)
    ax_ly.set_ylim(0, 1)
    ax_ly.axvspan(PEAK_START, PEAK_END, color="#fff3eb", zorder=0)
    if len(lynch_snps):
        antag = lynch_snps[lynch_snps["antagonistic"]]
        nonantag = lynch_snps[~lynch_snps["antagonistic"]]
        if len(nonantag):
            ax_ly.vlines(nonantag["POS"].values, 0.20, 0.80,
                          color="#bbbbbb", linewidth=1.0, alpha=0.85, zorder=2)
        if len(antag):
            ax_ly.vlines(antag["POS"].values, 0.10, 0.90,
                          color="#1a5e1a", linewidth=1.4, alpha=1.0, zorder=3)
        ax_ly.text(VIEW_END - (VIEW_END-VIEW_START)*0.005, 0.5,
                    f"{len(antag)}/{len(lynch_snps)} antag.",
                    ha="right", va="center", fontsize=5.5,
                    color="#666666", style="italic")
    ax_ly.text(VIEW_START - (VIEW_END-VIEW_START)*0.005, 0.5,
                "Lynch SNPs", ha="right", va="center", fontsize=7.5,
                color="#1a5e1a", fontweight="bold")
    ax_ly.set_xticks([])
    ax_ly.set_yticks([])
    plt.setp(ax_ly.get_xticklabels(), visible=False)
    plt.setp(ax_ly.get_yticklabels(), visible=False)
    ax_ly.tick_params(axis="both", which="both", length=0,
                       labelbottom=False, labelleft=False, labelright=False, labeltop=False)
    for sp in ("top","right","left","bottom"):
        ax_ly.spines[sp].set_visible(False)

    ax_pi = fig.add_subplot(gs[4], sharex=ax_p)
    ax_pi.set_xlim(VIEW_START, VIEW_END)
    ax_pi.axvspan(PEAK_START, PEAK_END, color="#fff3eb", zorder=0)
    pi_max = max(div["pi_F"].max() if len(div) else gw_pi, gw_pi) * 1.05
    pi_min = min(div["pi_F"].min() if len(div) else gw_pi, gw_pi) * 0.95
    if len(div):
        widths = (div["end"] - div["start"]).values
        ax_pi.bar(div["start"].values, div["pi_F"].values - pi_min,
                   bottom=pi_min, width=widths,
                   align="edge", color=C_PI, alpha=0.65, edgecolor="none", zorder=2)
    ax_pi.axhline(gw_pi, color="#222222", linestyle="--", linewidth=0.6, alpha=0.7, zorder=3)
    ax_pi.set_ylim(pi_min, pi_max)
    ax_pi.text(VIEW_END - (VIEW_END-VIEW_START)*0.005, gw_pi + (pi_max-pi_min)*0.02,
                f"genome-wide π = {gw_pi:.2f}", ha="right", va="bottom",
                fontsize=5.5, color="#666666", style="italic")
    ax_pi.text(VIEW_START - (VIEW_END-VIEW_START)*0.005, (pi_min+pi_max)/2,
                "wild π", ha="right", va="center", fontsize=7.5,
                color=C_PI, fontweight="bold")
    ax_pi.set_xticks([])
    ax_pi.tick_params(axis="x", length=0, labelbottom=False)
    ax_pi.tick_params(axis="y", labelsize=6)
    ax_pi.set_yticks([round(pi_min, 2), round(gw_pi, 2), round(pi_max, 2)])
    for sp in ("top","right","bottom"):
        ax_pi.spines[sp].set_visible(False)
    ax_pi.spines["left"].set_linewidth(0.5)

    ax_td = fig.add_subplot(gs[5], sharex=ax_p)
    ax_td.set_xlim(VIEW_START, VIEW_END)
    ax_td.axvspan(PEAK_START, PEAK_END, color="#fff3eb", zorder=0)
    if len(div):
        widths = (div["end"] - div["start"]).values
        ax_td.bar(div["start"].values, div["td_F"].values, width=widths,
                   align="edge", color=C_TD, alpha=0.65, edgecolor="none", zorder=2)
    ax_td.axhline(gw_td, color="#222222", linestyle="--", linewidth=0.6, alpha=0.7, zorder=3)
    td_max = max(div["td_F"].max() if len(div) else gw_td, gw_td) * 1.10
    td_min = 0
    ax_td.set_ylim(td_min, td_max)
    ax_td.text(VIEW_END - (VIEW_END-VIEW_START)*0.005, gw_td + (td_max-td_min)*0.02,
                f"genome-wide D = {gw_td:+.2f}", ha="right", va="bottom",
                fontsize=5.5, color="#666666", style="italic")
    ax_td.text(VIEW_START - (VIEW_END-VIEW_START)*0.005, td_max/2,
                "Tajima D", ha="right", va="center", fontsize=7.5,
                color=C_TD, fontweight="bold")
    ax_td.set_xticks([])
    ax_td.tick_params(axis="x", length=0, labelbottom=False)
    ax_td.tick_params(axis="y", labelsize=6)
    ax_td.set_yticks([0, round(gw_td, 1), round(td_max, 1)])
    for sp in ("top","right","bottom"):
        ax_td.spines[sp].set_visible(False)
    ax_td.spines["left"].set_linewidth(0.5)

    ax_g = fig.add_subplot(gs[6], sharex=ax_p)
    ax_g.set_xlim(VIEW_START, VIEW_END)
    ax_g.axvspan(PEAK_START, PEAK_END, color="#fff3eb", zorder=0)

    for _, gn in genes.iterrows():
        y = 0.6 if gn["strand"] == "+" else -0.6
        cat = gn["category"]
        col = CAT_COLOR[cat][0]
        is_bold = cat in CAT_BOLD
        in_peak = (gn["end"] >= PEAK_START) and (gn["start"] <= PEAK_END)
        h = 0.55 if is_bold else (0.50 if in_peak else 0.40)

        if not is_bold and in_peak:
            col = "#8a8a8a"
        ax_g.add_patch(Rectangle((gn["start"], y - h/2),
                                   gn["end"] - gn["start"], h,
                                   facecolor=col, edgecolor="none",
                                   alpha=0.95, zorder=4 if in_peak else (3 if is_bold else 2)))

    annotated = genes[genes["category"].isin(CAT_BOLD)].copy()
    annotated["label"] = annotated["dmel_symbol"].fillna("")
    annotated = annotated[annotated["label"] != ""].sort_values("start")
    lab_y_above = 1.30
    lab_y_below = -1.30
    last_x_above = -1e9
    last_x_below = -1e9
    MIN_GAP = (VIEW_END - VIEW_START) * 0.08
    for _, gn in annotated.iterrows():
        gx = (gn["start"] + gn["end"]) / 2
        col = CAT_COLOR[gn["category"]][0]
        if gn["strand"] == "+":
            tier_y = lab_y_above if gx - last_x_above > MIN_GAP else lab_y_below
            if tier_y == lab_y_above: last_x_above = gx
            else: last_x_below = gx
        else:
            tier_y = lab_y_below if gx - last_x_below > MIN_GAP else lab_y_above
            if tier_y == lab_y_below: last_x_below = gx
            else: last_x_above = gx
        gene_y = 0.6 if gn["strand"] == "+" else -0.6
        anchor_y = gene_y + (0.30 if tier_y > 0 else -0.30)
        ax_g.plot([gx, gx], [anchor_y, tier_y], color=col,
                   linewidth=0.6, alpha=0.6, zorder=2)
        ax_g.text(gx, tier_y, gn["label"],
                   ha="center",
                   va="bottom" if tier_y > 0 else "top",
                   fontsize=6.5, fontweight="bold", fontstyle="italic", color=col)

    peak_genes = genes[
        (genes["end"] >= PEAK_START) & (genes["start"] <= PEAK_END)
        & genes["dmel_symbol"].notna() & (genes["dmel_symbol"] != "")
    ].sort_values("start")
    for _, gn in peak_genes.iterrows():
        gx = (gn["start"] + gn["end"]) / 2
        gene_y = 0.6 if gn["strand"] == "+" else -0.6
        anchor_y = gene_y + (0.30 if gene_y > 0 else -0.30)
        tier_y = 1.65
        ax_g.plot([gx, gx], [anchor_y, tier_y - 0.05],
                   color="#555555", linewidth=0.5, alpha=0.6, zorder=2)
        ax_g.text(gx, tier_y, gn["dmel_symbol"],
                   ha="left", va="bottom", rotation=45,
                   fontsize=6, fontweight="bold", fontstyle="italic",
                   color="#333333", zorder=5)
    ax_g.set_ylim(-1.7, 2.5)
    ax_g.set_yticks([0.6, -0.6])
    ax_g.set_yticklabels(["+", "−"], fontsize=7)
    ax_g.set_xticks(np.linspace(VIEW_START, VIEW_END, 6))
    ax_g.set_xticklabels([f"{x/1e6:.2f}" for x in np.linspace(VIEW_START, VIEW_END, 6)],
                          fontsize=7)
    ax_g.set_xlabel(f"Position on {CHROM} (Mb)", fontsize=8.5)
    for sp in ("top","right","left"):
        ax_g.spines[sp].set_visible(False)
    ax_g.tick_params(axis="y", length=0)

    from matplotlib.patches import Patch
    present = [c for c in CAT_COLOR if (genes["category"] == c).any()]
    legend_handles = [Patch(facecolor=CAT_COLOR[c][0], edgecolor="none",
                             label=CAT_COLOR[c][1]) for c in present]
    leg = ax_g.legend(handles=legend_handles, loc="upper left",
                       bbox_to_anchor=(1.005, 1.4),
                       ncol=1, frameon=False, fontsize=7,
                       handletextpad=0.4)
    label_colors = [CAT_COLOR[c][0] for c in present]
    for tx, c in zip(leg.get_texts(), label_colors):
        tx.set_color(c)

    fig.savefig(f"{OUT_BASE}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.svg", bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.pdf", bbox_inches="tight")
    import re
    svg = Path(f"{OUT_BASE}.svg")
    txt = svg.read_text()
    txt = re.sub(r"<clipPath[^>]*>.*?</clipPath>", "", txt, flags=re.DOTALL)
    txt = re.sub(r'\s*clip-path="url\([^)]+\)"', "", txt)
    svg.write_text(txt)


if __name__ == "__main__":
    main()
