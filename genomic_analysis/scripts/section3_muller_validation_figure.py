#!/usr/bin/env python3
from __future__ import annotations
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

ROOT = Path(".")
WILD = ROOT / "final_plots/wild"

OUT_PREFIX = WILD / "section3_muller_validation"

MULLER_ORDER = ["A", "B", "C", "D", "E", "F"]
MULLER_COLORS = {
    "A": "#e41a1c",   # red
    "B": "#377eb8",   # blue
    "C": "#4daf4a",   # green
    "D": "#984ea3",   # purple
    "E": "#ff7f00",   # orange
    "F": "#a65628",   # brown
    "U": "#bdbdbd",   # unresolved
}

prot = pd.read_csv(WILD / "sfla_v2_proteome.tsv", sep="\t")
blast = pd.read_csv(
    WILD / "sfla_v2_dmel_blastp.tsv", sep="\t", header=None,
    names=["sfla_id", "dmel_uniprot", "dmel_desc",
           "pident", "length", "evalue", "bitscore"])
dmel = pd.read_csv(WILD / "dmel_gene_chrom.tsv", sep="\t")
scaffold_summary = pd.read_csv(WILD / "sfla_v2_scaffold_muller.tsv", sep="\t")


def parse_gn(desc: str) -> str | None:
    if not isinstance(desc, str):
        return None
    m = re.search(r"GN=(?:Dmel[\\\/])?(\S+)", desc)
    return m.group(1) if m else None


blast["dmel_symbol"] = blast["dmel_desc"].apply(parse_gn)
blast = blast.drop_duplicates("sfla_id", keep="first")

priority = {"canonical": 0, "canonical_bare": 1, "synonym": 2}
dmel_clean = dmel.dropna(subset=["symbol", "muller"]).copy()
dmel_clean["prio"] = dmel_clean["source"].map(priority).fillna(3).astype(int)
dmel_clean = dmel_clean.sort_values(["symbol", "prio"]).drop_duplicates(
    "symbol", keep="first")
chrom_lookup = dmel_clean[["symbol", "muller", "chrom"]]


genes = (prot.merge(blast[["sfla_id", "dmel_symbol"]],
                    left_on="gene_id", right_on="sfla_id", how="left")
              .merge(chrom_lookup, left_on="dmel_symbol", right_on="symbol",
                     how="left"))

genes = genes.rename(columns={"chrom_x": "scaffold", "chrom_y": "dmel_chrom"})
if "chrom" in genes.columns:
    genes = genes.rename(columns={"chrom": "scaffold"})

counts = (genes.dropna(subset=["muller"])
                .groupby(["scaffold", "muller"]).size().unstack(fill_value=0))
counts = counts.reindex(columns=MULLER_ORDER, fill_value=0)
totals = counts.sum(axis=1)

main_scaffs = totals[totals >= 30].sort_values(ascending=False).index.tolist()
panelA_data = counts.loc[main_scaffs]
panelA_data_pct = panelA_data.div(panelA_data.sum(axis=1), axis=0)
panelA_totals = panelA_data.sum(axis=1)


def short_name(scaff: str) -> str:
    m = re.match(r"chr_ScDA7r2_(\d+)(?:_HRSCAF_\d+)?(_unlocalized\.\d+)?", scaff)
    if m:
        suff = m.group(2) or ""
        return f"chr_{m.group(1)}{suff}"
    return scaff


def spatial_data(scaff: str) -> pd.DataFrame:
    sub = genes[genes["scaffold"] == scaff].copy()
    sub["mid_mb"] = (sub["start"] + sub["end"]) / 2 / 1e6
    return sub


CHR597 = scaffold_summary.iloc[0]["scaffold"]  # 4679 genes, fusion
CHR439 = "chr_ScDA7r2_439_HRSCAF_779"           # 1836 genes, Muller A

assert CHR597.startswith("chr_ScDA7r2_597"), CHR597
chr597_data = spatial_data(CHR597)
chr439_data = spatial_data(CHR439)

def _detect_break(df, n_bins=60):
    df = df.dropna(subset=["muller"]).sort_values("mid_mb")
    edges = np.linspace(df["mid_mb"].min(), df["mid_mb"].max(), n_bins + 1)
    cents = (edges[:-1] + edges[1:]) / 2
    fcs = []
    fds = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sub = df[(df["mid_mb"] >= lo) & (df["mid_mb"] < hi)]
        if len(sub) == 0:
            fcs.append(np.nan); fds.append(np.nan)
        else:
            fcs.append((sub["muller"] == "C").mean())
            fds.append((sub["muller"] == "D").mean())
    fcs = np.array(fcs); fds = np.array(fds)
    diff = fcs - fds
    sc = np.where(np.diff(np.sign(diff)) != 0)[0]
    return cents, fcs, fds, (cents[sc[len(sc)//2]] if len(sc) else None)

CHR597_BINS_X, CHR597_FRAC_C, CHR597_FRAC_D, CHR597_BREAK_MB = _detect_break(
    chr597_data, n_bins=60)
CHR597_SIZE_MB = chr597_data["end"].max() / 1e6

def rolling_fraction(df: pd.DataFrame, mullers: list[str], n_bins: int = 50):
    df = df.dropna(subset=["muller"]).sort_values("mid_mb")
    if len(df) == 0:
        return None
    edges = np.linspace(df["mid_mb"].min(), df["mid_mb"].max(), n_bins + 1)
    cents = (edges[:-1] + edges[1:]) / 2
    fracs = {}
    for m in mullers:
        f = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sub = df[(df["mid_mb"] >= lo) & (df["mid_mb"] < hi)]
            if len(sub) == 0:
                f.append(np.nan)
            else:
                f.append((sub["muller"] == m).mean())
        fracs[m] = np.array(f)
    return cents, fracs


plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 9,
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
})

fig = plt.figure(figsize=(7.8, 7.4))
gs_top = gridspec.GridSpec(2, 1, height_ratios=[1.4, 2.5],
                           hspace=0.42,
                           left=0.13, right=0.97, top=0.95, bottom=0.06)
gs_A = gs_top[0]
gs_BC = gridspec.GridSpecFromSubplotSpec(
    2, 2, subplot_spec=gs_top[1], height_ratios=[0.55, 1.0],
    width_ratios=[1.0, 1.15],
    hspace=0.05, wspace=0.32)

axA = fig.add_subplot(gs_A)
y = np.arange(len(panelA_data_pct))
left = np.zeros(len(panelA_data_pct))
for m in MULLER_ORDER:
    vals = panelA_data_pct[m].values
    axA.barh(y, vals, left=left, color=MULLER_COLORS[m],
             edgecolor="white", linewidth=0.4, label=f"Muller {m}")
    left += vals
labels = [short_name(s) for s in panelA_data_pct.index]
axA.set_yticks(y)
axA.set_yticklabels(labels, fontsize=8)
axA.invert_yaxis()
axA.set_xlim(0, 1)
axA.set_xlabel("Fraction of genes assigned to D. melanogaster Muller (top blastp hit)")
axA.set_title("A. Per-scaffold Muller composition", fontsize=10, loc="left",
              fontweight="bold", pad=6)


summary_idx = scaffold_summary.set_index("scaffold")
for i, scaff in enumerate(panelA_data_pct.index):
    n = int(summary_idx.loc[scaff, "n_blastp_resolved"])
    classif = summary_idx.loc[scaff, "classification"]
    assigned = summary_idx.loc[scaff, "muller_assigned"]
    if classif == "confident":
        tag = f"= {assigned}"
    elif classif == "ambiguous":
        tag = f"= {assigned}*"
    else:
        tag = "fusion"
    axA.text(1.02, i, f"n={n:,}  {tag}",
             ha="left", va="center", fontsize=7.5)

axA.legend(ncol=6, loc="upper center",
           bbox_to_anchor=(0.5, -0.18), frameon=False,
           fontsize=8, handlelength=1, handletextpad=0.5,
           columnspacing=1.0)


axB_top = fig.add_subplot(gs_BC[0, 0])
axB_bot = fig.add_subplot(gs_BC[1, 0], sharex=axB_top)

dB = chr597_data.dropna(subset=["muller"]).copy()
n_resolved = int(summary_idx.loc[CHR597, "n_blastp_resolved"])

rng = np.random.default_rng(42)
dB["yj"] = rng.uniform(-0.4, 0.4, size=len(dB))
for m in MULLER_ORDER:
    sub = dB[dB["muller"] == m]
    if len(sub) == 0:
        continue
    axB_top.scatter(sub["mid_mb"], sub["yj"], s=2.0,
                    c=MULLER_COLORS[m], alpha=0.45, linewidths=0,
                    label=f"{m} (n={len(sub):,})", rasterized=True)
axB_top.set_yticks([])
axB_top.set_ylim(-0.7, 0.7)
axB_top.set_title(f"B. chr_597 -- putative C+D fusion ({n_resolved:,} genes)",
                  fontsize=10, loc="left", fontweight="bold", pad=4)
axB_top.set_xlabel("")
axB_top.tick_params(labelbottom=False)


cents, fracs = rolling_fraction(dB, ["C", "D"], n_bins=60)
axB_bot.plot(cents, fracs["C"], color=MULLER_COLORS["C"], lw=1.5, label="Muller C")
axB_bot.plot(cents, fracs["D"], color=MULLER_COLORS["D"], lw=1.5, label="Muller D")
axB_bot.axhline(0.5, ls=":", lw=0.6, color="grey")


diff = fracs["C"] - fracs["D"]
sign_change = np.where(np.diff(np.sign(diff)) != 0)[0]
if len(sign_change) > 0:
    fb_idx = sign_change[len(sign_change) // 2]
    fb_mb = cents[fb_idx]
    for ax in (axB_top, axB_bot):
        ax.axvline(fb_mb, color="black", ls="--", lw=0.8, alpha=0.6)
    axB_top.text(fb_mb, 0.6, f"break ~{fb_mb:.0f} Mb",
                 ha="center", va="bottom", fontsize=7.5,
                 color="black")

axB_bot.set_ylim(0, 1)
axB_bot.set_ylabel("Fraction of genes\nin 60 bins")
axB_bot.set_xlabel("Position on chr_597 (Mb)")
axB_bot.legend(loc="center right", frameon=False, fontsize=7.5)


axC = fig.add_subplot(gs_BC[:, 1])

DMEL_SIZES_MB = {"X": 23.5, "2L": 23.5, "2R": 25.3,
                 "3L": 28.1, "3R": 32.1, "4": 1.35}
DMEL_TO_MULLER = {"X": "A", "2L": "B", "2R": "C",
                  "3L": "D", "3R": "E", "4": "F"}


SFLA_SIZES_MB = {
    s: prot[prot["chrom"] == s]["end"].max() / 1e6
    for s in [CHR439,
              "chr_ScDA7r2_126_HRSCAF_325",
              CHR597,
              "chr_ScDA7r2_110_HRSCAF_295",
              "chr_ScDA7r2_2_HRSCAF_23"]
}


ROWS = [
    ("Muller A", [("X",  DMEL_SIZES_MB["X"])],   CHR439,                       None),
    ("Muller B", [("2L", DMEL_SIZES_MB["2L"])], "chr_ScDA7r2_126_HRSCAF_325",   None),
    ("Muller D + C\n(2R + 3L fusion)",
                  [("3L", DMEL_SIZES_MB["3L"]),
                   ("2R", DMEL_SIZES_MB["2R"])], CHR597,            CHR597_BREAK_MB),
    ("Muller E", [("3R", DMEL_SIZES_MB["3R"])], "chr_ScDA7r2_110_HRSCAF_295",   None),
    ("Muller F", [("4",  DMEL_SIZES_MB["4"])],  "chr_ScDA7r2_2_HRSCAF_23",      None),
]

bar_h = 0.36
gap = 0.95   # vertical gap between rows
n_rows = len(ROWS)
y_max = n_rows * gap + 0.5

x_max = max(SFLA_SIZES_MB.values()) * 1.04  # leave a touch of right padding

# Title
axC.set_title("C. Karyotype mapping: Drosophila chromosomes paired with "
              "S. flava scaffolds",
              fontsize=10, loc="left", fontweight="bold", pad=6)

import matplotlib.patches as patches

for ri, (label, dmel_list, sfla, break_mb) in enumerate(ROWS):
    y_center = (n_rows - 1 - ri) * gap
    y_dmel_bottom = y_center + 0.35
    y_dmel_top    = y_dmel_bottom + bar_h
    y_sfla_top    = y_center - 0.05
    y_sfla_bottom = y_sfla_top - bar_h
    if break_mb is None:
        chrom, sz = dmel_list[0]
        muller = DMEL_TO_MULLER[chrom]
        sfla_size = SFLA_SIZES_MB[sfla]
        ortholog_lo, ortholog_hi = 0, sfla_size
        mid_orth = (ortholog_lo + ortholog_hi) / 2
        x_dmel_lo = mid_orth - sz / 2
        x_dmel_hi = mid_orth + sz / 2
        axC.add_patch(patches.Rectangle(
            (x_dmel_lo, y_dmel_bottom), sz, bar_h,
            facecolor=MULLER_COLORS[muller], edgecolor="black",
            linewidth=1.4, alpha=0.85))
        if sz / x_max > 0.04:
            axC.text(x_dmel_lo + sz / 2, y_dmel_bottom + bar_h / 2,
                     chrom, ha="center", va="center",
                     fontsize=8, color="white", fontweight="bold")
        else:
            axC.text(x_dmel_lo + sz / 2, y_dmel_top + 0.06,
                     chrom, ha="center", va="bottom",
                     fontsize=7, color="black", fontweight="bold")
        axC.text(x_dmel_hi + 0.6, y_dmel_bottom + bar_h / 2,
                 f"{sz:.1f} Mb", ha="left", va="center", fontsize=6.5,
                 color="dimgrey")
        axC.add_patch(patches.Rectangle(
            (0, y_sfla_bottom), sfla_size, bar_h,
            facecolor=MULLER_COLORS[muller], edgecolor="black",
            linewidth=0.6, alpha=0.85))
        if sfla_size / x_max > 0.06:
            axC.text(sfla_size / 2, y_sfla_bottom + bar_h / 2,
                     f"{short_name(sfla)}",
                     ha="center", va="center", fontsize=7.5,
                     color="white", fontweight="bold")
        else:
            axC.text(sfla_size + 2.5, y_sfla_bottom + bar_h / 2,
                     f"{short_name(sfla)}", ha="left", va="center",
                     fontsize=7, color="black", fontweight="bold")
        axC.text(sfla_size + 0.6, y_sfla_bottom - 0.18,
                 f"{sfla_size:.1f} Mb", ha="left", va="top",
                 fontsize=6.5, color="dimgrey") if sfla_size / x_max <= 0.06 \
            else axC.text(sfla_size + 0.6, y_sfla_bottom + bar_h / 2,
                          f"{sfla_size:.1f} Mb", ha="left", va="center",
                          fontsize=6.5, color="dimgrey")
        axC.add_patch(patches.Polygon(
            [(x_dmel_lo, y_dmel_bottom),
             (x_dmel_hi, y_dmel_bottom),
             (sfla_size, y_sfla_top),
             (0,         y_sfla_top)],
            closed=True,
            facecolor=MULLER_COLORS[muller], alpha=0.16,
            edgecolor="none"))
    else:
        sfla_size = SFLA_SIZES_MB[sfla]
        block_widths = [(0, break_mb, "3L", "D"),
                        (break_mb, sfla_size, "2R", "C")]
        for blo, bhi, chrom, muller in block_widths:
            sz = DMEL_SIZES_MB[chrom]
            mid = (blo + bhi) / 2
            x_dmel_lo = mid - sz / 2
            x_dmel_hi = mid + sz / 2
            axC.add_patch(patches.Rectangle(
                (x_dmel_lo, y_dmel_bottom), sz, bar_h,
                facecolor=MULLER_COLORS[muller], edgecolor="black",
                linewidth=1.4, alpha=0.85))
            axC.text(x_dmel_lo + sz / 2, y_dmel_bottom + bar_h / 2,
                     chrom, ha="center", va="center",
                     fontsize=8, color="white", fontweight="bold")
            axC.text(x_dmel_hi + 0.5, y_dmel_bottom + bar_h / 2,
                     f"{sz:.1f}", ha="left", va="center",
                     fontsize=6.0, color="dimgrey")
            axC.add_patch(patches.Polygon(
                [(x_dmel_lo, y_dmel_bottom),
                 (x_dmel_hi, y_dmel_bottom),
                 (bhi,       y_sfla_top),
                 (blo,       y_sfla_top)],
                closed=True,
                facecolor=MULLER_COLORS[muller], alpha=0.16,
                edgecolor="none"))
            axC.add_patch(patches.Rectangle(
                (blo, y_sfla_bottom), bhi - blo, bar_h,
                facecolor=MULLER_COLORS[muller], edgecolor="none",
                alpha=0.85))
        axC.add_patch(patches.Rectangle(
            (0, y_sfla_bottom), sfla_size, bar_h,
            facecolor="none", edgecolor="black", linewidth=0.6))
        axC.plot([break_mb, break_mb],
                 [y_sfla_bottom - 0.04, y_sfla_top + 0.04],
                 color="black", linewidth=0.9, linestyle="--",
                 alpha=0.75)
        axC.text(break_mb, y_sfla_bottom - 0.16,
                 f"break ~{break_mb:.0f} Mb",
                 ha="center", va="top", fontsize=6.5,
                 color="black")
        axC.text(sfla_size / 2, y_sfla_bottom + bar_h / 2,
                 f"{short_name(sfla)}",
                 ha="center", va="center", fontsize=7.5,
                 color="white", fontweight="bold")
        axC.text(sfla_size + 0.6, y_sfla_bottom + bar_h / 2,
                 f"{sfla_size:.1f} Mb", ha="left", va="center",
                 fontsize=6.5, color="dimgrey")

    axC.text(-x_max * 0.03, y_center,
             label, ha="right", va="center", fontsize=8,
             fontweight="bold")

axC.text(x_max * 0.5, y_max - 0.55,
         "top of each pair: D. melanogaster chromosomes   |   "
         "bottom: S. flava v2 scaffolds",
         ha="center", va="bottom", fontsize=8, fontstyle="italic",
         color="dimgrey")

axC.set_xlim(-x_max * 0.06, x_max + x_max * 0.05)
axC.set_ylim(-0.55, y_max - 0.15)
axC.set_xlabel("Length (Mb)")
axC.set_yticks([])
axC.spines["left"].set_visible(False)

for ext in ("png", "svg", "pdf"):
    out = OUT_PREFIX.with_suffix(f".{ext}")
    fig.savefig(out, dpi=300 if ext == "png" else None, bbox_inches="tight")
    print(f"wrote {out}")
plt.close(fig)
