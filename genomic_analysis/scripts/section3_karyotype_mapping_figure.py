#!/usr/bin/env python3
from __future__ import annotations
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

ROOT = Path(".")
WILD = ROOT / "final_plots/wild"
OUT_PREFIX = WILD / "section3_muller_karyotype"

MULLER_COLORS = {
    "A": "#e41a1c",  # red
    "B": "#377eb8",  # blue
    "C": "#4daf4a",  # green
    "D": "#984ea3",  # purple
    "E": "#ff7f00",  # orange
    "F": "#a65628",  # brown
}

DMEL_SIZES_MB = {"X": 23.5, "2L": 23.5, "2R": 25.3,
                 "3L": 28.1, "3R": 32.1, "4": 1.35}
DMEL_TO_MULLER = {"X": "A", "2L": "B", "2R": "C",
                  "3L": "D", "3R": "E", "4": "F"}

CHR597 = "chr_ScDA7r2_597_HRSCAF_953"
CHR439 = "chr_ScDA7r2_439_HRSCAF_779"


def short_name(scaff: str) -> str:
    m = re.match(r"chr_ScDA7r2_(\d+)(?:_HRSCAF_\d+)?", scaff)
    return f"chr_{m.group(1)}" if m else scaff


def parse_gn(desc: str) -> str | None:
    if not isinstance(desc, str):
        return None
    m = re.search(r"GN=(?:Dmel[\\\/])?(\S+)", desc)
    return m.group(1) if m else None


def load_chr597_break() -> tuple[float, float]:
    prot = pd.read_csv(WILD / "sfla_v2_proteome.tsv", sep="\t")
    blast = pd.read_csv(
        WILD / "sfla_v2_dmel_blastp.tsv", sep="\t", header=None,
        names=["sfla_id", "dmel_uniprot", "dmel_desc",
               "pident", "length", "evalue", "bitscore"])
    dmel = pd.read_csv(WILD / "dmel_gene_chrom.tsv", sep="\t")
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
                  .merge(chrom_lookup, left_on="dmel_symbol",
                         right_on="symbol", how="left"))
    if "chrom_x" in genes.columns:
        genes = genes.rename(columns={"chrom_x": "scaffold"})
    elif "chrom" in genes.columns:
        genes = genes.rename(columns={"chrom": "scaffold"})

    sub = genes[(genes["scaffold"] == CHR597) & (genes["muller"].notna())].copy()
    sub["mid_mb"] = (sub["start"] + sub["end"]) / 2 / 1e6
    sub = sub.sort_values("mid_mb")

    n_bins = 60
    edges = np.linspace(sub["mid_mb"].min(), sub["mid_mb"].max(), n_bins + 1)
    cents = (edges[:-1] + edges[1:]) / 2
    fcs, fds = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = sub[(sub["mid_mb"] >= lo) & (sub["mid_mb"] < hi)]
        if len(s) == 0:
            fcs.append(np.nan); fds.append(np.nan)
        else:
            fcs.append((s["muller"] == "C").mean())
            fds.append((s["muller"] == "D").mean())
    fcs, fds = np.array(fcs), np.array(fds)
    diff = fcs - fds
    sc = np.where(np.diff(np.sign(diff)) != 0)[0]
    break_mb = float(cents[sc[len(sc) // 2]]) if len(sc) else np.nan
    size_mb = float(genes[genes["scaffold"] == CHR597]["end"].max() / 1e6)
    return size_mb, break_mb


def load_sfla_sizes() -> dict[str, float]:
    prot = pd.read_csv(WILD / "sfla_v2_proteome.tsv", sep="\t")
    targets = [CHR439, "chr_ScDA7r2_126_HRSCAF_325", CHR597,
               "chr_ScDA7r2_110_HRSCAF_295", "chr_ScDA7r2_2_HRSCAF_23"]
    return {s: prot[prot["chrom"] == s]["end"].max() / 1e6 for s in targets}




def main():
    chr597_size, break_mb = load_chr597_break()
    sfla_sizes = load_sfla_sizes()


    rows = [
        ("Muller A", [("X",  DMEL_SIZES_MB["X"])],   CHR439,                       None),
        ("Muller B", [("2L", DMEL_SIZES_MB["2L"])], "chr_ScDA7r2_126_HRSCAF_325",   None),
        ("Muller D + C\n(2R + 3L fusion)",
                      [("3L", DMEL_SIZES_MB["3L"]),
                       ("2R", DMEL_SIZES_MB["2R"])], CHR597,                  break_mb),
        ("Muller E", [("3R", DMEL_SIZES_MB["3R"])], "chr_ScDA7r2_110_HRSCAF_295",   None),
        ("Muller F", [("4",  DMEL_SIZES_MB["4"])],  "chr_ScDA7r2_2_HRSCAF_23",      None),
    ]

    plt.rcParams.update({
        "font.family": "Helvetica",
        "font.size": 9,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.spines.left": False,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "svg.fonttype": "none",
    })

    fig, axC = plt.subplots(figsize=(7.0, 5.2))
    fig.subplots_adjust(left=0.16, right=0.92, top=0.92, bottom=0.10)
    axC.set_clip_on(False)

    bar_h = 0.36
    gap = 1.45
    n_rows = len(rows)
    y_max = n_rows * gap + 0.5
    x_max = max(sfla_sizes.values()) * 1.04

    dmel_artists, sfla_artists, ribbon_artists, anno_artists = [], [], [], []

    def _gid(art, gid):
        art.set_gid(gid)
        return art

    for ri, (label, dmel_list, sfla, fusion_break) in enumerate(rows):
        y_center = (n_rows - 1 - ri) * gap
        y_dmel_bottom = y_center + 0.35
        y_dmel_top    = y_dmel_bottom + bar_h
        y_sfla_top    = y_center - 0.05
        y_sfla_bottom = y_sfla_top - bar_h
        sfla_size = sfla_sizes[sfla]

        if fusion_break is None:
            chrom, sz = dmel_list[0]
            muller = DMEL_TO_MULLER[chrom]
            mid = sfla_size / 2
            x_dmel_lo = mid - sz / 2
            x_dmel_hi = mid + sz / 2

            poly = patches.Polygon(
                [(x_dmel_lo, y_dmel_bottom),
                 (x_dmel_hi, y_dmel_bottom),
                 (sfla_size, y_sfla_top),
                 (0,         y_sfla_top)],
                closed=True,
                facecolor=MULLER_COLORS[muller], alpha=0.16,
                edgecolor="none")
            poly.set_clip_on(False)
            axC.add_patch(_gid(poly, f"ribbon_{muller}"))
            ribbon_artists.append(poly)

            dmel_rect = patches.Rectangle(
                (x_dmel_lo, y_dmel_bottom), sz, bar_h,
                facecolor=MULLER_COLORS[muller], edgecolor="black",
                linewidth=1.4, alpha=0.85)
            dmel_rect.set_clip_on(False)
            axC.add_patch(_gid(dmel_rect, f"dmel_{chrom}_bar"))
            dmel_artists.append(dmel_rect)

            if sz / x_max > 0.05:
                t = axC.text(x_dmel_lo + sz / 2, y_dmel_bottom + bar_h / 2,
                             chrom, ha="center", va="center",
                             fontsize=8, color="white", fontweight="bold",
                             clip_on=False)
                _gid(t, f"dmel_{chrom}_label")
                dmel_artists.append(t)
                t = axC.text(x_dmel_hi + 0.6, y_dmel_bottom + bar_h / 2,
                             f"{sz:.1f} Mb", ha="left", va="center",
                             fontsize=6.5, color="dimgrey", clip_on=False)
                _gid(t, f"dmel_{chrom}_size")
                dmel_artists.append(t)
            else:
                t = axC.text(x_dmel_hi + 0.6, y_dmel_bottom + bar_h / 2,
                             f"{chrom}  ({sz:.1f} Mb)", ha="left",
                             va="center", fontsize=7, color="black",
                             fontweight="bold", clip_on=False)
                _gid(t, f"dmel_{chrom}_label")
                dmel_artists.append(t)

            sfla_rect = patches.Rectangle(
                (0, y_sfla_bottom), sfla_size, bar_h,
                facecolor=MULLER_COLORS[muller], edgecolor="black",
                linewidth=0.6, alpha=0.85)
            sfla_rect.set_clip_on(False)
            axC.add_patch(_gid(sfla_rect, f"sfla_{short_name(sfla)}_bar"))
            sfla_artists.append(sfla_rect)

            if sfla_size / x_max > 0.07:
                t = axC.text(sfla_size / 2, y_sfla_bottom + bar_h / 2,
                             short_name(sfla),
                             ha="center", va="center", fontsize=7.5,
                             color="white", fontweight="bold",
                             clip_on=False)
                _gid(t, f"sfla_{short_name(sfla)}_label")
                sfla_artists.append(t)
                t = axC.text(sfla_size + 0.6, y_sfla_bottom + bar_h / 2,
                             f"{sfla_size:.1f} Mb",
                             ha="left", va="center", fontsize=6.5,
                             color="dimgrey", clip_on=False)
                _gid(t, f"sfla_{short_name(sfla)}_size")
                sfla_artists.append(t)
            else:
                t = axC.text(sfla_size + 0.6, y_sfla_bottom + bar_h / 2,
                             f"{short_name(sfla)}  ({sfla_size:.1f} Mb)",
                             ha="left", va="center", fontsize=7,
                             color="black", fontweight="bold",
                             clip_on=False)
                _gid(t, f"sfla_{short_name(sfla)}_label")
                sfla_artists.append(t)

        else:
            block_widths = [(0, fusion_break, "3L", "D"),
                            (fusion_break, sfla_size, "2R", "C")]

            for blo, bhi, chrom, muller in block_widths:
                sz = DMEL_SIZES_MB[chrom]
                mid = (blo + bhi) / 2
                x_dmel_lo = mid - sz / 2
                x_dmel_hi = mid + sz / 2
                poly = patches.Polygon(
                    [(x_dmel_lo, y_dmel_bottom),
                     (x_dmel_hi, y_dmel_bottom),
                     (bhi,       y_sfla_top),
                     (blo,       y_sfla_top)],
                    closed=True,
                    facecolor=MULLER_COLORS[muller], alpha=0.16,
                    edgecolor="none")
                poly.set_clip_on(False)
                axC.add_patch(_gid(poly, f"ribbon_{muller}_fusion"))
                ribbon_artists.append(poly)

            for blo, bhi, chrom, muller in block_widths:
                seg = patches.Rectangle(
                    (blo, y_sfla_bottom), bhi - blo, bar_h,
                    facecolor=MULLER_COLORS[muller], edgecolor="none",
                    alpha=0.85)
                seg.set_clip_on(False)
                axC.add_patch(_gid(seg, f"sfla_chr_597_segment_{muller}"))
                sfla_artists.append(seg)
            outline = patches.Rectangle(
                (0, y_sfla_bottom), sfla_size, bar_h,
                facecolor="none", edgecolor="black", linewidth=0.6)
            outline.set_clip_on(False)
            axC.add_patch(_gid(outline, "sfla_chr_597_outline"))
            sfla_artists.append(outline)

            for blo, bhi, chrom, muller in block_widths:
                sz = DMEL_SIZES_MB[chrom]
                mid = (blo + bhi) / 2
                x_dmel_lo = mid - sz / 2
                x_dmel_hi = mid + sz / 2
                rect = patches.Rectangle(
                    (x_dmel_lo, y_dmel_bottom), sz, bar_h,
                    facecolor=MULLER_COLORS[muller], edgecolor="black",
                    linewidth=1.4, alpha=0.85)
                rect.set_clip_on(False)
                axC.add_patch(_gid(rect, f"dmel_{chrom}_bar"))
                dmel_artists.append(rect)

                t = axC.text(x_dmel_lo + sz / 2, y_dmel_bottom + bar_h / 2,
                             chrom, ha="center", va="center",
                             fontsize=8, color="white", fontweight="bold",
                             clip_on=False)
                _gid(t, f"dmel_{chrom}_label")
                dmel_artists.append(t)

                t = axC.text(x_dmel_hi + 0.5, y_dmel_bottom + bar_h / 2,
                             f"{sz:.1f}", ha="left", va="center",
                             fontsize=6.0, color="dimgrey",
                             clip_on=False)
                _gid(t, f"dmel_{chrom}_size")
                dmel_artists.append(t)

            line, = axC.plot(
                [fusion_break, fusion_break],
                [y_sfla_bottom - 0.04, y_sfla_top + 0.04],
                color="black", linewidth=0.9, linestyle="--",
                alpha=0.75, clip_on=False)
            _gid(line, "fusion_break_marker")
            anno_artists.append(line)

            t = axC.text(fusion_break, y_sfla_bottom - 0.16,
                         f"break ~{fusion_break:.0f} Mb",
                         ha="center", va="top", fontsize=6.5,
                         color="black", clip_on=False)
            _gid(t, "fusion_break_label")
            anno_artists.append(t)

            t = axC.text(sfla_size / 2, y_sfla_bottom + bar_h / 2,
                         short_name(sfla),
                         ha="center", va="center", fontsize=7.5,
                         color="white", fontweight="bold",
                         clip_on=False)
            _gid(t, "sfla_chr_597_label")
            sfla_artists.append(t)

            t = axC.text(sfla_size + 0.6, y_sfla_bottom + bar_h / 2,
                         f"{sfla_size:.1f} Mb",
                         ha="left", va="center", fontsize=6.5,
                         color="dimgrey", clip_on=False)
            _gid(t, "sfla_chr_597_size")
            sfla_artists.append(t)

        t = axC.text(-x_max * 0.03, y_center, label,
                     ha="right", va="center", fontsize=8,
                     fontweight="bold", clip_on=False)
        _gid(t, f"row_label_{ri}")
        anno_artists.append(t)

    t = axC.text(x_max * 0.5, y_max - 0.55,
                 "top of each pair: D. melanogaster chromosomes   |   "
                 "bottom: S. flava v2 scaffolds",
                 ha="center", va="bottom", fontsize=8, fontstyle="italic",
                 color="dimgrey", clip_on=False)
    _gid(t, "header_strip")
    anno_artists.append(t)

    axC.set_xlim(-x_max * 0.06, x_max + x_max * 0.05)
    axC.set_ylim(-0.55, y_max - 0.15)
    axC.set_xlabel("Length (Mb)")
    axC.set_yticks([])
    axC.spines["left"].set_visible(False)

    for ext in ("png", "pdf", "svg"):
        out = OUT_PREFIX.with_suffix(f".{ext}")
        fig.savefig(out, dpi=300 if ext == "png" else None,
                    bbox_inches="tight")
        print(f"wrote {out}")
    plt.close(fig)

    svg_path = OUT_PREFIX.with_suffix(".svg")
    cleaned = postprocess_svg(svg_path.read_text())
    svg_path.write_text(cleaned)
    print(f"post-processed {svg_path}")


def postprocess_svg(svg: str) -> str:
    svg = re.sub(r'\s*clip-path\s*=\s*"url\(#[^"]+\)"', "", svg)
    svg = re.sub(r"<clipPath[^>]*>.*?</clipPath>\s*", "", svg, flags=re.S)
    prefixes = {
        "ribbons":           ("ribbon_",),
        "dmel_chromosomes":  ("dmel_",),
        "sfla_scaffolds":    ("sfla_",),
        "annotations":       ("row_label_", "header_strip",
                              "fusion_break_marker", "fusion_break_label"),
    }

    def matches_prefix(gid: str, prefs: tuple[str, ...]) -> bool:
        return any(gid.startswith(p) for p in prefs)

    pat = re.compile(r'(<g [^>]*id="([^"]+)"[^>]*>.*?</g>)', re.S)
    chunks: list[tuple[str, str]] = []  # (gid, raw_chunk)
    last_end = 0
    out_parts: list[str] = []
    for m in pat.finditer(svg):
        chunk = m.group(1)
        gid = m.group(2)
        out_parts.append(svg[last_end:m.start()])
        out_parts.append(("__CHUNK__", gid, chunk))
        last_end = m.end()
    out_parts.append(svg[last_end:])

    rebuilt: list[str] = []
    current_group: str | None = None
    for part in out_parts:
        if isinstance(part, tuple) and part[0] == "__CHUNK__":
            _, gid, chunk = part
            target = None
            for grp, prefs in prefixes.items():
                if matches_prefix(gid, prefs):
                    target = grp
                    break
            if target != current_group:
                if current_group is not None:
                    rebuilt.append("</g>\n")
                if target is not None:
                    rebuilt.append(f'<g id="{target}">\n')
                current_group = target
            rebuilt.append(chunk)
        else:
            if current_group is not None:
                rebuilt.append("</g>\n")
                current_group = None
            rebuilt.append(part)
    if current_group is not None:
        rebuilt.append("</g>\n")
    return "".join(rebuilt)


if __name__ == "__main__":
    main()
