#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch
import numpy as np

ROOT = Path(".")
OUT_BASE = ROOT / "final_plots/wild/section3_fig4a_conceptual"

C_B       = "#4C86A8"   # Barbarea / B host
C_T       = "#E07B54"   # Turritis / T host
C_GREY    = "#888888"
C_LIGHT   = "#EAEAEA"
C_TEXT    = "#222222"
C_PARA    = "#9C7AA8"   # parasitoid (minor)


def rbox(ax, xy, w, h, label, color, fontsize=9, text_color="white",
         bold=False, alpha=1.0):
    x, y = xy
    box = FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.03",
        facecolor=color, edgecolor="black", linewidth=0.8, alpha=alpha, zorder=3)
    ax.add_patch(box)
    fw = "bold" if bold else "normal"
    ax.text(x + w/2, y + h/2, label, ha="center", va="center",
            fontsize=fontsize, fontweight=fw, color=text_color, zorder=4)


def arrow(ax, start, end, color="#444444", lw=1.0, style="-|>", rad=0):
    cs = f"arc3,rad={rad}"
    a = FancyArrowPatch(start, end, arrowstyle=style, color=color, lw=lw,
                         connectionstyle=cs, zorder=2, mutation_scale=12)
    ax.add_patch(a)


def draw_lab(ax, cx, cy):
    ax.text(cx, cy + 1.5, "Lab regime", ha="center", va="bottom",
             fontsize=12, fontweight="bold", color=C_TEXT)
    ax.text(cx, cy + 1.20, "(idealized 50:50)", ha="center", va="bottom",
             fontsize=8.5, fontstyle="italic", color=C_GREY)
    R = 0.45
    bx, by = cx - 0.55, cy
    tx, ty = cx + 0.55, cy
    ax.add_patch(Circle((bx, by), R, facecolor=C_B,
                         edgecolor="black", linewidth=1.0, zorder=3))
    ax.add_patch(Circle((tx, ty), R, facecolor=C_T,
                         edgecolor="black", linewidth=1.0, zorder=3))
    ax.text(bx, by, "Barbarea", ha="center", va="center",
             color="white", fontweight="bold", fontsize=8.5)
    ax.text(tx, ty, "Turritis", ha="center", va="center",
             color="white", fontweight="bold", fontsize=8.5)
    arrow(ax, (bx + R + 0.05, cy + 0.20), (tx - R - 0.05, cy + 0.20),
           color=C_GREY, lw=0.8, style="<->")
    arrow(ax, (bx + R + 0.05, cy - 0.20), (tx - R - 0.05, cy - 0.20),
           color=C_GREY, lw=0.8, style="<->")
    ax.text(bx, cy - R - 0.25, "50%", ha="center", va="top",
             fontsize=8, color=C_TEXT)
    ax.text(tx, cy - R - 0.25, "50%", ha="center", va="top",
             fontsize=8, color=C_TEXT)
    ax.text(cx, cy - 1.3,
             "Strict subdivided selection",
             ha="center", va="center", fontsize=9, color=C_TEXT)
    ax.text(cx, cy - 1.55,
             "(hard, Dempster 1955)",
             ha="center", va="center", fontsize=8, fontstyle="italic",
             color=C_GREY)


def draw_wild(ax, cx, cy):
    ax.text(cx, cy + 1.5, "Wild regime", ha="center", va="bottom",
             fontsize=12, fontweight="bold", color=C_TEXT)
    ax.text(cx, cy + 1.20, "(observed in nature)", ha="center", va="bottom",
             fontsize=8.5, fontstyle="italic", color=C_GREY)
    sites = [("AV", cx - 1.20, 0.85),
             ("PS", cx,         0.78),
             ("RM", cx + 1.20, 0.65)]   # site name, x, Barbarea fraction
    for name, sx, bf in sites:
        rB = 0.45 * np.sqrt(bf)
        ax.add_patch(Circle((sx - 0.05, cy), rB,
                              facecolor=C_B, edgecolor="black",
                              linewidth=0.9, zorder=3))
        rT = 0.45 * np.sqrt(1 - bf)
        ax.add_patch(Circle((sx + rB + rT * 0.7, cy - rB + rT),
                              rT, facecolor=C_T, edgecolor="black",
                              linewidth=0.9, zorder=4))
        ax.text(sx, cy - 0.85, name, ha="center", va="center",
                 fontsize=8.5, color=C_TEXT, fontweight="bold")
        ax.text(sx, cy - 1.10,
                 f"{int(bf*100)}% / {int((1-bf)*100)}%",
                 ha="center", va="center", fontsize=7, color=C_GREY)
    ax.text(cx, cy - 1.55,
             "Barbarea-dominant × geographic structure",
             ha="center", va="center", fontsize=9, color=C_TEXT)
    ax.text(cx, cy - 1.80,
             "(asymmetric host use, sites $\\approx$ hosts)",
             ha="center", va="center", fontsize=8, fontstyle="italic",
             color=C_GREY)


def draw_predictions(ax, x0, y0, w, lab=True):
    h = 1.55
    title = "Lab predictions" if lab else "Wild predictions"
    color = C_B if lab else C_T
    ax.add_patch(FancyBboxPatch((x0, y0 + h - 0.30), w, 0.30,
                                  boxstyle="round,pad=0.0",
                                  facecolor=color, edgecolor="black",
                                  linewidth=0.7, zorder=3))
    ax.text(x0 + w/2, y0 + h - 0.15, title, ha="center", va="center",
             color="white", fontsize=9.5, fontweight="bold", zorder=4)
    ax.add_patch(FancyBboxPatch((x0, y0), w, h - 0.30,
                                  boxstyle="round,pad=0.0",
                                  facecolor="white",
                                  edgecolor="black", linewidth=0.7, zorder=2))
    if lab:
        lines = [
            "• strong B$\\leftrightarrow$T antagonism",
            "• intermediate-AF excess",
            "• polymorphism excess",
            "  (HKA-positive)",
        ]
    else:
        lines = [
            "• partial intermediate-AF",
            "  excess  (p = 0.002, observed)",
            "• no polymorphism excess",
            "  (HKA-negative, observed)",
        ]
    for i, ln in enumerate(lines):
        ax.text(x0 + 0.10, y0 + h - 0.55 - i * 0.22, ln,
                 ha="left", va="center", fontsize=8.5, color=C_TEXT)


def main() -> None:
    plt.rcParams.update({
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    })

    fig, ax = plt.subplots(figsize=(10, 5.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5.6)
    ax.set_aspect("equal")
    ax.axis("off")

    LAB_CX = 2.0;  WILD_CX = 7.0;  REGIME_CY = 4.0
    draw_lab(ax, LAB_CX, REGIME_CY)
    draw_wild(ax, WILD_CX, REGIME_CY)

    arrow(ax, (LAB_CX + 1.40, REGIME_CY),
                (WILD_CX - 1.85, REGIME_CY),
                color="#555555", lw=1.4)
    ax.text((LAB_CX + WILD_CX) / 2, REGIME_CY + 0.20,
             "Real ecology", ha="center", va="bottom",
             fontsize=9, fontstyle="italic", color="#555555")

    mid_x = (LAB_CX + WILD_CX) / 2
    ax.text(mid_x, REGIME_CY - 0.20,
             "+ asymmetric host use\n+ geographic substructure",
             ha="center", va="top", fontsize=8, color=C_TEXT,
             linespacing=1.4)
    ax.text(mid_x, REGIME_CY - 0.95,
             "(+ tritrophic pressure;\n  not addressed here)",
             ha="center", va="top", fontsize=6.5, fontstyle="italic",
             color=C_GREY, linespacing=1.3)

    PRED_Y = 0.30
    draw_predictions(ax, 0.40, PRED_Y, 3.6, lab=True)
    draw_predictions(ax, 5.20, PRED_Y, 4.4, lab=False)

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
