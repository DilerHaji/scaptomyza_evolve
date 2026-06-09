#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.proj3d import proj_transform


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, **kwargs):
        super().__init__((0, 0), (0, 0), **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return float(np.min(zs))

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42


PALETTE = {
    "Barbarea_experimental": "#96CDFF",
    "Founders":              "#000000",
    "Mixed_experimental":    "#901442",
    "Turitus_experimental":  "#F3C43C",
    "Wild_Barbarea":         "#96CDFF",   # same host → same hue as lab Barbarea
    "Wild_Turitus":          "#F3C43C",   # same host → same hue as lab Turitus
}

WILD_TREATMENTS = {
    "AVB": "Wild_Barbarea", "PSB": "Wild_Barbarea", "RMB": "Wild_Barbarea",
    "AVT": "Wild_Turitus",  "PST": "Wild_Turitus",  "RMT": "Wild_Turitus",
}

WILD_CATEGORIES = {"Wild_Barbarea", "Wild_Turitus"}

TRAJECTORY_CATEGORIES = {"Barbarea_experimental",
                         "Turitus_experimental",
                         "Mixed_experimental"}

def classify(ind: str):
    m = re.match(r"^([BTM])([1-4])G(\d{2})$", ind)
    if m:
        trt, rp, gen = m.group(1), int(m.group(2)), int(m.group(3))
        trt_full = {"B": "Barbarea_experimental",
                    "T": "Turitus_experimental",
                    "M": "Mixed_experimental"}[trt]
        return trt_full, gen, rp
    m = re.match(r"^F([1-4])(G00)?$", ind)
    if m:
        return "Founders", 0, int(m.group(1))
    if ind in WILD_TREATMENTS:
        return WILD_TREATMENTS[ind], None, None
    return None, None, None


def draw_3d_axes(ax, df, pve, elev, azim, title=None, label_endpoints=False):
    for trt, color in PALETTE.items():
        if trt in TRAJECTORY_CATEGORIES:
            continue
        sub = df[df["treat"] == trt]
        if len(sub) == 0:
            continue
        edge = "none" if trt in WILD_CATEGORIES else "#222222"
        size = 36 if trt in WILD_CATEGORIES else 18
        ax.scatter(sub["PC1"], sub["PC2"], sub["PC3"],
                   c=color, s=size, alpha=0.95,
                   marker="o", edgecolors=edge, linewidths=0.4,
                   depthshade=False)

    for (rp, treat), grp in df.groupby(["rp", "treat"]):
        if pd.isna(rp) or len(grp) < 2:
            continue
        grp = grp.sort_values("gen")
        xs, ys, zs = grp["PC1"].values, grp["PC2"].values, grp["PC3"].values
        color = PALETTE[treat]
        if len(xs) >= 2:
            ax.plot(xs[:-1], ys[:-1], zs[:-1],
                    color=color, linewidth=1.8, alpha=0.9, solid_capstyle="round")
            arrow = Arrow3D([xs[-2], xs[-1]], [ys[-2], ys[-1]], [zs[-2], zs[-1]],
                            arrowstyle="-|>,head_length=1.0,head_width=0.7",
                            mutation_scale=42,
                            lw=2.2, color=color, alpha=1.0, zorder=10)
            ax.add_artist(arrow)
            if label_endpoints:
                tip_name = grp.iloc[-1]["ind"]
                ax.text(xs[-1], ys[-1], zs[-1], f"  {tip_name}",
                        color=color, fontsize=7, alpha=0.95,
                        ha="left", va="center", zorder=11)
        else:
            ax.plot(xs, ys, zs, color=color, linewidth=1.8, alpha=0.9)

    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.set_facecolor("white")
        a.pane.set_edgecolor("#BBBBBB")
        a.pane.set_alpha(1.0)
        a.set_major_locator(plt.MaxNLocator(4))
    ax.grid(False)

    ax.set_xlabel(f"PC1 ({pve[0]:.1f}%)", labelpad=4, fontsize=9)
    ax.set_ylabel(f"PC2 ({pve[1]:.1f}%)", labelpad=4, fontsize=9)
    ax.set_zlabel(f"PC3 ({pve[2]:.1f}%)", labelpad=4, fontsize=9)
    ax.set_box_aspect((float(pve[0]), float(pve[1]), float(pve[2])))
    ax.view_init(elev=elev, azim=azim)
    if title:
        ax.set_title(title, fontsize=9)
    ax.tick_params(axis="both", labelsize=7)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores",    required=True)
    ap.add_argument("--eigenvals", required=True)
    ap.add_argument("--out",       required=True,
                    help="Output path. A PNG and an SVG are written (extension is appended if missing).")
    ap.add_argument("--angle",     type=float, default=None, help="Azimuth (single view)")
    ap.add_argument("--elev",      type=float, default=None, help="Elevation (single view)")
    ap.add_argument("--title",     default=None, help="Override plot title (single view only)")
    ap.add_argument("--label-endpoints", action="store_true",
                    help="Label each trajectory's G10 endpoint with the sample name.")
    args = ap.parse_args()

    scores = pd.read_csv(args.scores)
    pve_df = pd.read_csv(args.eigenvals)
    pve = pve_df["pve"].values

    scores = scores.drop(columns=[c for c in ("treat", "gen", "rp", "spp")
                                   if c in scores.columns])

    meta = scores["ind"].apply(lambda x: pd.Series(classify(x), index=["treat", "gen", "rp"]))
    df = pd.concat([scores, meta], axis=1)
    df = df.dropna(subset=["treat"]).reset_index(drop=True)

    legend_handles = []
    for name, c in PALETTE.items():
        if name in TRAJECTORY_CATEGORIES:
            legend_handles.append(
                Line2D([0], [0], color=c, linewidth=2.2, label=name))
        else:
            edge = "none" if name in WILD_CATEGORIES else "#222222"
            legend_handles.append(
                Line2D([0], [0], marker="o", linestyle="",
                       markerfacecolor=c, markeredgecolor=edge,
                       markersize=6, label=name))

    out_base = re.sub(r"\.(png|svg|pdf)$", "", args.out, flags=re.IGNORECASE)

    if args.angle is not None and args.elev is not None:
        fig = plt.figure(figsize=(10, 9))
        ax = fig.add_subplot(111, projection="3d")
        title = args.title if args.title else f"elev={args.elev}, azim={args.angle}"
        draw_3d_axes(ax, df, pve, args.elev, args.angle, title=title,
                     label_endpoints=args.label_endpoints)
        fig.legend(handles=legend_handles, loc="upper right",
                   bbox_to_anchor=(0.98, 0.95), frameon=False, fontsize=9)
        plt.tight_layout()
    else:
        view_spec = [
            (20,  30), (20, 120), (20, 210),
            (40,  60), (40, 150), (40, 240),
            (60,  45), (60, 135), (60, 225),
        ]
        views = [(e, a, f"elev={e}, azim={a}") for (e, a) in view_spec]
        nrow, ncol = 3, 3
        fig = plt.figure(figsize=(4.8 * ncol, 4.6 * nrow))
        for i, (elev, azim, title) in enumerate(views, 1):
            ax = fig.add_subplot(nrow, ncol, i, projection="3d")
            draw_3d_axes(ax, df, pve, elev, azim, title=title,
                         label_endpoints=args.label_endpoints)
        fig.legend(handles=legend_handles, loc="lower center",
                   ncol=len(legend_handles), frameon=False, fontsize=10,
                   bbox_to_anchor=(0.5, 0.0))
        fig.subplots_adjust(left=0.04, right=0.96, top=0.96, bottom=0.06,
                            wspace=0.30, hspace=0.30)

    save_kw = {} if (args.angle is None and args.elev is None) \
              else dict(bbox_inches="tight")
    fig.savefig(f"{out_base}.png", dpi=200, pad_inches=0.3, **save_kw)
    fig.savefig(f"{out_base}.svg",            pad_inches=0.3, **save_kw)

if __name__ == "__main__":
    main()
