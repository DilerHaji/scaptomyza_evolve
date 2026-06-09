#!/usr/bin/env python3
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", required=True, help="PCA scores CSV")
    ap.add_argument("--n-pcs", type=int, default=10,
                    help="Number of PCs to use (default 10).")
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--seed", type=int, default=42)
    return ap.parse_args()


def classify(ind: str):
    m = re.match(r"^([BTM])([1-4])G10$", ind)
    if m:
        return m.group(1), int(m.group(2))
    m = re.match(r"^F([1-4])(G00)?$", ind)
    if m:
        return "F", int(m.group(1))
    return None, None


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    df = pd.read_csv(args.scores)
    pc_cols = [f"PC{i+1}" for i in range(args.n_pcs)]
    assert all(c in df.columns for c in pc_cols), \
        f"need {pc_cols} in {args.scores}"

    info = df["ind"].apply(lambda s: pd.Series(classify(s), index=["trt", "rep"]))
    df = pd.concat([df, info], axis=1).dropna(subset=["trt"])

    founders = df[df["trt"] == "F"]
    t0 = founders[pc_cols].mean(axis=0).values   # shape (n_pcs,)

    g10 = df[df["trt"].isin({"B", "T", "M"}) & df["ind"].str.endswith("G10")].copy()
    g10_vec = g10[pc_cols].values - t0           # shape (n_reps, n_pcs)
    g10_trt = g10["trt"].values
    g10_ind = g10["ind"].values

    assert len(g10) == 12, f"expected 12 G10 replicates, got {len(g10)}"

    def stat_from_vectors(vecs, trts):
        vB = vecs[trts == "B"].mean(axis=0)
        vT = vecs[trts == "T"].mean(axis=0)
        vM = vecs[trts == "M"].mean(axis=0)
        u  = 0.5 * (vB + vT)

        d_cent = float(np.linalg.norm(vM - u))
        denom = np.linalg.norm(vM) * np.linalg.norm(u)
        cos = float(np.dot(vM, u) / denom) if denom > 0 else 0.0
        cos = max(-1.0, min(1.0, cos))
        theta_deg = float(np.degrees(np.arccos(cos)))

        n_B = int((trts == "B").sum())
        n_T = int((trts == "T").sum())
        n_M = int((trts == "M").sum())
        n_tot = n_B + n_T + n_M
        inv_n_eff = (1.0 / n_M) + (0.25 / n_B) + (0.25 / n_T)
        n_eff = 1.0 / inv_n_eff
        L = vM - u
        ss_contrast = n_eff * float(np.sum(L * L))
        ss_within = 0.0
        for (trt, vbar) in (("B", vB), ("T", vT), ("M", vM)):
            reps = vecs[trts == trt]
            ss_within += float(np.sum((reps - vbar) ** 2))
        df_within = n_tot - 3
        F = ss_contrast / (ss_within / df_within) if ss_within > 0 else np.inf

        axis = vB - vT
        axis_norm = np.linalg.norm(axis)
        if axis_norm > 0:
            w = axis / axis_norm
            sB = float(vecs[trts == "B"].dot(w).mean())
            sT = float(vecs[trts == "T"].dot(w).mean())
            sM = float(vecs[trts == "M"].dot(w).mean())
            midpoint_BT = 0.5 * (sB + sT)
            delta1d = sM - midpoint_BT

            def orth_norm(v):
                return float(np.linalg.norm(v - (v.dot(w)) * w))
            orthB = orth_norm(vB)
            orthT = orth_norm(vT)
            orthM = orth_norm(vM)
            orth_contrast = orthM - 0.5 * (orthB + orthT)
        else:
            sB = sT = sM = 0.0
            delta1d = 0.0
            orthB = orthT = orthM = 0.0
            orth_contrast = 0.0

        return {
            "d_cent":    d_cent,
            "cos":       cos,
            "theta":     theta_deg,
            "F":         F,
            "delta_1d":  float(delta1d),
            "sB_mean":   sB,
            "sT_mean":   sT,
            "sM_mean":   sM,
            "orth_B":    orthB,
            "orth_T":    orthT,
            "orth_M":    orthM,
            "orth_contrast": float(orth_contrast),
            "ss_contrast": ss_contrast,
            "ss_within":   ss_within,
            "vB": vB, "vT": vT, "vM": vM, "u": u,
        }

    obs = stat_from_vectors(g10_vec, g10_trt)
    vB, vT, vM, u = obs["vB"], obs["vT"], obs["vM"], obs["u"]
    d_obs, theta_obs, cos_obs, F_obs = obs["d_cent"], obs["theta"], obs["cos"], obs["F"]
    d_norm = d_obs / np.linalg.norm(vM) if np.linalg.norm(vM) > 0 else np.nan

    n = len(g10_vec)
    d_null        = np.empty(args.n_perm)
    theta_null    = np.empty(args.n_perm)
    F_null        = np.empty(args.n_perm)
    delta1d_null  = np.empty(args.n_perm)
    orthcon_null  = np.empty(args.n_perm)
    for k in range(args.n_perm):
        perm = rng.permutation(n)
        trts_perm = g10_trt[perm]
        st = stat_from_vectors(g10_vec, trts_perm)
        d_null[k]       = st["d_cent"]
        theta_null[k]   = st["theta"]
        F_null[k]       = st["F"]
        delta1d_null[k] = st["delta_1d"]
        orthcon_null[k] = st["orth_contrast"]

    p_d     = (np.sum(d_null     >= d_obs    ) + 1) / (args.n_perm + 1)
    p_theta = (np.sum(theta_null >= theta_obs) + 1) / (args.n_perm + 1)
    p_F     = (np.sum(F_null     >= F_obs    ) + 1) / (args.n_perm + 1)

    delta1d_obs = obs["delta_1d"]
    p_delta1d   = (np.sum(np.abs(delta1d_null) >= abs(delta1d_obs)) + 1) / (args.n_perm + 1)

    orth_obs = obs["orth_contrast"]
    p_orth   = (np.sum(orthcon_null >= orth_obs) + 1) / (args.n_perm + 1)

    fig, axes = plt.subplots(1, 5, figsize=(22, 3.8))

    axes[0].hist(d_null, bins=50, color="#BBBBBB", edgecolor="white", linewidth=0.3)
    axes[0].axvline(d_obs, color="#901442", linewidth=2.0,
                    label=f"obs = {d_obs:.1f}")
    axes[0].set_xlabel("d = ‖v̄_M − (v̄_B + v̄_T)/2‖")
    axes[0].set_title(f"10-D centroid distance\np = {p_d:.4f}", fontsize=10)

    axes[1].hist(theta_null, bins=50, color="#BBBBBB", edgecolor="white", linewidth=0.3)
    axes[1].axvline(theta_obs, color="#901442", linewidth=2.0,
                    label=f"obs = {theta_obs:.1f}°")
    axes[1].set_xlabel("angle(v̄_M, (v̄_B + v̄_T)/2)  [degrees]")
    axes[1].set_title(f"Angle only\np = {p_theta:.4f}", fontsize=10)

    axes[2].hist(F_null, bins=50, color="#BBBBBB", edgecolor="white", linewidth=0.3)
    axes[2].axvline(F_obs, color="#901442", linewidth=2.0,
                    label=f"obs = {F_obs:.2f}")
    axes[2].set_xlabel("pseudo-F  (replicate-level)")
    axes[2].set_title(f"PERMANOVA-style\np = {p_F:.4f}", fontsize=10)

    axes[3].hist(delta1d_null, bins=50, color="#BBBBBB", edgecolor="white", linewidth=0.3)
    axes[3].axvline(delta1d_obs, color="#901442", linewidth=2.0,
                    label=f"obs = {delta1d_obs:+.2f}")
    axes[3].axvline(-abs(delta1d_obs), color="#901442", linewidth=2.0, linestyle="--")
    axes[3].set_xlabel("Δ = s_M − midpoint(B,T)   (on B-T axis)")
    axes[3].set_title(f"1D projection (on-axis)\np = {p_delta1d:.4f}", fontsize=10)

    axes[4].hist(orthcon_null, bins=50, color="#BBBBBB", edgecolor="white", linewidth=0.3)
    axes[4].axvline(orth_obs, color="#901442", linewidth=2.0,
                    label=f"obs = {orth_obs:+.2f}")
    axes[4].set_xlabel("‖v_M_orth‖ − mean(‖v_B_orth‖, ‖v_T_orth‖)")
    axes[4].set_title(f"Off-axis excess (M)\np = {p_orth:.4f}", fontsize=10)

    for ax in axes:
        ax.set_ylabel("count")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.legend(frameon=False, loc="upper right", fontsize=9)

    fig.suptitle(f"Subdivided- vs uniform-selection tests  (n_perm = {args.n_perm})",
                 fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    out = args.out_prefix
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out}_null.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{out}_null.svg",            bbox_inches="tight")

    summary = pd.DataFrame([{
        "d_observed":          round(float(d_obs), 3),
        "d_over_vM_norm":      round(float(d_norm), 3),
        "angle_deg":           round(float(theta_obs), 2),
        "cos_similarity":      round(float(cos_obs), 4),
        "pseudo_F_observed":   round(float(F_obs), 3),
        "ss_contrast":         round(float(obs["ss_contrast"]), 3),
        "ss_within":           round(float(obs["ss_within"]),   3),
        "p_d_permutation":     float(p_d),
        "p_theta_permutation": float(p_theta),
        "p_F_permutation":     float(p_F),
        "delta_1d_observed":   round(float(delta1d_obs), 3),
        "p_delta1d_permutation": float(p_delta1d),
        "orth_B":              round(float(obs["orth_B"]), 3),
        "orth_T":              round(float(obs["orth_T"]), 3),
        "orth_M":              round(float(obs["orth_M"]), 3),
        "orth_contrast_observed": round(float(orth_obs), 3),
        "p_orth_permutation":  float(p_orth),
        "n_perm":              args.n_perm,
        "norm_vB":             round(float(np.linalg.norm(vB)), 3),
        "norm_vT":             round(float(np.linalg.norm(vT)), 3),
        "norm_vM":             round(float(np.linalg.norm(vM)), 3),
        "norm_u":              round(float(np.linalg.norm(u)),  3),
    }])
    summary.to_csv(f"{out}_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
