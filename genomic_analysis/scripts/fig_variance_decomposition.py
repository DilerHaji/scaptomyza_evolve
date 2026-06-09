#!/usr/bin/env python3

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

PHI = 2.11
N_POOL = 80
N_HAP = 2 * N_POOL
R_TYP = 40
NE = 500
K_REP = 4         # replicates per treatment
GENS = np.arange(1, 11)

C_SEL   = "#2ca02c"   # green — selection signal
C_DRIFT = "#D64545"   # red — drift
C_POOL  = "#2E86AB"   # blue — pool sampling
C_SEQ   = "#E8CC5B"   # yellow — sequencing
C_MEAS  = "#7B68AE"   # purple — total measurement
C_NOISE = "#888888"   # grey — total noise
C_BG    = "#F5F5F5"


def var_drift(p, t, Ne=NE):
    return p * (1 - p) * (1 - (1 - 1/(2*Ne))**t)

def var_pool(p, n_hap=N_HAP):
    return p * (1 - p) * PHI / n_hap

def var_seq(p, R=R_TYP):
    return p * (1 - p) / R

def var_meas(p, n_hap=N_HAP, R=R_TYP):
    return var_pool(p, n_hap) + var_seq(p, R)

def var_delta_single(p, t):
    return var_drift(p, t) + 2 * var_meas(p)

def var_delta_replicate_mean(p, t, k=K_REP):
    return var_delta_single(p, t) / k

def expected_delta_p(p, s, t):
    return s * p * (1 - p) * t

def snr_single(p, s, t):
    signal = expected_delta_p(p, s, t)
    noise = np.sqrt(var_delta_single(p, t))
    return np.abs(signal) / noise

def snr_replicated(p, s, t, k=K_REP):
    signal = expected_delta_p(p, s, t)
    noise = np.sqrt(var_delta_replicate_mean(p, t, k))
    return np.abs(signal) / noise

def min_detectable_s(p, t, k=K_REP, z=2.0):
    noise_sd = np.sqrt(var_delta_replicate_mean(p, t, k))
    # z * noise_sd = s * p * (1-p) * t
    return z * noise_sd / (p * (1 - p) * t)


def main():
    fig = plt.figure(figsize=(14, 11))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.32,
                           left=0.08, right=0.96, top=0.92, bottom=0.07)

    p_mid = 0.25  # representative frequency

    ax_a = fig.add_subplot(gs[0, 0])

    p_range = np.linspace(0.02, 0.98, 200)
    t = 10

    v_d = np.array([var_drift(p, t) for p in p_range])
    v_p = np.array([var_pool(p) for p in p_range])
    v_s = np.array([var_seq(p) for p in p_range])

    ax_a.fill_between(p_range, 0, v_s, alpha=0.4, color=C_SEQ,
                      label=f"Sequencing (R={R_TYP})")
    ax_a.fill_between(p_range, v_s, v_s + v_p, alpha=0.4, color=C_POOL,
                      label=f"Pool sampling (n={N_POOL}, φ={PHI})")
    ax_a.fill_between(p_range, v_s + v_p, v_s + v_p + v_d, alpha=0.4,
                      color=C_DRIFT, label=f"Drift ({t} gen, Ne={NE})")

    for s, ls in [(0.05, ":"), (0.1, "--"), (0.2, "-")]:
        delta = np.array([expected_delta_p(p, s, t) for p in p_range])
        ax_a.plot(p_range, delta**2, color=C_SEL, lw=1.5, ls=ls,
                  label=f"Selection E[Δp]² (s={s})")

    ax_a.set_xlabel("True allele frequency (p)")
    ax_a.set_ylabel("Variance / Expected Δp²")
    ax_a.set_title("A. Noise components vs. selection signal (F→G10)",
                   fontweight="bold", loc="left", fontsize=9)
    ax_a.legend(fontsize=6, loc="upper center", framealpha=0.9, ncol=2)
    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(bottom=0)

    ax_b = fig.add_subplot(gs[0, 1])

    s_vals = [0.02, 0.05, 0.1, 0.2, 0.5]
    colors_s = plt.cm.Greens(np.linspace(0.3, 0.9, len(s_vals)))

    noise_1rep = [np.sqrt(var_delta_single(p_mid, t)) for t in GENS]
    noise_2sd_1rep = [2 * n for n in noise_1rep]

    noise_4rep = [np.sqrt(var_delta_replicate_mean(p_mid, t)) for t in GENS]
    noise_2sd_4rep = [2 * n for n in noise_4rep]

    ax_b.fill_between(GENS, 0, noise_2sd_1rep, alpha=0.12, color=C_NOISE,
                      label="2 SD noise (1 rep)")
    ax_b.fill_between(GENS, 0, noise_2sd_4rep, alpha=0.25, color=C_MEAS,
                      label="2 SD noise (4 rep avg)")

    for s, c in zip(s_vals, colors_s):
        delta = [np.abs(expected_delta_p(p_mid, s, t)) for t in GENS]
        ax_b.plot(GENS, delta, color=c, lw=2, marker="o", ms=3,
                  label=f"s = {s}")

    ax_b.set_xlabel("Generations (t)")
    ax_b.set_ylabel("|E[Δp]| or noise SD")
    ax_b.set_title(f"B. Selection signal vs. noise envelope (p={p_mid})",
                   fontweight="bold", loc="left", fontsize=9)
    ax_b.legend(fontsize=6, loc="upper left", framealpha=0.9, ncol=2)
    ax_b.set_xlim(0.5, 10.5)
    ax_b.set_xticks(GENS)
    ax_b.set_ylim(bottom=0)

    ax_c = fig.add_subplot(gs[1, 0])

    s_show = [0.05, 0.1, 0.2]
    ls_map = {"1 rep": "--", "4 rep": "-"}

    for s, c in zip(s_show, [C_SEL, "#1a7a1a", "#0d4d0d"]):
        snr_1 = [snr_single(p_mid, s, t) for t in GENS]
        snr_4 = [snr_replicated(p_mid, s, t, K_REP) for t in GENS]
        ax_c.plot(GENS, snr_1, color=c, lw=1.5, ls="--", alpha=0.5,
                  marker="^", ms=3)
        ax_c.plot(GENS, snr_4, color=c, lw=2.5, ls="-",
                  marker="o", ms=4, label=f"s={s} (4 rep)")

    ax_c.axhline(2.0, color="black", lw=0.8, ls=":", alpha=0.4)
    ax_c.text(0.7, 2.1, "z = 2 (p ≈ 0.05)", fontsize=6, color="black", alpha=0.6)
    ax_c.axhline(3.0, color="black", lw=0.8, ls=":", alpha=0.3)
    ax_c.text(0.7, 3.1, "z = 3 (p ≈ 0.003)", fontsize=6, color="black", alpha=0.5)

    ax_c.plot([], [], color="grey", lw=1.5, ls="--", label="1 replicate (dashed)")
    ax_c.plot([], [], color="grey", lw=2.5, ls="-", label="4 replicates (solid)")

    ax_c.set_xlabel("Generations (t)")
    ax_c.set_ylabel("Signal-to-noise ratio (|E[Δp]| / SD[Δp])")
    ax_c.set_title(f"C. Detection power: 1 vs. 4 replicates (p={p_mid})",
                   fontweight="bold", loc="left", fontsize=9)
    ax_c.legend(fontsize=6, loc="upper left", framealpha=0.9)
    ax_c.set_xlim(0.5, 10.5)
    ax_c.set_xticks(GENS)
    ax_c.set_ylim(bottom=0)

    ax_d = fig.add_subplot(gs[1, 1])

    p_vals = [0.1, 0.25, 0.5]
    p_colors = [C_POOL, C_MEAS, C_DRIFT]
    rep_vals = [1, 4]

    for p_val, pc in zip(p_vals, p_colors):
        for k, ls, alpha in zip(rep_vals, ["--", "-"], [0.4, 1.0]):
            min_s = [min_detectable_s(p_val, t, k=k, z=2.0) for t in GENS]
            label = f"p={p_val}, {k} rep" if p_val == 0.25 or k == 4 else None
            ax_d.plot(GENS, min_s, color=pc, lw=2, ls=ls, alpha=alpha,
                      marker="s" if k == 4 else "^", ms=3, label=label)

    for s_ref, label in [(0.1, "s = 0.1"), (0.05, "s = 0.05")]:
        ax_d.axhline(s_ref, color=C_SEL, lw=0.8, ls=":", alpha=0.5)
        ax_d.text(10.6, s_ref, label, fontsize=6, color=C_SEL, va="center")

    ax_d.set_xlabel("Generations (t)")
    ax_d.set_ylabel("Minimum detectable s (z = 2)")
    ax_d.set_title(f"D. Detection limit vs. duration",
                   fontweight="bold", loc="left", fontsize=9)
    ax_d.legend(fontsize=6, loc="upper right", framealpha=0.9)
    ax_d.set_xlim(0.5, 10.5)
    ax_d.set_xticks(GENS)
    ax_d.set_ylim(0, 0.6)

    ax_d.text(5.5, 0.55,
              "Solid = 4 replicates,  Dashed = 1 replicate\n"
              "Intermediate AF (p≈0.25–0.5) most detectable",
              fontsize=6.5, ha="center", va="top", color="#555555",
              bbox=dict(boxstyle="round,pad=0.3", facecolor=C_BG,
                        edgecolor="#CCCCCC"))

    fig.suptitle(
        "Allele frequency variance decomposition: selection detectability\n"
        f"n_pool = {N_POOL}  |  φ = {PHI}  |  R = {R_TYP}×  |  "
        f"Ne = {NE}  |  k = {K_REP} replicates",
        fontsize=11, fontweight="bold", y=0.98
    )

    out = ("./"
           "final_plots/fig_variance_decomposition")
    fig.savefig(out + ".png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(out + ".pdf", bbox_inches="tight", facecolor="white")

    for s in [0.02, 0.05, 0.1, 0.2, 0.5]:
        dp = expected_delta_p(p_mid, s, 10)
        snr1 = snr_single(p_mid, s, 10)
        snr4 = snr_replicated(p_mid, s, 10)

    for t in [1, 5, 10]:
        s1 = min_detectable_s(p_mid, t, k=1)
        s4 = min_detectable_s(p_mid, t, k=4)

if __name__ == "__main__":
    main()
