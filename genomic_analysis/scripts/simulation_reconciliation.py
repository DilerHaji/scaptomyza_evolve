#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd
from scipy import stats as sp_stats
from itertools import combinations

N_SIMS = 50
N_SITES = 50_000
NE_TRUE = 250
N_REPS = 4
TIMEPOINTS = [0, 1, 2, 6, 7, 8, 9]  # G00 through G09 (4-rep gens only)
N_EFF = 29     # effective pool diploids
DEPTH_DEFAULT = 83
SEED_BASE = 12345

DEPTH_TRAJECTORIES = {
    "B": [179, 87, 88, 75, 84, 68, 73],
    "T": [179, 72, 73, 84, 76, 78, 86],
    "M": [179, 85, 103, 76, 70, 62, 66],
}


def simulate_wf(ne, n_sites, timepoints, s_matrix, rng):
    max_gen = max(timepoints)
    n_tp = len(timepoints)

    p0 = rng.uniform(0.1, 0.9, n_sites)
    current_p = p0.copy()

    freq_true = np.zeros((n_tp, n_sites))
    tp_idx = 0

    for g in range(max_gen + 1):
        if tp_idx < n_tp and g == timepoints[tp_idx]:
            freq_true[tp_idx, :] = current_p
            tp_idx += 1

        if g < max_gen:
            s = s_matrix[g, :]
            if np.any(s != 0):
                delta = s * current_p * (1 - current_p)
                current_p = current_p + delta

            n_hap = 2 * ne
            new_counts = rng.binomial(n_hap, np.clip(current_p, 1e-4, 1-1e-4))
            current_p = new_counts / n_hap
            current_p = np.clip(current_p, 1e-4, 1-1e-4)

    return freq_true


def observe_poolseq(freq_true, n_eff, depths, rng):
    n_tp, n_sites = freq_true.shape
    if np.isscalar(depths):
        depths = np.full(n_tp, depths)

    freq_obs = np.zeros_like(freq_true)
    for t in range(n_tp):
        n_hap_pool = 2 * n_eff
        pool_counts = rng.binomial(n_hap_pool, freq_true[t, :])
        p_pool = pool_counts / n_hap_pool
        d = int(depths[t])
        seq_counts = rng.binomial(d, p_pool)
        freq_obs[t, :] = seq_counts / d

    return freq_obs


def compute_F_regression(freq_obs, timepoints, n_eff, depths):
    n_reps, n_tp, n_sites = freq_obs.shape
    k = n_reps

    gens, F_raw, F_bio = [], [], []
    for ti in range(n_tp):
        t = timepoints[ti]
        if t == 0:
            continue  # skip founder
        p_reps = freq_obs[:, ti, :]  # k × n_sites
        pbar = p_reps.mean(axis=0)
        maf = np.minimum(pbar, 1 - pbar)
        ok = maf >= 0.05

        if ok.sum() < 100:
            continue

        ss = np.sum((p_reps[:, ok] - pbar[ok][None, :])**2, axis=0) / (k - 1)
        het = pbar[ok] * (1 - pbar[ok])
        valid = het > 1e-10
        F = np.mean(ss[valid]) / np.mean(het[valid])
        F_raw.append(F)

        if np.isscalar(depths):
            d = depths
        else:
            d = depths[ti]
        f_noise = 1/(2*n_eff) + 1/d
        F_bio.append(F - f_noise)
        gens.append(t)

    if len(gens) < 3:
        return np.nan, np.nan

    sr, *_ = sp_stats.linregress(gens, F_raw)
    sb, *_ = sp_stats.linregress(gens, F_bio)
    ne_raw = 1/(2*sr) if sr > 0 else np.inf
    ne_bio = 1/(2*sb) if sb > 0 else np.inf
    return ne_raw, ne_bio


def compute_poolseq_ne(freq_obs, timepoints, n_eff, depths):
    n_reps, n_tp, n_sites = freq_obs.shape
    t0_idx, t1_idx = 1, n_tp - 1
    t_diff = timepoints[t1_idx] - timepoints[t0_idx]

    ne_per_rep = []
    for r in range(n_reps):
        p0 = freq_obs[r, t0_idx, :]
        pt = freq_obs[r, t1_idx, :]
        pbar = (p0 + pt) / 2
        ok = (np.minimum(pbar, 1-pbar) >= 0.05) & np.isfinite(p0) & np.isfinite(pt)

        if ok.sum() < 100:
            continue

        het = pbar[ok] * (1 - pbar[ok])
        valid = het > 1e-10
        Fc = (p0[ok][valid] - pt[ok][valid])**2 / het[valid]

        if np.isscalar(depths):
            d0, dt = depths, depths
        else:
            d0, dt = depths[t0_idx], depths[t1_idx]
        inv_neff0 = 1/(2*n_eff) + 1/d0 - 1/(2*n_eff*d0)
        inv_neff1 = 1/(2*n_eff) + 1/dt - 1/(2*n_eff*dt)
        Fc_corrected = Fc - inv_neff0 - inv_neff1

        phi_mean = np.mean(Fc_corrected)
        if phi_mean > 0:
            ne_r = t_diff / (2 * phi_mean)
        else:
            ne_r = np.inf
        if np.isfinite(ne_r) and 0 < ne_r < 1e7:
            ne_per_rep.append(ne_r)

    return float(np.mean(ne_per_rep)) if ne_per_rep else np.nan


def compute_G_k2(freq_obs, timepoints):
    n_reps, n_tp, n_sites = freq_obs.shape
    tp_idx = list(range(1, n_tp))
    T = len(tp_idx)

    G_per_rep = []
    for r in range(n_reps):
        deltas = np.diff(freq_obs[r, tp_idx, :], axis=0)  # (T-1) × n_sites
        pbar = np.mean(freq_obs[r, tp_idx, :], axis=0)
        ok = (np.minimum(pbar, 1-pbar) >= 0.05)
        deltas = deltas[:, ok]
        n_int = deltas.shape[0]

        tcov = np.cov(deltas, bias=True)  # n_int × n_int

        offdiag_k2 = np.tril(tcov, -2) + np.triu(tcov, 2)
        total_cov = np.sum(offdiag_k2)
        total_var = np.sum(tcov)  # total variance including diagonal
        G = total_cov / total_var if total_var > 0 else np.nan
        G_per_rep.append(G)

    return float(np.mean(G_per_rep)), G_per_rep


def compute_cc(freq_obs, timepoints, mode="diagonal"):
    n_reps, n_tp, n_sites = freq_obs.shape

    if mode == "total":
        delta = freq_obs[:, -1, :] - freq_obs[:, 1, :]  # n_reps × n_sites
        num_vals, denom_vals = [], []
        for a, b in combinations(range(n_reps), 2):
            xa = delta[a] - delta[a].mean()  # center across sites
            xb = delta[b] - delta[b].mean()
            cov_ab = np.sum(xa * xb) / (n_sites - 1)
            var_a = np.sum(xa**2) / (n_sites - 1)
            var_b = np.sum(xb**2) / (n_sites - 1)
            num_vals.append(cov_ab)
            denom_vals.append(np.sqrt(max(var_a, 1e-12) * max(var_b, 1e-12)))
        return float(np.mean(num_vals)) / float(np.mean(denom_vals))

    elif mode == "diagonal":
        tp_idx = list(range(1, n_tp))  # skip founder
        deltas = np.diff(freq_obs[:, tp_idx, :], axis=1)  # n_reps × (n_int) × n_sites
        n_int = deltas.shape[1]

        cc_diag = []
        for ti in range(n_int):
            d = deltas[:, ti, :]  # n_reps × n_sites
            num_vals, denom_vals = [], []
            for a, b in combinations(range(n_reps), 2):
                xa = d[a] - d[a].mean()  # center across sites
                xb = d[b] - d[b].mean()
                cov_ab = np.sum(xa * xb) / (n_sites - 1)
                var_a = np.sum(xa**2) / (n_sites - 1)
                var_b = np.sum(xb**2) / (n_sites - 1)
                num_vals.append(cov_ab)
                denom_vals.append(np.sqrt(max(var_a, 1e-12) * max(var_b, 1e-12)))
            cc_ti = float(np.mean(num_vals)) / float(np.mean(denom_vals))
            cc_diag.append(cc_ti)
        return float(np.mean(cc_diag))


def build_s_matrix(scenario, n_sites, max_gen, rng, **kwargs):
    s_matrix = np.zeros((max_gen + 1, n_sites))

    if scenario == "drift":
        pass  # all zeros

    elif scenario == "sustained_parallel":
        s = kwargs.get("s", 0.05)
        frac = kwargs.get("frac", 0.05)
        n_sel = int(n_sites * frac)
        sel_idx = rng.choice(n_sites, n_sel, replace=False)
        signs = rng.choice([-1, 1], n_sel)
        for g in range(max_gen + 1):
            s_matrix[g, sel_idx] = s * signs

    elif scenario == "episodic_parallel":
        s = kwargs.get("s", 0.10)
        frac = kwargs.get("frac", 0.02)
        n_sel = int(n_sites * frac)
        interval_starts = TIMEPOINTS[:-1]
        interval_ends = TIMEPOINTS[1:]
        for gs, ge in zip(interval_starts, interval_ends):
            sel_idx = rng.choice(n_sites, n_sel, replace=False)
            signs = np.where(rng.random(n_sel) < 0.8, 1, -1)
            for g in range(gs, ge):
                s_matrix[g, sel_idx] = s * signs

    elif scenario == "sustained_weak":
        s = kwargs.get("s", 0.02)
        frac = kwargs.get("frac", 0.10)
        n_sel = int(n_sites * frac)
        sel_idx = rng.choice(n_sites, n_sel, replace=False)
        signs = rng.choice([-1, 1], n_sel)
        for g in range(max_gen + 1):
            s_matrix[g, sel_idx] = s * signs

    return s_matrix


def run_one_simulation(scenario, sim_id, depths, **sel_kwargs):
    rng = np.random.RandomState(SEED_BASE + sim_id * 1000 + hash(scenario) % 10000)
    max_gen = max(TIMEPOINTS)
    n_tp = len(TIMEPOINTS)

    s_matrix = build_s_matrix(scenario, N_SITES, max_gen, rng, **sel_kwargs)

    freq_obs_all = np.zeros((N_REPS, n_tp, N_SITES))
    for r in range(N_REPS):
        rep_rng = np.random.RandomState(rng.randint(0, 2**31))
        freq_true = simulate_wf(NE_TRUE, N_SITES, TIMEPOINTS, s_matrix, rep_rng)

        if np.isscalar(depths):
            d = depths
        else:
            d = np.array(depths)
        freq_obs_all[r] = observe_poolseq(freq_true, N_EFF, d, rep_rng)

    ne_raw, ne_bio = compute_F_regression(freq_obs_all, TIMEPOINTS, N_EFF, depths)
    ne_poolseq = compute_poolseq_ne(freq_obs_all, TIMEPOINTS, N_EFF, depths)
    G_k2, G_per_rep = compute_G_k2(freq_obs_all, TIMEPOINTS)
    cc_diag = compute_cc(freq_obs_all, TIMEPOINTS, mode="diagonal")
    cc_total = compute_cc(freq_obs_all, TIMEPOINTS, mode="total")

    return {
        "scenario": scenario,
        "sim_id": sim_id,
        "Ne_true": NE_TRUE,
        "Ne_F_raw": ne_raw,
        "Ne_F_corrected": ne_bio,
        "Ne_poolseq": ne_poolseq,
        "F_pool_ratio": ne_bio / ne_poolseq if (ne_poolseq > 0 and np.isfinite(ne_bio) and np.isfinite(ne_poolseq)) else np.nan,
        "G_k2": G_k2,
        "G_all_positive": all(g > 0 for g in G_per_rep),
        "cc_diagonal": cc_diag,
        "cc_total": cc_total,
    }


def main():
    results = []

    for i in range(N_SIMS):
        if i % 10 == 0: print(f"  sim {i}/{N_SIMS}", file=sys.stderr)
        results.append(run_one_simulation("drift", i, DEPTH_DEFAULT))

    for i in range(N_SIMS):
        if i % 10 == 0: print(f"  sim {i}/{N_SIMS}", file=sys.stderr)
        results.append(run_one_simulation("sustained_parallel", i, DEPTH_DEFAULT,
                                          s=0.05, frac=0.05))

    for i in range(N_SIMS):
        if i % 10 == 0: print(f"  sim {i}/{N_SIMS}", file=sys.stderr)
        results.append(run_one_simulation("episodic_parallel", i, DEPTH_DEFAULT,
                                          s=0.10, frac=0.02))

    for i in range(N_SIMS):
        if i % 10 == 0: print(f"  sim {i}/{N_SIMS}", file=sys.stderr)
        results.append(run_one_simulation("drift_T_depths", i, DEPTH_TRAJECTORIES["T"]))

    for i in range(N_SIMS):
        if i % 10 == 0: print(f"  sim {i}/{N_SIMS}", file=sys.stderr)
        results.append(run_one_simulation("drift_M_depths", i, DEPTH_TRAJECTORIES["M"]))

    for i in range(N_SIMS):
        if i % 10 == 0: print(f"  sim {i}/{N_SIMS}", file=sys.stderr)
        results.append(run_one_simulation("sustained_weak", i, DEPTH_DEFAULT,
                                          s=0.02, frac=0.10))

    df = pd.DataFrame(results)
    df.to_csv("variance_analysis/section1_rigorous/simulation_reconciliation.tsv",
              sep="\t", index=False)

    for scenario in df["scenario"].unique():
        sub = df[df["scenario"] == scenario]
        def med_iqr(col):
            vals = sub[col].dropna()
            vals = vals[np.isfinite(vals)]
            if len(vals) == 0:
                return "   n/a"
            return f"{np.median(vals):8.0f}" if col.startswith("Ne") else f"{np.median(vals):8.3f}"

        pct_g_pos = 100 * sub["G_all_positive"].mean()

if __name__ == "__main__":
    main()
