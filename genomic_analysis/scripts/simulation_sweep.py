#!/usr/bin/env python3

import sys
import os
import numpy as np
import pandas as pd
from scipy import stats as sp_stats
from itertools import combinations, product

NE_TRUE = 250
N_REPS = 4
TIMEPOINTS = [0, 1, 2, 6, 7, 8, 9]
N_EFF = 29
DEPTH = 83
N_SITES = 50_000
N_SIMS = 20
SEED_BASE = 31415


def simulate_wf_clustered(ne, n_sites, timepoints, n_targets, ld_size, s,
                           rng, sustained=True):

    max_gen = max(timepoints)
    n_tp = len(timepoints)
    block_size = 1 + ld_size if ld_size > 0 else 1


    if n_targets == 0:
        target_positions = np.array([], dtype=int)
    elif ld_size == 0:
        target_positions = rng.choice(n_sites, min(n_targets, n_sites), replace=False)
    else:
        spacing = n_sites // max(n_targets, 1)
        target_positions = []
        pos = spacing // 2
        for _ in range(n_targets):
            if pos + block_size <= n_sites:
                target_positions.append(pos)
            pos += spacing
        target_positions = np.array(target_positions[:n_targets], dtype=int)

    n_targets_actual = len(target_positions)

    target_signs = np.zeros(n_sites)
    for tp in target_positions:
        target_signs[tp] = 1 if rng.random() < 0.5 else -1

    block_leader = np.full(n_sites, -1, dtype=int)  # -1 = independent
    if ld_size > 0:
        for tp in target_positions:
            end = min(tp + block_size, n_sites)
            block_leader[tp:end] = tp

    p0 = rng.uniform(0.1, 0.9, n_sites)
    if ld_size > 0:
        for tp in target_positions:
            end = min(tp + block_size, n_sites)
            p0[tp:end] = p0[tp]

    current_p = p0.copy()
    freq_true = np.zeros((n_tp, n_sites))
    tp_idx = 0

    for g in range(max_gen + 1):
        if tp_idx < n_tp and g == timepoints[tp_idx]:
            freq_true[tp_idx, :] = current_p.copy()
            tp_idx += 1

        if g < max_gen:
            n_hap = 2 * ne

            if ld_size == 0:
                for tp in target_positions:
                    delta = s * target_signs[tp] * current_p[tp] * (1 - current_p[tp])
                    current_p[tp] = np.clip(current_p[tp] + delta, 1e-4, 1-1e-4)
                counts = rng.binomial(n_hap, np.clip(current_p, 1e-4, 1-1e-4))
                current_p = counts / n_hap
            else:
                processed = np.zeros(n_sites, dtype=bool)

                for tp in target_positions:
                    end = min(tp + block_size, n_sites)
                    delta = s * target_signs[tp] * current_p[tp] * (1 - current_p[tp])
                    current_p[tp] = np.clip(current_p[tp] + delta, 1e-4, 1-1e-4)
                    new_count = rng.binomial(n_hap, np.clip(current_p[tp], 1e-4, 1-1e-4))
                    new_p = new_count / n_hap
                    current_p[tp:end] = np.clip(new_p, 1e-4, 1-1e-4)
                    processed[tp:end] = True

                indep = ~processed
                if indep.any():
                    counts = rng.binomial(n_hap, np.clip(current_p[indep], 1e-4, 1-1e-4))
                    current_p[indep] = counts / n_hap

            current_p = np.clip(current_p, 1e-4, 1-1e-4)

    return freq_true


def observe_poolseq(freq_true, n_eff, depth, rng):
    n_tp, n_sites = freq_true.shape
    freq_obs = np.zeros_like(freq_true)
    for t in range(n_tp):
        pool = rng.binomial(2 * n_eff, np.clip(freq_true[t], 1e-6, 1-1e-6))
        p_pool = pool / (2 * n_eff)
        seq = rng.binomial(depth, p_pool)
        freq_obs[t] = seq / depth
    return freq_obs



def compute_stats(freq_obs_all, timepoints, n_eff, depth):
    n_reps, n_tp, n_sites = freq_obs_all.shape
    gens, F_bio_list = [], []
    for ti in range(n_tp):
        t = timepoints[ti]
        if t == 0: continue
        p_reps = freq_obs_all[:, ti, :]
        pbar = p_reps.mean(axis=0)
        ok = np.minimum(pbar, 1-pbar) >= 0.05
        if ok.sum() < 100: continue
        ss = np.sum((p_reps[:, ok] - pbar[ok][None, :])**2, axis=0) / (n_reps-1)
        het = pbar[ok] * (1-pbar[ok]); valid = het > 1e-10
        F = np.mean(ss[valid]) / np.mean(het[valid])
        F_bio_list.append(F - 1/(2*n_eff) - 1/depth)
        gens.append(t)
    if len(gens) >= 3:
        sb, *_ = sp_stats.linregress(gens, F_bio_list)
        ne_f = 1/(2*sb) if sb > 0 else np.inf
    else:
        ne_f = np.nan

    t0_idx, t1_idx = 1, n_tp - 1
    t_diff = timepoints[t1_idx] - timepoints[t0_idx]
    ne_pool_reps = []
    for r in range(n_reps):
        p0 = freq_obs_all[r, t0_idx, :]; pt = freq_obs_all[r, t1_idx, :]
        pbar = (p0+pt)/2
        ok = (np.minimum(pbar, 1-pbar) >= 0.05) & np.isfinite(p0) & np.isfinite(pt)
        if ok.sum() < 1000: continue
        het = pbar[ok]*(1-pbar[ok]); valid = het > 1e-10
        Fc = (p0[ok][valid]-pt[ok][valid])**2 / het[valid]
        inv_n = 1/(2*n_eff) + 1/depth - 1/(2*n_eff*depth)
        Fc_c = Fc - 2*inv_n
        n_v = len(Fc_c); wnd = 1000; n_wnd = n_v // wnd
        if n_wnd < 5: continue
        wnd_ne = []
        for w in range(n_wnd):
            phi = np.mean(Fc_c[w*wnd:(w+1)*wnd])
            if phi > 0:
                ne_w = t_diff/(2*phi)
                if 0 < ne_w < 1e7: wnd_ne.append(ne_w)
        if wnd_ne:
            ne_pool_reps.append(np.median(wnd_ne))
    ne_pool = np.mean(ne_pool_reps) if ne_pool_reps else np.nan

    G_per_rep = []
    for r in range(n_reps):
        tp_idx = list(range(1, n_tp))
        deltas = np.diff(freq_obs_all[r, tp_idx, :], axis=0)
        pbar_r = np.mean(freq_obs_all[r, tp_idx, :], axis=0)
        ok = np.minimum(pbar_r, 1-pbar_r) >= 0.05
        d = deltas[:, ok]
        if d.shape[1] < 100: continue
        tcov = np.cov(d, bias=True)
        offdiag_k2 = np.tril(tcov, -2) + np.triu(tcov, 2)
        tv = np.sum(tcov)
        G_per_rep.append(np.sum(offdiag_k2) / tv if tv > 0 else np.nan)
    G_k2 = float(np.mean(G_per_rep)) if G_per_rep else np.nan
    G_all_pos = all(g > 0 for g in G_per_rep) if G_per_rep else False

    tp_idx = list(range(1, n_tp))
    deltas_all = np.diff(freq_obs_all[:, tp_idx, :], axis=1)
    n_int = deltas_all.shape[1]
    cc_diag_vals = []
    for ti in range(n_int):
        d = deltas_all[:, ti, :]
        num_v, den_v = [], []
        for a, b in combinations(range(n_reps), 2):
            xa = d[a] - d[a].mean(); xb = d[b] - d[b].mean()
            num_v.append(np.sum(xa*xb)/(n_sites-1))
            va = np.sum(xa**2)/(n_sites-1); vb = np.sum(xb**2)/(n_sites-1)
            den_v.append(np.sqrt(max(va,1e-12)*max(vb,1e-12)))
        cc_diag_vals.append(np.mean(num_v) / np.mean(den_v))
    cc_diag = float(np.mean(cc_diag_vals))

    dt = freq_obs_all[:, -1, :] - freq_obs_all[:, 1, :]
    num_v, den_v = [], []
    for a, b in combinations(range(n_reps), 2):
        xa = dt[a]-dt[a].mean(); xb = dt[b]-dt[b].mean()
        num_v.append(np.sum(xa*xb)/(n_sites-1))
        va = np.sum(xa**2)/(n_sites-1); vb = np.sum(xb**2)/(n_sites-1)
        den_v.append(np.sqrt(max(va,1e-12)*max(vb,1e-12)))
    cc_total = float(np.mean(num_v)/np.mean(den_v))

    ratio = ne_f / ne_pool if (ne_pool > 0 and np.isfinite(ne_f) and np.isfinite(ne_pool)) else np.nan

    return dict(Ne_F=ne_f, Ne_pool=ne_pool, F_pool_ratio=ratio,
                G_k2=G_k2, G_all_pos=G_all_pos, cc_diag=cc_diag, cc_total=cc_total)


def run_one(s, frac, ld_size, sim_id):
    rng = np.random.RandomState(SEED_BASE + sim_id * 997 +
                                 int(s*100) * 31 + int(frac*100) * 7 + ld_size)
    n_targets = int(N_SITES * frac)
    freq_obs_all = np.zeros((N_REPS, len(TIMEPOINTS), N_SITES))
    for r in range(N_REPS):
        rr = np.random.RandomState(rng.randint(0, 2**31))
        ft = simulate_wf_clustered(NE_TRUE, N_SITES, TIMEPOINTS,
                                    n_targets, ld_size, s, rr)
        freq_obs_all[r] = observe_poolseq(ft, N_EFF, DEPTH, rr)
    return compute_stats(freq_obs_all, TIMEPOINTS, N_EFF, DEPTH)


def main():
    S_VALS = [0.0, 0.05, 0.10, 0.20, 0.30]
    FRAC_VALS = [0.01, 0.05, 0.10, 0.20]
    LD_VALS = [0, 10, 50, 100]

    results = []
    total = len(S_VALS) * len(FRAC_VALS) * len(LD_VALS)
    count = 0

    for s_val in S_VALS:
        for frac_val in FRAC_VALS:
            for ld_val in LD_VALS:
                count += 1
                if s_val == 0 and (frac_val != FRAC_VALS[0] or ld_val != LD_VALS[0]):
                    continue
                label = f"s={s_val:.2f}_f={frac_val:.2f}_ld={ld_val}"
   
                sim_stats = []
                for i in range(N_SIMS):
                    st = run_one(s_val, frac_val, ld_val, i)
                    sim_stats.append(st)

                for stat_key in ["Ne_F", "Ne_pool", "F_pool_ratio",
                                 "G_k2", "cc_diag", "cc_total"]:
                    vals = [d[stat_key] for d in sim_stats
                            if np.isfinite(d[stat_key])]
                    med = np.median(vals) if vals else np.nan
                    results.append(dict(
                        s=s_val, frac=frac_val, ld_size=ld_val,
                        statistic=stat_key, median=med,
                        n_sims=len(vals),
                    ))

                gpos = [d["G_all_pos"] for d in sim_stats]
                results.append(dict(
                    s=s_val, frac=frac_val, ld_size=ld_val,
                    statistic="pct_G_allpos", median=100*np.mean(gpos),
                    n_sims=N_SIMS,
                ))

    df = pd.DataFrame(results)
    outpath = os.path.join("variance_analysis/section1_rigorous",
                           "simulation_sweep.tsv")
    df.to_csv(outpath, sep="\t", index=False)

    for stat in ["Ne_F", "Ne_pool", "G_k2", "cc_diag", "cc_total"]:
        sub = df[df["statistic"] == stat]
        for ld in LD_VALS:
            ld_sub = sub[sub["ld_size"] == ld]
            if ld_sub.empty: continue
            piv = ld_sub.pivot(index="s", columns="frac", values="median")


if __name__ == "__main__":
    main()
