#!/usr/bin/env python3

import argparse
import glob
import os
import re
import sys
import numpy as np
import pandas as pd

try:
    from scipy import stats as sstats
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pcangsd-dir", required=True,
                   help="Dir containing founders_K{N}.cov and founders_K{N}.admix.{N}.Q")
    p.add_argument("--logs-dir",    required=True,
                   help="Dir containing pcangsd_founders_K{N}.log")
    p.add_argument("--metadata",    required=True, help="Founder metadata TSV")
    p.add_argument("--bamlist",     required=True, help="Founder bamlist")
    p.add_argument("--output",      required=True, help="QC report path (.txt)")
    p.add_argument("--strict",      action="store_true",
                   help="Fail exit code if any WARN message")
    return p.parse_args()


class Reporter:
    def __init__(self):
        self.lines = []
        self.fail = False
        self.warn = False

    def ok(self, msg):    self.lines.append(f"[ OK ]   {msg}")
    def info(self, msg):  self.lines.append(f"[info]   {msg}")
    def warn_(self, msg): self.lines.append(f"[WARN]   {msg}"); self.warn = True
    def fail_(self, msg): self.lines.append(f"[FAIL]   {msg}"); self.fail = True

    def header(self, s):
        self.lines.append("")
        self.lines.append(f"== {s} ==")


def count_lines(path):
    with open(path) as fh:
        return sum(1 for _ in fh if _.strip())


def list_k_from_dir(pcangsd_dir, prefix="founders_K"):
    ks = set()
    for p in glob.glob(os.path.join(pcangsd_dir, f"{prefix}*.cov")):
        m = re.match(rf".*{prefix}(\d+)\.cov$", p)
        if m:
            ks.add(int(m.group(1)))
    return sorted(ks)


def main():
    args = parse_args()
    rep = Reporter()


    rep.header("1. Sample counts")
    n_bam = count_lines(args.bamlist)
    meta  = pd.read_csv(args.metadata, sep="\t")
    n_meta = len(meta)
    if n_bam == n_meta:
        rep.ok(f"bamlist and metadata agree: n={n_bam}")
    else:
        rep.fail_(f"bamlist n={n_bam} != metadata n={n_meta}")


    rep.header("2. PCAngsd .Q matrices")
    ks = list_k_from_dir(args.pcangsd_dir)
    if not ks:
        rep.fail_(f"No founders_K*.cov files in {args.pcangsd_dir}")
    else:
        rep.info(f"K values found: {ks}")

    q_by_k = {}
    for k in ks:
        qpath = os.path.join(args.pcangsd_dir, f"founders_K{k}.admix.{k}.Q")
        covpath = os.path.join(args.pcangsd_dir, f"founders_K{k}.cov")
        if not os.path.exists(qpath):
            rep.fail_(f"K={k}: missing {qpath}")
            continue
        q = np.loadtxt(qpath)
        if q.ndim == 1:
            q = q[:, None]
        q_by_k[k] = q

        # Shape
        if q.shape != (n_meta, k):
            rep.fail_(f"K={k}: Q shape {q.shape}, expected ({n_meta},{k})")
        else:
            rep.ok(f"K={k}: Q shape correct ({q.shape[0]} x {q.shape[1]})")

        # NaN / negatives
        if np.isnan(q).any():
            rep.fail_(f"K={k}: Q contains NaN")
        if (q < -1e-6).any():
            rep.fail_(f"K={k}: Q contains negative values")

        # Row sums
        rs = q.sum(axis=1)
        if np.allclose(rs, 1.0, atol=1e-3):
            rep.ok(f"K={k}: all Q rows sum to 1.0")
        else:
            worst = np.abs(rs - 1.0).max()
            rep.fail_(f"K={k}: Q row sums deviate from 1 (max |err|={worst:.2e})")

        # Cov shape + PC variance
        if not os.path.exists(covpath):
            rep.warn_(f"K={k}: missing {covpath}")
        else:
            cov = np.loadtxt(covpath)
            if cov.shape != (n_meta, n_meta):
                rep.fail_(f"K={k}: cov shape {cov.shape}, expected ({n_meta},{n_meta})")
            else:
                eigvals = np.linalg.eigvalsh(cov)
                eigvals = np.sort(np.abs(eigvals))[::-1]
                pct = 100 * eigvals / eigvals.sum()
                rep.info(f"K={k}: PC1={pct[0]:.2f}%  PC2={pct[1]:.2f}%  PC3={pct[2]:.2f}%")
                if pct[0] <= 0 or not np.isfinite(pct[0]):
                    rep.fail_(f"K={k}: PC1 variance is non-positive/NaN")


    rep.header("3. Convergence")
    for k in ks:
        logp = os.path.join(args.logs_dir, f"pcangsd_founders_K{k}.log")
        if not os.path.exists(logp):
            rep.warn_(f"K={k}: log missing ({logp})")
            continue
        with open(logp) as fh:
            text = fh.read()
        if "did not converge" in text.lower():
            rep.warn_(f"K={k}: log reports non-convergence "
                      f"(pipeline still produced Q, but treat with caution)")
        else:
            rep.ok(f"K={k}: no convergence warning in log")


    rep.header("4. Component stability across K")
    prior_k = None
    for k in sorted(q_by_k):
        if prior_k is None:
            prior_k = k
            continue
        q_a = q_by_k[prior_k]
        q_b = q_by_k[k]

        best = []
        for ca in range(q_a.shape[1]):
            corrs = [np.corrcoef(q_a[:, ca], q_b[:, cb])[0, 1]
                     for cb in range(q_b.shape[1])]
            best.append(max(corrs))
        rep.info(f"K={prior_k}→K={k}: per-component max correlation = "
                 f"{[round(x, 3) for x in best]}")
        prior_k = k


    rep.header("5. B-vs-T homogeneity test at K=2")
    if 2 not in q_by_k:
        rep.warn_("K=2 not run; skipping B-vs-T test")
    elif not HAVE_SCIPY:
        rep.warn_("scipy not available; skipping t-test")
    else:
        q2 = q_by_k[2]
        b_mask = (meta["host_plant"] == "B").values
        t_mask = (meta["host_plant"] == "T").values
        if b_mask.sum() < 2 or t_mask.sum() < 2:
            rep.warn_(f"Not enough B ({b_mask.sum()}) / T ({t_mask.sum()}) samples")
        else:
            q1_b = q2[b_mask, 0]
            q1_t = q2[t_mask, 0]
            tstat, pval = sstats.ttest_ind(q1_b, q1_t, equal_var=False)
            rep.info(f"B: mean Q1={q1_b.mean():.3f} (n={b_mask.sum()}); "
                     f"T: mean Q1={q1_t.mean():.3f} (n={t_mask.sum()})")
            rep.info(f"Welch t = {tstat:.3f}, p = {pval:.3g}")
            if pval < 0.01:
                rep.warn_("B and T differ significantly at K=2 — founder "
                          "colony may retain host-source structure "
                          "(user's prior was: no difference expected)")
            else:
                rep.ok("B and T not significantly different at K=2 "
                       "(consistent with homogeneous admixed colony)")





    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as fh:
        status = "FAIL" if rep.fail else ("WARN" if rep.warn else "PASS")
        fh.write(f"popstructure QC report: {status}\n")
        fh.write("=" * 50 + "\n")
        fh.write("\n".join(rep.lines) + "\n")

    if rep.fail:
        sys.exit(2)
    elif rep.warn:
        sys.exit(1 if args.strict else 0)
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()
