#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_baypass_wild import read_vcf_positions, parse_mpileup  

SAMPLES_22 = [
    ('AVB',   'wild',    50), ('AVT',   'wild',    50),
    ('PSB',   'wild',    50), ('PST',   'wild',    50),
    ('RMB',   'wild',    50), ('RMT',   'wild',    50),
    ('F1G00', 'founder', 80), ('F2G00', 'founder', 80),
    ('F3G00', 'founder', 80), ('F4G00', 'founder', 80),
    ('B1G10', 'g10',    160), ('B2G10', 'g10',    160),
    ('B3G10', 'g10',    160), ('B4G10', 'g10',    160),
    ('T1G10', 'g10',    160), ('T2G10', 'g10',    160),
    ('T3G10', 'g10',    160), ('T4G10', 'g10',    160),
    ('M1G10', 'g10',    160), ('M2G10', 'g10',    160),
    ('M3G10', 'g10',    160), ('M4G10', 'g10',    160),
]


def group_of(sample: str) -> str:
    if sample in ('AVB', 'PSB', 'RMB'):
        return 'Wild B'
    if sample in ('AVT', 'PST', 'RMT'):
        return 'Wild T'
    if sample.startswith('F'):
        return 'Founders'
    if sample.startswith('B') and 'G10' in sample:
        return 'G10 B'
    if sample.startswith('T') and 'G10' in sample:
        return 'G10 T'
    if sample.startswith('M') and 'G10' in sample:
        return 'G10 M'
    raise ValueError(f'Unknown sample: {sample}')



GROUP_COLOR = {
    'Wild B':   '#D55E00',
    'Wild T':   '#0072B2',
    'Founders': '#7f7f7f',
    'G10 B':    '#D55E00',
    'G10 T':    '#0072B2',
    'G10 M':    '#009E73',
}

GROUP_MARKER = {
    'Wild B':   'o', 'Wild T':   'o',
    'Founders': 's',
    'G10 B':    '^', 'G10 T':    '^', 'G10 M':    'D',
}


def build_af_matrix(vcf_path: str, dir_map: dict, samples: list,
                    min_cov: int, max_cov: int) -> pd.DataFrame:

    positions = read_vcf_positions(vcf_path)
    ordered = sorted(positions.keys())
    n = len(ordered)


    af = np.full((n, len(samples)), np.nan, dtype=np.float32)
    for j, (name, dirkey, _psz) in enumerate(samples):
        mp = Path(dir_map[dirkey]) / f'{name}.mpileup'
        if not mp.exists():
            sys.exit(f'error')
        counts = parse_mpileup(str(mp), positions)
        for i, key in enumerate(ordered):
            r, a, d = counts.get(key, (0, 0, 0))
            if min_cov <= d <= max_cov and (r + a) > 0:
                af[i, j] = a / (r + a)
        n_ok = int(np.sum(np.isfinite(af[:, j])))


    chrpos = [f'{c}:{p}' for c, p in ordered]
    out = pd.DataFrame(af, columns=[s[0] for s in samples])
    out.insert(0, 'chrom_pos', chrpos)
    return out


def filter_to_sites(af_df: pd.DataFrame, sites_file: str) -> pd.DataFrame:
    raw = pd.read_csv(sites_file, sep=None, engine='python')
    if 'chrom' in raw.columns and 'pos' in raw.columns:
        keep = set(f'{c}:{int(p)}' for c, p in zip(raw['chrom'], raw['pos']))
    elif set(raw.columns[:3]) >= {'chrom', 'start', 'end'} or len(raw.columns) >= 3:
        c, _, e = raw.columns[0], raw.columns[1], raw.columns[2]
        keep = set(f'{cc}:{int(pp)}' for cc, pp in zip(raw[c], raw[e]))
    else:
        sys.exit(f'error')
    before = len(af_df)
    af_df = af_df[af_df['chrom_pos'].isin(keep)].reset_index(drop=True)
    return af_df


def run_pca(af_df: pd.DataFrame, samples: list) -> tuple:
    names = [s[0] for s in samples]
    mat = af_df[names].to_numpy(dtype=np.float32)

    keep = np.all(np.isfinite(mat), axis=1)
    mat = mat[keep]

    p = mat.mean(axis=1)
    var_mask = (p > 0.01) & (p < 0.99)
    mat = mat[var_mask]
    p = p[var_mask]

    scale = np.sqrt(p * (1 - p) + 1e-12)
    Z = (mat - p[:, None]) / scale[:, None]

    pca = PCA(n_components=min(10, Z.shape[1]))
    coords = pca.fit_transform(Z.T)
    pv = pca.explained_variance_ratio_ * 100

    return coords, pv, mat.shape[0]


def plot_pca(coords: np.ndarray, pv: np.ndarray, samples: list,
             outdir: Path, n_sites: int, suffix: str = ''):
    names = [s[0] for s in samples]
    groups = [group_of(n) for n in names]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), dpi=200)
    for ax, (xi, yi) in zip(axes, [(0, 1), (0, 2)]):
        for i, (name, g) in enumerate(zip(names, groups)):
            ax.scatter(coords[i, xi], coords[i, yi],
                       c=GROUP_COLOR[g], marker=GROUP_MARKER[g],
                       s=130, edgecolors='k', linewidths=0.5, zorder=3)
            ax.annotate('  ' + name, (coords[i, xi], coords[i, yi]),
                        fontsize=7, va='center', zorder=4)
        ax.set_xlabel(f'PC{xi+1} ({pv[xi]:.1f}%)')
        ax.set_ylabel(f'PC{yi+1} ({pv[yi]:.1f}%)')
        ax.axhline(0, color='k', lw=0.3, alpha=0.3)
        ax.axvline(0, color='k', lw=0.3, alpha=0.3)

    handles = [Line2D([0], [0], marker=GROUP_MARKER[g], color='w',
                      markerfacecolor=GROUP_COLOR[g], markeredgecolor='k',
                      markersize=10, label=g)
               for g in ['Wild B', 'Wild T', 'Founders', 'G10 B', 'G10 T', 'G10 M']]
    axes[1].legend(handles=handles, fontsize=8, loc='best')
    title = f'AF-based PCA of 22 pools — {n_sites:,} SNPs, Patterson-scaled'
    if suffix:
        title += f'\n(subset: {suffix.lstrip("_")})'
    fig.suptitle(title, fontsize=10, y=1.02)
    fig.tight_layout()

    for ext in ('png', 'svg'):
        fig.savefig(outdir / f'af_pca_22pools{suffix}.{ext}', bbox_inches='tight')
    print(f'  wrote {outdir/f"af_pca_22pools{suffix}.png"} / .svg')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--wild-dir',    default='baypass_wild/pileups')
    ap.add_argument('--founder-dir', default='grenfst/diversity_combined/pileups')
    ap.add_argument('--g10-dir',     default='grenfst/diversity_attrition/pileups')
    ap.add_argument('--outdir', default='final_plots/wild')
    ap.add_argument('--min-cov', type=int, default=10)
    ap.add_argument('--max-cov', type=int, default=500)
    ap.add_argument('--af-cache', default=None,
                    help='Path to cached AF matrix CSV. If exists, re-use; else write there after parsing.')
    ap.add_argument('--force-rebuild', action='store_true')
    ap.add_argument('--sites-file', default=None,
                    help='Optional TSV/CSV with chrom+pos columns. If given, the PCA is restricted to these sites.')
    ap.add_argument('--suffix', default='',
                    help='Suffix appended to output filenames (e.g. "_baypass_top1pct").')
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cache = Path(args.af_cache) if args.af_cache else (outdir / 'af_matrix_22pools.csv')

    dir_map = {'wild': args.wild_dir, 'founder': args.founder_dir, 'g10': args.g10_dir}

    if cache.exists() and not args.force_rebuild:
        print(f'Loading cached AF matrix: {cache}')
        af_df = pd.read_csv(cache)
    else:
        af_df = build_af_matrix(args.vcf, dir_map, SAMPLES_22,
                                args.min_cov, args.max_cov)
        print(f'Caching AF matrix → {cache}')
        af_df.to_csv(cache, index=False)

    if args.sites_file:
        af_df = filter_to_sites(af_df, args.sites_file)

    coords, pv, n_sites = run_pca(af_df, SAMPLES_22)

    suffix = args.suffix
    scores = pd.DataFrame(
        coords, columns=[f'PC{i+1}' for i in range(coords.shape[1])]
    )
    scores.insert(0, 'sample', [s[0] for s in SAMPLES_22])
    scores.insert(1, 'group', [group_of(s[0]) for s in SAMPLES_22])
    scores.to_csv(outdir / f'af_pca_22pools{suffix}_scores.tsv', sep='\t',
                  index=False, float_format='%.4f')

    pd.DataFrame({'pc': [f'PC{i+1}' for i in range(len(pv))],
                  'variance_explained_pct': pv}).to_csv(
        outdir / f'af_pca_22pools{suffix}_variance.tsv', sep='\t',
        index=False, float_format='%.3f')

    plot_pca(coords, pv, SAMPLES_22, outdir, n_sites, suffix=suffix)
    print(f'Done. Scores: {outdir/f"af_pca_22pools{suffix}_scores.tsv"}')


if __name__ == '__main__':
    main()
