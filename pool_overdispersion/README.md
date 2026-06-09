# pool_overdispersion — pool technical-replicate variance

Estimates pool-sequencing overdispersion (φ) and effective sample size from G10 **technical-replicate
pool pairs** and the founder pools, used in the variance-decomposition / Ne sections (Supp. Fig S7).

**This is a [grenepipe](https://github.com/moiexpositoalonsolab/grenepipe) overlay** — only the
study-specific files are here. Install grenepipe and add these files in (see [INSTALL.md](../INSTALL.md)).

## Entry points
- `Snakefile_pool_comparison` (+ `config_pool_comparison.yml`) — compute per-site AF across technical-replicate pools
- `make_qc_table.py` — per-pool sequencing/mapping QC summary (Supp. Table S3)
- `scripts/` — `compute_overdispersion.py`, `extract_af_from_pileup.py`, `merge_pool_af.py`, `plot_pool_variation.py`
- `submit_pool_comparison.sh` — SLURM submission (set `${SLURM_*}` first)
- `profiles/slurm/`, `envs/` — Snakemake SLURM profile + conda specs

See [../SCRIPTS.md](../SCRIPTS.md) for a one-line description of every script.
