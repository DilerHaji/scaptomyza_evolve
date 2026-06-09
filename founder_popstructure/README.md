# founder_popstructure — individual founder sequencing

Population-structure and overdispersion analyses from low-coverage **individual** sequencing of the
128 founder flies (vs. the pooled time series). Tests for residual *B*- vs *T*-host structure in the
starting colony (Supp. Fig S4) and estimates pool-sequencing overdispersion (Supp. Fig S7).

**This is a [grenepipe](https://github.com/moiexpositoalonsolab/grenepipe) overlay** — only the
study-specific files are here. Install grenepipe and add these files in (see [INSTALL.md](../INSTALL.md)).

## Entry points
- `Snakefile_popstructure` (+ `config/popstructure.yaml`) — ANGSD genotype likelihoods → PCAngsd PCA + admixture (K=2..8)
- `Snakefile_founders_analyses` — downstream founder analyses
- `scripts/` — `plot_pcangsd.py`, `plot_popstructure_pca.py`, `plot_admix_ksweep.py`, `summarize_admix_by_source.py`, `estimate_overdispersion.py`, `plot_af_comparison.py`
- `submit_popstructure.sh`, `submit_founders_analyses.sh` — SLURM submission (set `${SLURM_*}` first)
- `envs/` — conda specs (angsd, pcangsd, bcftools, plotting)

See [../SCRIPTS.md](../SCRIPTS.md) for a one-line description of every script.
