# genomic_analysis — main pool-seq pipeline

The core Snakemake pipeline and analysis/figure scripts for the *S. flava* experimental-evolution
pool-seq time series. See the top-level [README](../README.md) for the full figure→code map and
configuration notes, and [../SCRIPTS.md](../SCRIPTS.md) for a one-line description of every script.

## Entry points

| Workflow | File | Produces |
|----------|------|----------|
| Mapping → SNP calling → AF tables | `Snakefile` (+ `config.yml`) | sync / allele-frequency tables |
| Genetic diversity | `Snakefile_diversity` | grenedalf π, Watterson θ, Tajima's D (windows) |
| Founder LD | `Snakefile_ngsld` (+ `submit_ngsld.sh`) | ngsLD r² blocks (scaffold 439) |
| Haplotype blocks | `Snakefile_haplovalidate` (+ `scripts/run_haplovalidate.sh`) | temporal haplotype blocks |

Selection scans, the 5-test vote scan (Fig 3), and the lab-to-wild convergence analysis (Fig 4)
are driven by the `scripts/section1_*.py`, `scripts/section2_*.py`, and `scripts/section3_*.py`
families plus `rules/*.smk`. Run the per-analysis `run_*.sh` helpers after editing the
`PATH/TO/...` and `${SLURM_*}` placeholders.

## Layout

- `rules/` — Snakemake rule files (`*.smk`)
- `scripts/` — Python/R analysis + figure scripts
- `maps/` — small sample/locus lookup tables
- `config.yml` — all parameters (paths are placeholders — edit before use)
- `environments/`, `envs/` — conda environment specs

> Note: `config.yml` retains commented rule blocks and parameters for analyses explored during the
> study (e.g. freebayes, simulation sweeps) that are **not** part of the final paper; the active
> includes in each `Snakefile` define what was used.
