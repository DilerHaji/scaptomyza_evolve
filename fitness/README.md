# fitness — performance, viability, heritability & selection theory

R analyses of the reciprocal performance assays, population
viability across treatments (Fig 2D, Fig 3A), realized heritability of host specialization
(modified breeder's equation), and fitness-set theory connecting the empirical
results to the multiple-niche model.

Run the scripts from inside the `fitness/` directory (the `*.csv` inputs are referenced by
relative path). Raw phenotype data are included as the small CSVs in this directory.

- `sflava_fitness.R`, `sflava_fitness_svg.R` — main driver / figure export
- `response_to_selection_final.R` — realized heritability + viability (Fig 2D / 3A)
- `Levins.R`, `levene_test.R`, `levene_test_grid.R`, `levene_dempster_comparison.R`,
  `dempster_test_grid.R`, `dominance_reversal.R`, `reach_equilibrium.R` — fitness-set / equilibrium theory
- `fig_specialization_with_absolute.R` — specialization figure
- `scripts/` — shared helpers (`data_setup.R`, `contrasts.R`, `colors.R`, `dev_time.R`, `helper.R`, `viz_data.R`)
- `*.csv` — raw experimental data (fly counts, selection-experiment generations)
