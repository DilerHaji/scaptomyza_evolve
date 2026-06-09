This repository contains the genomic, population-genetic, and fitness-analysis code for an
evolve-and-resequence (E&R) experiment in the leaf-mining drosophilid *Scaptomyza flava*, an
endoparasite of mustard plants (Brassicales). Replicate laboratory populations were evolved for
ten generations on two host plants –– *Barbarea vulgaris* (**B**) and *Turritis glabra* (**T**) —
and on a combination of both (**B+T**) and then resequenced as pools over the time series. The
analyses test whether multiple-niche (Levene Model) selection maintains host-associated genomic polymorphism, characterize the genomic targets of antagonistic selection, and compare the laboratory signatures to wild host-associated differentiation.

This is analysis source code was developed to run on an HPC cluster (Berkeley Savio) using Snakemake. Absolute cluster/local paths and cluster-specific SLURM settings are replaced with placeholders. You must set these for your environment before running anything, see
Configuration. Third-party tools are documented in INSTALL.md.

Repository Organization
```
sflava-expevol-release/
├── genomic_analysis/        Main pool-seq Snakemake pipeline (the bulk of the analyses)
│   ├── Snakefile            Mapping, PoolSNP calling, allele-frequency tables
│   ├── Snakefile_diversity  grenedalf diversity (π, Watterson θ, Tajima's D)
│   ├── Snakefile_ngsld      Founder linkage-disequilibrium (ngsLD)
│   ├── Snakefile_haplovalidate  Temporal haplotype-block detection
│   ├── config.yml           Pipeline parameters (paths are placeholders)
│   ├── rules/               Snakemake rule files (*.smk)
│   ├── scripts/             Python/R analysis & figure scripts
│   ├── maps/                Analysis-specific configurations
│   └── environments/, envs/ Environments for Snakemake
├── founder_popstructure/    Individual founder sequencing
│                            PCAngsd population structure and admixture
├── pool_overdispersion/     Poolseq overdispersion and variance decomposition
└── fitness/                 Reciprocal performance, viability, heritability, theory
```

#### Figure to code map
| Figure | What it shows                                                              | Key script(s)                                                                                                                                                                                                                                                                                             |
| :----- | :------------------------------------------------------------------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **1**  | Experimental design, host-plant occurrence map                             | occurrence map drawn from iNaturalist records                                                                                                                                                                                                                                 |
| **2A** | Allele-frequency-trajectory PCA (uniform B/T axis vs subdivided B+T)       | build_af_pca_bt.py,plot_af_pca_3d.py, plot_af_pca_trajectory.R,test_subdivided_vs_uniform.py                                                                                                                                                                                                              |
| **2B** | Per-pool diversity (π, θ_W, Tajima's D) across wild/founder/G10            | process_diversity.py, plot_s3_main_diversity_panels.py (grenedalf via Snakefile_diversity)                                                                                                                                                                                                                |
| **2C** | Idiosyncratic rare-variant retention (scaffold 110)                        | section2_wild_qstar_diversity_filter_v3.py, fig_section1_all.py                                                                                                                                                                                                                                           |
| **2D** | Population viability across treatments                                     | fitness/response_to_selection_final.R                                                                                                                                                                                                                                                                     |
| **3A** | Reciprocal performance / realized heritability                             | fitness/response_to_selection_final.R, fitness/Levins.R                                                                                                                                                                                                                                                   |
| **3B** | 5-test genome scan for candidate regions (scaffold 439); chr439 gene track | section2_candidate_identification_v2.py, section2_vote_sweep.py, section2_per_test_perm_null.py, section2_fig3b_chr439_clean_panel.py, glmm.r, fst_slope_peaks.py/fst_lmm.py, rules/baypass_wild.smk, extract_chr439_region_proteins.py, parse_chr439_dmel_annot.py, section2_chr439_snp_feature_annot.py |
| **3C** | Polarized AF trajectory and diversity retention at the signal block        | section2_fig3c_diversity_strengthening.py, section2_fig3c_pointcloud_diversity.py, run_af_trajectory_*.py                                                                                                                                                                                                 |
| **3E** | Equilibrium-frequency (q*) landscapes and selection intensity              | section2_fig3e_* , recompute_qstar_binomial_glm.py                                                                                                                                                                                                                                                        |
| **4**  | Lab-to-wild  concordance at scaffold 439 (D. mel X homolog)                | section3_main_panel_clean.py, section3_fig4_three_panels_v2.py, section3_fig4_concordance_delta_dp.py, section3_x_enrichment_analysis.py, section3_muller_validation_figure.py                                                                                                                            |Supplementary analyses (representative scripts):


### Supplementary
| Supp.   | Analysis                                                                      | Key script(s)                                                                                                                                                                                                          |
| :------ | :---------------------------------------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| S5–S6   | Per-pool & per-replicate diversity trajectories                               | plot_s3_fig1_wild_diversity.py, plot_s3_trajectory_diversity.py, plot_s3_fig1B_thetaW_trajectory_G1_G10.py, s3_tables_summary_and_wilcoxon.py                                                                          |
| S7      | Pool-sequencing overdispersion                                                | founder_popstructure/scripts/estimate_overdispersion.py, fig_S7_technical_noise.py, pool_overdispersion/                                                                                                               |
| S8      | Between-replicate F-regression Ne (noise-corrected)                           | f_regression_noise_corrected.py, fig_f_regression_noise_correction.py                                                                                                                                                  |
| S9–S10  | Within-replicate poolSeq Ne; method-of-moments Ne; 3-estimator reconciliation | section1_rigorous_analysis.py, reconcile_ne_methods.py, fig_SYY_ne_summary.py, run_lynch_s.sh                                                                                                                          |
| S11–S16 | Temporal-replicate covariance, G statistic, convergence correlation, power    | run_cvtkpy_final.py, fig_temporal_cov_matrices.py, fig_temporal_cov_by_rep.py, make_supp_figures.py, make_supp_svgs_v2.py, cvtkpy_power_analysis.py, cvtkpy_power_grid.py, simulation_sweep.py/fig_simulation_sweep.py |
| S17–S18 | Haplotype-block detection (haplovalidate) + CMH (ACER)                        | Snakefile_haplovalidate, run_haplovalidate.sh, rules/cmh.smk                                                                                                                                                           |
| S19–S20 | Founder linkage disequilibrium (ngsLD) and scaffold-439 LD block              | Snakefile_ngsld, submit_ngsld.sh                                                                                                                                                                                       |
| —       | Founder population structure (ANGSD + PCAngsd)                                | founder_popstructure/Snakefile_popstructure, plot_popstructure_pca.py, plot_admix_ksweep.py                                                                                                                            |
| —       | Wild host-associated differentiation (BayPass)                                | rules/baypass_wild.smk, run_baypass_wild.sh                                                                                                                                                                            |


### Third-party tools
| Tool                             | Upstream                                            | License  | Used for                                 |
| :------------------------------- | :-------------------------------------------------- | :------- | :--------------------------------------- |
| grenedalf                        | github.com/lczech/grenedalf                         | GPL-3    | diversity (π/θ_W/D) & FST                |
| ANGSD / PCAngsd                  | github.com/ANGSD/angsd, github.com/Rosemeis/pcangsd | GPL      | founder genotype likelihoods & structure |
| ngsLD                            | github.com/fgvieira/ngsLD                           | GPL-2    | founder linkage disequilibrium           |
| BayPass                          | www1.montpellier.inrae.fr/CBGP/software/baypass     | CeCILL-B | wild host-association scan               |
| haplovalidate / haploReconstruct | github.com/popgenvienna/haplovalidate; CRAN         | GPL-3    | temporal haplotype blocks                |
| ACER                             | github.com/MartaPelizzola/ACER                      | GPL      | CMH test for AF change                   |
| poolSeq                          | github.com/ThomasTaus/poolSeq                       | GPL      | within-replicate Ne                      |
| cvtkpy                           | github.com/vsbuffalo/cvtk                           | GPL-3    | temporal covariance / G statistic        |
| PoolSNP                          | github.com/capoony/PoolSNP                          | MIT      | pooled SNP calling                       |
| popoolation2                     | sourceforge.net/p/popoolation2                      | GPL-3    | mpileup→sync                             |
| diamond                          | github.com/bbuchfink/diamond                        | GPL-3    | Muller-element assignment                |
| grenepipe                        | github.com/moiexpositoalonsolab/grenepipe           | GPL-3    | base workflow for the two overlay dirs   |


## AI assistance

Portions of this code and documentation were developed with the assistance of AI coding tools (Anthropic's Claude). All analyses, scientific decisions, and final code were designed, reviewed, and validated by the authors. 


#### License

This project's code is released under the **GNU General Public License v3.0**.

