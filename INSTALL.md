Installation & dependencies
This code uses a Linux HPC environment with conda/mamba, Snakemake (≥8), and SLURM.

- Python ≥ 3.10 (numpy, pandas, polars, scipy, scikit-learn, matplotlib)
- R ≥ 4.2 (lme4, data.table, tidyverse)
- Snakemake ≥ 8.15 with the SLURM executor plugin (snakemake-executor-plugin-slurm)
- conda or mamba

Conda environments are provided and referenced by the Snakemake rules:
genomic_analysis/environments/   # ngsld.yaml, haplovalidate.yaml
genomic_analysis/envs/           # relernn.yml (and others)
founder_popstructure/envs/       # angsd.yaml, pcangsd.yaml, bcftools.yaml, plotting_founders.yaml
pool_overdispersion/envs/        # samtools.yaml, qualimap.yaml, pool_comparison.yml

Read processing, mapping & SNP calling (Table S4)
| Tool                  | Version | Install                                                                                        |
| :-------------------- | :------ | :--------------------------------------------------------------------------------------------- |
| fastp                 | 0.20.1  | conda install -c bioconda fastp=0.20.1                                                         |
| BWA-MEM               | 0.7.17  | conda install -c bioconda bwa=0.7.17                                                           |
| samtools              | 1.x     | conda install -c bioconda samtools                                                             |
| Picard MarkDuplicates | 2.27.5  | conda install -c bioconda picard=2.27.5                                                        |
| GATK (IndelRealigner) | 3.x     | legacy GATK 3 (RealignerTargetCreator/IndelRealigner)                                          |
| PoolSNP               | latest  | https://github.com/capoony/PoolSNP (MIT)                                                       |
| popoolation2          | 1.201   | https://sourceforge.net/p/popoolation2 — provides mpileup2sync.jar → set PATH/TO/popoolation2/ |Population gemomics
| Tool             | Version | Upstream                                                |
| :--------------- | :------ | :------------------------------------------------------ |
| grenedalf        | 0.6.2   | https://github.com/lczech/grenedalf                     |
| ANGSD            | 0.940   | https://github.com/ANGSD/angsd                          |
| PCAngsd          | latest  | https://github.com/Rosemeis/pcangsd                     |
| ngsLD            | latest  | https://github.com/fgvieira/ngsLD                       |
| BayPass          | latest  | https://www1.montpellier.inrae.fr/CBGP/software/baypass |
| haplovalidate    | latest  | https://github.com/popgenvienna/haplovalidate           |
| haploReconstruct | latest  | CRAN / https://github.com/popgenvienna/haploReconstruct |
| ACER             | latest  | https://github.com/MartaPelizzola/ACER (R)              |
| poolSeq          | latest  | https://github.com/ThomasTaus/poolSeq (R)               |
| diamond          | 2.1.10  | https://github.com/bbuchfink/diamond                    |cvtkpy (temporal covariance / G statistic)
git clone https://github.com/vsbuffalo/cvtk
pip install -e ./cvtk 

grenepipe
git clone https://github.com/moiexpositoalonsolab/grenepipe
