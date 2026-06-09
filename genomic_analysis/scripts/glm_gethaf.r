#!/usr/bin/env Rscript

if (!require("data.table", quietly = TRUE)) {
    install.packages("data.table", repos = "http://cran.us.r-project.org")
    library(data.table)
}

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 5) {
  stop("Usage: Rscript glm_gethaf.r <neff> <af> <sites> <samps> <output_prefix>")
}

neff_file <- args[1]
af_file <- args[2]
sites_file <- args[3]
samps_file <- args[4]
output_prefix <- args[5]

neff_dt <- fread(neff_file, header = TRUE, showProgress = FALSE)
neff <- as.matrix(neff_dt)
rm(neff_dt); gc()

af_dt <- fread(af_file, header = TRUE, showProgress = FALSE)
afmat <- as.matrix(af_dt)
rm(af_dt); gc()

sites <- fread(sites_file, header = TRUE, showProgress = FALSE)
colnames(sites) <- c("chrom", "pos")

samps <- fread(samps_file, header = TRUE, showProgress = FALSE)
samps <- as.data.frame(samps)

if(nrow(afmat) != nrow(neff)) stop("Error: Neff and AF matrices have different row counts")
if(ncol(afmat) != ncol(neff)) stop("Error: Neff and AF matrices have different column counts")
if(nrow(afmat) != nrow(sites)) stop(paste("Error: AF matrix rows (", nrow(afmat), ") do not match sites (", nrow(sites), ")"))

mat_samples <- colnames(afmat)
meta_samples <- samps$SourceFile

if(!identical(mat_samples, meta_samples)) {
    stop("Error: Sample order in AF matrix columns does not match 'SourceFile' column in samps CSV.")
}

out_rds <- paste0(output_prefix, "_neff.RDS")
out_rdata <- paste0(output_prefix, ".Rdata")

saveRDS(neff, file = out_rds)

save(afmat, sites, samps, file = out_rdata)