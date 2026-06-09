#!/usr/bin/env Rscript

# glm_merge_haf.r
# Merges multiple split .Rdata and .RDS files into single files.

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 4) {
  stop("Not enough arguments. Need outputs and at least one input pair.")
}

output_rdata <- args[1]
output_rds   <- args[2]

remaining_args <- args[3:length(args)]
num_inputs <- length(remaining_args)

if(num_inputs %% 2 != 0) {
  stop("Uneven number of input files. Rdata and RDS counts must match.")
}

half_point <- num_inputs / 2
rdata_files <- remaining_args[1:half_point]
rds_files   <- remaining_args[(half_point + 1):num_inputs]

list_afmat <- list()
list_sites <- list()
list_neff  <- list()
final_samps <- NULL

valid_chunks <- 0

for(i in seq_along(rdata_files)) {
  f_rdata <- rdata_files[i]
  f_rds   <- rds_files[i]

  if(file.size(f_rdata) > 0 && file.size(f_rds) > 0) {
    
    chunk_neff <- readRDS(f_rds)

    e <- new.env()
    load(f_rdata, envir = e)
    
    if(is.null(final_samps)) {
      final_samps <- e$samps
    } else {
    }
    
    valid_chunks <- valid_chunks + 1
    list_afmat[[valid_chunks]] <- e$afmat
    list_sites[[valid_chunks]] <- e$sites
    list_neff[[valid_chunks]]  <- chunk_neff
    
    rm(e, chunk_neff) # cleanup
  } 
}

if(valid_chunks == 0) {
  file.create(output_rdata)
  file.create(output_rds)
  quit(save="no")
}


full_afmat <- do.call(rbind, list_afmat)
full_sites <- do.call(rbind, list_sites)
full_neff  <- do.call(rbind, list_neff)

if(nrow(full_afmat) != nrow(full_neff) || nrow(full_afmat) != nrow(full_sites)) {
  stop("Dimension mismatch after merge!")
}

saveRDS(full_neff, file = output_rds)

afmat <- full_afmat
sites <- full_sites
samps <- final_samps

save(afmat, sites, samps, file = output_rdata)
