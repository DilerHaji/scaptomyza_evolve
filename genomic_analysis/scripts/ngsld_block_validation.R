#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

RESULTS_DIR <- "hv_results_cluster"
OUT_DIR <- "ngsld/block_validation"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SCAFFOLDS <- c("chr_ScDA7r2_110_HRSCAF_295", "chr_ScDA7r2_126_HRSCAF_325",
               "chr_ScDA7r2_439_HRSCAF_779", "chr_ScDA7r2_597_HRSCAF_953")
PARAMSET <- "B_effective"

all_blocks <- list()
for (trt in c("B", "T", "M")) {
  f <- file.path(RESULTS_DIR, trt, PARAMSET, "dominant_blocks.tsv")
  if (file.exists(f)) {
    dt <- fread(f)
    if (nrow(dt) > 0) {
      dt[, treatment := trt]
      all_blocks[[trt]] <- dt
    }
  }
}
blocks <- rbindlist(all_blocks)

ld_list <- list()
for (cohort in c("founders", "T2G07")) {
  for (scaff in SCAFFOLDS) {
    f <- file.path("ngsld/ld", cohort, paste0(scaff, ".ld.gz"))
    if (!file.exists(f)) {
      cat("  MISSING:", f, "\n")
      next
    }
    d <- fread(f, header = FALSE, col.names = c("site1", "site2", "dist",
                                                  "r2_ExpG", "D", "Dp", "r2"))
    d[, pos1 := as.integer(sub(".*:", "", site1))]
    d[, pos2 := as.integer(sub(".*:", "", site2))]
    d[, chr := scaff]
    d[, cohort := cohort]
    ld_list[[paste(cohort, scaff, sep = "_")]] <- d
  }
}
ld <- rbindlist(ld_list, fill = TRUE)


block_ld <- list()
for (i in 1:nrow(blocks)) {
  blk <- blocks[i]
  for (cohort in c("founders", "T2G07")) {
    pairs_inside <- ld[cohort == cohort & chr == blk$chr &
                        pos1 >= blk$start & pos1 <= blk$end &
                        pos2 >= blk$start & pos2 <= blk$end]
    if (nrow(pairs_inside) < 5) next
    block_ld[[length(block_ld) + 1]] <- data.table(
      block_tag = blk$tag, treatment = blk$treatment,
      chr = blk$chr, start = blk$start, end = blk$end,
      span_Mb = blk$span_Mb,
      cohort = cohort,
      n_pairs = nrow(pairs_inside),
      mean_r2 = mean(pairs_inside$r2, na.rm = TRUE),
      median_r2 = median(pairs_inside$r2, na.rm = TRUE),
      mean_r2_expG = mean(pairs_inside$r2_ExpG, na.rm = TRUE),
      mean_dist = mean(pairs_inside$dist, na.rm = TRUE)
    )
  }
  if (i %% 50 == 0) cat("  Block", i, "/", nrow(blocks), "\n")
}

block_summary <- rbindlist(block_ld)
fwrite(block_summary, file.path(OUT_DIR, "block_LD_summary.tsv"), sep = "\t")

block_ranges <- blocks[, .(starts = list(start), ends = list(end)), by = chr]

ld[, in_block := FALSE]
for (i in 1:nrow(block_ranges)) {
  starts <- unlist(block_ranges$starts[i])
  ends <- unlist(block_ranges$ends[i])
  for (j in seq_along(starts)) {
    ld[chr == block_ranges$chr[i] &
       ((pos1 >= starts[j] & pos1 <= ends[j]) |
        (pos2 >= starts[j] & pos2 <= ends[j])),
       in_block := TRUE]
  }
}

bg_summary <- ld[in_block == FALSE,
                  .(n_pairs = .N,
                    mean_r2 = mean(r2, na.rm = TRUE),
                    median_r2 = median(r2, na.rm = TRUE)),
                  by = cohort]

print(block_summary[, .(n_blocks = .N,
                         mean_mean_r2 = mean(mean_r2),
                         median_mean_r2 = median(mean_r2)),
                     by = cohort])