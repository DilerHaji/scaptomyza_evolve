#!/usr/bin/env Rscript

suppressPackageStartupMessages({ library(data.table) })

setwd(".")

dyn <- fread("hv_results_cluster/block_dynamics/block_dynamics_summary.tsv")

bed <- unique(dyn[, .(chr, start, end, block_tag)])

fwrite(bed, "ngsld_inputs/blocks.bed", sep = "\t", col.names = FALSE)

pbs <- fread("variance_analysis/cvtkpy_final/pbs_vs_bt_aligned.tsv")
dyn_pbs <- copy(dyn)
dyn_pbs[, `:=`(rz_pbs_B = NA_real_, rz_pbs_T = NA_real_)]
for (i in 1:nrow(dyn_pbs)) {
  sub_pbs <- pbs[chrom == dyn_pbs$chr[i] &
                   start <= dyn_pbs$end[i] & end >= dyn_pbs$start[i]]
  if (nrow(sub_pbs) > 0) {
    dyn_pbs[i, rz_pbs_B := mean(sub_pbs$rz_pbs_B, na.rm = TRUE)]
    dyn_pbs[i, rz_pbs_T := mean(sub_pbs$rz_pbs_T, na.rm = TRUE)]
  }
}
dyn_pbs <- dyn_pbs[!is.na(rz_pbs_B) & !is.na(rz_pbs_T)]
dyn_pbs[, combined_z := (rz_pbs_B + rz_pbs_T) / 2]
dyn_pbs[, direction := ifelse(combined_z > 0,
                                (rz_pbs_T - rz_pbs_B) / (rz_pbs_T + rz_pbs_B), 0)]
dyn_pbs[, sym_score := ifelse(combined_z > 0,
                                combined_z * (1 - pmin(abs(direction), 1)), 0)]

by_tag <- dyn_pbs[, .(chr = chr[1], start = min(start), end = max(end),
                       sym_score = max(sym_score), span_Mb = max(span_Mb)),
                   by = block_tag]

by_tag <- by_tag[span_Mb >= 1]
setorder(by_tag, -sym_score)

TOP_N <- 15
top <- head(by_tag, TOP_N)

fwrite(top[, .(block_tag, chr, start, end, sym = sym_score)],
       "ngsld_inputs/top_blocks.tsv", sep = "\t")