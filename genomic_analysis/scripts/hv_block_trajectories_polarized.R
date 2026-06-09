#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

setwd(".")

RESULTS_DIR <- "hv_results_cluster"
AD_FILE <- "variance_analysis/merged_ad.tsv"
SAMPLE_FILE <- "variance_analysis/sample_list.txt"
OUT_DIR <- "hv_results_cluster/trajectories_polarized"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PARAMSET <- "B_effective"
GENS <- c(0, 1, 2, 6, 7, 8, 9)
SCAFFOLDS <- c("chr_ScDA7r2_110_HRSCAF_295", "chr_ScDA7r2_126_HRSCAF_325",
               "chr_ScDA7r2_439_HRSCAF_779", "chr_ScDA7r2_597_HRSCAF_953")

COL_B <- "#4EA2FF"; COL_T <- "#EDB72F"; COL_M <- "#9BAB96"; COL_ANTAG <- "#D76161"
TRT_COL <- c(B = COL_B, T = COL_T, M = COL_M, BT = "#8B4DAF", BM = "#2E8B57", TM = "#FF8C00")


all_blocks <- list()
for (trt in c("B", "T", "M", "BT", "BM", "TM")) {
  f <- file.path(RESULTS_DIR, trt, PARAMSET, "dominant_blocks.tsv")
  if (file.exists(f)) {
    dt <- fread(f)
    if (nrow(dt) > 0) { dt[, treatment := trt]; all_blocks[[trt]] <- dt }
  }
}
blocks <- rbindlist(all_blocks)


all_samples <- readLines(SAMPLE_FILE)
all_samples <- all_samples[nzchar(all_samples)]

get_sample_idx <- function(trt_letter, gen, reps = 1:4) {
  names <- if (gen == 0) paste0("F", reps, "G00") else sprintf("%s%dG%02d", trt_letter, reps, gen)
  idx <- match(names, all_samples)
  idx[!is.na(idx)]
}


t0 <- proc.time()
ad <- fread(AD_FILE, header = FALSE, sep = "\t", colClasses = "character")
ad <- ad[V1 %in% SCAFFOLDS]

n_samples <- length(all_samples)
af_mat <- matrix(NA_real_, nrow = nrow(ad), ncol = n_samples)
for (j in 1:n_samples) {
  parts <- strsplit(ad[[j + 4]], ",", fixed = TRUE)
  rc <- as.integer(sapply(parts, `[`, 1))
  ac <- as.integer(sapply(parts, `[`, 2))
  tot <- rc + ac
  af_mat[, j] <- ifelse(tot > 0, ac / tot, NA_real_)
}
positions <- data.table(chr = ad$V1, pos = as.integer(ad$V2))
rm(ad); gc()
cat("  AF matrix:", nrow(af_mat), "sites,", round((proc.time() - t0)[3]), "sec\n")


get_polarized_trajectory <- function(snp_rows, trt_letter, polarize_using = trt_letter) {
  g0_idx <- get_sample_idx(polarize_using, 0)
  g9_idx <- get_sample_idx(polarize_using, 9)

  af_g0 <- rowMeans(af_mat[snp_rows, g0_idx, drop = FALSE], na.rm = TRUE)
  af_g9 <- rowMeans(af_mat[snp_rows, g9_idx, drop = FALSE], na.rm = TRUE)

  flip <- ifelse(af_g9 < af_g0, TRUE, FALSE)

  rows <- list()
  for (g in GENS) {
    s_idx <- get_sample_idx(trt_letter, g)
    for (r in seq_along(s_idx)) {
      af_raw <- af_mat[snp_rows, s_idx[r]]
      af_pol <- ifelse(flip, 1 - af_raw, af_raw)
      rows[[length(rows) + 1]] <- data.table(
        rep = r, gen = g,
        polarized_af = mean(af_pol, na.rm = TRUE),
        n_snps_used = sum(!is.na(af_pol))
      )
    }
  }
  rbindlist(rows)
}

traj_list <- list()

for (b_idx in seq_len(nrow(blocks))) {
  blk <- blocks[b_idx]
  snp_rows <- which(positions$chr == blk$chr & positions$pos >= blk$start & positions$pos <= blk$end)
  if (length(snp_rows) == 0) next

  trt <- blk$treatment

  if (nchar(trt) == 1) {
    t_data <- get_polarized_trajectory(snp_rows, trt, polarize_using = trt)
    t_data[, trt_letter := trt]
    t_data[, polarized_using := trt]
    t_data[, block_tag := blk$tag]
    t_data[, block_treatment := trt]
    t_data[, chr := blk$chr]
    t_data[, block_start := blk$start]
    t_data[, block_end := blk$end]
    t_data[, span_Mb := blk$span_Mb]
    t_data[, n_snps := blk$n_snps]
    t_data[, noThreshold := blk$noThreshold]
    traj_list[[length(traj_list) + 1]] <- t_data
  } else if (nchar(trt) == 2) {
    t1 <- substr(trt, 1, 1); t2 <- substr(trt, 2, 2)
    for (pl in c(t1, t2)) {
      for (use_trt in c(t1, t2)) {
        t_data <- get_polarized_trajectory(snp_rows, use_trt, polarize_using = pl)
        t_data[, trt_letter := use_trt]
        t_data[, polarized_using := pl]
        t_data[, block_tag := blk$tag]
        t_data[, block_treatment := trt]
        t_data[, chr := blk$chr]
        t_data[, block_start := blk$start]
        t_data[, block_end := blk$end]
        t_data[, span_Mb := blk$span_Mb]
        t_data[, n_snps := blk$n_snps]
        t_data[, noThreshold := blk$noThreshold]
        traj_list[[length(traj_list) + 1]] <- t_data
      }
    }
  }

  if (b_idx %% 50 == 0) { cat("  ", b_idx, "/", nrow(blocks), "\n"); flush.console() }
}

traj <- rbindlist(traj_list)
fwrite(traj, file.path(OUT_DIR, "block_trajectories_polarized.tsv"), sep = "\t")


single_summary <- traj[block_treatment == trt_letter & polarized_using == trt_letter]
blk_summary <- single_summary[, .(
  af_G0 = mean(polarized_af[gen == 0], na.rm = TRUE),
  af_G9 = mean(polarized_af[gen == 9], na.rm = TRUE),
  deltaAF = mean(polarized_af[gen == 9], na.rm = TRUE) - mean(polarized_af[gen == 0], na.rm = TRUE)
), by = .(block_tag, block_treatment, chr, block_start, block_end, span_Mb, n_snps, noThreshold)]

fwrite(blk_summary, file.path(OUT_DIR, "block_deltaAF_polarized.tsv"), sep = "\t")

png(file.path(OUT_DIR, "span_vs_deltaAF_polarized.png"), width = 1400, height = 500, res = 120)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

ymax <- max(blk_summary$deltaAF, na.rm = TRUE) * 1.05
for (trt in c("B", "T", "M")) {
  d <- blk_summary[block_treatment == trt]
  plot(d$span_Mb, d$deltaAF,
       pch = ifelse(d$noThreshold, 1, 19),
       col = adjustcolor(TRT_COL[trt], 0.7),
       xlab = "Block span (Mb)", ylab = "ΔAF polarized (G9 - G0)",
       main = paste0(trt, " (", nrow(d), " blocks)"),
       log = "x", xlim = c(0.1, 100), ylim = c(0, ymax))
  d <- d[order(-deltaAF)]
  top <- head(d, 5)
  for (k in 1:nrow(top)) {
    text(top$span_Mb[k], top$deltaAF[k],
         labels = paste0(gsub("chr_ScDA7r2_|_HRSCAF_.*", "", top$chr[k]), ":",
                         round(top$block_start[k] / 1e6, 1)),
         pos = 4, cex = 0.55, col = TRT_COL[trt])
  }
}
dev.off()

png(file.path(OUT_DIR, "top_trajectories_polarized.png"), width = 1800, height = 900, res = 120)
par(mfrow = c(3, 6), mar = c(3, 3.5, 2.5, 0.5), oma = c(2, 2, 2, 0))

for (trt in c("B", "T", "M")) {
  top_blocks <- blk_summary[block_treatment == trt][order(-deltaAF)][1:6]
  for (k in 1:nrow(top_blocks)) {
    blk_tag <- top_blocks$block_tag[k]
    t_data <- traj[block_tag == blk_tag & trt_letter == trt & polarized_using == trt]

    plot(NULL, xlim = c(0, 9), ylim = c(0, 1),
         xlab = "", ylab = "",
         main = paste0(gsub("chr_ScDA7r2_|_HRSCAF_.*", "", top_blocks$chr[k]), ":",
                       round(top_blocks$block_start[k] / 1e6, 1), "-",
                       round(top_blocks$block_end[k] / 1e6, 1), "Mb\n",
                       round(top_blocks$span_Mb[k], 1), "Mb, ΔAF=",
                       sprintf("%.2f", top_blocks$deltaAF[k])),
         cex.main = 0.7, las = 1)

    for (r in 1:4) {
      rd <- t_data[rep == r][order(gen)]
      if (nrow(rd) > 0) {
        lines(rd$gen, rd$polarized_af, col = adjustcolor(TRT_COL[trt], 0.4), lwd = 1)
        points(rd$gen, rd$polarized_af, col = adjustcolor(TRT_COL[trt], 0.7), pch = 19, cex = 0.5)
      }
    }
    mt <- t_data[, .(m = mean(polarized_af, na.rm = TRUE)), by = gen][order(gen)]
    lines(mt$gen, mt$m, col = TRT_COL[trt], lwd = 2.5)
  }
}
mtext("Generation", side = 1, outer = TRUE, line = 0.5, cex = 0.9)
mtext("Polarized allele frequency (rising allele)", side = 2, outer = TRUE, line = 0.5, cex = 0.9)
mtext("Top blocks by polarized ΔAF (B_effective)", outer = TRUE, cex = 1, font = 2)
dev.off()



bt_traj <- traj[block_treatment == "BT"]
if (nrow(bt_traj) > 0) {
  bt_b_pol <- bt_traj[polarized_using == "B"]
  bt_b_summary <- bt_b_pol[, .(
    af_G0 = mean(polarized_af[gen == 0], na.rm = TRUE),
    af_G9 = mean(polarized_af[gen == 9], na.rm = TRUE)
  ), by = .(block_tag, trt_letter, chr, block_start, block_end, span_Mb)]
  bt_b_summary[, deltaAF := af_G9 - af_G0]
  bt_wide <- dcast(bt_b_summary, block_tag + chr + block_start + block_end + span_Mb ~ trt_letter,
                   value.var = "deltaAF")
  setnames(bt_wide, c("B", "T"), c("deltaAF_B", "deltaAF_T"))

  bt_wide[, antagonistic := !is.na(deltaAF_T) & deltaAF_T < -0.01]  # T falls by at least 0.01 when B rises

  fwrite(bt_wide, file.path(OUT_DIR, "BT_block_divergence_polarized.tsv"), sep = "\t")
  png(file.path(OUT_DIR, "BT_antagonism_polarized.png"), width = 900, height = 700, res = 120)
  par(mar = c(4, 4, 3, 1))
  plot(bt_wide$deltaAF_B, bt_wide$deltaAF_T,
       pch = 19,
       col = ifelse(bt_wide$antagonistic, adjustcolor(COL_ANTAG, 0.8), adjustcolor("grey40", 0.5)),
       cex = pmin(2.5, bt_wide$span_Mb / 5 + 0.5),
       xlab = "ΔAF in B (polarized by B's rising allele)",
       ylab = "ΔAF in T (polarized by B's rising allele)",
       main = paste0("BT blocks: B vs T trajectories\n",
                     sum(bt_wide$antagonistic, na.rm = TRUE), " antagonistic / ",
                     nrow(bt_wide), " total"),
       cex.main = 0.9, las = 1)
  abline(h = 0, v = 0, col = "grey70")
  abline(a = 0, b = 1, col = "grey70", lty = 2)
  legend("topleft", legend = c("Concordant", "Antagonistic (T falls)"),
         col = c(adjustcolor("grey40", 0.7), adjustcolor(COL_ANTAG, 0.8)),
         pch = 19, cex = 0.75, bg = "white")
  dev.off()
}
