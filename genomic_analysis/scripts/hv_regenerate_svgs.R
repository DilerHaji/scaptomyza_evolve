#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(svglite)
})

setwd(".")

blend_alpha <- function(col, alpha, bg = "white") {
  rgb_col <- col2rgb(col) / 255
  rgb_bg <- col2rgb(bg) / 255
  rgb_out <- rgb_col * alpha + rgb_bg * (1 - alpha)
  rgb(rgb_out[1,1], rgb_out[2,1], rgb_out[3,1])
}

COL_B <- "#4EA2FF"; COL_T <- "#EDB72F"; COL_M <- "#9BAB96"
COL_HIGHLIGHT <- "#D76161"
COL_BT <- "#8B4DAF"; COL_BM <- "#2E8B57"; COL_TM <- "#FF8C00"
TRT_COL <- c(B = COL_B, T = COL_T, M = COL_M, BT = COL_BT, BM = COL_BM, TM = COL_TM)

COL_B_35  <- blend_alpha(COL_B, 0.35)
COL_T_35  <- blend_alpha(COL_T, 0.35)
COL_M_35  <- blend_alpha(COL_M, 0.35)
COL_HL_50 <- blend_alpha(COL_HIGHLIGHT, 0.5)
COL_HL_30 <- blend_alpha(COL_HIGHLIGHT, 0.3)
COL_HL_08 <- blend_alpha(COL_HIGHLIGHT, 0.08)

col_with_alpha <- function(trt, alpha) blend_alpha(TRT_COL[trt], alpha)

PARAMSET <- "B_effective"
OUT_SVG <- "hv_results_cluster/svg"
dir.create(OUT_SVG, recursive = TRUE, showWarnings = FALSE)

blocks_list <- list()
for (trt in c("B","T","M","BT","BM","TM")) {
  f <- file.path("hv_results_cluster", trt, PARAMSET, "dominant_blocks.tsv")
  if (file.exists(f)) {
    dt <- fread(f)
    if (nrow(dt) > 0) { dt[, treatment := trt]; blocks_list[[trt]] <- dt }
  }
}
blocks <- rbindlist(blocks_list)

dyn <- fread("hv_results_cluster/block_dynamics/block_dynamics_summary.tsv")

scaff_order <- c("chr_ScDA7r2_110_HRSCAF_295", "chr_ScDA7r2_126_HRSCAF_325",
                 "chr_ScDA7r2_439_HRSCAF_779", "chr_ScDA7r2_597_HRSCAF_953")
scaff_sizes <- c(chr_ScDA7r2_110_HRSCAF_295 = 48877389,
                 chr_ScDA7r2_126_HRSCAF_325 = 30545709,
                 chr_ScDA7r2_439_HRSCAF_779 = 31830347,
                 chr_ScDA7r2_597_HRSCAF_953 = 92653032)
scaff_labels <- c("110", "126", "439", "597")
gap <- 3e6
offsets <- cumsum(c(0, head(scaff_sizes[scaff_order], -1) + gap))
mids <- offsets + scaff_sizes[scaff_order] / 2
genome_len <- max(offsets + scaff_sizes[scaff_order])

dyn[, scaff_idx := match(chr, scaff_order)]
dyn[, g_start := start + offsets[scaff_idx]]
dyn[, g_end := end + offsets[scaff_idx]]

m_giants <- dyn[treatment == "M" & span_Mb > 10][order(-span_Mb)]
m_giants[, g_start := start + offsets[match(chr, scaff_order)]]
m_giants[, g_end := end + offsets[match(chr, scaff_order)]]
m_giants[, label := sprintf("M%d", 1:.N)]

boot_median_ci <- function(x, n = 2000, conf = 0.95) {
  b <- replicate(n, median(sample(x, replace = TRUE)))
  c(median(x), quantile(b, (1 - conf) / 2), quantile(b, 1 - (1 - conf) / 2))
}
ci <- blocks[treatment %in% c("B","T","M"),
             {res <- boot_median_ci(span_Mb)
              list(median = res[1], ci_lo = res[2], ci_hi = res[3], n = .N)},
             by = treatment]
ci[, treatment := factor(treatment, levels = c("B","T","M"))]
setorder(ci, treatment)


svglite(file.path(OUT_SVG, "block_size_final.svg"), width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), mgp = c(2.8, 0.8, 0))

x_pos <- c(B = 1, T = 2, M = 3)
plot(NULL, xlim = c(0.5, 3.5), ylim = c(0.1, 30), log = "y",
     xaxt = "n", xlab = "Treatment", ylab = "Block span (Mb)",
     main = "A. Block size per treatment", las = 1, cex.main = 1)
abline(h = c(0.45, 0.9, 1.8, 4.5), col = "grey90", lty = 3)
text(3.47, c(0.45, 0.9, 1.8, 4.5), labels = c("s≈0.005","s≈0.01","s≈0.02","s≈0.05"),
     pos = 4, offset = 0.2, cex = 0.55, col = "grey50", xpd = TRUE)

set.seed(42)
for (trt in c("B","T","M")) {
  xs <- blocks[treatment == trt]
  xj <- x_pos[trt] + runif(nrow(xs), -0.25, 0.25)
  is_giant <- trt == "M" & xs$span_Mb > 10
  points(xj[!is_giant], xs$span_Mb[!is_giant],
         pch = ifelse(xs$noThreshold[!is_giant], 1, 19),
         col = col_with_alpha(trt, 0.35), cex = 0.8)
  if (any(is_giant)) {
    points(xj[is_giant], xs$span_Mb[is_giant],
           pch = 21, col = COL_HIGHLIGHT, bg = COL_HL_50,
           cex = 1.5, lwd = 1.5)
  }
}

for (trt in c("B","T","M")) {
  cr <- ci[treatment == trt]
  segments(x_pos[trt] - 0.15, cr$median, x_pos[trt] + 0.15, cr$median,
           col = TRT_COL[trt], lwd = 3)
  segments(x_pos[trt], cr$ci_lo, x_pos[trt], cr$ci_hi, col = TRT_COL[trt], lwd = 2)
  segments(x_pos[trt] - 0.08, cr$ci_lo, x_pos[trt] + 0.08, cr$ci_lo, col = TRT_COL[trt], lwd = 2)
  segments(x_pos[trt] - 0.08, cr$ci_hi, x_pos[trt] + 0.08, cr$ci_hi, col = TRT_COL[trt], lwd = 2)
}
axis(1, at = 1:3, labels = paste0(c("B","T","M"), "\n(n=", ci$n, ")"),
     tick = FALSE, cex.axis = 0.95, padj = 0.3)
segments(x_pos["B"], 18, x_pos["M"], 18, col = "black")
text(mean(x_pos[c("B","M")]), 22, labels = "*  (p = 0.022)", cex = 0.8)

legend("bottomleft",
       legend = c("Dominant block", "noThreshold block", "M >10 Mb (n=4)"),
       pch = c(19, 1, 21), col = c("grey50", "grey50", COL_HIGHLIGHT),
       pt.bg = c(NA, NA, COL_HL_50),
       pt.cex = c(0.8, 0.8, 1.5), cex = 0.7, bg = "white")

size_classes <- c(0, 0.5, 1, 2, 5, 10, 100)
size_labels <- c("<0.5","0.5-1","1-2","2-5","5-10",">10")
blocks_abt <- blocks[treatment %in% c("B","T","M")]
blocks_abt[, size_class := cut(span_Mb, breaks = size_classes, labels = size_labels, include.lowest = TRUE)]
tab <- dcast(blocks_abt[, .N, by = .(treatment, size_class)],
             size_class ~ treatment, value.var = "N", fill = 0)
setorder(tab, size_class)
mat <- as.matrix(tab[, c("B","T","M")])
rownames(mat) <- tab$size_class

par(mar = c(4.5, 4.5, 3, 1))
bp <- barplot(t(mat), beside = TRUE, col = TRT_COL[colnames(mat)],
              xlab = "Block span (Mb)", ylab = "Number of blocks",
              main = "B. Blocks stratified by span", las = 1, cex.names = 0.9,
              cex.main = 1, ylim = c(0, max(mat) * 1.25))

idx_10 <- which(rownames(mat) == ">10")
m_col <- which(colnames(mat) == "M")
rect_x <- bp[m_col, idx_10] - 0.4
rect(rect_x, 0, rect_x + 0.8, mat[idx_10, "M"],
     border = COL_HIGHLIGHT, lwd = 2.5)
arrows(bp[m_col, idx_10], mat[idx_10, "M"] + 4,
       bp[m_col, idx_10], mat[idx_10, "M"] + 0.5,
       length = 0.08, col = COL_HIGHLIGHT, lwd = 2)
text(bp[m_col, idx_10], mat[idx_10, "M"] + 6,
     sprintf("M has %d\nblocks >10 Mb", mat[idx_10, "M"]),
     cex = 0.7, col = COL_HIGHLIGHT, font = 2)
legend("topright", legend = colnames(mat), fill = TRT_COL[colnames(mat)],
       cex = 0.9, bg = "white", inset = c(0.02, 0.02))
dev.off()


svglite(file.path(OUT_SVG, "M_giants_manhattan.svg"), width = 14, height = 7)
par(mfrow = c(3, 1), mar = c(1.5, 5, 1.5, 1), oma = c(3, 0, 3, 0))
ymax <- max(dyn$deltaAF, na.rm = TRUE) * 1.05

for (trt in c("B","T","M")) {
  s <- dyn[treatment == trt]
  plot(NULL, xlim = c(0, genome_len), ylim = c(0, ymax),
       xaxt = "n", xlab = "", ylab = "ΔAF", las = 1,
       main = trt, cex.main = 1, font.main = 2, col.main = TRT_COL[trt])

  for (j in seq_along(scaff_order)) {
    rect(offsets[j], -1, offsets[j] + scaff_sizes[scaff_order[j]], 2,
         col = ifelse(j %% 2 == 0, "grey96", "grey92"), border = NA)
  }

  for (i in 1:nrow(m_giants)) {
    rect(m_giants$g_start[i], -1, m_giants$g_end[i], 2,
         col = COL_HL_08, border = COL_HL_30, lty = 2, lwd = 1)
  }

  if (nrow(s) > 0) {
    s_thresh <- s[noThreshold == FALSE]
    s_nt <- s[noThreshold == TRUE]
    col_t_90 <- col_with_alpha(trt, 0.90)
    col_t_50 <- col_with_alpha(trt, 0.50)
    if (nrow(s_thresh) > 0) {
      segments(s_thresh$g_start, s_thresh$deltaAF, s_thresh$g_end, s_thresh$deltaAF,
               col = col_t_90, lwd = 2.2)
    }
    if (nrow(s_nt) > 0) {
      segments(s_nt$g_start, s_nt$deltaAF, s_nt$g_end, s_nt$deltaAF,
               col = col_t_50, lwd = 1.2)
    }
  }

  if (trt == "M") {
    for (i in 1:nrow(m_giants)) {
      segments(m_giants$g_start[i], m_giants$deltaAF[i],
               m_giants$g_end[i], m_giants$deltaAF[i],
               col = COL_HIGHLIGHT, lwd = 4)
      text((m_giants$g_start[i] + m_giants$g_end[i]) / 2,
           ymax * 0.92, labels = m_giants$label[i],
           col = COL_HIGHLIGHT, font = 2, cex = 0.9)
    }
    axis(1, at = mids, labels = scaff_labels, tick = FALSE)
    mtext("Scaffold", side = 1, line = 2, cex = 0.85)
  }

  mtext(sprintf("n=%d blocks", nrow(s)), side = 3, line = -1.2,
        adj = 0.99, cex = 0.7, col = TRT_COL[trt])

  if (trt %in% c("B","T")) {
    for (i in 1:nrow(m_giants)) {
      n_over <- nrow(s[chr == m_giants$chr[i] &
                        start <= m_giants$end[i] & end >= m_giants$start[i]])
      text((m_giants$g_start[i] + m_giants$g_end[i]) / 2,
           ymax * 0.92, labels = sprintf("%s: %d blocks", m_giants$label[i], n_over),
           col = COL_HIGHLIGHT, font = 2, cex = 0.7)
    }
  }
}
mtext("M giant regions (>10 Mb blocks) — fragmentation across treatments",
      outer = TRUE, cex = 1.1, font = 2, line = 0.5)
dev.off()

svglite(file.path(OUT_SVG, "deltaAF_distribution.svg"), width = 9, height = 4.5)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1), mgp = c(2.8, 0.8, 0))

breaks <- seq(0, 0.10, by = 0.005)
h_B <- hist(dyn[treatment == "B"]$deltaAF, breaks = breaks, plot = FALSE)
h_T <- hist(dyn[treatment == "T"]$deltaAF, breaks = breaks, plot = FALSE)
h_M <- hist(dyn[treatment == "M"]$deltaAF, breaks = breaks, plot = FALSE)
ymax <- max(c(h_B$counts, h_T$counts, h_M$counts)) * 1.1

plot(h_B$mids, h_B$counts, type = "l", lwd = 2.5, col = COL_B,
     xlim = c(0, 0.10), ylim = c(0, ymax),
     xlab = "ΔAF (polarized, G9 - G0)", ylab = "Number of blocks",
     main = "A. Block ΔAF distribution", las = 1)
lines(h_T$mids, h_T$counts, lwd = 2.5, col = COL_T)
lines(h_M$mids, h_M$counts, lwd = 2.5, col = COL_M)
legend("topright", legend = c("B","T","M"), col = c(COL_B, COL_T, COL_M),
       lwd = 2.5, bg = "white", cex = 0.9)

ci_dyn <- dyn[treatment %in% c("B","T","M"),
              {res <- boot_median_ci(deltaAF)
               list(median = res[1], ci_lo = res[2], ci_hi = res[3])},
              by = treatment]
ci_dyn[, treatment := factor(treatment, levels = c("B","T","M"))]
setorder(ci_dyn, treatment)
x_pos <- c(B=1, T=2, M=3)

plot(NULL, xlim = c(0.5, 3.5), ylim = c(0, 0.08),
     xaxt = "n", xlab = "Treatment", ylab = "ΔAF (polarized)",
     main = "B. Block ΔAF per treatment", las = 1, cex.main = 1)
set.seed(42)
for (trt in c("B","T","M")) {
  d <- dyn[treatment == trt]
  xj <- x_pos[trt] + runif(nrow(d), -0.25, 0.25)
  points(xj, d$deltaAF,
         pch = ifelse(d$noThreshold, 1, 19),
         col = col_with_alpha(trt, 0.35), cex = 0.8)
}
for (trt in c("B","T","M")) {
  cr <- ci_dyn[treatment == trt]
  segments(x_pos[trt]-0.15, cr$median, x_pos[trt]+0.15, cr$median, col=TRT_COL[trt], lwd=3)
  segments(x_pos[trt], cr$ci_lo, x_pos[trt], cr$ci_hi, col=TRT_COL[trt], lwd=2)
  segments(x_pos[trt]-0.08, cr$ci_lo, x_pos[trt]+0.08, cr$ci_lo, col=TRT_COL[trt], lwd=2)
  segments(x_pos[trt]-0.08, cr$ci_hi, x_pos[trt]+0.08, cr$ci_hi, col=TRT_COL[trt], lwd=2)
}
axis(1, at = 1:3, labels = c("B","T","M"), tick = FALSE, cex.axis = 0.95)
dev.off()
cat("  ", file.path(OUT_SVG, "deltaAF_distribution.svg"), "\n")

svglite(file.path(OUT_SVG, "trajectory_shape.svg"), width = 12, height = 4.5)
par(mfrow = c(1, 3), mar = c(4, 4.5, 3, 1))
for (trt in c("B","T","M")) {
  s <- dyn[treatment == trt]
  plot(s$early_slope, s$late_slope,
       pch = 19, col = col_with_alpha(trt, 0.5),
       cex = pmin(2.5, s$span_Mb / 3 + 0.5),
       xlab = "Early slope (G0→G2)", ylab = "Late slope (G6→G9)",
       main = paste0(trt, " (n=", nrow(s), ")"),
       xlim = c(-0.01, 0.03), ylim = c(-0.02, 0.03),
       las = 1, cex.main = 1)
  abline(a = 0, b = 1, col = "grey60", lty = 2)
  abline(h = 0, v = 0, col = "grey80")
  mtext("Accelerating", side = 3, line = -1.5, adj = 0.98, cex = 0.6, col = "grey40")
  mtext("Decelerating", side = 3, line = -1.5, adj = 0.02, cex = 0.6, col = "grey40")
}
dev.off()

svglite(file.path(OUT_SVG, "span_vs_deltaAF.svg"), width = 12, height = 4.5)
par(mfrow = c(1, 3), mar = c(4, 4.5, 3, 1))
ymax <- max(dyn$deltaAF, na.rm = TRUE) * 1.05
for (trt in c("B","T","M")) {
  s <- dyn[treatment == trt]
  plot(s$span_Mb, s$deltaAF,
       pch = ifelse(s$noThreshold, 1, 19),
       col = col_with_alpha(trt, 0.7),
       xlab = "Block span (Mb)", ylab = "ΔAF polarized",
       main = paste0(trt, " (n=", nrow(s), ")"),
       log = "x", xlim = c(0.1, 100), ylim = c(0, ymax), las = 1)
  if (nrow(s) > 3) {
    r <- cor(log10(s$span_Mb), s$deltaAF, method = "spearman", use = "pairwise.complete.obs")
    mtext(sprintf("Spearman r = %.2f", r), side = 3, line = -1.5, adj = 0.02, cex = 0.8)
  }
}
dev.off()

svglite(file.path(OUT_SVG, "block_manhattan_all.svg"), width = 14, height = 6)
par(mfrow = c(3, 1), mar = c(1.5, 5, 1.5, 1), oma = c(3, 0, 2, 0))
ymax <- max(dyn$deltaAF, na.rm = TRUE) * 1.05

for (trt in c("B","T","M")) {
  s <- dyn[treatment == trt]
  plot(NULL, xlim = c(0, genome_len), ylim = c(0, ymax),
       xaxt = "n", xlab = "", ylab = "ΔAF polarized", las = 1,
       main = trt, cex.main = 1, col.main = TRT_COL[trt])
  for (j in seq_along(scaff_order)) {
    rect(offsets[j], -1, offsets[j] + scaff_sizes[scaff_order[j]], 2,
         col = ifelse(j %% 2 == 0, "grey96", "grey92"), border = NA)
  }
  if (nrow(s) > 0) {
    s_th <- s[noThreshold == FALSE]
    s_nt <- s[noThreshold == TRUE]
    c_th <- col_with_alpha(trt, 0.8)
    c_nt <- col_with_alpha(trt, 0.5)
    if (nrow(s_th) > 0) segments(s_th$g_start, s_th$deltaAF, s_th$g_end, s_th$deltaAF, col = c_th, lwd = 2)
    if (nrow(s_nt) > 0) segments(s_nt$g_start, s_nt$deltaAF, s_nt$g_end, s_nt$deltaAF, col = c_nt, lwd = 1.2)
    top5 <- s[order(-deltaAF)][1:5]
    points((top5$g_start + top5$g_end) / 2, top5$deltaAF,
           pch = 19, col = TRT_COL[trt], cex = 1.2)
  }
  if (trt == "M") axis(1, at = mids, labels = scaff_labels, tick = FALSE)
}
mtext("Scaffold", side = 1, outer = TRUE, line = 1.5, cex = 0.9)
mtext("Polarized block ΔAF along the genome (B_effective, segments = block spans)",
      outer = TRUE, cex = 1, font = 2, line = 0.5)
dev.off()

traj <- fread("hv_results_cluster/block_dynamics/rep_trajectories.tsv")
svglite(file.path(OUT_SVG, "top_trajectories.svg"), width = 16, height = 8)
par(mfrow = c(3, 6), mar = c(3, 3.5, 2.5, 0.5), oma = c(2, 2, 2, 0))

for (trt in c("B","T","M")) {
  top_blocks <- dyn[treatment == trt][order(-deltaAF)][1:6]
  for (k in 1:nrow(top_blocks)) {
    blk_tag <- top_blocks$block_tag[k]
    t_data <- traj[block_tag == blk_tag & treatment == trt]
    plot(NULL, xlim = c(0, 9), ylim = c(0, 1), xlab = "", ylab = "",
         main = paste0(gsub("chr_ScDA7r2_|_HRSCAF_.*", "", top_blocks$chr[k]), ":",
                       round(top_blocks$start[k] / 1e6, 1), "-",
                       round(top_blocks$end[k] / 1e6, 1), "Mb\n",
                       round(top_blocks$span_Mb[k], 1), "Mb, ΔAF=",
                       sprintf("%.2f", top_blocks$deltaAF[k])),
         cex.main = 0.7, las = 1)
    c_rep <- col_with_alpha(trt, 0.4)
    c_pt <- col_with_alpha(trt, 0.7)
    for (r in 1:4) {
      rd <- t_data[rep == r][order(gen)]
      if (nrow(rd) > 0) {
        lines(rd$gen, rd$polarized_af, col = c_rep, lwd = 1)
        points(rd$gen, rd$polarized_af, col = c_pt, pch = 19, cex = 0.5)
      }
    }
    mt <- t_data[, .(m = mean(polarized_af, na.rm = TRUE)), by = gen][order(gen)]
    lines(mt$gen, mt$m, col = TRT_COL[trt], lwd = 2.5)
  }
}
mtext("Generation", side = 1, outer = TRUE, line = 0.5, cex = 0.9)
mtext("Polarized allele frequency", side = 2, outer = TRUE, line = 0.5, cex = 0.9)
mtext("Top blocks by polarized ΔAF", outer = TRUE, cex = 1, font = 2)
dev.off()

