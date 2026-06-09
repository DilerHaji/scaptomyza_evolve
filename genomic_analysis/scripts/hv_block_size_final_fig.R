#!/usr/bin/env Rscript

suppressPackageStartupMessages({ library(data.table) })

setwd(".")
PARAMSET <- "B_effective"
OUT_DIR <- "hv_results_cluster/block_sizes"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

COL_B <- "#4EA2FF"; COL_T <- "#EDB72F"; COL_M <- "#9BAB96"
TRT_COL <- c(B = COL_B, T = COL_T, M = COL_M)

blocks_list <- list()
for (trt in c("B","T","M")) {
  dt <- fread(file.path("hv_results_cluster", trt, PARAMSET, "dominant_blocks.tsv"))
  if (nrow(dt) > 0) { dt[, treatment := trt]; blocks_list[[trt]] <- dt }
}
blocks <- rbindlist(blocks_list)
blocks[, treatment := factor(treatment, levels = c("B","T","M"))]

boot_median_ci <- function(x, n = 2000, conf = 0.95) {
  if (length(x) < 3) return(c(NA, NA, NA))
  b <- replicate(n, median(sample(x, replace = TRUE)))
  c(median(x), quantile(b, (1 - conf) / 2), quantile(b, 1 - (1 - conf) / 2))
}

ci <- blocks[, {
  res <- boot_median_ci(span_Mb)
  list(median = res[1], ci_lo = res[2], ci_hi = res[3], n = .N)
}, by = treatment]
print(ci)

pdf(file.path(OUT_DIR, "block_size_final.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), mgp = c(2.8, 0.8, 0))

x_pos <- c(B = 1, T = 2, M = 3)
set.seed(42)
plot(NULL, xlim = c(0.5, 3.5), ylim = c(0.1, 30), log = "y",
     xaxt = "n", xlab = "Treatment", ylab = "Block span (Mb)",
     main = "A. Block size per treatment", las = 1, cex.main = 1)

abline(h = c(0.45, 0.9, 1.8, 4.5), col = "grey90", lty = 3)
text(3.47, c(0.45, 0.9, 1.8, 4.5), labels = c("s≈0.005", "s≈0.01", "s≈0.02", "s≈0.05"),
     pos = 4, offset = 0.2, cex = 0.55, col = "grey50", xpd = TRUE)

for (trt in c("B","T","M")) {
  xs <- blocks[treatment == trt]
  ys <- xs$span_Mb
  xj <- x_pos[trt] + runif(nrow(xs), -0.25, 0.25)
  points(xj, ys,
         pch = ifelse(xs$noThreshold, 1, 19),
         col = adjustcolor(TRT_COL[trt], 0.35),
         cex = 0.8)
}

for (trt in c("B","T","M")) {
  c_row <- ci[treatment == trt]
  segments(x_pos[trt] - 0.15, c_row$median, x_pos[trt] + 0.15, c_row$median,
           col = TRT_COL[trt], lwd = 3)
  segments(x_pos[trt], c_row$ci_lo, x_pos[trt], c_row$ci_hi,
           col = TRT_COL[trt], lwd = 2)
  segments(x_pos[trt] - 0.08, c_row$ci_lo, x_pos[trt] + 0.08, c_row$ci_lo,
           col = TRT_COL[trt], lwd = 2)
  segments(x_pos[trt] - 0.08, c_row$ci_hi, x_pos[trt] + 0.08, c_row$ci_hi,
           col = TRT_COL[trt], lwd = 2)
}

axis(1, at = 1:3, labels = paste0(c("B","T","M"), "\n(n=", ci$n, ")"),
     tick = FALSE, cex.axis = 0.95, padj = 0.3)

y_star <- 18
segments(x_pos["B"], y_star, x_pos["M"], y_star, col = "black")
text(mean(x_pos[c("B","M")]), y_star * 1.1, labels = "*  (p = 0.022)", cex = 0.8)

size_classes <- c(0, 0.5, 1, 2, 5, 10, 100)
size_labels <- c("<0.5", "0.5-1", "1-2", "2-5", "5-10", ">10")
blocks[, size_class := cut(span_Mb, breaks = size_classes, labels = size_labels, include.lowest = TRUE)]
tab <- dcast(blocks[, .N, by = .(treatment, size_class)],
             size_class ~ treatment, value.var = "N", fill = 0)
setorder(tab, size_class)
mat <- as.matrix(tab[, -1])
rownames(mat) <- tab$size_class

par(mar = c(4.5, 4.5, 2.5, 1))
bp <- barplot(t(mat), beside = TRUE, col = TRT_COL[colnames(mat)],
              xlab = "Block span (Mb)", ylab = "Number of blocks",
              main = "B. Blocks stratified by span", las = 1, cex.names = 0.9,
              cex.main = 1, ylim = c(0, max(mat) * 1.25))

idx_25 <- which(rownames(mat) == "2-5")
x_mid <- mean(bp[, idx_25])
x_span <- diff(range(bp[, idx_25])) + 0.5
rect(x_mid - x_span/2 - 0.2, 0, x_mid + x_span/2 + 0.2, max(mat[idx_25,]) * 1.15,
     col = adjustcolor("#D76161", 0.08), border = adjustcolor("#D76161", 0.3), lty = 2)
text(x_mid, max(mat[idx_25,]) * 1.22,
     "2-5 Mb: s ≈ 0.01-0.02\nsingle-host selection",
     cex = 0.65, col = "#D76161", font = 3)

legend("topright", legend = colnames(mat), fill = TRT_COL[colnames(mat)],
       cex = 0.9, bg = "white", inset = c(0.02, 0.02))

mtext(paste0("Ne≈250, r≈2 cM/Mb, 9 generations. Dotted lines in A show hitchhiking-predicted block sizes for given selection coefficients."),
      side = 1, outer = TRUE, line = -1.5, cex = 0.65, col = "grey40")

dev.off()

png(file.path(OUT_DIR, "block_size_final.png"), width = 2000, height = 1000, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), mgp = c(2.8, 0.8, 0))

plot(NULL, xlim = c(0.5, 3.5), ylim = c(0.1, 30), log = "y",
     xaxt = "n", xlab = "Treatment", ylab = "Block span (Mb)",
     main = "A. Block size per treatment", las = 1, cex.main = 1)
abline(h = c(0.45, 0.9, 1.8, 4.5), col = "grey90", lty = 3)
text(3.47, c(0.45, 0.9, 1.8, 4.5), labels = c("s≈0.005", "s≈0.01", "s≈0.02", "s≈0.05"),
     pos = 4, offset = 0.2, cex = 0.55, col = "grey50", xpd = TRUE)

set.seed(42)
for (trt in c("B","T","M")) {
  xs <- blocks[treatment == trt]
  xj <- x_pos[trt] + runif(nrow(xs), -0.25, 0.25)
  points(xj, xs$span_Mb,
         pch = ifelse(xs$noThreshold, 1, 19),
         col = adjustcolor(TRT_COL[trt], 0.35), cex = 0.8)
}
for (trt in c("B","T","M")) {
  c_row <- ci[treatment == trt]
  segments(x_pos[trt] - 0.15, c_row$median, x_pos[trt] + 0.15, c_row$median,
           col = TRT_COL[trt], lwd = 3)
  segments(x_pos[trt], c_row$ci_lo, x_pos[trt], c_row$ci_hi, col = TRT_COL[trt], lwd = 2)
  segments(x_pos[trt] - 0.08, c_row$ci_lo, x_pos[trt] + 0.08, c_row$ci_lo, col = TRT_COL[trt], lwd = 2)
  segments(x_pos[trt] - 0.08, c_row$ci_hi, x_pos[trt] + 0.08, c_row$ci_hi, col = TRT_COL[trt], lwd = 2)
}
axis(1, at = 1:3, labels = paste0(c("B","T","M"), "\n(n=", ci$n, ")"),
     tick = FALSE, cex.axis = 0.95, padj = 0.3)
segments(x_pos["B"], 18, x_pos["M"], 18, col = "black")
text(mean(x_pos[c("B","M")]), 20, labels = "*  (p = 0.022)", cex = 0.8)

par(mar = c(4.5, 4.5, 2.5, 1))
bp <- barplot(t(mat), beside = TRUE, col = TRT_COL[colnames(mat)],
              xlab = "Block span (Mb)", ylab = "Number of blocks",
              main = "B. Blocks stratified by span", las = 1, cex.names = 0.9,
              cex.main = 1, ylim = c(0, max(mat) * 1.25))
idx_25 <- which(rownames(mat) == "2-5")
x_mid <- mean(bp[, idx_25])
x_span <- diff(range(bp[, idx_25])) + 0.5
rect(x_mid - x_span/2 - 0.2, 0, x_mid + x_span/2 + 0.2, max(mat[idx_25,]) * 1.15,
     col = adjustcolor("#D76161", 0.08), border = adjustcolor("#D76161", 0.3), lty = 2)
text(x_mid, max(mat[idx_25,]) * 1.22,
     "2-5 Mb: s ≈ 0.01-0.02\nsingle-host selection",
     cex = 0.65, col = "#D76161", font = 3)
legend("topright", legend = colnames(mat), fill = TRT_COL[colnames(mat)],
       cex = 0.9, bg = "white", inset = c(0.02, 0.02))

dev.off()

for (cmp in list(c("B","T"), c("B","M"), c("T","M"))) {
  x <- blocks[treatment == cmp[1]]$span_Mb
  y <- blocks[treatment == cmp[2]]$span_Mb
  wt <- suppressWarnings(wilcox.test(x, y))
  cat(sprintf("  %s vs %s: Wilcoxon p=%.3g\n", cmp[1], cmp[2], wt$p.value))
}

bt_25 <- sum(mat["2-5", c("B","T")])
m_25 <- mat["2-5", "M"]
bt_total <- sum(mat[, c("B","T")])
m_total <- sum(mat[, "M"])
ft <- fisher.test(matrix(c(bt_25, bt_total - bt_25, m_25, m_total - m_25), nrow = 2))