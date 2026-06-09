#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(svglite)
})

setwd(".")

COL_B    <- "#2B8A3E"
COL_T    <- "#AD1A57"
COL_M    <- "#1F497D"
COL_HIGH <- "#000000"
COL_LOW  <- "#CCCCCC"
COL_REF  <- "grey55"

blend_alpha <- function(col, alpha, bg = "white") {
  rc <- col2rgb(col) / 255; bc <- col2rgb(bg) / 255
  out <- rc * alpha + bc * (1 - alpha)
  rgb(out[1,1], out[2,1], out[3,1])
}

OUT_SVG <- "hv_results_cluster/svg"
OUT_PNG <- "hv_results_cluster/block_dynamics"

r2 <- fread("ngsld_inputs/block_r2_summary.tsv")
dyn <- fread("hv_results_cluster/block_dynamics/block_dynamics_summary.tsv")

pbs <- fread("variance_analysis/cvtkpy_final/pbs_vs_bt_aligned.tsv")
dyn[, `:=`(rz_pbs_B = NA_real_, rz_pbs_T = NA_real_)]
for (i in 1:nrow(dyn)) {
  sub_pbs <- pbs[chrom == dyn$chr[i] & start <= dyn$end[i] & end >= dyn$start[i]]
  if (nrow(sub_pbs) > 0) {
    dyn[i, rz_pbs_B := mean(sub_pbs$rz_pbs_B, na.rm = TRUE)]
    dyn[i, rz_pbs_T := mean(sub_pbs$rz_pbs_T, na.rm = TRUE)]
  }
}
dyn[, combined_z := (rz_pbs_B + rz_pbs_T) / 2]
dyn[, direction := ifelse(!is.na(combined_z) & combined_z > 0,
                           (rz_pbs_T - rz_pbs_B) / (rz_pbs_T + rz_pbs_B), 0)]
dyn[, sym_score := ifelse(!is.na(combined_z) & combined_z > 0,
                           combined_z * (1 - pmin(abs(direction), 1)), 0)]

by_tag <- dyn[, .(chr = chr[1], start = min(start), end = max(end),
                   span_Mb = max(span_Mb), deltaAF = max(deltaAF),
                   sym_score = max(sym_score, na.rm = TRUE),
                   treatment = paste(unique(treatment), collapse = "+")),
               by = block_tag]
by_tag[is.infinite(sym_score), sym_score := 0]

r2_wide <- dcast(r2, block_tag + dist_bin ~ cohort + region,
                  value.var = "mean_r2")

r2_wide <- merge(r2_wide, by_tag, by = "block_tag")

r2_wide[, dr2_in := T2G07_inside   - founders_inside]
r2_wide[, dr2_fl := T2G07_flanking - founders_flanking]
r2_wide[, dr2_excess := dr2_in - dr2_fl]   # block-specific LD accumulation

r2_wide <- r2_wide[!is.na(dr2_in) & !is.na(dr2_fl)]

stats_by_bin <- r2_wide[, .(
    n = .N,
    median_dr2_in  = median(dr2_in),
    median_dr2_fl  = median(dr2_fl),
    median_excess  = median(dr2_excess),
    wilcox_p = wilcox.test(dr2_in, dr2_fl, paired = TRUE, alternative = "greater")$p.value,
    spearman_r_sym = cor(dr2_excess, sym_score, method = "spearman", use = "complete.obs")
  ), by = dist_bin]
setorder(stats_by_bin, dist_bin)

hisym <- r2_wide[sym_score > 0.5]
if (nrow(hisym) > 10) {
  print(hisym[, .(n = .N, median_excess = median(dr2_excess),
                   wilcox_p = wilcox.test(dr2_in, dr2_fl, paired = TRUE,
                                          alternative = "greater")$p.value),
               by = dist_bin])
}

make_fig <- function() {
  layout(matrix(c(1,2,3), nrow = 1), widths = c(1.1, 1.0, 1.0))
  par(mar = c(4.5, 4.8, 3.2, 1), mgp = c(2.8, 0.75, 0))
  sub10 <- r2_wide[dist_bin == "3_10-100kb"]
  max_sym <- max(r2_wide$sym_score, na.rm = TRUE)
  sym_to_col <- function(s) {
    s_norm <- pmin(s / max(0.1, max_sym), 1)
    rgb(1 - s_norm, 1 - s_norm, 1 - s_norm)
  }
  rng <- range(c(sub10$dr2_in, sub10$dr2_fl), na.rm = TRUE) * 1.05

  plot(sub10$dr2_fl, sub10$dr2_in,
       xlim = rng, ylim = rng,
       xlab = expression("Flanking " * Delta * r^2 * " (T2G07 − founders)"),
       ylab = expression("Within-block " * Delta * r^2 * " (T2G07 − founders)"),
       main = "A. 10-100 kb pairs",
       pch = 21, cex = 0.9 + log10(pmax(sub10$span_Mb, 0.1)) * 0.4,
       bg = sapply(sub10$sym_score, sym_to_col),
       col = "grey30", lwd = 0.6,
       las = 1, cex.main = 1.05, cex.lab = 1, cex.axis = 0.95)
  abline(h = 0, v = 0, col = COL_REF, lty = 3)
  abline(0, 1, col = "black", lty = 2, lwd = 1.3)

  s_row <- stats_by_bin[dist_bin == "3_10-100kb"]
  text(rng[1] + 0.02 * diff(rng), rng[2] * 0.98,
       bquote("Wilcoxon paired (inside > flank): " * italic(p) *
              " = " * .(formatC(s_row$wilcox_p, format = "e", digits = 1))),
       pos = 4, cex = 0.85)
  text(rng[1] + 0.02 * diff(rng), rng[2] * 0.92,
       bquote("Spearman " * italic(r) * "(excess " * Delta * r^2 *
              ", sym) = " * .(sprintf("%.2f", s_row$spearman_r_sym))),
       pos = 4, cex = 0.85)

  par(mar = c(4.5, 4.8, 3.2, 1), mgp = c(2.8, 0.75, 0))
  levels_bin <- c("1_<1kb","2_1-10kb","3_10-100kb")
  labels_bin <- c("<1 kb","1-10 kb","10-100 kb")
  bp_data <- lapply(levels_bin, function(b) r2_wide[dist_bin == b]$dr2_excess)
  max_y <- max(sapply(bp_data, function(v) max(abs(v), na.rm = TRUE))) * 1.1

  boxplot(bp_data, names = labels_bin,
          col = c(blend_alpha("grey", 0.3), blend_alpha("grey", 0.4),
                  blend_alpha("grey", 0.5)),
          ylab = expression("Excess " * Delta * r^2 * " (inside − flanking)"),
          xlab = "Pair distance bin",
          main = "B. Block-specific LD accumulation",
          las = 1, ylim = c(-max_y, max_y), outline = FALSE,
          cex.main = 1.05, cex.lab = 1, cex.axis = 0.95)
  stripchart(bp_data, vertical = TRUE, method = "jitter", jitter = 0.15,
             pch = 19, cex = 0.4, add = TRUE,
             col = c(blend_alpha("black", 0.4), blend_alpha("black", 0.5),
                     blend_alpha("black", 0.6)))
  abline(h = 0, col = COL_REF, lty = 3, lwd = 1)

  for (i in seq_along(levels_bin)) {
    p <- stats_by_bin[dist_bin == levels_bin[i]]$wilcox_p
    text(i, max_y * 0.92,
         sprintf("p = %s", formatC(p, format = "e", digits = 1)),
         cex = 0.75, font = 2)
  }

  par(mar = c(4.5, 4.8, 3.2, 1), mgp = c(2.8, 0.75, 0))
  sub10 <- r2_wide[dist_bin == "3_10-100kb"]
  plot(sub10$sym_score, sub10$dr2_excess,
       pch = 19, cex = 0.9 + log10(pmax(sub10$span_Mb, 0.1)) * 0.4,
       col = blend_alpha("black", 0.5),
       xlab = "Block symmetry score",
       ylab = expression("Excess " * Delta * r^2 * " (10-100 kb)"),
       main = "C. Symmetry predicts LD accumulation",
       las = 1, cex.main = 1.05, cex.lab = 1, cex.axis = 0.95)
  abline(h = 0, col = COL_REF, lty = 3)
  fit <- lm(dr2_excess ~ sym_score, data = sub10)
  abline(fit, col = "black", lwd = 1.8)

  r_sym <- cor(sub10$sym_score, sub10$dr2_excess, method = "spearman", use = "c")
  p_sym <- cor.test(sub10$sym_score, sub10$dr2_excess, method = "spearman", use = "c")$p.value
  text(min(sub10$sym_score) + 0.02 * diff(range(sub10$sym_score)),
       max(sub10$dr2_excess, na.rm = TRUE) * 0.95,
       bquote("Spearman " * italic(r) * " = " * .(sprintf("%.2f", r_sym))),
       pos = 4, cex = 0.9, font = 2)
  text(min(sub10$sym_score) + 0.02 * diff(range(sub10$sym_score)),
       max(sub10$dr2_excess, na.rm = TRUE) * 0.85,
       bquote(italic(p) * " = " * .(formatC(p_sym, format = "e", digits = 1))),
       pos = 4, cex = 0.85)
}

svglite(file.path(OUT_SVG, "figure6B_r2_enrichment.svg"), width = 15, height = 5)
make_fig(); dev.off()

png(file.path(OUT_PNG, "figure6B_r2_enrichment.png"), width = 3000, height = 1000, res = 200)
make_fig(); dev.off()

fwrite(r2_wide, "ngsld_inputs/block_r2_wide.tsv", sep = "\t")
