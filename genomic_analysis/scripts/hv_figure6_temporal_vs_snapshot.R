#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(svglite)
})

setwd(".")

COL_FND  <- "#000000"
COL_T2   <- "#D55E00"
COL_REF  <- "grey70"

blend_alpha <- function(col, alpha, bg = "white") {
  rc <- col2rgb(col) / 255; bc <- col2rgb(bg) / 255
  out <- rc * alpha + bc * (1 - alpha)
  rgb(out[1,1], out[2,1], out[3,1])
}

sym_to_col <- function(s, max_s = 1.3) {
  s_norm <- pmin(pmax(s, 0), max_s) / max_s
  bg <- col2rgb("#DDDDDD") / 255
  fg <- col2rgb("#000000") / 255
  sapply(s_norm, function(x) {
    rgb_out <- fg * x + bg * (1 - x)
    rgb(rgb_out[1], rgb_out[2], rgb_out[3])
  })
}

OUT_SVG <- "hv_results_cluster/svg"
OUT_PNG <- "hv_results_cluster/block_dynamics"

dyn <- fread("hv_results_cluster/block_dynamics/block_dynamics_summary.tsv")
dyn <- dyn[treatment %in% c("B","T","M")]

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

dyn <- dyn[!is.na(rz_pbs_B) & !is.na(rz_pbs_T)]

winr2 <- fread("ngsld_inputs/windowed_r2.tsv")
winr2[, win_mid := (win_start + win_end) / 2]

scaff_order <- c("chr_ScDA7r2_110_HRSCAF_295",
                  "chr_ScDA7r2_126_HRSCAF_325",
                  "chr_ScDA7r2_439_HRSCAF_779",
                  "chr_ScDA7r2_597_HRSCAF_953")
scaff_labels <- c("110", "126", "439", "597")
scaff_sizes <- winr2[, .(max_pos = max(win_end)), by = scaffold]
setkey(scaff_sizes, scaffold)

gap <- 2e6
cum_offset <- numeric(4)
for (i in seq_along(scaff_order)) {
  if (i == 1) cum_offset[i] <- 0
  else cum_offset[i] <- cum_offset[i-1] + scaff_sizes[scaff_order[i-1]]$max_pos + gap
}
names(cum_offset) <- scaff_order

to_genome_x <- function(scaff, pos) cum_offset[scaff] + pos
genome_max <- cum_offset[4] + scaff_sizes[scaff_order[4]]$max_pos

winr2_wide <- dcast(winr2, scaffold + win_start + win_end + win_mid ~ cohort,
                     value.var = c("mean_r2", "n_pairs"))
winr2_wide[, gx := to_genome_x(scaffold, win_mid)]
winr2_wide[, avg_r2 := rowMeans(.SD, na.rm = TRUE),
            .SDcols = c("mean_r2_founders", "mean_r2_T2G07")]
winr2_wide[, min_pairs := pmin(n_pairs_founders, n_pairs_T2G07, na.rm = TRUE)]
winr2_wide[is.na(min_pairs), min_pairs := 0]
winr2_wide[, tier := fifelse(min_pairs < 50000, "lo",
                      fifelse(min_pairs < 100000, "med", "hi"))]

make_fig <- function() {
  layout(matrix(1:3, nrow = 3), heights = c(1, 1, 1.2))

  par(mar = c(1, 5, 3, 1.5), mgp = c(3.0, 0.75, 0))

  trt_y <- c(B = 3, T = 2, M = 1)
  trt_labels <- c("B", "T", "B+T")

  plot(NA, xlim = c(0, genome_max), ylim = c(0.3, 3.8),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "",
       main = "Haplovalidate blocks (pool-seq: 10 generations × 4 replicates)",
       cex.main = 1.1)
  axis(2, at = c(1, 2, 3), labels = rev(trt_labels), las = 1, cex.axis = 0.9, tick = FALSE)

  for (i in seq_along(scaff_order)) {
    x0 <- cum_offset[i]
    x1 <- x0 + scaff_sizes[scaff_order[i]]$max_pos
    if (i > 1) abline(v = x0 - gap/2, col = "grey85", lty = 1)
    text((x0 + x1) / 2, 3.7, scaff_labels[i], cex = 0.9, col = "grey40")
  }

  for (i in seq_len(nrow(dyn))) {
    d <- dyn[i]
    trt <- d$treatment
    if (!(trt %in% names(trt_y))) next
    y_center <- trt_y[trt]
    x0 <- to_genome_x(d$chr, d$start)
    x1 <- to_genome_x(d$chr, d$end)
    col <- sym_to_col(d$sym_score)
    rect(x0, y_center - 0.35, x1, y_center + 0.35, col = col, border = NA)
    rect(x0, y_center - 0.35, x1, y_center + 0.35,
         col = NA, border = blend_alpha("black", 0.3), lwd = 0.3)
  }

  for (trt in c("B","T","M")) {
    n_t <- sum(dyn$treatment == trt)
    text(genome_max * 1.0, trt_y[trt],
         sprintf("n=%d", n_t), pos = 2, cex = 0.8, col = "grey40")
  }

  legend("topright", bty = "n", cex = 0.8,
         legend = c("High symmetry", "Low symmetry"),
         fill = c(sym_to_col(1.3), sym_to_col(0)),
         border = "grey50")

  draw_r2_track <- function(cohort_label, cohort_name, is_bottom) {
    par(mar = c(if (is_bottom) 4.5 else 1.5, 5, 1.5, 1.5), mgp = c(3.0, 0.75, 0))

    r2_range <- range(winr2$mean_r2, na.rm = TRUE)
    r2_range <- c(r2_range[1] * 0.95, r2_range[2] * 1.05)

    plot(NA, xlim = c(0, genome_max), ylim = r2_range,
         xaxt = "n",
         xlab = if (is_bottom) "Genomic position (Mb)" else "",
         ylab = expression("Mean " * r^2 * " (10-100 kb)"),
         main = "", las = 1, cex.lab = 1, cex.axis = 0.9)

    for (i in seq_along(scaff_order)) {
      x0 <- cum_offset[i]
      x1 <- x0 + scaff_sizes[scaff_order[i]]$max_pos
      if (i > 1) abline(v = x0 - gap/2, col = "grey85")
      if (is_bottom) {
        ticks_mb <- seq(0, scaff_sizes[scaff_order[i]]$max_pos / 1e6, by = 20)
        axis(1, at = x0 + ticks_mb * 1e6, labels = ticks_mb, cex.axis = 0.7)
      }
      text((x0 + x1) / 2, r2_range[2] * 0.98, scaff_labels[i],
           cex = 0.9, col = "grey40")
    }

    sub <- winr2[cohort == cohort_name]
    sub[, gx := to_genome_x(scaffold, win_mid)]
    sub[, tier := fifelse(n_pairs < 50000, "lo",
                   fifelse(n_pairs < 100000, "med", "hi"))]

    COL_LO  <- "grey85"
    COL_MED <- "grey55"
    COL_HI  <- "grey25"

    mean_hi <- mean(sub[tier == "hi"]$mean_r2, na.rm = TRUE)
    abline(h = mean_hi, col = "grey40", lty = 3, lwd = 0.8)

    for (t in c("lo", "med", "hi")) {
      s <- sub[tier == t]
      col <- switch(t, lo = COL_LO, med = COL_MED, hi = COL_HI)
      sz  <- switch(t, lo = 0.3, med = 0.5, hi = 0.6)
      ptype <- switch(t, lo = 1, med = 19, hi = 19)
      points(s$gx, s$mean_r2, pch = ptype, cex = sz, col = col)
    }

    mtext(cohort_label, side = 3, line = 0.2, cex = 0.85, font = 2, adj = 0)

    if (!is_bottom) {
      legend("topleft", bty = "n", cex = 0.75, inset = c(0.12, 0),
             legend = c(expression("">=100*"k pairs"),
                        "50-100k pairs",
                        "<50k pairs"),
             pch = c(19, 19, 1),
             col = c(COL_HI, COL_MED, COL_LO),
             pt.cex = c(0.8, 0.6, 0.4))
    }
  }

  draw_r2_track("Founders (n=128, 0.5\u00d7)", "founders", FALSE)
  draw_r2_track("T2G07 (n=64, 0.5\u00d7)", "T2G07", TRUE)
}

svglite(file.path(OUT_SVG, "figure6_temporal_vs_snapshot.svg"),
        width = 14, height = 9)
make_fig(); dev.off()

png(file.path(OUT_PNG, "figure6_temporal_vs_snapshot.png"),
    width = 2800, height = 1800, res = 200)
make_fig(); dev.off()