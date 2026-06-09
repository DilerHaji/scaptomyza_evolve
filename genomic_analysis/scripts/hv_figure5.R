#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(svglite)
})

setwd(".")

COL_FND <- "#000000"     # founders  = black
COL_T2  <- "#D55E00"     # T2G07     = vermillion
COL_REF <- "grey55"

blend_alpha <- function(col, alpha, bg = "white") {
  rc <- col2rgb(col) / 255; bc <- col2rgb(bg) / 255
  out <- rc * alpha + bc * (1 - alpha)
  rgb(out[1,1], out[2,1], out[3,1])
}

OUT_SVG <- "hv_results_cluster/svg"
OUT_PNG <- "hv_results_cluster/block_dynamics"
dir.create(OUT_SVG, recursive = TRUE, showWarnings = FALSE)

dat <- data.table(
  cohort   = rep(c("founders","T2G07"), each = 6),
  region   = rep(rep(c("hotspot","flanking"), each = 3), 2),
  dist_bin = rep(c("<1kb","1-10kb","10-100kb"), 4),
  dist_mid = rep(c(500, 3162, 31623), 4),   # log-midpoints
  n_pairs  = rep(c(240780, 871735, 7955579, 163750, 501205, 4112455), 2),
  r2       = c(
    0.083910, 0.018525, 0.017385,   # fnd hotspot
    0.093847, 0.019834, 0.018614,   # fnd flanking
    0.079088, 0.019639, 0.018798,   # T2 hotspot
    0.086885, 0.019828, 0.018822    # T2 flanking
  )
)
dat[, dist_bin := factor(dist_bin, levels = c("<1kb","1-10kb","10-100kb"))]

panelA <- function() {
  par(mar = c(4.8, 4.8, 3, 1.2), mgp = c(2.9, 0.75, 0))
  x_at <- 1:3
  plot(NA, NA,
       xlim = c(0.7, 3.3), ylim = c(0, 0.11),
       xaxt = "n", xlab = "Pair distance",
       ylab = expression("Mean " * italic(r)^2 * " (ngsLD, GL-aware)"),
       main = "A. LD decay, scaffold 439",
       las = 1, cex.main = 1.1, cex.lab = 1, cex.axis = 0.95)
  axis(1, at = x_at, labels = levels(dat$dist_bin), cex.axis = 0.95)

  series <- list(
    list(cohort = "founders", region = "hotspot",  col = COL_FND, lty = 1, pch = 16, cex = 1.5),
    list(cohort = "founders", region = "flanking", col = COL_FND, lty = 2, pch = 1,  cex = 1.5),
    list(cohort = "T2G07",    region = "hotspot",  col = COL_T2,  lty = 1, pch = 17, cex = 1.5),
    list(cohort = "T2G07",    region = "flanking", col = COL_T2,  lty = 2, pch = 2,  cex = 1.5)
  )
  for (s in series) {
    d <- dat[cohort == s$cohort & region == s$region]
    setorder(d, dist_bin)
    lines(x_at, d$r2, col = s$col, lty = s$lty, lwd = 2.2)
    points(x_at, d$r2, col = s$col, pch = s$pch, cex = s$cex, lwd = 2,
           bg = if (s$pch %in% c(1, 2)) "white" else s$col)
  }

  legend("topright", bty = "n", cex = 0.88,
         legend = c("Founders · hotspot", "Founders · flanking",
                    "T2G07 · hotspot",    "T2G07 · flanking"),
         col   = c(COL_FND, COL_FND, COL_T2, COL_T2),
         lty   = c(1, 2, 1, 2),
         pch   = c(16, 1, 17, 2),
         pt.cex = 1.2, lwd = 2)
}

ratios <- dcast(dat[, .(cohort, region, dist_bin, r2)],
                cohort + dist_bin ~ region, value.var = "r2")
ratios[, ratio := hotspot / flanking]
setorder(ratios, cohort, dist_bin)

panelB <- function() {
  par(mar = c(4.8, 4.8, 3, 1.2), mgp = c(2.9, 0.75, 0))
  bin_levels <- levels(dat$dist_bin)
  n_bins <- length(bin_levels)
  bar_w <- 0.38
  gap   <- 0.08
  positions <- cbind(
    (1:n_bins) - (bar_w/2 + gap/2),
    (1:n_bins) + (bar_w/2 + gap/2)
  )
  ymin <- 0.85; ymax <- 1.05

  plot(NA, NA, xlim = c(0.5, n_bins + 0.5), ylim = c(ymin, ymax),
       xaxt = "n", xlab = "Pair distance",
       ylab = "r² ratio (hotspot / flanking)",
       main = "B. Hotspot LD deficit closes in T2G07",
       las = 1, cex.main = 1.1, cex.lab = 1, cex.axis = 0.95)
  axis(1, at = 1:n_bins, labels = bin_levels, cex.axis = 0.95)

  abline(h = 1, col = COL_REF, lty = 3, lwd = 1.3)
  text(n_bins + 0.35, 1.005, "no deficit",
       cex = 0.72, col = COL_REF, pos = 2, font = 3)

  for (i in 1:n_bins) {
    rf <- ratios[cohort == "founders" & dist_bin == bin_levels[i]]$ratio
    rt <- ratios[cohort == "T2G07"    & dist_bin == bin_levels[i]]$ratio

    rect(positions[i,1] - bar_w/2, ymin,
         positions[i,1] + bar_w/2, rf,
         col = COL_FND, border = "grey20")

    rect(positions[i,2] - bar_w/2, ymin,
         positions[i,2] + bar_w/2, rt,
         col = COL_T2, border = "grey20")

    text(positions[i,1], rf + 0.008, sprintf("%.3f", rf),
         cex = 0.78, col = COL_FND)
    text(positions[i,2], rt + 0.008, sprintf("%.3f", rt),
         cex = 0.78, col = COL_T2)
  }

  legend("bottomright", bty = "n", cex = 0.88,
         legend = c("Founders", "T2G07"),
         fill = c(COL_FND, COL_T2), border = "grey20")
}

make_fig <- function() {
  layout(matrix(1:2, nrow = 1), widths = c(1, 1))
  panelA()
  panelB()
}

svglite(file.path(OUT_SVG, "figure5_ngsld_validation.svg"),
        width = 10, height = 5)
make_fig()
dev.off()

png(file.path(OUT_PNG, "figure5_ngsld_validation.png"),
    width = 2000, height = 1000, res = 200)
make_fig()
dev.off()

print(ratios)
wide <- dcast(dat[region == "hotspot"], dist_bin ~ cohort, value.var = "r2")
wide[, delta := T2G07 - founders]
wide[, pct := 100 * delta / founders]
print(wide)
