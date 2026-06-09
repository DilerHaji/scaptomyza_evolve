#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({
  library(readr)
  library(dplyr)
  if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    stop("Package 'scatterplot3d' is required. Install with: install.packages('scatterplot3d')")
  }
  library(scatterplot3d)
}))

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); opt[[key]] <- args[i + 1]; i <- i + 2
}
if (is.null(opt[["scores"]]) || is.null(opt[["eigenvals"]]) || is.null(opt[["out-prefix"]])) {
  stop("Usage: Rscript plot_af_pca_3d.R --scores SCORES.csv --eigenvals EIG.csv --out-prefix PREFIX")
}

pca <- read_csv(opt[["scores"]], show_col_types = FALSE)
pve <- read_csv(opt[["eigenvals"]], show_col_types = FALSE)

pca$spp <- NA_character_
pca$spp[grepl("^B[1-4]G",       pca$ind)] <- "Barbarea_experimental"
pca$spp[grepl("^T[1-4]G",       pca$ind)] <- "Turitus_experimental"
pca$spp[grepl("^M[1-4]G",       pca$ind)] <- "Mixed_experimental"
pca$spp[grepl("^F[1-4](G00)?$", pca$ind)] <- "Founders"

pca$gen <- NA_integer_
for (g in 0:10) pca$gen[grepl(sprintf("G%02d", g), pca$ind)] <- g
pca$gen[grepl("^F[1-4]$", pca$ind)] <- 0

pca$rp <- NA_integer_
for (p in 1:4) {
  pca$rp[grepl(sprintf("^[BTM]%d(G|$)", p), pca$ind)] <- p
  pca$rp[grepl(sprintf("^F%d(G00)?$",  p), pca$ind)] <- p
}

pca <- pca[!is.na(pca$spp), ]
pca$treat <- pca$spp
pca <- pca[order(pca$rp, pca$gen), ]

palette_trt <- c(
  Barbarea_experimental = "#96CDFF",
  Founders              = "#B1B1B1",
  Mixed_experimental    = "#901442",
  Turitus_experimental  = "#F3C43C"
)
pca$pt_color <- palette_trt[pca$spp]

make_3d <- function(angle, out_path, pch = 19, cex = 0.9) {
  png(out_path, width = 2200, height = 2000, res = 220)
  par(mar = c(3, 3, 2, 3))
  s3d <- scatterplot3d(
    x = pca$PC1, y = pca$PC2, z = pca$PC3,
    color = pca$pt_color, pch = pch,
    cex.symbols = cex,
    angle = angle, scale.y = 0.8, box = TRUE,
    xlab = sprintf("PC1 (%.1f%%)", pve$pve[1]),
    ylab = sprintf("PC2 (%.1f%%)", pve$pve[2]),
    zlab = sprintf("PC3 (%.1f%%)", pve$pve[3]),
    main = sprintf("AF-PCA (3D), angle = %d deg", angle)
  )
  for (grp in split(pca, list(pca$rp, pca$treat), drop = TRUE)) {
    if (nrow(grp) < 2) next
    grp <- grp[order(grp$gen), ]
    xy <- s3d$xyz.convert(grp$PC1, grp$PC2, grp$PC3)
    lines(xy$x, xy$y, col = grp$pt_color[1], lwd = 1.2)
    n <- length(xy$x)
    if (n >= 2) {
      arrows(xy$x[n-1], xy$y[n-1], xy$x[n], xy$y[n],
             col = grp$pt_color[1], length = 0.07, angle = 25, lwd = 1.2)
    }
  }
  legend("topright", legend = names(palette_trt),
         col = palette_trt, pch = pch, pt.cex = 1, cex = 0.7,
         bty = "n", inset = c(-0.02, 0.02), xpd = NA)
  dev.off()
}

for (a in c(20, 40, 60, 80)) {
  out_path <- sprintf("%s_angle%d.png", opt[["out-prefix"]], a)
  make_3d(a, out_path)
}

cat("done.\n")
