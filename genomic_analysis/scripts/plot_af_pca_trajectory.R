#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(dplyr)
}))

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); opt[[key]] <- args[i + 1]; i <- i + 2
}
if (is.null(opt[["scores"]]) || is.null(opt[["eigenvals"]]) || is.null(opt[["out-prefix"]])) {
  stop("Usage: Rscript plot_af_pca_trajectory.R --scores SCORES.csv --eigenvals EIG.csv --out-prefix PREFIX")
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

safe_colorblind_palette <- c(
  Barbarea_experimental = "#96CDFF",
  Founders              = "#B1B1B1",
  Mixed_experimental    = "#901442",
  Turitus_experimental  = "#F3C43C"
)

make_plot <- function(xcol, ycol, xlab, ylab) {
  ggplot(pca, aes(.data[[xcol]], .data[[ycol]], col = spp, label = ind)) +
    geom_point(size = 1.6) +
    scale_color_manual(values = safe_colorblind_palette) +
    coord_equal() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab(xlab) + ylab(ylab) +
    geom_path(aes(group = interaction(rp, treat)),
              arrow = arrow(type = "closed", angle = 30, length = unit(0.10, "inches")),
              linewidth = 0.8, linetype = "solid", alpha = 0.6) +
    geom_text_repel(max.overlaps = 20, force = 6,
                    segment.linetype = 2, box.padding = unit(0.8, "lines"),
                    size = 2.5)
}

g1 <- make_plot("PC1", "PC2",
                sprintf("PC1 (%.1f%%)", pve$pve[1]),
                sprintf("PC2 (%.1f%%)", pve$pve[2]))
g2 <- make_plot("PC1", "PC3",
                sprintf("PC1 (%.1f%%)", pve$pve[1]),
                sprintf("PC3 (%.1f%%)", pve$pve[3]))

ggsave(sprintf("%s_PC1_PC2.png", opt[["out-prefix"]]), g1, width = 9, height = 7, dpi = 200)
ggsave(sprintf("%s_PC1_PC3.png", opt[["out-prefix"]]), g2, width = 9, height = 7, dpi = 200)