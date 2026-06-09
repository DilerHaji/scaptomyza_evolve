#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

ld <- fread(cmd = "zcat ngsld/ld/founders/chr_ScDA7r2_439_HRSCAF_779.ld.gz")
ld[, pos1 := as.integer(sub(".*:", "", site1))]
ld[, pos2 := as.integer(sub(".*:", "", site2))]

ld[, in_hotspot := pos1 >= 2e6 & pos1 <= 12e6 &
                    pos2 >= 2e6 & pos2 <= 12e6]

ld[, dist_bin := cut(dist,
                      breaks = c(0, 1e3, 1e4, 1e5),
                      labels = c("<1kb", "1-10kb", "10-100kb"),
                      include.lowest = TRUE)]

summ <- ld[, .(n = .N, mean_r2 = mean(r2_ExpG, na.rm = TRUE)),
            by = .(in_hotspot, dist_bin)]