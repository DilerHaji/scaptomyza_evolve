#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(haploReconstruct)
  library(haplovalidate)
  library(ACER)
})

args <- commandArgs(trailingOnly = TRUE)
pa <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  args[idx + 1]
}

sync_file   <- pa("--sync")
ad_file     <- pa("--ad", "variance_analysis/merged_ad.tsv")
depth_file  <- pa("--depth", "variance_analysis/merged_depth.tsv")
sample_file <- pa("--sample_list", "variance_analysis/sample_list.txt")
treatment   <- pa("--treatment", "B")
Ne_val      <- as.integer(pa("--ne", "250"))
ps_val      <- as.integer(pa("--poolsize", "160"))
n_repl      <- as.integer(pa("--n_repl", "4"))
outdir      <- pa("--outdir", paste0("hv_results/hv_", treatment))
scaff_str   <- pa("--scaffolds", "ALL")
max_cands   <- as.integer(pa("--max_cands_per_scaffold", "50000"))
maf_min     <- as.numeric(pa("--maf_min", "0.05"))
min_dp      <- as.integer(pa("--min_depth", "20"))
mncs_val    <- as.numeric(pa("--mncs", "0.01"))
takerandom  <- as.integer(pa("--takerandom", "2000"))
filterrange <- as.integer(pa("--filterrange", "5000"))
score_thresh <- as.numeric(pa("--score_threshold", "1.3"))
use_bh       <- !is.null(pa("--use_bh"))  # if flag present, filter on BH q < 0.05 instead of raw score

acer_samples_str <- pa("--acer_samples")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

gens <- c(0, 1, 2, 6, 7, 8, 9)
n_gens <- length(gens)
n_samples <- n_repl * n_gens

base.pops <- c(rep(TRUE, n_repl), rep(FALSE, n_repl * (n_gens - 1)))
pop.ident <- rep(1:n_repl, n_gens)
pop.generation <- rep(gens, each = n_repl)
compare <- rep(TRUE, n_samples)
last_gen_start <- n_repl * (n_gens - 1) + 1
polaRise <- lapply(1:n_repl, function(r) c(r, last_gen_start + r - 1))

t0 <- proc.time()
cands.all <- sync_to_frequencies(sync_file, base.pops = base.pops,
                                  header = FALSE, mincov = 15, polaRise = polaRise)
cands.all[, pos := as.integer(pos)]


all_samples <- readLines(sample_file)
all_samples <- all_samples[nzchar(all_samples)]


if (!is.null(acer_samples_str)) {
  acer_names <- strsplit(acer_samples_str, ",")[[1]]
} else {
  trt_letters <- strsplit(treatment, "")[[1]]
  if (length(trt_letters) == 1) {
    acer_names <- character(0)
    for (r in 1:n_repl) {
      acer_names <- c(acer_names, sprintf("F%dG00", r))
      for (g in gens[-1]) acer_names <- c(acer_names, sprintf("%s%dG%02d", treatment, r, g))
    }
  } else {
    stop("For multi-treatment runs, provide --acer_samples explicitly")
  }
}
si <- match(acer_names, all_samples)
if (any(is.na(si))) stop("Missing samples: ", paste(acer_names[is.na(si)], collapse = ", "))

founder_si <- match(c("F1G00", "F2G00", "F3G00", "F4G00"), all_samples)

if (scaff_str != "ALL") {
  target_scaffolds <- strsplit(scaff_str, ",")[[1]]
} else {
  target_scaffolds <- NULL  # use all
}

t0 <- proc.time()
depth_raw <- fread(depth_file, header = FALSE, sep = "\t")
if (!is.null(target_scaffolds)) depth_raw <- depth_raw[V1 %in% target_scaffolds]
cm <- as.matrix(depth_raw[, .SD, .SDcols = si + 2])


t0 <- proc.time()
ad_raw <- fread(ad_file, header = FALSE, sep = "\t", colClasses = "character")
if (!is.null(target_scaffolds)) ad_raw <- ad_raw[V1 %in% target_scaffolds]

fm <- matrix(NA_real_, nrow(ad_raw), length(si))
for (j in seq_along(si)) {
  parts <- strsplit(ad_raw[[si[j] + 4]], ",", fixed = TRUE)
  rc <- as.integer(sapply(parts, `[`, 1))
  ac <- as.integer(sapply(parts, `[`, 2))
  tot <- rc + ac
  fm[, j] <- ifelse(tot > 0, ac / tot, NA_real_)
}

fr <- fa <- matrix(0L, nrow(ad_raw), 4)
for (j in 1:4) {
  parts <- strsplit(ad_raw[[founder_si[j] + 4]], ",", fixed = TRUE)
  fr[, j] <- as.integer(sapply(parts, `[`, 1))
  fa[, j] <- as.integer(sapply(parts, `[`, 2))
}
founder_maf <- pmin(rowSums(fa), rowSums(fr)) / (rowSums(fr) + rowSums(fa))
min_depth <- apply(cm, 1, min)

keep <- founder_maf >= maf_min & min_depth >= min_dp


fm_f <- fm[keep, ]
cm_f <- cm[keep, ]
pos_f <- data.table(chr = ad_raw$V1[keep], pos = as.integer(ad_raw$V2[keep]))
cands.all <- merge(cands.all, pos_f, by = c("chr", "pos"))

t0 <- proc.time()
acer_result <- suppressWarnings(adapted.cmh.test(
  freq = fm_f, coverage = cm_f,
  Ne = rep(Ne_val, n_repl), gen = gens, repl = 1:n_repl,
  poolSize = rep(ps_val, length(acer_names)),
  mincov = min_dp, MeanStart = TRUE, IntGen = TRUE, TA = FALSE,
  order = 0, RetVal = 2
))
pv <- acer_result[, 2]
qv <- p.adjust(pv, "BH")

cmh <- data.table(chr = pos_f$chr, pos = pos_f$pos, score = -log10(pv), qval = qv)
cmh <- cmh[!is.na(score) & is.finite(score)]
fwrite(cmh, file.path(outdir, "acer_scores.tsv"), sep = "\t")

if (use_bh) {
  cmh_sig <- cmh[qval < 0.05]
} else {
  cmh_sig <- cmh[score > score_thresh]
}

cmh_capped <- cmh_sig[, {
  if (.N > max_cands) .SD[order(-score)][1:max_cands] else .SD
}, by = chr]

cands <- merge(cands.all, cmh_capped[, .(chr, pos)], by = c("chr", "pos"))

parameters <- tryCatch(
  get.mncs.win(cands, cmh, wins = seq(0.1, 10, 0.05), mncs = mncs_val),
  error = function(e) {
    if (mncs_val < 0.03) {
      tryCatch(get.mncs.win(cands, cmh, wins = seq(0.1, 10, 0.05), mncs = 0.03),
               error = function(e2) { cat("  Fallback also failed\n"); NULL })
    } else { NULL }
  }
)

t0 <- proc.time()
happy <- tryCatch({
  invisible(capture.output(
    res <- haplovalidate(
      cands = cands, cmh = cmh, parameters = parameters,
      pop.ident = pop.ident, pop.generation = pop.generation,
      base.pops = base.pops, compare = compare,
      takerandom = takerandom, filterrange = filterrange
    ), type = "output"))
  res
}, error = function(e) {
  NULL
})
elapsed <- round((proc.time() - t0)[3])

if (!is.null(happy) && length(happy) > 0) {
  dh <- happy$dominant_haplotypes
  ah <- happy$all_haplotypes
  n_blocks <- uniqueN(dh$tag)
  n_nt <- sum(grepl("noThreshold", unique(dh$tag)))
  dh[, pos := as.integer(pos)]
  summ <- dh[, .(blocks = uniqueN(tag), snps = .N,
                  noThresh = sum(grepl("noThreshold", unique(tag))),
                  min_pos = min(pos), max_pos = max(pos)), by = chr]
  setorder(summ, chr)

  block_table <- dh[, .(chr = chr[1],
                         start = min(as.integer(pos)),
                         end = max(as.integer(pos)),
                         n_snps = .N,
                         span_Mb = round((max(as.integer(pos)) - min(as.integer(pos))) / 1e6, 2),
                         noThreshold = any(grepl("noThreshold", tag))), by = tag]
  setorder(block_table, chr, start)
  fwrite(block_table, file.path(outdir, "dominant_blocks.tsv"), sep = "\t")
} else {
  empty_blocks <- data.table(tag = character(), chr = character(), start = integer(),
                              end = integer(), n_snps = integer(), span_Mb = numeric(),
                              noThreshold = logical())
  fwrite(empty_blocks, file.path(outdir, "dominant_blocks.tsv"), sep = "\t")
}

saveRDS(happy, file.path(outdir, "haplovalidate_result.rds"))
saveRDS(cands, file.path(outdir, "cands.rds"))

sink(file.path(outdir, "summary.txt"))
cat("Treatment:", treatment, "\n")
cat("Ne:", Ne_val, "\n")
cat("poolSize:", ps_val, "\n")
cat("MAF_min:", maf_min, "\n")
cat("min_depth:", min_dp, "\n")
cat("max_cands_per_scaffold:", max_cands, "\n")
cat("mncs:", mncs_val, "\n")
cat("takerandom:", takerandom, "\n")
cat("filterrange:", filterrange, "\n")
cat("score_threshold:", score_thresh, "\n")
cat("total_sites:", nrow(cmh), "\n")
cat("p05:", sum(pv < 0.05, na.rm = TRUE), "\n")
cat("bh05:", sum(qv < 0.05, na.rm = TRUE), "\n")
cat("candidates:", nrow(cmh_capped), "\n")
cat("dominant_blocks:", if (!is.null(happy$dominant_haplotypes)) uniqueN(happy$dominant_haplotypes$tag) else 0, "\n")
cat("noThreshold_blocks:", if (!is.null(happy$dominant_haplotypes)) sum(grepl("noThreshold", unique(happy$dominant_haplotypes$tag))) else 0, "\n")
cat("runtime_sec:", elapsed, "\n")
sink()