#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(poolSeq)
  library(data.table)
})

samples <- readLines("variance_analysis/sample_list.txt")
samples <- samples[nchar(samples) > 0]

ad <- fread("variance_analysis/merged_ad.tsv", header = FALSE, sep = "\t",
            colClasses = c(rep("character", 4), rep("character", 106)))
n_sites   <- nrow(ad)
n_samples <- length(samples)


parse_ad_col <- function(col) {
  col[is.na(col)] <- "."
  result <- matrix(0L, nrow = length(col), ncol = 2)
  good   <- col != "." & col != ".,."
  if (any(good)) {
    parts <- strsplit(col[good], ",", fixed = TRUE)
    refs  <- vapply(parts, function(p) as.integer(p[1]), integer(1))
    alts  <- vapply(parts, function(p) as.integer(p[2]), integer(1))
    result[good, 1] <- refs
    result[good, 2] <- alts
  }
  result
}


freq_mat  <- matrix(NA_real_, nrow = n_sites, ncol = n_samples)
depth_mat <- matrix(0L,       nrow = n_sites, ncol = n_samples)
for (i in seq_len(n_samples)) {
  parsed <- parse_ad_col(as.character(ad[[4 + i]]))
  total  <- parsed[, 1] + parsed[, 2]
  freq_mat[, i]  <- ifelse(total > 0, parsed[, 2] / total, NA_real_)
  depth_mat[, i] <- as.integer(total)
}
colnames(freq_mat)  <- samples
colnames(depth_mat) <- samples

chroms    <- as.character(ad[[1]])
positions <- as.integer(ad[[2]])

POOL_SIZES   <- c(80, 29, 24)
METHODS      <- c("P.planI", "P.planII")
TIME_WINDOWS <- list(
  list(name = "G01-G09", t1 = 1, t2 = 9),
  list(name = "G01-G06", t1 = 1, t2 = 6),
  list(name = "G02-G09", t1 = 2, t2 = 9)
)
WND_SIZE <- 1000   # SNPs per window (Barghi 2019)


ne_col_name <- function(method) {
  paste0("N", tolower(substr(method, 1, 1)), substring(method, 2))
}

run_one <- function(trt, rep, t1, t2, pool_size, method) {
  s1 <- sprintf("%s%dG%02d", trt, rep, t1)
  s2 <- sprintf("%s%dG%02d", trt, rep, t2)
  i1 <- which(samples == s1); i2 <- which(samples == s2)
  if (length(i1) == 0 || length(i2) == 0) return(NULL)

  p0   <- freq_mat[, i1];  pt <- freq_mat[, i2]
  cov0 <- depth_mat[, i1]; covt <- depth_mat[, i2]

  ok <- !is.na(p0) & !is.na(pt) & cov0 >= 10 & covt >= 10
  pbar <- (p0 + pt) / 2
  ok <- ok & pmin(pbar, 1 - pbar) >= 0.05
  if (sum(ok) < 5000) return(NULL)

  res <- tryCatch(
    estimateWndNe(
      chr      = chroms[ok],
      pos      = positions[ok],
      wndSize  = WND_SIZE,
      unit     = "SNP",
      p0       = p0[ok],   pt   = pt[ok],
      cov0     = as.integer(cov0[ok]),
      covt     = as.integer(covt[ok]),
      t        = t2 - t1,
      ploidy   = 2,
      method   = method,
      Ncensus  = 500,
      poolSize = c(pool_size, pool_size),
      truncAF  = 0.05
    ),
    error = function(e) { cat("    error:", e$message, "\n"); NULL }
  )
  if (is.null(res)) return(NULL)

  ne_col <- ne_col_name(method)
  if (!ne_col %in% names(res)) {
    cat("    expected column", ne_col, "not found; have:", paste(names(res), collapse=","), "\n")
    return(NULL)
  }
  ne_vals <- res[[ne_col]]
  ne_vals <- ne_vals[is.finite(ne_vals) & ne_vals > 0]
  if (length(ne_vals) < 10) return(NULL)

  data.table(
    treatment = trt, replicate = rep,
    time_window = sprintf("G%02d-G%02d", t1, t2),
    pool_size = pool_size,
    method = method,
    n_windows = length(ne_vals),
    Ne_median = median(ne_vals),
    Ne_mean = mean(ne_vals),
    Ne_q25 = quantile(ne_vals, 0.25, names = FALSE),
    Ne_q75 = quantile(ne_vals, 0.75, names = FALSE)
  )
}

results <- list()
total_configs <- length(POOL_SIZES) * length(METHODS) * length(TIME_WINDOWS)
cfg_n <- 0
for (ps in POOL_SIZES) {
  for (mth in METHODS) {
    for (tw in TIME_WINDOWS) {
      cfg_n <- cfg_n + 1
      cat(sprintf("\n[%d/%d] poolSize=%d  method=%s  window=%s\n",
                  cfg_n, total_configs, ps, mth, tw$name))
      for (trt in c("B", "T", "M")) {
        for (rep in 1:4) {
          r <- run_one(trt, rep, tw$t1, tw$t2, ps, mth)
          if (!is.null(r)) results[[length(results) + 1]] <- r
        }
      }
    }
  }
}

if (length(results) == 0) {
  stop("No results collected — every estimateWndNe call failed. Check column name / params.")
}
dt <- rbindlist(results)
fwrite(dt, "variance_analysis/section1_rigorous/poolseq_sweep_per_rep.tsv", sep = "\t")

summary_dt <- dt[, .(
  Ne_mean_of_medians = mean(Ne_median),
  Ne_min_rep = min(Ne_median),
  Ne_max_rep = max(Ne_median),
  n_reps = .N
), by = .(treatment, pool_size, method, time_window)]

fwrite(summary_dt, "variance_analysis/section1_rigorous/poolseq_sweep_summary.tsv",
       sep = "\t")


for (mth in METHODS) {
  for (tw_name in sapply(TIME_WINDOWS, function(x) x$name)) {
    sub <- summary_dt[method == mth & time_window == tw_name]
    if (nrow(sub) == 0) { cat("  (no results)\n\n"); next }
    wide <- dcast(sub, pool_size ~ treatment, value.var = "Ne_mean_of_medians")
    print(wide)
    cat("\n")
  }
}
