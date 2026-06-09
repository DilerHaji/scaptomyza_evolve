#!/usr/bin/env Rscript
suppressMessages(library(data.table))
suppressMessages(library(arrow))
suppressMessages(library(argparser))
suppressMessages(library(dplyr))

p <- arg_parser("Convert HAF data to Arrow for Julia")
p <- add_argument(p, "--haf", help="Path to merged_gethaf.Rdata")
p <- add_argument(p, "--rd", help="Path to merged_gethaf_neff.RDS")
p <- add_argument(p, "--chunk", help="Chunk TSV (to filter sites for test)")
p <- add_argument(p, "--repName", help="Random effect column (e.g. pop)")
p <- add_argument(p, "--selectTrts", help="Treatments (e.g. B,T)")
p <- add_argument(p, "--out", help="Output .arrow file")

p <- add_argument(p, "--min-coverage", type="integer", default=10,
                  help="Minimum coverage threshold (default: 10)")
p <- add_argument(p, "--min-samples", type="integer", default=4,
                  help="Minimum samples with variant present (default: 4)")
p <- add_argument(p, "--min-distinct-freqs", type="integer", default=2,
                  help="Minimum distinct frequency values (default: 2)")

opts <- parse_args(p)

load(opts$haf) # afmat, samps, sites
rd_full <- readRDS(opts$rd)
if (!file.exists(opts$chunk)) stop("Chunk file not found: ", opts$chunk)
cands <- fread(opts$chunk)

colnames(cands)[grep("^chrom$", colnames(cands), ignore.case=TRUE)] <- "chrom"
colnames(cands)[grep("^pos$", colnames(cands), ignore.case=TRUE)] <- "pos"

if (!all(c("chrom", "pos") %in% colnames(cands))) {
  stop("Chunk file missing 'chrom' or 'pos' columns after standardization.")
}

cands$chrom <- trimws(as.character(cands$chrom))
cands$pos   <- as.integer(cands$pos)
sites$chrom <- trimws(as.character(sites$chrom))
sites$pos   <- as.integer(sites$pos)

if (!(opts$repName %in% colnames(samps))) stop("repName not found in samps.")
sampIX <- which(!is.na(samps[[opts$repName]]))

if (!is.na(opts$selectTrts)) {
  trts <- trimws(unlist(strsplit(opts$selectTrts, ",")))
  sampIX <- intersect(sampIX, which(samps$trt %in% trts))
}
samps_sub <- samps[sampIX, ]

site_keys <- paste(sites$chrom, sites$pos, sep="_")
cand_keys <- paste(cands$chrom, cands$pos, sep="_")
match_indices <- match(cand_keys, site_keys)

found_indices <- match_indices[!is.na(match_indices)]

if (length(found_indices) == 0) {
  stop("No matching sites found! Check chromosome naming conventions.")
}

data_list <- lapply(found_indices, function(i) {
  af <- afmat[i, sampIX]
  rd <- rd_full[i, sampIX]
  succ <- round(rd * af)
  fail <- round(rd * (1 - af))
  df <- samps_sub[, c(opts$repName, "gen", "trt")]
  colnames(df)[1] <- "pop" # Standardization for Julia
  df$chrom <- as.character(sites$chrom[i])
  df$pos   <- as.integer(sites$pos[i])
  df$success <- as.integer(succ)
  df$total   <- as.integer(succ + fail)
  return(df)
})

long_df <- bind_rows(data_list)

initial_sites <- length(unique(paste(long_df$chrom, long_df$pos)))
initial_rows <- nrow(long_df)

long_df <- long_df %>%
  group_by(chrom, pos) %>%
  filter(
    all(total >= opts$min_coverage),
    sum(success > 0) >= opts$min_samples,
    n_distinct(success / total) >= opts$min_distinct_freqs
  ) %>%
  ungroup()

final_sites <- length(unique(paste(long_df$chrom, long_df$pos)))
final_rows <- nrow(long_df)
removed_sites <- initial_sites - final_sites
removal_pct <- 100 * removed_sites / initial_sites

long_df$pop <- as.factor(long_df$pop)
long_df$trt <- as.factor(long_df$trt)
long_df$chrom <- as.factor(long_df$chrom)

write_feather(long_df, opts$out)
