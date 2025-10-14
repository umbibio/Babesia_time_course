#!/usr/bin/env Rscript
# LOESS smoothing → one curve per gene (50 evaluation points)
# Robust to missing samples after noise removal:
# - averages tech reps (BE1/BE2, CK1/CK2) per time using available columns
# - averages biological reps (BE vs CK) per time using available values
# - for each gene, uses only times with non-NA to fit LOESS
# - skips genes with <3 usable timepoints and logs them
# - writes BdyValues1..50 and BdNormyValues1..50 + metrics

suppressPackageStartupMessages({
  library(openxlsx); library(dplyr); library(stringr); library(readr)
})

in_counts_clean <- "./tables/raw_counts_normal_growth_clean.xlsx"
out_loess_table <- "./tables/BdAllData.xlsx"
out_skip_log    <- "./tables/loess_skipped_genes.tsv"

dir.create(dirname(out_loess_table), showWarnings = FALSE, recursive = TRUE)

# ---- load cleaned counts ----
bdata <- read.xlsx(in_counts_clean)
stopifnot("GeneName" %in% names(bdata))
names(bdata)[names(bdata) == "GeneName"] <- "Gene"

# ---- helpers ----
.time_from_name <- function(x) suppressWarnings(as.numeric(stringr::str_extract(x, "(?<=-)(\\d+)$")))
times_target    <- c(0, 2, 4, 6, 8, 10, 12)

# pull columns by block label
sel_block <- function(df, tag) df %>% dplyr::select(dplyr::starts_with(tag))

# average the two tech blocks (A,B) per time *using available columns*
avg_two_blocks_by_time <- function(A, B, times) {
  rn <- rownames(A)
  out <- matrix(NA_real_, nrow = nrow(A), ncol = length(times),
                dimnames = list(rn, as.character(times)))
  # time extraction for each block
  tA <- .time_from_name(colnames(A)); tB <- .time_from_name(colnames(B))
  for (j in seq_along(times)) {
    t <- times[j]
    ia <- which(!is.na(tA) & tA == t)
    ib <- which(!is.na(tB) & tB == t)
    if (length(ia) + length(ib) == 0) next
    if (length(ia) > 0 && length(ib) > 0) {
      out[, j] <- rowMeans(cbind(A[, ia[1], drop = FALSE], B[, ib[1], drop = FALSE]), na.rm = TRUE)
    } else if (length(ia) > 0) {
      out[, j] <- as.numeric(A[, ia[1]])
    } else if (length(ib) > 0) {
      out[, j] <- as.numeric(B[, ib[1]])
    }
  }
  out
}

# ---- split blocks (whatever remains after cleaning) ----
BE1 <- sel_block(bdata, "BdC9-BE1")
BE2 <- sel_block(bdata, "BdC9-BE2")
CK1 <- sel_block(bdata, "BdC9-CK1")
CK2 <- sel_block(bdata, "BdC9-CK2")
rownames(BE1) <- rownames(BE2) <- rownames(CK1) <- rownames(CK2) <- bdata$Gene

# average tech reps per time (tolerant to missing columns)
BE_avg <- avg_two_blocks_by_time(BE1, BE2, times_target)   # genes x times
CK_avg <- avg_two_blocks_by_time(CK1, CK2, times_target)

# average biological reps per time using available values
# (if only one of BE/CK is present at a time, use it; if both present, mean)
CombinedAvg <- matrix(NA_real_, nrow = nrow(BE_avg), ncol = ncol(BE_avg),
                      dimnames = list(rownames(BE_avg), colnames(BE_avg)))
for (j in seq_len(ncol(BE_avg))) {
  CombinedAvg[, j] <- rowMeans(cbind(BE_avg[, j], CK_avg[, j]), na.rm = TRUE)
  # if both NA, rowMeans(..., na.rm=TRUE) returns NA — expected
}

# ---- per-gene LOESS with variable time availability ----
eval_n <- 50
n_genes <- nrow(CombinedAvg)

BdAllData <- matrix(NA_real_, nrow = n_genes, ncol = 6 + 2*eval_n) |> as.data.frame()
colnames(BdAllData)[1:6] <- c("Gene","maxTime","minTime","maxValue","minValue","FoldChange")
colnames(BdAllData)[7:(6+eval_n)]            <- paste0("BdyValues", seq_len(eval_n))
colnames(BdAllData)[(7+eval_n):(6+2*eval_n)] <- paste0("BdNormyValues", seq_len(eval_n))
BdAllData$Gene <- rownames(CombinedAvg)

skipped <- list()

for (i in seq_len(n_genes)) {
  y <- as.numeric(CombinedAvg[i, ])
  keep <- !is.na(y)
  t_avail <- times_target[keep]
  y_avail <- y[keep]
  
  # require at least 3 timepoints to fit loess
  if (length(y_avail) < 3) {
    skipped[[length(skipped) + 1]] <- data.frame(
      Gene = BdAllData$Gene[i],
      reason = sprintf("only %d usable timepoints", length(y_avail)),
      used_times = paste(t_avail, collapse = ", ")
    )
    next
  }
  
  # loess on available times; evaluate on an even grid between min/max available
  fit2 <- stats::loess.smooth(x = t_avail, y = y_avail,
                              span = 0.5, degree = 2, family = "gaussian",
                              evaluation = eval_n)
  
  maxTime    <- fit2$x[which.max(fit2$y)]
  minTime    <- fit2$x[which.min(fit2$y)]
  maxValue   <- max(fit2$y)
  minValue   <- min(fit2$y)
  FoldChange <- maxValue / max(minValue, .Machine$double.eps)
  
  yValues     <- fit2$y
  SumYvalues  <- sum(yValues)
  NormYvalues <- if (SumYvalues > 0) (yValues / SumYvalues) * 100 else rep(NA_real_, eval_n)
  
  BdAllData$maxTime[i]     <- maxTime
  BdAllData$minTime[i]     <- minTime
  BdAllData$maxValue[i]    <- maxValue
  BdAllData$minValue[i]    <- minValue
  BdAllData$FoldChange[i]  <- FoldChange
  BdAllData[i, 7:(6+eval_n)]            <- yValues
  BdAllData[i, (7+eval_n):(6+2*eval_n)] <- NormYvalues
}

# write outputs + skip log
openxlsx::write.xlsx(BdAllData, out_loess_table)

if (length(skipped)) {
  skipped_df <- dplyr::bind_rows(skipped)
  readr::write_tsv(skipped_df, out_skip_log)
  message("LOESS complete. Wrote: ", out_loess_table, 
          " | Skipped genes: ", nrow(skipped_df),
          " (logged to ", out_skip_log, ")")
} else {
  message("LOESS complete. Wrote: ", out_loess_table, " | No genes skipped.")
}
