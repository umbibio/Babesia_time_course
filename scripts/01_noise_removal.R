#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(openxlsx); library(edgeR); library(dplyr); library(stringr); library(readr)
})

in_counts <- "./input/raw_counts_normal_growth.xlsx"
in_design <- "./input/experiment_design.xlsx"
out_counts_clean <- "./tables/raw_counts_normal_growth_clean.xlsx"
out_design_clean <- "./tables/experiment_design_clean.xlsx"

dir.create(dirname(out_counts_clean), showWarnings = FALSE, recursive = TRUE)

detect_gene_col <- function(df) {
  if ("GeneName" %in% names(df)) "GeneName" else if ("Gene" %in% names(df)) "Gene" else names(df)[1]
}

set.seed(1234)
count  <- read.xlsx(in_counts)
design <- read.xlsx(in_design)

gene_col <- detect_gene_col(count)
sample_cols <- setdiff(names(count), gene_col)
stopifnot(all(c("Sample","treatment") %in% names(design)))

design <- design %>% filter(Sample %in% sample_cols) %>% arrange(match(Sample, sample_cols))
stopifnot(identical(design$Sample, sample_cols))

CPM  <- edgeR::cpm(as.matrix(count[, sample_cols, drop = FALSE]))
keep <- rowSums(CPM > 2) >= 3
x    <- count[keep, sample_cols, drop = FALSE]
rownames(x) <- count[[gene_col]][keep]

y <- edgeR::DGEList(counts = x, group = factor(design$treatment, levels = unique(design$treatment)))
y <- edgeR::calcNormFactors(y)
logCPM <- edgeR::cpm(y, log = TRUE, prior.count = 3, normalized.lib.sizes = TRUE)

hc <- hclust(dist(t(logCPM)))
clusterCut <- cutree(hc, k = 4)
cltab <- tibble(Sample = names(clusterCut), Batch = as.integer(clusterCut))
design2 <- left_join(design, cltab, by = "Sample")
design_clean <- design2 %>% filter(Batch != 3 & Batch != 4)

keep_samps <- design_clean$Sample
counts_clean <- count %>% select(all_of(c(gene_col, keep_samps)))
names(counts_clean)[names(counts_clean) == gene_col] <- "GeneName"

write.xlsx(counts_clean, out_counts_clean)
write.xlsx(design_clean, out_design_clean)
message("Wrote: ", out_counts_clean, " and ", out_design_clean)
