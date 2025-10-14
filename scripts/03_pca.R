#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(openxlsx); library(edgeR); library(dplyr); library(stringr)
  library(tidyr); library(ggplot2); library(ggrepel); library(readr)
})

in_counts_clean <- "./tables/raw_counts_normal_growth_clean.xlsx"
in_design_clean <- "./tables/experiment_design_clean.xlsx"
out_scores <- "./tables/pca_scores.tsv"
out_plot   <- "./figures/pca_scores.pdf"
dir.create(dirname(out_scores), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_plot),   showWarnings = FALSE, recursive = TRUE)

min_fold_change <- NULL  # set to 1.5 to apply FC≥1.5 filter across time

time_from_name <- function(x) suppressWarnings(as.numeric(stringr::str_extract(x, "(?<=-)(\\d+)$")))

count  <- read.xlsx(in_counts_clean)
design <- read.xlsx(in_design_clean)
stopifnot("GeneName" %in% names(count))
sample_cols <- setdiff(names(count), "GeneName")

design <- design %>% filter(Sample %in% sample_cols) %>% arrange(match(Sample, sample_cols))
stopifnot(identical(design$Sample, sample_cols))

x <- as.matrix(count[, sample_cols, drop = FALSE])
rownames(x) <- count$GeneName
y <- edgeR::DGEList(counts = x, group = factor(design$treatment, levels = unique(design$treatment)))
y <- edgeR::calcNormFactors(y)
logCPM <- edgeR::cpm(y, log = TRUE, prior.count = 3, normalized.lib.sizes = TRUE)

logCPM_df <- as.data.frame(logCPM) %>% tibble::rownames_to_column("GeneName")
long <- logCPM_df %>% pivot_longer(-GeneName, names_to = "Sample", values_to = "expn") %>%
  mutate(time = time_from_name(Sample)) %>%
  filter(!is.na(time))

avg_by_time <- long %>% group_by(GeneName, time) %>%
  summarise(mean = mean(expn), .groups = "drop") %>%
  pivot_wider(names_from = time, values_from = mean, names_sort = TRUE)

if (!is.null(min_fold_change)) {
  log2_thresh <- log2(min_fold_change)
  time_cols <- setdiff(names(avg_by_time), "GeneName")
  avg_by_time <- avg_by_time %>%
    mutate(log2_range = apply(across(all_of(time_cols)), 1, function(v) max(v, na.rm = TRUE) - min(v, na.rm = TRUE))) %>%
    filter(log2_range >= log2_thresh) %>%
    select(-log2_range)
}

if (ncol(avg_by_time) < 3) stop("Not enough timepoints for PCA.")
if (nrow(avg_by_time) < 2) stop("Not enough genes for PCA.")

mat <- as.data.frame(avg_by_time[, -1, drop = FALSE]); rownames(mat) <- avg_by_time$GeneName
df <- t(mat) %>% as.data.frame(); df$Sample <- rownames(df)

PCA <- prcomp(df[, setdiff(names(df), "Sample"), drop = FALSE], scale. = TRUE, center = TRUE)
pvar <- round((PCA$sdev^2) / sum(PCA$sdev^2) * 100, 1)

pca.data <- tibble(Sample = paste0(df$Sample, "h"), PC1 = PCA$x[,1], PC2 = PCA$x[,2])
pca.data$Sample <- factor(pca.data$Sample, levels = unique(pca.data$Sample))

p <- ggplot(pca.data, aes(PC1, PC2, color = Sample)) +
  geom_point(size = 4, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = Sample), size = 4, fontface = "bold",
                           max.overlaps = Inf, box.padding = 0.4, show.legend = FALSE) +
  xlab(paste0("PC1 (", pvar[1], "%)")) + ylab(paste0("PC2 (", pvar[2], "%)")) +
  ggtitle("PCA Projection of Babesia Growth Time Course") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_color_brewer(palette = "Set2")

write_tsv(pca.data, out_scores)
ggsave(out_plot, p, width = 5, height = 4, dpi = 300)
message("Wrote: ", out_scores, " and ", out_plot)
