# Babesia Time Course — Reproducible Analysis

This repository provides a **minimal, reviewer-ready** pipeline to reproduce:

1) **Noise removal / batch correction** (edgeR TMM; outlier clusters removed)  
2) **LOESS smoothing** → **one curve per gene** (average technical → average biological)  
3) **PCA** on cleaned counts (optional FC ≥ 1.5 filter; off by default)

The code is RStudio-friendly and runnable from the command line.

---

## Repository Structure

- input/ # <- Put the two Excel inputs here (not tracked by git)
- scripts/ # 01_noise_removal.R, 02_loess.R, 03_pca.R
- tables/ # outputs: cleaned tables, LOESS tables, PCA scores
- figures/ # outputs: PCA plot
- legacy/ # archived older scripts (not used by the pipeline)

## Requirements

- R (≥ 4.1 recommended)  
- Install these packages once:

```r
install.packages(c(
  "openxlsx","edgeR","dplyr","stringr","tidyr","ggplot2","ggrepel","readr", "tidyverse"
))
```

## How to Run (Command Line )

In a terminal:

```
git clone https://github.com/umbibio/Babesia_time_course.git
```

From the repository root (```Babesia_time_course```):

## STEP 1: Noise Removal / Batch Correction

```
Rscript scripts/01_noise_removal.R
```


### Outputs

- tables/raw_counts_normal_growth_clean.xlsx — cleaned counts (first col GeneName)

- tables/experiment_design_clean.xlsx — cleaned design with Batch (clusters 3 & 4 removed)

### Summary

- Filters low expression (CPM > 2 in ≥ 3 samples)

- edgeR TMM normalization → logCPM

- Hierarchical clustering on samples → cut into 4 clusters → drop clusters 3 & 4

- Standardizes gene column to GeneName



## STEP 2: LOESS Smoothing (One Curve per Gene)

```
Rscript scripts/02_loess.R
```

### Outputs

- tables/BdAllData.xlsx with columns:

- Gene, maxTime, minTime, maxValue, minValue, FoldChange

- BdyValues1..50 (smoothed curve values)

- BdNormyValues1..50 (same curve, normalized to sum = 100)

- If some genes had < 3 usable timepoints after cleaning:

### Summary

- Averages technical replicates within each biological replicate (BE1+BE2 → BE, CK1+CK2 → CK) per timepoint, using available columns

- Averages biological replicates (BE vs CK) per timepoint, using available values

- Fits one curve per gene with adaptive LOESS (more stable with few points);
- falls back to linear interpolation when needed

- Evaluates each curve at 50 points; keeps legacy column names for compatibility

- tables/loess_skipped_genes.tsv (gene ID, reason, times used)

## STEP 3: PCA on Cleaned Counts

```
Rscript scripts/03_pca.R
```

### Outputs

- tables/pca_scores.tsv — PC1/PC2 for each timepoint (after averaging reps)

- figures/pca_scores.pdf — publication-ready PCA plot

### Summary

- edgeR TMM → logCPM on cleaned counts

- Averages all replicates at each timepoint (BE/CK; 1/2)

- Optional amplitude filter across time (set min_fold_change <- 1.5 in the script to enable)

- PCA on timepoints (rows) × genes (columns)


## Running in RStudio (Optional)

- Open the repo in RStudio.

- Source & run, in order:

- scripts/01_noise_removal.R

- scripts/02_loess.R

- scripts/03_pca.R

- Inspect outputs under tables/ and figures/.


## Naming Conventions & Assumptions

- Sample names: BdC9-<BE|CK><rep>-<time> (e.g., BdC9-BE1-6)

- Time is parsed from the trailing integer (hours).

- LOESS uses only available timepoints per gene after cleaning.

- A minimum of 3 usable timepoints is required to fit a curve.

- Seeds are set where applicable for reproducibility.


## License
This repository is licensed under the [MIT License](./LICENSE).
