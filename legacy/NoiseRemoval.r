
#' removeBatch
#' @description removes the noisy samples from raw counts table. 
#' 
#' @usage removeBatch(counts.file, expmnt.file)
#' 
#' @param counts.file An excel file containing counts of each gene per sample.
#' 
#' @param expmnt.file An excel file containing information about each sample of experiment.
#' 
#' @return Two data frames of raw counts and experiment design without noisy samples.
#'
#' @example
#'
#' counts <- removeBatch(count.file, expmnt.file)

removeBatch <- function(count.file, expmnt.file){
  
  count <- read.xlsx(count.file)
  expmnt <- read.xlsx(expmnt.file)
  
  # cpm normalization
  CPM <- cpm(count[2:ncol(count)]) 
  
  # filter low expressed genes
  keep <- rowSums(CPM > 2) >= 3     
  x <- count[keep, 2:ncol(count)]
  rownames(x) <- count$Gene[keep]
  treatment <- expmnt$treatment
  treatment <- factor(treatment, levels = unique(treatment))
  
  # Creating DGElist object 
  y <- DGEList(counts = x, group = treatment)  
  y <- calcNormFactors(y)
  # logarithmic scale
  logCPM <- cpm(y, log=TRUE, prior.count=3, normalized.lib.sizes=TRUE)  
  
  # Detecting batches/noisy samples
  clusters <- hclust(dist(t(logCPM)))
  clusterCut <- cutree(clusters, 4)
  CC <- data.frame(Sample = names(clusterCut), Batch = clusterCut, stringsAsFactors = F)
  expmnt <- left_join(expmnt, CC, by = "Sample")
  expmnt.cleaned <- expmnt %>% dplyr::filter(Batch != 3 & Batch !=4)
  
  count.cleaned <- count %>% 
    select(c("GeneName",names(count)[colnames(count) %in% expmnt$Sample]))
  
  
  return(list(counts = count.cleaned, expmnt = expmnt.cleaned))
  
}



abs.path <- "~/work/babesia/output/paper_material /"
count.file <-  paste(abs.path, "raw_counts_normal_growth.xlsx", sep = "")
expmnt.file <- paste(abs.path, "experiment_design.xlsx", sep = "")
gtf.file <- paste(abs.path, "Bdiv.gff", sep = "")
Gene.list.file <- paste(abs.path, "BdALL 1.5 fold change.xlsx", sep = "")

example <- removeBatch(count.file, expmnt.file)



