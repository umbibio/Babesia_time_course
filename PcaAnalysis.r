
#' pcaProjection
#' @description Dimension reduction using Principal Component Analysis.  
#' 
#' @usage pcaProjection(count.file, expmnt.file, Gene.list.file)
#' 
#' @param counts.file An excel file containing counts of each gene per sample.
#' 
#' @param expmnt.file An excel file containing information about each sample of experiment.
#' 
#' @param Gene.list.file An excel file containg genes of interest to run pca on.
#' 
#' @return A data frame containing the first 2 principal compents and visualization in the new 
#' coordinate system. 
#' 
#' @examples
#'
#' pca.example <- pcaProjection(count.file = count.file, expmnt.file = expmnt.file, Gene.list.file = Gene.list.file)
#'
#' note: you need to pass the clean version of counts and expmnt design file to this function 


abs.path <- "~/work/babesia/output/paper_material /"
count.file.clean <-  paste(abs.path, "raw_counts_normal_growth_clean_data.xlsx", sep = "")
expmnt.Info.clean <- paste(abs.path, "expmnt_Info_clean.xlsx", sep = "")
gtf.file <- paste(abs.path, "Bdiv.gff", sep = "")
Gene.list.file <- paste(abs.path, "BdALL 1.5 fold change.xlsx", sep = "")

count.file <- count.file.clean
expmnt.file <- expmnt.Info.clean
Gene.list.file <- Gene.list.file
gtf.file <- gtf.file


pcaProjection <- function(count.file, expmnt.file, Gene.list.file){
  
  count <- read.xlsx(count.file)
  expmnt <- read.xlsx(expmnt.file)
  genes <- read.xlsx(Gene.list.file)
  
  #GenesCount <- right_join(count, genes, by = "GeneName")
  
  y <- count 
  CPM <- cpm(y[2:ncol(y)])
  keep <- rowSums(CPM > 2) >= 3
  x <- y[keep, 2:ncol(y)]
  rownames(x) <- y$Gene[keep]
  
  treatment <- expmnt$treatment
  treatment <- factor(treatment, levels = unique(treatment))
  
  y <- DGEList(counts = x, group = treatment)
  y <- calcNormFactors(y)
  logCPM <- cpm(y, log=TRUE, prior.count=3, normalized.lib.sizes=TRUE)
  
  # average the normalized values of biological replicates per sample
  logCPM <- logCPM %>% data.frame() %>% rownames_to_column(var = "GeneName")
  
  logCPM.mean.sd <- logCPM %>% 
    gather(-GeneName, key = samples, value = expn) %>%
    mutate(bio_rep = ifelse(str_detect(samples, 'CK1'), 'CK1', "BE1"),
           time_point = paste("BdC9",gsub(".*\\.","", samples), sep = ".")) %>%
    group_by(GeneName, time_point) %>%
    summarise(mean = mean(expn) , sd = sd(expn)) 
  
  
  logCPM.mean <- logCPM.mean.sd %>% ungroup() %>%  mutate(Name = time_point) %>%
    select(-c('time_point', 'sd')) %>% spread(key=Name, value=mean)
  colnames(logCPM.mean)[2:ncol(logCPM.mean)] <- 
    paste(colnames(logCPM.mean)[2:ncol(logCPM.mean)], 'mean', sep = '_')
  
  logCPM.sd <- logCPM.mean.sd %>% ungroup() %>% mutate(Name = time_point) %>%
    select(-c('time_point', 'mean')) %>% spread(key=Name, value=sd)
  colnames(logCPM.sd)[2:ncol(logCPM.sd)] <- 
    paste(colnames(logCPM.sd)[2:ncol(logCPM.sd)], 'sd', sep = '_')
  
  logCPM.mean.sd <- left_join(logCPM.mean, logCPM.sd, by = 'GeneName')
  final.logCPM.mean.sd <- left_join(logCPM, logCPM.mean.sd, by = 'GeneName')
  
  
  # selecting the mean values
  tab <- final.logCPM.mean.sd %>% select(contains("mean")) %>% 
    mutate(GeneName = final.logCPM.mean.sd$GeneName) %>% 
    select(GeneName, everything())
  
  # intersection of all genes with genes with FC > 1.5
  tab <- left_join(genes, tab, by = "GeneName") %>% na.omit()
  
  # pca projection
  df <- data.frame(t(tab[, 2:ncol(tab)]))
  colnames(df) <- tab$GeneName
  df <- df %>% data.frame() %>% rownames_to_column("Sample")
  PCA <- prcomp(df[,2:ncol(df)], scale. = T, center = TRUE)
  
  pca.var <- PCA$sdev ^ 2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  
  pca.data <- data.frame(Sample = paste(gsub(".*\\.","", gsub("_.*", "", df$Sample)), "h", sep = ""),
                         PC1 = PCA$x[,1],
                         PC2 = PCA$x[,2])
  
  pca.data$Sample <- factor(pca.data$Sample, levels = unique(pca.data$Sample))
  
  pca.plot <- ggplot(data = pca.data, mapping = aes(x = PC1, y = PC2, 
                                                    color = Sample)) +
    geom_point(size =3) +
    geom_text(aes(label=Sample),hjust=0.5, vjust=-0.75, show.legend = FALSE)+
    xlab(paste("PC1", paste(pca.var.per[1], "%", sep = ""), sep = "_"))+
    ylab(paste("PC2", paste(pca.var.per[2], "%", sep = ""), sep = "_"))+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    guides(col=guide_legend("Time Points"))+
    theme(plot.title = element_text(size = 10),
          legend.title = element_blank(),
          legend.position = "none",
          axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"))
   
  pca.plot <-  pca.plot +  ggtitle("PCA Projection")
  
  
  return(list(pca.data = pca.data, pca.plot = pca.plot))
  
}

pca.example <- pcaProjection(count.file = count.file, expmnt.file = expmnt.file, Gene.list.file = Gene.list.file)
write.xlsx(pca.example$pca.data, "~/work/babesia/output/paper_material /pca_table.xlsx")
