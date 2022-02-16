#' orthologousGenes
#' @description performing gene orthology inference using the reciprocal best hit (RBH) method.  
#' 
#' @usage orthologousGenes(query.fasta, subject.fasta)
#' 
#' @param query.file Path to the fasta file of interest (query organism)
#' 
#' @param subject.file Path to the fasta file of interest (subject organism).
#' 
#' @return A data frame containing the gene ids of orthologous genes identified by RBH method.
#' 
#' @examples
#'
#' rec_B.div_P.berg <- prthologousGenes(quary.fasta = Bdivergens.fasta, subject.fasta = Pberghei.fasta)
#' 


orthologousGenes <- function(quary.fasta, subject.fasta){
  
  reciprocal_best_hit <- blast_rec(query_file   = quary.fasta,
                                   subject_file = subject.fasta,
                                   delete_corrupt_cds = T, seq_type = "protein",
                                   format = "fasta", blast_algorithm = "blastp",
                                   eval = 0.0001)
  
  return(reciprocal_best_hit)
  
} 

