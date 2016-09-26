devtools::use_package('dplyr')

get_frequent_ppis <- function(pca_results,min_frequency,max_frequency=Inf){
  freq_table <- table(pca_results$PPI.short)
  passing_ppis <- names(which(freq_table >= min_frequency & freq_table <= max_frequency))
  return(passing_ppis)
}

get_genes_from_ppis <- function(pca_file,ppis){
  return(as.vector(unique(unlist(dplyr::filter(pca_file,PPI.short %in% ppis) %>% select(ORF.1,ORF.2)))))
}