library(knitr)

freq_perturb_output_path <- paste(c(output_path,'freq_peturb'),collapse='/')
dir.create(freq_perturb_output_path, showWarnings = FALSE)

pca_univ <- read.table(pca_universe,head=T,sep='\t')
pca_enh <- read.table(pca_enhanced_calls,head=T,sep='\t')
pca_depl <- read.table(pca_depleted_calls,head=T,sep='\t')

pca_perturbed <- rbind(pca_enh,pca_depl)

frequent_ppis <- get_frequent_ppis(pca_perturbed,4)
frequent_genes <- get_genes_from_ppis(pca_perturbed,frequent_ppis)

go_lines <- readLines(system.file('go_map_sgd.txt',package='bcPcaAnalysis'))
membrane_go_line <- go_lines[grep('GO:0005886',go_lines)]
membrane_genes <- map_gene_names(strsplit(strsplit(membrane_go_line,split='\t')[[1]][3],split=' ')[[1]])
frequent_membrane_ppis <- c()
for(i in 1:length(frequent_ppis)){
  ppi <- frequent_ppis[i]
  genes <- strsplit(ppi,split=':')[[1]]
  if(genes[1] %in% membrane_genes & genes[2] %in% membrane_genes){
    frequent_membrane_ppis <- c(frequent_membrane_ppis,ppi)
  }
}

frequent_membrane_genes <- get_genes_from_ppis(pca_perturbed,frequent_membrane_ppis)

ppis_with_frequent_genes <- dplyr::filter(pca_perturbed,ORF.1 %in% frequent_membrane_genes & ORF.2 %in% frequent_membrane_genes)
all_ppis_with_frequent_genes <- dplyr::filter(pca_univ,ORF.1 %in% frequent_membrane_genes & ORF.2 %in% frequent_membrane_genes)
#
ppi_freq_table <- table(ppis_with_frequent_genes$PPI.short)

#ppi_freq_table_2

output_matrix <- t(sapply(as.vector(unique(all_ppis_with_frequent_genes$PPI.short)),function(ppi){
  names <- strsplit(ppi,split=':')[[1]]
  if(is.na(ppi_freq_table[ppi])){
    return(c(names,0))
  }else{
    return(c(names,ppi_freq_table[ppi]))
  }
}))
colnames(output_matrix) <- c('Protein1','Protein2','Frequency')
write.table(output_matrix,row.names=F,quote=F,sep='\t',file=paste(c(freq_perturb_output_path,'ppis_with_frequeny_genes.tsv'),collapse='/'))

all_ppis <- get_frequent_ppis(pca_univ,0)
all_genes <- get_genes_from_ppis(pca_univ,all_ppis)

go_enrichment <- funcassociate(frequent_genes,all_genes,order_mode='unordered')
outfile <- paste(c(freq_perturb_output_path,'go_table_frequent_perturbations.html'),collapse='/')
write(kable(format_funcassociate(go_enrichment),format='html',caption='Go Enrichment amongst proteins participating in frequently dynamic complexes'),outfile)