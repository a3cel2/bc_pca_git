library(utils)
library(seqinr)
library(org.Sc.sgd.db)
library(dplyr)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

download_unzip_genome <- function(genome_http_dir,genome_dest_dir){
  gunzip_command <- paste(c('gunzip',genome_dest_dir),collapse=' ')
  download.file(genome_http_dir,genome_dest_dir)
  system(gunzip_command)
}

get_genes_matching_regex_motif <- function(
  motif_regexp = 'c(g|a|v|l|i){2,2}(m|s|q|a|c|l|e)\\*',
  genome_http_dir = 'http://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz',
  genome_dest_dir = '../data/yeast_seq.gz'
){
  
  if(!file.exists(strsplit(genome_dest_dir,split='.gz')[[1]])){
    download_unzip_genome(genome_http_dir = genome_http_dir,
                          genome_dest_dir = genome_dest_dir)
  }
  
  #Parse_genome
  fasta_file <- strsplit(genome_dest_dir,'.gz')[[1]][1]
  yeast_proteome <- read.fasta(fasta_file,as.string=T,set.attributes = F)
  motif_search <- sapply(yeast_proteome,function(x){grep(motif_regexp,x,value=T)})
  
  return(names(unlist(motif_search)))
}

bpc_notation_to_yeast_orf_pair <- function(bpc_list){
  name_matr <- t(as.matrix(sapply(bpc_list,function(bpc){
    split_name <- strsplit(bpc,'::')[[1]]
    parsed_names <- sapply(split_name,function(unparsed_name){
      strsplit(unparsed_name,'_')[[1]][1]
    })
    sapply(parsed_names,function(name){
      test_name <- org.Sc.sgdCOMMON2ORF[[name]]
      if(!is.null(test_name)){
        return(test_name)
      }
      return(name)
    })
  })))
  colnames(name_matr) <- c('Gene1','Gene2')
  return(name_matr)
}

given_motif_enrichment <- function(universe_file,
                                   pca_call_file,
                                   motif_regexp,
                                   drug,
                                   direction){
  orf_universe <- read.csv(universe_file,sep='\t',stringsAsFactors = F)
  all_orfs <- unique(unlist(select(orf_universe,ORF.1,ORF.2)))
  
  #Found in condition
  pca_calls <- read.csv(pca_call_file,sep='\t',stringsAsFactors = F)
  
  filtered_pca_calls <- filter(pca_calls, Condition == drug, Effect.on.growth == direction)
  
  atorvastatin_enhanced_orfs <- unique(bpc_notation_to_yeast_orf_pair(as.vector(filtered_pca_calls$BPC)))
  
  all_regex_matching_orfs <- get_genes_matching_regex_motif(motif_regexp=motif_regexp)
  all_regex_matching_orfs <- intersect(all_orfs,all_regex_matching_orfs )
  
  atorvastatin_prenylated_orfs <- intersect(atorvastatin_enhanced_orfs,all_regex_matching_orfs)
  
  stat_table <- rbind(c(length(all_regex_matching_orfs),length(all_orfs)),
                      c(length(atorvastatin_prenylated_orfs),length(atorvastatin_enhanced_orfs)))
  
  return(fisher.test(stat_table))
}

