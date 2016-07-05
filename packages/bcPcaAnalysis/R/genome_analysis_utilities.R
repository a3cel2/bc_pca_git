devtools::use_package('utils')
devtools::use_package('seqinr')
devtools::use_package('org.Sc.sgd.db')
devtools::use_package('dplyr')


#' Download a gzipped file and unzip it (used for genome downloading)
#'
#' @param genome_http_dir Web address of genome
#' @param genome_dest_dir Where you want it downloaded
#'
download_unzip_genome <- function(genome_http_dir,genome_dest_dir){
  gunzip_command <- paste(c('gunzip',genome_dest_dir),collapse=' ')
  download.file(genome_http_dir,genome_dest_dir)
  system(gunzip_command)
}


#' Get yeast genes matching a regular expression
#' @description Downloads yeast genome to a given directory
#' if it does not already exist, and checks the genome
#' (a gzipped multi-fasta file) for genes matching a given
#' regular expression
#'
#' @param motif_regexp Desired regular expression
#' @param genome_http_dir Web address of genome, a gzipped multi-fasta file
#' @param genome_dest_dir Where you want genome downloaded
#'
#' @return names of yeast genes (ORF name - e.g. YKR103W) matching desired expression
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

#' Converts a list of BPC notations to pairs of yeast ORFS
#' @description given a list of BPC IDs such as 'MET3_F1::MET3_F3', convert to a pair of yeast ORFs
#'
#' @param bpc_list list of BPC ids
#' @return two column table of Yeast ORFs corresponding to BPC Ids
bpc_notation_to_yeast_orf_pair <- function(bpc_list){
  name_matr <- t(as.matrix(sapply(bpc_list,function(bpc){
    split_name <- strsplit(bpc,'::')[[1]]
    parsed_names <- sapply(split_name,function(unparsed_name){
      strsplit(unparsed_name,'_')[[1]][1]
    })
    sapply(parsed_names,function(name){
      test_name <- org.Sc.sgd.db::org.Sc.sgdCOMMON2ORF[[name]]
      if(!is.null(test_name)){
        return(test_name)
      }
      return(name)
    })
  })))
  colnames(name_matr) <- c('Gene1','Gene2')
  return(name_matr)
}



#' Test for conditional motif enrichment
#' @description Calculates possible enrichment for a given motif
#' under a given condition
#' @param universe_file A data frame which contains all the pairs
#' in the PCA experiment under the columns 'ORF.1','ORF.2'
#' @param pca_call_file A data frame which contains all signficant calls
#' in the PCA experiment (either enhanced or depleted), 
#' with the condition of interest under the
#' 'Condition' column and the orf pairs under the columns 'ORF.1','ORF.2'
#' @param motif_regexp The regular expression to search for
#' @param condition the name of the condition, as it appears in pca_call_file
given_motif_enrichment <- function(universe_file,
                                   pca_call_file,
                                   motif_regexp,
                                   condition){
  orf_universe <- read.csv(universe_file,sep='\t',stringsAsFactors = F)
  orf_universe <- tbl_df(orf_universe)
  
  all_orfs <- unique(unlist(dplyr::select(orf_universe,ORF.1,ORF.2)))
  
  
  #Found in condition
  pca_calls <- read.csv(pca_call_file,sep='\t',stringsAsFactors = F)
  pca_calls <- tbl_df(pca_calls)
  
  
  filtered_pca_calls <- dplyr::filter(pca_calls, Condition == condition)#, Effect.on.growth == direction)
  
  condition_orfs <- unique(unlist(dplyr::select(filtered_pca_calls,ORF.1,ORF.2)))
  
  
  
  all_regex_matching_orfs <- get_genes_matching_regex_motif(motif_regexp=motif_regexp)
  all_regex_matching_orfs <- intersect(all_orfs,all_regex_matching_orfs )
  
  atorvastatin_prenylated_orfs <- intersect(condition_orfs,all_regex_matching_orfs)
  
  stat_table <- rbind(c(length(all_regex_matching_orfs),length(all_orfs) - length(all_regex_matching_orfs)),
                      c(length(atorvastatin_prenylated_orfs),length(condition_orfs) - length(atorvastatin_prenylated_orfs)))
  
  colnames(stat_table) <- c('Matching Motif','Not Matching Motif')
  rownames(stat_table) <- c('Universe',condition)
  
  return(list(
    fisher_table = stat_table,
    fisher_test = fisher.test(stat_table)
  ))
}