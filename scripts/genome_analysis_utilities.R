library(utils)
library(seqinr)
library(org.Sc.sgd.db)

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

bpc_notation_to_yeast_orf <- function(bpc_list){
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

