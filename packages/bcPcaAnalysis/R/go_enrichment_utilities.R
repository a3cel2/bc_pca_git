require(rPython)

#' Converts a vector to a literal python list representation (all items in strings)
#'
#' @param contents a vector of string
#'
#' @return a list which contains one item per line in the file name
#'
#' @examples 
#' vector_to_python_list(c('A','B'))
#' "['A','B']"
vector_to_python_list <- function(contents){
  contents <- sapply(contents,function(x){paste(c('"',x,'"'),collapse='')})
  contents <- paste(contents,collapse=',')
  file_str <- paste(c('[',contents,']'),collapse='')
  return(file_str)
}


#' Submit a query to funcassociate and get the results
#'
#' @param query a vector of gene identifiers, or if mode is 'edgewise', a two column matrix
#' of gene names
#' @param universe a vector of all possible gene identifiers in the experiment, in the
#' same format as query.  If 'edgewise', PPIs must have the same order as in the query
#' @param associations_file the file used to build the GO associations, see examples
#' @param order_mode either ordered or unordered
#' @param go_type either default_funcassociate or slim, each has their own
#' ways to parse the file, see examples for the appropriate format
#' @param nametype either 'common' or 'ORF'.  'common' is not compatible with
#' the default_funcassociate mapping mode
#' @param network_mode either 'nodewise' or 'edgewise', 'nodewise' is the classic
#' GO analysis, in 'edgewise' 
#' @param p_value_cutoff p value cutoff for the enrichment analysis, passed directly to the server
#'
#' @return The funcassociate server output parsed into a list
#'
funcassociate <- function(query,
                          universe,
                          associations_file,
                          order_mode='ordered',
                          go_type='default_funcassociate',
                          nametype='ORF',
                          network_mode='nodewise',
                          p_value_cutoff=0.05){
  if(network_mode == 'edgewise'){
    query <- apply(query,1,function(x){paste(x,collapse='\t')})
    universe <- apply(universe,1,function(x){paste(x,collapse='\t')})
  }

  query <- vector_to_python_list(query)
  universe <- vector_to_python_list(universe)
  
  rPython::python.load(system.file('python/test.py',package='bcPcaAnalysis'))
  return(fromJSON(rPython:::python.call("funcassociate",query,
                        universe,
                        associations_file,
                        order_mode,
                        go_type,
                        nametype,
                        network_mode,
                        p_value_cutoff)))
}


#Helper function for bcpca_funcassociate_analysis
edgewise_to_nodewise <- function(edge_list){
  node_list <- c()
  for(i in 1:nrow(edge_list)){
    ppi_pair <- edge_list[i,]
    for(protein in ppi_pair){
      if(!(protein %in% node_list)){
        node_list <- c(node_list,protein)
      }
    }
  }
  return(node_list)
}


bcpca_funcassociate_analysis <- function(universe_file,
                                   enhanced_calls_file,
                                   depleted_calls_file,
                                   conditions,
                                   directions=c('enhanced','depleted'),
                                   order_mode="ordered",
                                   network_modes=c('nodewise','edgewise')){
  funcassociate_output_list <- list()
  
  pca_universe_table <- read.csv(universe_file,sep='\t',stringsAsFactors=F)
  pca_enhanced_calls_table <- read.csv(enhanced_calls_file,sep='\t',stringsAsFactors=F)
  pca_depleted_calls_table <- read.csv(depleted_calls_file,sep='\t',stringsAsFactors=F)
  
  
  if(sum(directions %in% c('enhanced','depleted')) != length(directions)){
    stop("directions must be one or both of 'enhanced' or 'depleted'")
  }
  
  if(sum(network_modes %in% c('nodewise','edgewise')) != length(network_modes)){
    stop("Network modes must be one or both of 'nodewise' or 'edgewise'")
  }
  
  for(condition in conditions){
    pca_universe_table_cond <- dplyr::filter(pca_universe_table, Condition == condition)
    pca_enhanced_calls_table_cond <- dplyr::filter(pca_enhanced_calls_table, Condition == condition)
    pca_depleted_calls_table_cond <- dplyr::filter(pca_depleted_calls_table, Condition == condition)
    for(direction in directions){
      for(network_mode in network_modes){
        if(direction == 'enhanced'){
          pca_calls = pca_enhanced_calls_table_cond
        }
        if(direction == 'depleted'){
          pca_calls = pca_depleted_calls_table_cond
        }
        if(order_mode == 'ordered'){
          effect_size_sort <- sort(abs(pca_calls$FC..UP.tag.+pca_calls$FC..DN.tag.),index.return=T,decreasing=T)$ix
          query <- pca_calls[effect_size_sort,c('ORF.1','ORF.2')]
        }
        if(order_mode == 'unordered'){
          query <- pca_calls[,c('ORF.1','ORF.2')]
        }
        universe <- pca_universe_table_cond[,c('ORF.1','ORF.2')]
        if(network_mode == 'nodewise'){
          query <- edgewise_to_nodewise(query)
          universe <- edgewise_to_nodewise(universe)
        }
        funcassociate_output_list[[condition]][[direction]][[network_mode]] <-
          funcassociate(query,universe,go_association_file,network_mode=network_mode,order_mode = order_mode)
      }
    }
  }
  return(funcassociate_output_list)
}