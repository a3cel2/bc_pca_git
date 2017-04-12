devtools::use_package('rPython')
devtools::use_package('jsonlite')
devtools::use_package('stargazer')
devtools::use_package('R.matlab') 
devtools::use_package('stringr') 
devtools::use_package('MASS') 

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
                          associations_file=system.file('go_map.txt',package='bcPcaAnalysis'),
                          order_mode='ordered',
                          go_type='default_funcassociate',
                          nametype='ORF',
                          network_mode='nodewise',
                          p_value_cutoff=0.05,
                          reps=1000){
  if(network_mode == 'edgewise'){
    query <- apply(query,1,function(x){paste(x,collapse='\t')})
    universe <- apply(universe,1,function(x){paste(x,collapse='\t')})
  }
  query <- vector_to_python_list(query)
  universe <- vector_to_python_list(universe)

  go_map <- read.csv(associations_file,head=F,sep='\t')
  
  rPython::python.load(system.file('python/funcassociate_submit.py',package='bcPcaAnalysis'))
  return(jsonlite::fromJSON(rPython:::python.call("funcassociate",query,
                        universe,
                        associations_file,
                        order_mode,
                        go_type,
                        nametype,
                        network_mode,
                        p_value_cutoff,
                        reps)))
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
                                   network_modes=c('nodewise','edgewise'),
                                   reps=1000){
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
        }else if(network_mode == 'edgewise'){
          #Remove homodimers
          query <- query[query[,1] != query[,2], ]
          universe <- universe[universe[,1] != universe[,2], ]
        }
        funcassociate_output_list[[condition]][[direction]][[network_mode]] <-
          funcassociate(query,universe,network_mode=network_mode,order_mode = order_mode,reps=reps)
      }
    }
  }
  return(funcassociate_output_list)
}

#Formats funcassociate results to table
format_funcassociate <- function(funcassociate_json){
  result_table <- funcassociate_json$result$over
  colnames(result_table) <- c('N','M','X','LOD','P','P_adj','Gene-Ontology-ID','Gene-Ontology-Attribute')
  
  result_table[,'LOD'] <- sapply(result_table[,'LOD'],function(x){format(as.numeric(x),digits=2)})
  result_table[,'P'] <- sapply(result_table[,'P'],function(x){format(as.numeric(x),digits=2)})
  result_table[,'P_adj'] <- sapply(result_table[,'P_adj'],function(x){format(as.numeric(x),digits=2)})

  return(result_table)
}

#Generalizes above for multiple conditions and strategies
format_bc_pca_funcassociate <- function(bcpca_funcassociate_result,conditions,directions=c("enhanced","depleted"),modes=c("nodewise","edgewise")){
  resulting_table <- c()
  for(condition in conditions){
    for(mode in modes){
      for(direction in directions){
        object <- bcpca_funcassociate_result[[condition]][[direction]][[mode]]
        if(length(object$result$over) > 0){
          primitive_df <- format_funcassociate(object)
          primitive_df <- cbind(cbind(rep(condition,nrow(primitive_df)),
                                  rep(mode,nrow(primitive_df)),
                                  rep(direction,nrow(primitive_df))
                                  ),primitive_df)
          colnames(primitive_df)[1:3] <- c('Condition','Mode','Direction')
          resulting_table <- rbind(resulting_table,primitive_df)
        }
      }
    }
  }
  colnames(resulting_table)[1:3] <- c('Condition','Mode','Direction')
  return(resulting_table)
}

convert_costanzo_matlab_data <- function(matlab_file,output_file){
  if(!file.exists(output_file)){
    mat_go <- R.matlab::readMat(matlab_file)
    term_ids <- as.vector(mat_go$go[[1]])
    term_names <- unlist(mat_go$go[[2]])
    orfs <- unlist(mat_go$go[[3]])
    matrix <- mat_go$go[[4]]
    output_df <- c()
    for(i in 1:length(term_ids)){
      term_id <- term_ids[i]
      term_id <- paste(c('GO:',stringr::str_pad(term_id,7,pad="0")),collapse='')
      term_name <- term_names[i]
      genes <- paste(orfs[which(matrix[i,] > 0)],collapse=' ')
      output_df <- rbind(output_df,c(term_id,term_name,genes))
    }
    write.table(output_df,sep='\t',row.names=F,col.names=F,quote=F,file=output_file)
  }
}

convert_sgd_slim_go_data <- function(slim_file,output_file){
  if(!file.exists(output_file)){
    slim_file <- read.csv(slim_file, head=F, sep='\t',stringsAsFactors=F)
    output_list <- list()
    for(i in 1:nrow(slim_file)){
      term <- slim_file[i,6]
      if(is.null(output_list[[term]])){
        output_list[[term]] <- list('name'='','genes'=c())
      }
      output_list[[term]][['name']] <- slim_file[i,5]
      output_list[[term]][['genes']] <- c(output_list[[term]][['genes']],slim_file[i,1])
    }
    
    output_df <- c()
    go_names <- names(output_list)
    go_names <- go_names[sapply(go_names,nchar) > 3]
    for(term in go_names){
      name <- output_list[[term]][['name']]
      genes <- paste(output_list[[term]][['genes']],collapse=' ')
      output_df <- rbind(output_df,c(term,name,genes))
    }
    write.table(output_df,sep='\t',row.names=F,col.names=F,quote=F,file=output_file)
  }
}

convert_sgd_full_go_data <- function(go_file,term_map,output_file){
  if(!exists(output_file)){
    term_file <- read.delim(term_map,sep='\t',head=F,stringsAsFactors=F)
    term_mapping_list <- list()
    for(i in 1:nrow(term_file)){
      term_id <- term_file[i,1]
      term_name <- term_file[i,2]
      term_id <- paste(c('GO:',stringr::str_pad(term_id,7,pad="0")),collapse='')
      term_mapping_list[[term_id]] <- term_name
    }
    
    go_map <- read.delim(go_file,comment.char='!',sep='\t',head=F,stringsAsFactors=F,quote = "")
    output_list <- list()
    for(i in 1:nrow(go_map)){
      term <- go_map[i,5]
      if(is.null(output_list[[term]])){
        output_list[[term]] <- list('name'='','genes'=c())
      }
      gene <- strsplit(go_map[i,11],split='\\|')[[1]][1]
      output_list[[term]][['genes']] <- c(output_list[[term]][['genes']],gene)
      output_list[[term]][['name']] <- term_mapping_list[[term]]
    }
    output_df <- c()
    go_names <- names(output_list)
    go_names <- go_names[sapply(go_names,nchar) > 3]
    print(length(go_names))
    for(term in go_names){
      name <- output_list[[term]][['name']]
      genes <- paste(output_list[[term]][['genes']],collapse=' ')
      output_df <- rbind(output_df,c(term,name,genes))
    }
    write.table(output_df,sep='\t',row.names=F,col.names=F,quote=F,file=output_file)
  }
}