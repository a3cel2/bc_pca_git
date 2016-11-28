modify_edges <- function(universe_matrix,
                      matched_empty_matrix,
                      node=NULL,
                      edge,
                      n_edge_indeces,
                      sample_mode,
                      modification_mode='addition'){
  if(sample_mode == 'node'){
    if(modification_mode == 'addition'){
      v1 <- universe_matrix[node,]
      v2 <- universe_matrix[,node]
      matched_empty_matrix[node,] <- v1
      matched_empty_matrix[,node] <- v2
    } else if(modification_mode == 'subtraction'){
      matched_empty_matrix[node,] <- 0
      matched_empty_matrix[,node] <- 0
    }
  } else if(sample_mode == 'edge'){
    if(modification_mode == 'addition'){
      
      matched_empty_matrix[edge[1],edge[2]] <- universe_matrix[edge[1],edge[2]]
    } else if(modification_mode == 'subtraction'){
      matched_empty_matrix[edge[1],edge[2]] <- 0
    }
  }
  return(matched_empty_matrix)
}

nodewise_simulation <- function(universe_matrix,
                                matched_empty_matrix,
                                n_edges_enhanced,
                                n_edges_depleted,
                                nodelist,
                                prob_node=1){
  
  remaining_edges_enhanced <- n_edges_enhanced
  remaining_edges_depleted <- n_edges_depleted
  
  edge_indeces <- Matrix::which(as.matrix(universe_matrix) != 0,arr.ind=T)
  n_edge_indeces <- nrow(edge_indeces)
  
  n_edges_edgewise_enhanced <- round((1 - prob_node)*n_edges_enhanced)
  n_edges_nodewise_enhanced <- n_edges_enhanced - n_edges_edgewise_enhanced
  
  n_edges_edgewise_depleted <- round((1 - prob_node)*n_edges_depleted)
  n_edges_nodewise_depleted <- n_edges_depleted - n_edges_edgewise_depleted
  
  remaining_edges_nodewise_enhanced <- n_edges_nodewise_enhanced
  remaining_edges_nodewise_depleted <- n_edges_nodewise_depleted
  
  #Two separate matrices here make for more efficient calculations
  matched_empty_matrix_enhanced <- matched_empty_matrix
  matched_empty_matrix_depleted <- matched_empty_matrix
  
  #Alternate between enhanced and depleted sampling
  i <- 0
  while(remaining_edges_enhanced > 0 | remaining_edges_depleted > 0){
    #Enhanced sampling mode
    if(i %% 2 == 0){
      matched_empty_matrix <- matched_empty_matrix_enhanced
      remaining_edges_nodewise <- remaining_edges_nodewise_enhanced
      remaining_edges <- remaining_edges_enhanced
      n_edges <- n_edges_enhanced
      n_edges_nodewise <- n_edges_nodewise_enhanced
    #Depleted sampling mode
    }else if(i %% 2 == 1){
      matched_empty_matrix <- matched_empty_matrix_depleted
      remaining_edges_nodewise <- remaining_edges_nodewise_depleted
      remaining_edges <- remaining_edges_depleted
      n_edges <- n_edges_depleted
      n_edges_nodewise <- n_edges_nodewise_depleted
    }
    if(remaining_edges > 0){
      #Add 'node-centric' interaction changes until you're done
      #then start adding edge-centric interaction changes afterwards
      if(remaining_edges_nodewise > 0){
        sample_mode <- 'node'
        node <- sample(nodelist,1)
        edge <- NULL
      } else{
        sample_mode <- 'edge'
        node <- NULL
        edge <- edge_indeces[sample(n_edge_indeces,1),]
      }
      
      #print('we are sampling in this mode:')
      #print(sample_mode)
      
      test_matched_empty_matrix <- modify_edges(universe_matrix=universe_matrix,
                                                matched_empty_matrix,
                                                node,
                                                edge,
                                                n_edge_indeces,
                                                sample_mode=sample_mode)
      if(!identical(test_matched_empty_matrix,matched_empty_matrix)){
        sum_nonzero <- Matrix::nnzero(test_matched_empty_matrix)
        if(sample_mode == 'node'){
          remaining_edges <- n_edges - sum_nonzero
          remaining_edges_nodewise <- n_edges_nodewise - sum_nonzero
        } else {
          remaining_edges <- n_edges - sum_nonzero
        }
        #print('i took a sample and need the this made more edges now:')
        #print(remaining_edges_nodewise)
        
        #Don't oversample edges if doing node-based sampling - keep sampling a node until you don't over-sample
        #(e.g. if 1 edge needed, you can only sample a 'node' with interaction degree 1)
        #print('the amount of nodewise edges is:')
        #print(remaining_edges_nodewise)
        while(remaining_edges_nodewise < 0 & sample_mode == 'node'){
        #  print('I oversampled, attempting to correct')
          node <- sample(nodelist,1)
          test_matched_empty_matrix <- modify_edges(universe_matrix,
                                                    matched_empty_matrix,
                                                    node,
                                                    edge,
                                                    n_edge_indeces,
                                                    sample_mode=sample_mode)
          sum_nonzero <- Matrix::nnzero(test_matched_empty_matrix)
          remaining_edges_nodewise <- n_edges_nodewise - sum_nonzero
          #print(remaining_edges_nodewise)
        }
      }
      
      #Update, if there is a conflict in edges, remove it from the other matrix
      #to ensure exclusivity of edges
      if(i %% 2 == 0){
        matched_empty_matrix_enhanced <- test_matched_empty_matrix
        matched_empty_matrix_depleted <- modify_edges(universe_matrix,
                                                      matched_empty_matrix_depleted,
                                                      node,
                                                      edge,
                                                      n_edge_indeces,
                                                      sample_mode=sample_mode,
                                                      modification_mode='subtraction')
        remaining_edges_enhanced <- remaining_edges
        nzero_alt <- Matrix::nnzero(matched_empty_matrix_depleted)
        remaining_edges_depleted <- n_edges_depleted - nzero_alt
        if(sample_mode == 'node'){
          remaining_edges_nodewise_enhanced <- remaining_edges_nodewise
          remaining_edges_nodewise_depleted <- n_edges_nodewise_depleted - nzero_alt
        }
      }
      if(i %% 2 == 1){
        matched_empty_matrix_depleted <- test_matched_empty_matrix
        matched_empty_matrix_enhanced <- modify_edges(universe_matrix,
                                                      matched_empty_matrix_enhanced,
                                                      node,
                                                      edge,
                                                      n_edge_indeces,
                                                      sample_mode=sample_mode,
                                                      modification_mode='subtraction')
        remaining_edges_depleted <- remaining_edges
        nzero_alt <- Matrix::nnzero(matched_empty_matrix_enhanced)
        remaining_edges_enhanced <- n_edges_enhanced - nzero_alt
        if(sample_mode == 'node'){
          remaining_edges_nodewise_depleted <- remaining_edges_nodewise
          remaining_edges_nodewise_enhanced <- n_edges_nodewise_enhanced - nzero_alt
        }
      }
    }
    i <- i + 1
    #print(i)
  }
  edgelist_enh <- as.matrix(Matrix::summary(matched_empty_matrix_enhanced))
  edgelist_enh <- cbind(rownames(matched_empty_matrix_enhanced)[edgelist_enh[,1]],colnames(matched_empty_matrix)[edgelist_enh[,2]])
  
  edgelist_depl <- as.matrix(Matrix::summary(matched_empty_matrix_depleted))
  edgelist_depl <- cbind(rownames(matched_empty_matrix_depleted)[edgelist_depl[,1]],colnames(matched_empty_matrix_depleted)[edgelist_depl[,2]])
  
  
  return(list('enhanced'=edgelist_enh,
              'depleted'=edgelist_depl))
  #return(cbind(rownames(matched_empty_matrix)[matr_sum[,1]],colnames(matched_empty_matrix)[matr_sum[,2]]))
}

make_network_iterations <- function(pca_universe,
                                    pca_enhanced_file,
                                    pca_depleted_file,
                                    condition,
                                    n_iters,
                                    mode = 'edgewise',
                                    node_sensitivity = 1,
                                    prob_node = 1,
                                    metric = get_largest_component_from_edgelist,
                                    n_parallel = 4,
                                    cluster = NULL) {
  
  
  real_edgelist_enhanced <- as.matrix(dplyr::filter(pca_enhanced_file,Condition==condition) 
                             %>% dplyr::select(ORF.1,ORF.2))
  real_edgelist_depleted <- as.matrix(dplyr::filter(pca_depleted_file,Condition==condition) 
                                      %>% dplyr::select(ORF.1,ORF.2))
  
  
  n_edges_enhanced <- nrow(real_edgelist_enhanced)
  n_nodes_enhanced <- length(unique(as.vector(real_edgelist_enhanced)))
  
  n_edges_depleted <- nrow(real_edgelist_depleted)
  n_nodes_depleted <- length(unique(as.vector(real_edgelist_depleted)))
  
  n_edges_total <- n_edges_enhanced + n_edges_depleted
  
  
  potential_edgelist <- as.matrix(dplyr::filter(pca_universe,Condition==condition) 
                                  %>% dplyr::select(ORF.1,ORF.2))
  potential_nrow <- nrow(potential_edgelist)
  
  if(is.null(cluster)){
    cl <- snow::makeCluster(n_parallel)
  }
  else{
    cl <- cluster
  }
  if(mode == 'edgewise'){
    snow::clusterExport(cl,c('metric'),environment())
    retval <- snow::parLapply(cl,1:n_iters,function(n){
      edge_samples <- sample(1:potential_nrow,n_edges_total,replace=F)
      enhanced_sample <- edge_samples[1:n_edges_enhanced]
      depleted_sample <- edge_samples[(n_edges_enhanced+1):n_edges_total]
      #depleted_sample <- sample(1:potential_nrow[!(1:potential_nrow)],n_edges,replace=F)
      edgelist_enhanced <- potential_edgelist[enhanced_sample,,drop=F]
      edgelist_depleted <- potential_edgelist[depleted_sample,,drop=F]
      
      if(!is.null(metric)){
        return(list('enhanced'=metric(edgelist_enhanced),
                    'depleted'=metric(edgelist_depleted)))
      }
      return(edgelist)
    })
  } else if( mode == 'nodewise'){
    nodelist <- unique(as.vector(unlist(potential_edgelist)))
    universe_matrix <- igraph::as_adjacency_matrix(igraph::graph_from_edgelist(as.matrix(potential_edgelist)))
    universe_matrix[universe_matrix > 1] <- 1
    
    matched_empty_matrix <- universe_matrix
    matched_empty_matrix[matched_empty_matrix > 0 ] <- 0
    
    snow::clusterExport(cl,c('metric',
                             'nodewise_simulation',
                             'n_edges_enhanced',
                             'n_edges_depleted',
                             'modify_edges',
                             'prob_node',
                             'nodelist',
                             'universe_matrix',
                             'matched_empty_matrix'),
                        environment())
    snow::clusterEvalQ(cl, {library(Matrix)})
    #nodewise_simulation(universe_matrix,matched_empty_matrix,n_edges,nodelist,prob_node=prob_node)
    
    retval <- snow::parLapply(cl,1:n_iters,function(n){
      edgelist <- nodewise_simulation(universe_matrix,
                                      matched_empty_matrix,
                                      n_edges_enhanced,
                                      n_edges_depleted,
                                      nodelist,
                                      prob_node=prob_node)
      if(!is.null(metric)){
        return(list('enhanced'=metric(edgelist$enhanced),
                    'depleted'=metric(edgelist$depleted)))
      }
      return(edgelist)
    })
    
    
  }
  if(is.null(cluster)){
    snow::stopCluster(cl)
  }
  return(retval)
}

network_simulation_significance <- function(pca_universe,
                                            pca_enhanced,
                                            pca_depleted,
                                            excluded_condition_grep = 'CRISPR',
                                            iterations = 1000,
                                            seed=123,
                                            sigval=0.05,
                                            mode='edgewise',
                                            node_sensitivity = 1,
                                            prob_node=1,
                                            metric='component_size',
                                            n_parallel=8,
                                            cluster=NULL){
  set.seed(seed)
  conditions <- as.vector(unique(c(levels(pca_enhanced$Condition),levels(pca_depleted$Condition))))
  conditions <- conditions[!(grepl(excluded_condition_grep,conditions))]
  record_list <- list()
  iteration_list <- list()
  return_matrix <- matrix(nrow=2,ncol=length(conditions),data=1)
  rownames(return_matrix) <- c('enhanced','depleted')
  colnames(return_matrix) <- conditions
  
  #for(direction in c('accumulated','depleted')){
  #  record_list[[direction]] <- list()
  #  if(direction == 'accumulated'){
  #    pca_file <- pca_enhanced
  #  }
  #  else{
  #    pca_file <- pca_depleted
  #  }
    for(condition in conditions){
      write(condition,stderr())
      if(metric == 'component_size'){
        real_metric_enhanced <- get_largest_component_size(pca_enhanced,condition)
        real_metric_depleted <- get_largest_component_size(pca_depleted,condition)
        metric_function <- get_largest_component_from_edgelist
      }
      else if (metric == 'density'){
        real_metric_enhanced <- get_density_from_edgelist(pca_enhanced,condition)
        real_metric_depleted <- get_density_from_edgelist(pca_depleted,condition)
        metric_function <- get_density_from_edgelist
      }
      if(mode == 'edgewise'){
        results <- unlist(make_network_iterations(pca_universe,
                                                  pca_enhanced,
                                                  pca_depleted,
                                                  condition=condition,
                                                  mode=mode,
                                                  n_iters=iterations,
                                                  metric=metric_function,
                                                  n_parallel=n_parallel,
                                                  cluster=cluster))
        
      iteration_list[['enhanced']][[condition]] <- results[names(results) == 'enhanced']
      iteration_list[['depleted']][[condition]] <- results[names(results) == 'depleted']
      
      }
      else if(mode =='nodewise'){
        results <- unlist(make_network_iterations(pca_universe,
                                                  pca_enhanced,
                                                  pca_depleted,
                                                  condition,
                                                  n_iters=iterations,
                                                  mode = 'nodewise',
                                                  node_sensitivity,
                                                  prob_node=prob_node,
                                                  metric=metric_function,
                                                  n_parallel=n_parallel,
                                                  cluster=cluster))
        
        iteration_list[['enhanced']][[condition]] <- results[names(results) == 'enhanced']
        iteration_list[['depleted']][[condition]] <- results[names(results) == 'depleted']
        
      }
      else{
        stop('Invalid mode specified')
      }
      record_list[['enhanced']][[condition]] <- real_metric_enhanced
      record_list[['depleted']][[condition]] <- real_metric_depleted
      for(direction in c('enhanced','depleted')){
      return_matrix[direction,condition] <- (empirical_two_tailed_p(record_list[[direction]][[condition]],iteration_list[[direction]][[condition]])+1e-04)*
        sign(record_list[[direction]][[condition]] - median(iteration_list[[direction]][[condition]]))#sum(iteration_list[[direction]][[condition]] <= record_list[[direction]][[condition]])/iterations
      }
    }
  
  #print(return_matrix)
  return_matrix[which(abs(return_matrix) <= sigval + 1e-04 & sign(return_matrix) == -1)] <- -1.01
  return_matrix[which(abs(return_matrix) <= sigval + 1e-04 & sign(return_matrix) == 1)] <- 1.01
  return_matrix[which(!(return_matrix %in% c(-1.01,1.01)))] <- 0
  
  return_matrix[return_matrix == -1.01] <- -1
  return_matrix[return_matrix == 1.01] <- 1
  return(return_matrix)
}


network_simulation_significance_node_edge_search_matrix <- function(pca_universe,
                                                                    pca_enhanced,
                                                                    pca_depleted,
                                                                    iterations=100,
                                                                    node_probs=c(0:10)/10,
                                                                    metric='component size',
                                                                    conditions=NULL,
                                                                    excluded_condition_grep='CRISPR',
                                                                    n_parallel=8,
                                                                    load_saved=F,
                                                                    save_output=F,
                                                                    save_directory=NULL,
                                                                    save_filename=NULL,
                                                                    seed=99){
  
  saved_file_path <- paste(c(save_directory,save_filename),collapse='/')
  if(file.exists(saved_file_path) & load_saved == T){
    output_table <- read.table(saved_file_path)
    colnames(output_table) <- sapply(colnames(output_table),function(x){strsplit(x,split='X')[[1]][2]})
    return(as.matrix(output_table))
  }
  
  set.seed(seed)
  if(is.null(conditions)){
    conditions <- levels(pca_universe$Condition)
  }
  conditions <- conditions[!(grepl(excluded_condition_grep,conditions))]
  result_matrix <- c()
  cl <- snow::makeCluster(n_parallel)
  for(condition in conditions){
   # for(direction in c('accumulated','depleted')){
      write(paste(c("working on:",condition),collapse=' '), stderr())
      #if(direction == 'accumulated'){
      #  pca_file <- pca_enhanced
      #}
      #if(direction == 'depleted'){
      #  pca_file <- pca_depleted
      #}
      sub_pca_universe <- dplyr::filter(pca_universe,Condition==condition) 
      val_list_enhanced <- c()
      val_list_depleted <- c()
      if(metric == 'density'){
        observed_val_enhanced <- get_density(pca_enhanced,condition)
        observed_val_depleted <- get_density(pca_depleted,condition)
      }
      else if(metric == 'component size'){
        observed_val_enhanced <- get_largest_component_size(pca_enhanced,condition)
        observed_val_depleted <- get_largest_component_size(pca_depleted,condition)
      }
      for(i in node_probs){
        write(paste(c("current node proportion:",i),collapse=' '), stderr())
        my_net <- make_network_iterations(sub_pca_universe,
                                          pca_enhanced,
                                          pca_depleted,
                                          condition,
                                          n_iters=iterations,
                                          mode='nodewise',
                                          node_sensitivity=1,
                                          prob_node=i,
                                          metric=NULL,
                                          n_parallel=8,
                                          cluster=cl)
        #print('done')
        #print(my_net)
        #stop()
        if(metric == 'component size'){
          shuffled_vals_enhanced <- sapply(my_net,function(x){get_largest_component_from_edgelist(x$enhanced)})
          shuffled_vals_depleted <- sapply(my_net,function(x){get_largest_component_from_edgelist(x$depleted)})
        }
        if(metric == 'density'){
          shuffled_vals_enhanced <- sapply(my_net,function(x){get_density_from_edgelist(x$enhanced)})
          shuffled_vals_depleted <- sapply(my_net,function(x){get_density_from_edgelist(x$depleted)})
        }
        val_list_enhanced <- c(val_list_enhanced,(empirical_two_tailed_p(observed_val_enhanced,shuffled_vals_enhanced)+1e-04)*sign(observed_val_enhanced-mean(shuffled_vals_enhanced)))
        val_list_depleted <- c(val_list_depleted,(empirical_two_tailed_p(observed_val_depleted,shuffled_vals_depleted)+1e-04)*sign(observed_val_enhanced-mean(shuffled_vals_enhanced)))
      }
      result_matrix <- rbind(result_matrix,val_list_enhanced)
      rownames(result_matrix)[nrow(result_matrix)] <- paste(c(condition,'enhanced'),collapse=' ')
      result_matrix <- rbind(result_matrix,val_list_depleted)
      rownames(result_matrix)[nrow(result_matrix)] <- paste(c(condition,'depleted'),collapse=' ')
    #}
  }
  colnames(result_matrix) <- node_probs
  snow::stopCluster(cl)
  if(save_output == T){
    write.table(result_matrix,file=saved_file_path)
  }
  return(result_matrix)
}
