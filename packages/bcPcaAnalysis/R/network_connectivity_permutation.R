devtools::use_package('reshape2')
devtools::use_package('ggplot2')
devtools::use_package('dplyr')
devtools::use_package('Hmisc')
devtools::use_package('igraph')
devtools::use_package('snow')
devtools::use_package('Matrix')
library(ggplot2)

#my_color_list <- c(
#  rgb(1,0.45,0.25),
#  rgb(0.8,0.25,0.25),
#  rgb(0,0,0),
#  rgb(0.25,0.45,0.8),
# rgb(0.25,0.75,1)
#)


#' Return largest connected component
#'
#' @param graph and igraph graph object
#'
#' @return the subgraph consisting of the largest connected component in the graph
keep_only_largest_connected_component <- function(graph){
  components <- igraph::clusters(graph)
  largest_component <- which.max(components$csize)
  vertices_in_other_components <- names(which(components$membership != largest_component))
  new_graph <- igraph::delete.vertices(graph,vertices_in_other_components)
  return(new_graph)
}



#' Sets graph aesthetics for an igraph object - nodes are given a singe colour, edges are given multiple
#'
#' @param graph igraph object
#' @param node_size size to plot nodes
#' @param node_colour colour for all nodes
#' @param node_text_colour text colour for all nodes
#' @param edge_width width to plot edges
#' @param edge_attribute attribute by which to set the edge colour scale by
#' @param edge_colour_scale the colour scale object for the edge - a list of colours e.g. c('red','green')
#' @param edge_colour_min the minimum value for the colour scale
#' @param edge_colour_max the maximum value for the colour scale
#' @param edge_ncolours how many colours in the gradient
#' @param text_size_constant used to set character size based on number of characters, higher = smaller text
#'
#' @return an igraph graph with the given graphics parameters modified
set_network_aesthethics <- function(graph,
                                    node_size,
                                    node_colour,
                                    node_text_colour,
                                    edge_width,
                                    edge_attribute,
                                    edge_colour_scale,
                                    edge_colour_min,
                                    edge_colour_max,
                                    edge_ncolours=100,
                                    text_size_constant=25
){
  igraph::V(graph)$color <- node_colour
  igraph::V(graph)$label.cex <- ((node_size/text_size_constant)*4)/sapply(igraph::V(graph)$name,nchar)
  igraph::V(graph)$label.family <- "Arial Black"
  igraph::V(graph)$label.color <- node_text_colour
  igraph::E(graph)$color <- set_colours(edge_attribute,edge_colour_scale,edge_ncolours,edge_colour_min,edge_colour_max)
  igraph::E(graph)$width <- edge_width
  return(graph) 
}


#' Collapse multiple measurements for PCA file - returns mean of numeric columns in duplicate rows of a file,
#' with Gene.1 and Gene.2 columns as keys
#'
#' @param pca_file a pca file, standardly formatted - object will be grouped by identity in the Gene.1 and Gene.2 columns
#'
#' @return a pca file with multiple measurements per ORF pair averaged
collapse_multigenes <- function(pca_file){
  pca_file <- dplyr::group_by(pca_file, Gene.1,Gene.2)
  sum_fun <- function(x){
    if(typeof(x) == "double"){
      return(mean(x))
    } else {
      return(x[1])
    }
  }
  return(as.data.frame(dplyr::summarise_each(pca_file,funs(sum_fun))))
}

#' Plots a graph corresponding to the largest connected components in the enhanced and depleted
#' PPIs for a given PCA object
#'
#' @param pca_universe A PCA file showing all possible interactions for any given condition
#' @param pca_enhanced A PCA file showing enhanced interactions for each condition
#' @param pca_depleted A PCA file showing depleted interactions for each condition
#' @param condition the condition the scale should be plotted for
#' @param color_scale a list of colours to be made into a gradient
#' @param default_node_color what to colour the nodes if they are not to be made into a colour scale
#' @param node_size 
#' @param edge_width 
#' @param min_fc minimum value for colour scale
#' @param max_fc maximum "
#' @param node_text_colour tex value for nodes 
#' @param seed random seed
#' @param edit Should graph be edited prior to plotting? Opens a tkplot window to allow editing of plot
#' @param load_saved Should a prior layout be loaded?
#' @param layout_save_dir Location of prior layour (absolute path)
#' @param layout_file Name of layout file
#' @param my_title what to title the plot?
#'
#' @return A plot of the largest connected component in a given condition
network_connectivity_graph <- function(pca_universe,
                                       pca_enhanced,
                                       pca_depleted,
                                       condition,
                                       color_scale,
                                       default_node_color=rgb(0.3,0.3,0.3),
                                       node_size=10,
                                       edge_width=7,
                                       min_fc=-1,
                                       max_fc=1,
                                       node_text_colour='white',
                                       seed=104,
                                       edit=F,
                                       load_saved=T,
                                       to_plot='both',
                                       layout_save_dir='',
                                       layout_file='default.layout',
                                       network_export_file='default.tsv',
                                       my_title='\n\nDoxorubicin'){
  condition_enhanced <- collapse_multigenes(dplyr::filter(pca_enhanced,Condition==condition))
  condition_depleted <- collapse_multigenes(dplyr::filter(pca_depleted,Condition==condition))
  
  
  #Filter by largest connected component separately
  g1 <- igraph::graph_from_edgelist(as.matrix(condition_enhanced[,c('Gene.1','Gene.2')]),directed=T)
  g2 <- igraph::graph_from_edgelist(as.matrix(condition_depleted[,c('Gene.1','Gene.2')]),directed=T)
  g1_n <- keep_only_largest_connected_component(g1)
  g2_n <- keep_only_largest_connected_component(g2)
  
  if(to_plot == 'both'){
    g <- igraph::graph.union(g1_n,g2_n)
  } else if(to_plot == 'enhanced'){
    g <- g1_n
  } else if(to_plot == 'depleted'){
    g <- g2_n
  }
  
  #Get changes for new graph
  edges <- igraph::get.edgelist(g)
  colnames(edges) <- c('Gene.1','Gene.2')
  combined_vals <- rbind(condition_enhanced,condition_depleted)
  combined_vals <- merge(edges,combined_vals,by=c('Gene.1','Gene.2'),sort=F,all.x=F,all.y=F)
  
  expr <- apply(combined_vals[,c('FC..UP.tag.','FC..DN.tag.')],1,mean)
  g <- igraph::as.undirected(g,mode='each')
  g <- set_network_aesthethics(g,
                                node_size,
                                default_node_color,
                                node_text_colour,
                                edge_width,
                                edge_attribute = expr,
                                edge_colour_scale = color_scale,
                                edge_colour_min = min_fc,
                                edge_colour_max = max_fc
  )
  par(oma=c(0,0,0,0),mar=c(0,0,0,0))
  set.seed(seed)
  layout_file_path <- paste(c(layout_save_dir,layout_file),collapse='/')
  net_export_file_path <- paste(c(layout_save_dir,network_export_file),collapse='/')
  if(file.exists(layout_file_path) & load_saved == T){
   l <- read.table(layout_file_path,head=F) 
   l <- as.matrix(l)
  } else{
    l <- igraph::layout.lgl(g)
  }
  
  if(edit == T){
    editplot <- igraph::tkplot(g,layout=l)
    n <- readline(prompt="Press Enter when done editing layout, then close the window")
    l <- igraph::tkplot.getcoords(editplot)
  }
  write.table(l,file = layout_file_path, col.names=F,quote=F,row.names=F)
  write.table(cbind(igraph::get.edgelist(g),expr),quote=F,col.names=F,row.names=F,file=net_export_file_path)
  plot(g,layout=l,vertex.size=node_size,edge.width=edge_width)
  hub_legend_draw(min_fc,max_fc,my_color_list,main_font_size = 0.7,side_font_size = 0.7 ,width=0.2,new_plot=F,x_adjust=-1,y_adjust=-1)
  
  par(mar=c(0,0,1,0))
  par(oma=c(0,0,0,0))
  
  par(cex=1)
  title(my_title)
}


#' Returns the size of the largest component from a list of edges
#'
#' @param edgelist a two-column matrix
#'
#' @return the size of the largest connected component  (numeric)
get_largest_component_from_edgelist <- function(edgelist){
  graph <- igraph::graph_from_edgelist(edgelist)
  return(max(igraph::clusters(graph)$csize))
}

#' Given a PCA file, returns the size of the largest connected component for a given condition
#'
#' @param pca_file a list of interactions, split by Condition in the Condition columm
#' and each pair indexed by ORF.1 and ORF.2
#' @param condition in what condition the largest connected component is needed (e.g. FK506)
#'
#' @return the size of the largest connected component (numeric)
get_largest_component_size <- function(pca_file,
                                       condition){
  edgelist <- as.matrix(dplyr::filter(pca_file,Condition==condition) 
                        %>% dplyr::select(ORF.1,ORF.2))
  return(get_largest_component_from_edgelist(edgelist))
}

#' Returns the size of the largest component from a list of edges
#'
#' @param edgelist a two-column matrix
#'
#' @return graph density  (numeric)
get_density_from_edgelist <- function(edgelist){
  n_vertices <- length(unique(as.vector(edgelist)))
  return(2*nrow(edgelist)/(n_vertices*(n_vertices-1)))
}

#' Given a PCA file, returns the graph density for a given condition
#'
#' @param pca_file a list of interactions, split by Condition in the Condition columm
#' and each pair indexed by ORF.1 and ORF.2
#' @param condition in what condition the largest connected component is needed (e.g. FK506)
#'
#' @return the size of the largest connected component (numeric)
get_density <- function(pca_file,
                        condition){
  edgelist <- as.matrix(dplyr::filter(pca_file,Condition==condition) 
                        %>% dplyr::select(ORF.1,ORF.2))
  return(get_density_from_edgelist(edgelist))
}

#Indexes an edgelist to be accessible by a given node, used to speed up some simulations, now defunct
edge_index_universe <- function(universe_edgelist){
  node_index <- list()
  for(i in 1:nrow(universe_edgelist)){
    node_index[[universe_edgelist[i,1]]] <- rbind(node_index[[universe_edgelist[i,1]]],universe_edgelist[i,])
  }
}


#' Give two-tailed empirical p value from a bootstrap
#'
#' @param observed_value value observed in the experiment 
#' @param iterated_value_list values obtained by re-sampling
#'
#' @return proportion of values as extreme as one observed or more
empirical_two_tailed_p <- function(observed_value,iterated_value_list){
  lower_p <- sum(observed_value <= iterated_value_list)/length(iterated_value_list)
  higher_p <- sum(observed_value >= iterated_value_list)/length(iterated_value_list)
  return(min(c(min(c(lower_p,higher_p))*2,1)))
}

#' Modify edges from a protein-protein interaction matrix into a sampled subnetwork
#'
#' @param universe_matrix a sparse matrix with all possible edges
#' @param matched_empty_matrix a sparse matric that is a subset of universe_matrix, to which edges are to be added
#' @param nodelist a list of possible nodes
#' @param edge_indeces an index of all entries in universe_matrix corresponding to an edge, in row,column format
#' @param n_edge_indeces number of rows in edge_indeces
#' @param mode either 'node' or 'edge', whether to sample one itneraction, or all interactions involving one protein
#'
#' @return a modified version of matched_empty_matrix
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

#' Simulate a subnetwork from a larger protein-protein interaction network by randomly sampling either nodes or edges
#'
#' @param universe_matrix all possible protein-protein interactions
#' @param matched_empty_matrix an empty version of universe_matrix to store the values
#' @param n_edges number of edges needed for the simulation
#' @param nodelist nodes in universe_matrix
#' @param prob_node proportion of edges which are to be obtained by node-wise sampling rather than direct edge-based sampling
#'
#' @return a sparse matrix representing a sample of the larger network
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

#' Simulate a network that is the same size as one obtained from an experimental file
#'
#' @param pca_universe A PCA file corresponding to all protein-protein interactions
#' available in a given condition
#' @param pca_file A PCA file corresponding to all significantly changed
#' protein-protein interactions
#' @param condition the condition of interest, a string present in both pca_universe and pca_file
#' @param n_iters number of simulations to make
#' @param mode if 'edgewise', samples edges directly (much faster), if 'nodewise'
#' samples a proportion of nodes and edges given by prob_node
#' @param node_sensitivity in nodewise sampling, what is the probability that a given PPI will be
#' captured by the assay? (numeric between 0 and 1)
#' @param prob_node if mode is 'nodewise', what proportion of the sampled interactions
#' should be from a sampled node compared to a sampled edge
#' @param metric summary statistic for the returned network - should work on an edgelist
#' (a two colum matrix), if not specified, returns a list of two-column matrices
#' @param n_parallel number of parallel computing sessions to start
#' @param cluster the identity of a cluster created by the SNOW package, if not specified, it creates ones
#'
#' @return a list of two-column matrices if metric is not specified, or a list of the output of metric
#' on the two column-matrix if it is
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
  
  
  n_edges_enhanced <- nrow(unique(real_edgelist_enhanced))
  n_nodes_enhanced <- length(unique(as.vector(real_edgelist_enhanced)))
  
  n_edges_depleted <- nrow(unique(real_edgelist_depleted))
  n_nodes_depleted <- length(unique(as.vector(real_edgelist_depleted)))
  
  n_edges_total <- n_edges_enhanced + n_edges_depleted
  
  
  potential_edgelist <- unique(as.matrix(dplyr::filter(pca_universe,Condition==condition) 
                                  %>% dplyr::select(ORF.1,ORF.2)))
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


#' Returns a matrix corresponding to significant deviations from an observed network, given a defined property
#'
#' @param pca_universe A PCA file corresponding to all protein-protein interactions
#' available in a given condition
#' @param pca_enhanced A PCA file corresponding to all significantly enhanced
#' protein-protein interactions 
#' @param pca_depleted A PCA file corresponding to all significantly depleted
#' protein-protein interactions 
#' @param excluded_condition_grep a regular expression which matches all conditions to be excluded from the matrix
#' @param iterations number of sampled iterations to simulate
#' @param seed seed for random number generator
#' @param sigval significance cutoff (default 0.05)
#' @param mode either 'edgewise' or 'nodewise', whether the simulations should just sample edges from the main
#' network or sample a combination of edges, and nodes by which to obtain edges
#' @param node_sensitivity if a node is sampled, what proportion of the edges are recovered on average?
#' @param prob_node what proportion of edges in a given simulation should be obtained by sampling PPIs from nodes?
#' @param metric what summary statistic of the real and simulated network should be used to evaluate significance?
#' @param n_parallel number of parallel processes to spawn for simulation
#' @param cluster a pre-defined SNOW cluster, default is to create a new one
#'
#' @return a matrix where the rows correspond to directions and the columns are conditions,
#' values are either 1 (significant, real value is greater than simulations), 0 (nonsigificant), -1
#' (significant, real value is less than simulations)
network_simulation_significance <- function(pca_universe,
                                            pca_enhanced,
                                            pca_depleted,
                                            excluded_condition_grep = 'CRISPR',
                                            conditions = NULL,
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
  if(is.null(conditions)){
    conditions <- as.vector(unique(c(levels(pca_enhanced$Condition),levels(pca_depleted$Condition))))
    conditions <- conditions[!(grepl(excluded_condition_grep,conditions))]
  }
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
      edgelist_enhanced <- as.matrix(dplyr::filter(pca_enhanced,Condition==condition) 
                                     %>% dplyr::select(ORF.1,ORF.2))
      edgelist_depleted <- as.matrix(dplyr::filter(pca_depleted,Condition==condition) 
                                     %>% dplyr::select(ORF.1,ORF.2))
      
      real_metric_enhanced <- get_density_from_edgelist(edgelist_enhanced)
      real_metric_depleted <- get_density_from_edgelist(edgelist_depleted)
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
        sign(record_list[[direction]][[condition]] - mean(iteration_list[[direction]][[condition]]))#sum(iteration_list[[direction]][[condition]] <= record_list[[direction]][[condition]])/iterations
    }
  }
  
  print(return_matrix)
  return_matrix[which(abs(return_matrix) <= sigval + 1e-04 & sign(return_matrix) == -1)] <- -1.01
  return_matrix[which(abs(return_matrix) <= sigval + 1e-04 & sign(return_matrix) == 1)] <- 1.01
  return_matrix[which(!(return_matrix %in% c(-1.01,1.01)))] <- 0
  
  return_matrix[return_matrix == -1.01] <- -1
  return_matrix[return_matrix == 1.01] <- 1
  return(return_matrix)
}




#' Given an object from network_simulation_significance, makes a heatmap-like table representing the significance matrix
#'
#' @param my_matr a matrix from network_simulation_significance, where the rows correspond to directions and the columns are conditions,
#' values are either 1 (significant, real value is greater than simulations), 0 (nonsigificant), -1
#' (significant, real value is less than simulations)
#' @param my_color_list a colour list used to construct the legend - a gradient will be made and the first colour
#' will be used for -1, the middle will be 0, and the last colour will be 1
#' @param border_colour the border colour to draw around each grid in the heatmap
#' @param border_size the width of the border
#' @param legend_labels vector of form c('-1'=label1,'0'=label2,'1'=label3)
#'
#' @return a plot
connectivity_graph <- function(my_matr,
                               my_color_list,
                               border_colour="grey30",
                               border_size=0.8,
                               legend_labels=c('-1'="Decreased component size",
                                               '0'="Expected component size",
                                               '1'="Increased component size")
                               ){
  vals <- sort(apply(my_matr,2,function(x){sum(x[2] != 0)+2*sum(x[1] != 0)}),decreasing=T)
  level_order <- names(vals)
  cols <- grDevices::colorRampPalette(my_color_list)
  cols <- cols(3)
  my_matr.m <- reshape2::melt(my_matr)
  colnames(my_matr.m) <- c('Direction','Condition','value')
  my_matr.m$Condition <- factor(my_matr.m$Condition,
                                levels=level_order)
  my_matr.m$Direction <- factor(my_matr.m$Direction,
                                levels=c('depleted','enhanced'))
  my_matr.m[,3] <- as.character(my_matr.m[,3])
  
  myplot <- ggplot2::ggplot(data = my_matr.m, ggplot2::aes(x=Condition, y=Direction, fill=factor(value))) + 
    ggplot2::geom_tile(color=border_colour,size=border_size) +
    ggplot2::scale_fill_manual(values=c('-1'=cols[1],'0'=cols[2],'1'=cols[3]),
                      labels = legend_labels,
                      guide = ggplot2::guide_legend(reverse = TRUE, title='')) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.background = ggplot2::element_rect(fill = "white"),
          panel.grid = ggplot2::element_line(size = 0),
          axis.title = ggplot2::element_text(size=ggplot2::rel(0)),
          axis.text = ggplot2::element_text(size=ggplot2::rel(1.5)),
          legend.text = ggplot2::element_text(size=ggplot2::rel(1.1)),
          legend.key.height=ggplot2::unit(2,"line"),
          legend.key.width=ggplot2::unit(2,"line"))
  plot(myplot)
}


#' Performs a search over different proportions of node-based (indirect) protein complex changes in
#' a subnetwork simulation, and returns a matrix in which the rows are all condition x direction combinations,
#' and the columns correspond to the given proportions of node-based (indirect) protein complex changes
#' @param pca_universe A PCA file corresponding to all protein-protein interactions
#' available in a given condition
#' @param pca_enhanced A PCA file corresponding to all significantly enhanced
#' protein-protein interactions 
#' @param pca_depleted A PCA file corresponding to all significantly depleted
#' protein-protein interactions 
#' @param iterations number of iterations
#' @param node_probs a vector indicating the proportions of node-based interaction changes to be simulated
#' @param metric the metric to be tested for significance, e.g. largest component size or graph density
#' @param conditions conditions to be tested, defaults to all in pca_universe
#' @param excluded_condition_grep a regular expression which matches all conditions to be excluded from the matrix
#' @param n_parallel number of parallel processes to spawn for the simulations
#' @param load_saved if True, attempts to load a previous run instead of executing from the path specified by save_directory and save_filename
#' @param save_output if True, will save the output to the path specified by save_directory and save_filename
#' @param save_directory the directory of where the files will be loaded/saved
#' @param save_filename the filename of where the files will be loaded/saved
#' @param seed seed for the random number generator
#'
#' @return a matrix where the rows correspond to conditions x directions and the columns are given proportions of node-based changes,
#' values are either 1 (significant, real value is greater than simulations), 0 (nonsigificant), -1
#' (significant, real value is less than simulations)
network_simulation_significance_node_edge_search_matrix <- function(pca_universe,
                                                                    pca_enhanced,
                                                                    pca_depleted,
                                                                    iterations=1000,
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
      val_list_depleted <- c(val_list_depleted,(empirical_two_tailed_p(observed_val_depleted,shuffled_vals_depleted)+1e-04)*sign(observed_val_depleted-mean(shuffled_vals_depleted)))
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

#' A histogram of edgewise-shuffled largest component size iterations
#'
#' @param pca_universe A PCA file corresponding to all protein-protein interactions
#' available in a given condition
#' @param pca_enhanced A PCA file corresponding to all significantly enhanced
#' protein-protein interactions
#' @param pca_depleted A PCA file corresponding to all significantly depleted
#' protein-protein interactions
#' @param condition condition, a string
#' @param iterations number of iterations to plot in the histogram
#' @param metric the metric to be tested, defaults to largest connected component size
#' @param seed seed for the random number generator
#'
#' @return plots a histogram
connectivity_histogram <- function(pca_universe,
                                   pca_enhanced,
                                   pca_depleted,
                                   condition,
                                   iterations=10000,
                                   metric = get_largest_component_from_edgelist,
                                   seed=123,
                                   text_size=1.5,
                                   xlim=c(0,25),
                                   breakstep=1,
                                   xlab='Simulated largest component',
                                   sampling_mode='edgewise',
                                   prob_node=0){
  set.seed(seed)
  par(mfrow=c(1,2),
      mar=c(6,4.5,3,0),
      oma=c(0,0,0,0))
  
  edgelist_enhanced <- as.matrix(dplyr::filter(pca_enhanced,Condition==condition) 
                        %>% dplyr::select(ORF.1,ORF.2))
  edgelist_depleted <- as.matrix(dplyr::filter(pca_depleted,Condition==condition) 
                                 %>% dplyr::select(ORF.1,ORF.2))
  
  
  largest_component_enhanced <- metric(edgelist_enhanced)
  largest_component_depleted <- metric(edgelist_depleted)
  shuffled_iters <- unlist(make_network_iterations(pca_universe,
                                                   pca_enhanced,
                                                   pca_depleted,
                                                   condition,
                                                   n_iters=iterations,
                                                   metric = metric,
                                                   prob_node=prob_node,
                                                   mode=sampling_mode))
  print(shuffled_iters[names(shuffled_iters)=='enhanced'])
  #plot(density(shuffled_iters,from=0,to=largest_component+largest_component*0.1),xlim=c(min(shuffled_iters),largest_component+largest_component*0.05))
  my_hist_enh <- hist(shuffled_iters[names(shuffled_iters)=='enhanced'],
                      breaks=seq(0,max(shuffled_iters) + 0.1*max(shuffled_iters),by=breakstep),
                      cex.lab=text_size,
                      plot=F)
  
  my_hist_depl <- hist(shuffled_iters[names(shuffled_iters)=='depleted'],
                       breaks=seq(0,max(shuffled_iters) + 0.1*max(shuffled_iters),by=breakstep),
                       cex.lab=text_size,
                       plot=F)
  histlist <- list('accumulated'=list('histvar'=my_hist_enh,
                                      'largest_component'=largest_component_enhanced),
                   'depleted'=list('histvar'=my_hist_depl,
                                   'largest_component'=largest_component_depleted))
  for(i in 1:length(histlist)){
    hist <- histlist[[i]]
    my_hist <- hist$histvar
    direction <- names(histlist)[i]
    plot(my_hist,
         xlim=xlim,
         xlab=xlab,
         ylim=c(0,max(my_hist$counts)*1.33),
         main='',
         ylab='Frequency',
         col='black',
         cex.lab=1.3)
    
    
    abline(v=hist$largest_component,lwd=6,col=rgb(1,0,0,0.7))
    text(hist$largest_component,mean(c(max(my_hist$counts),max(my_hist$counts)*1.3)),'Observed',srt=90,adj=c(0.5,1.5),col='gray60',cex=1.3,font=2)
    mtext(paste(c(Hmisc::capitalize(condition),direction,'\ncomplexes'),collapse=' '),side=3,cex=1.5,line=-1)
  }
  #}
}



#' Returns a heatmap which plots the output from network_simulation_significance_node_edge_search_matrix
#' @param node_edge_sig_matrix the output from network_simulation_significance_node_edge_search_matrix,
#' where the rows correspond to conditions x directions and the columns are given proportions of node-based changes,
#' values are either 1 (significant, real value is greater than simulations), 0 (nonsigificant), -1
#' (significant, real value is less than simulations) 
#' @param my_color_function a gradient generating function, the first value will be used for -1 values in matrix,
#' the middle for 0 values, the last for 1 values
#' @param legend_labels vector of form c('-1'=label1,'0'=label2,'1'=label3,2=label4)
#' @param x_label the title to be plotted on the X axis
#' @param y_label the title to be plotted on the Y axis
#' @param border_colour border outline colour
#' @param border_size border outline width
#' @param conditions_plotted conditions to include, if NULL plots all conditions
#'
#' @return plots a heatmap
node_edge_search_heatmap <- function(node_edge_sig_matrix,
                                     condition_delimiter=' ',
                                     conditions_plotted=NULL,
                                     collapse_enh_depl=T,
                                     my_color_function,
                                     legend_labels=c('-1'="Large component sizes (p < 0.05)",
                                                     '0'="Expected component sizes",
                                                     '1'="Small component sizes (p < 0.05)",
                                                     '2'="Most consistent simulation"),
                                     x_label='',
                                     y_label='Proportion of Protein-\nCentered Complex Changes',
                                     border_colour='grey40',
                                     zero_colour='grey90',
                                     border_size=0.25,
                                     width=0.7,
                                     top_label='Protein-Centric\nModel',
                                     bottom_label='Interaction-Specific\nModel',
                                     right_margin=11){
  
  
  
    
  my_matr <- node_edge_sig_matrix
  print('here')
  print(my_matr)
  my_matr.m <- reshape2::melt(my_matr)
  colnames(my_matr.m) <- c('Condition','PercentNodes','value')
  my_matr.m[,3] <- as.character(my_matr.m[,3])
  my_matr.m[,2] <- as.character(my_matr.m[,2])
  
  ##Ggvis2 version
#   my_matr.m %>%
#     ggvis(~PercentNodes,~Condition,fill=~value) %>%
#     layer_rects(width = band(), height = band()) %>%
#     scale_nominal("x", padding = 0, points = FALSE) %>%
#     scale_nominal("y", padding = 0, points = FALSE) %>%
#     add_axis("x",
#              offset=-10,
#              orient = "top",
#              title = 'Proportion of Node-Based Changes',
#              title_offset = 50,
#              properties=axis_props(
#                axis=list(strokeWidth=0),
#                grid=list(strokeWidth=0),
#                ticks=list(strokeWidth=0),
#                labels=list(angle=90,fontSize=15,align="right"),
#                title=list(fontSize=20)
#              )) %>%
#     add_axis("y",
#              title_offset = 170,
#              properties=axis_props(
#                axis=list(strokeWidth=0),
#                ticks=list(strokeWidth=0),
#                labels=list(fontSize = 15),
#                title=list(fontSize=20))
#              )
  
  cols <- my_color_function(10)
  if(is.null(zero_colour)){
    zero_colour = cols[3]
  }
  myplot <- ggplot2::ggplot(data = my_matr.m, ggplot2::aes(x=Condition, y=PercentNodes, fill=factor(value))) + 
   ggplot2::geom_tile(color=border_colour,size=border_size,width=width) +
   ggplot2::scale_fill_manual(values=c('-1'=cols[9],'0'=zero_colour,'1'=cols[2],'2'='black'),
                              labels = legend_labels,
                              guide = ggplot2::guide_legend(reverse = TRUE, title='')) +
   ggplot2::xlab(x_label) +
   ggplot2::ylab(y_label) + 
   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                  panel.background = ggplot2::element_rect(fill = "white"),
                  panel.grid = ggplot2::element_line(size = 0),
                  axis.title = ggplot2::element_text(size=ggplot2::rel(1.5)),
                  axis.text = ggplot2::element_text(size=ggplot2::rel(1.1)),
                  legend.text = ggplot2::element_text(size=ggplot2::rel(1.1)),
                  legend.key.height=ggplot2::unit(1.5,"line"),
                  legend.key.width=ggplot2::unit(1.5,"line"),
                  legend.position='right',
                  plot.margin = ggplot2::unit(c(1,right_margin,1,1), "lines")) +
  ggplot2::coord_equal()
  #ggplot2::coord_flip()
  myplot <- myplot + ggplot2::annotation_custom(grob = grid::textGrob(label=bottom_label,just='left'),
                                               ymin=1,
                                               ymax=1,
                                               xmin=nrow(my_matr)+1,
                                               xmax=nrow(my_matr)+1)
  
  myplot <- myplot + ggplot2::annotation_custom(grob = grid::textGrob(label=top_label,just='left'),
                                       ymin=ncol(my_matr),
                                       ymax=ncol(my_matr),
                                       xmin=nrow(my_matr)+1,
                                       xmax=nrow(my_matr)+1)
  
  p <- myplot
  gt <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid::grid.draw(gt)
  #return(myplot)
}

filter_matching_conditions <- function(node_edge_sig_matrix,conditions_plotted){
  condnames <- rownames(node_edge_sig_matrix)
  match_matr <- sapply(condnames,function(x){
    sapply(conditions_plotted,function(condition){
      sum(grepl(condition,x))
    })
  })
  return(node_edge_sig_matrix[apply(match_matr,2,sum) > 0,])
}

process_node_edge_sig_matrix <- function(node_edge_sig_matrix,conditions_plotted,collapse_enh_depl=T){
  #Filter out given conditions e.g. if they have only a few changes
  if(!is.null(conditions_plotted)){
    node_edge_sig_matrix <- filter_matching_conditions(node_edge_sig_matrix,conditions_plotted)
  }
  if(collapse_enh_depl == T){
    node_edge_sig_matrix <- t(sapply(conditions_plotted,function(condition){
      matching_rows <- node_edge_sig_matrix[grep(condition,rownames(node_edge_sig_matrix)),]
      apply(matching_rows,2,function(x){x[which.min(abs(x))]})
    }))
  }
  
  #For colour encoding
  #eyo <<- node_edge_sig_matrix
  node_edge_sig_matrix[t(apply(node_edge_sig_matrix,1,function(x){abs(x)==max(abs(x))}))] <- 4
  node_edge_sig_matrix[node_edge_sig_matrix < 0.05 & node_edge_sig_matrix > 0] <- 2
  node_edge_sig_matrix[node_edge_sig_matrix > -0.05 & node_edge_sig_matrix < 0] <- 3
  node_edge_sig_matrix[!(node_edge_sig_matrix %in% c(2,3,4))] <- 0
  node_edge_sig_matrix[node_edge_sig_matrix == 2] <- 1
  node_edge_sig_matrix[node_edge_sig_matrix == 3] <- -1
  node_edge_sig_matrix[node_edge_sig_matrix == 4] <- 2
  
  return(node_edge_sig_matrix)
}


hub_bias_heatmap <- function(hub_df,
                             color_function,
                             legend_labels=c('-1'="Complex\ndepletion bias",
                                             '0'="Non-significant\nbias",
                                             '1'="Complex\naccumulation bias"),
                             border_colour= 'grey40',
                             border_size = 0.25,
                             legend_position=c(1,1),
                             legend_justification=c(0,1),
                             nonsig_colour='grey10',
                             q_cutoff = 0.05
                             
){
  
  #Format hub dataframe for plotting
  new_hub_df <- select(hub_df,Hub,Condition,q.value.BH.)
  new_hub_df[,3] <- as.numeric(new_hub_df[,3] < q_cutoff)
  new_hub_df[,3] <- new_hub_df[,3]*sign(hub_df$Delta)
  new_hub_df[,3][is.na(hub_df[,3])] <- 0
  colnames(new_hub_df)[3] <- 'value'
  
  
  
  #Number of nonzero conditions
  sorted_new_hub_df <- as.data.frame(new_hub_df %>% dplyr::group_by(Hub) %>% dplyr::summarize(nzero=sum(abs(value))))
  hub_count_list <- unlist(apply(sorted_new_hub_df,1,function(x){
    retval <- list()
    retval[x[1]] <- as.numeric(x[2])
    #names(retval) <- x[1]
    return(retval)
  }))
  
  sorted_hubs <- sorted_new_hub_df[,'Hub'][sort(sorted_new_hub_df[,'nzero'],index.return=T,decreasing=F)$ix]
  
  sorted_new_hub_df <- as.data.frame(new_hub_df %>% dplyr::group_by(Condition) %>% dplyr::summarize(nzero=sum(abs(value))))
  sorted_conditions <- sorted_new_hub_df[,'Condition'][sort(sorted_new_hub_df[,'nzero'],index.return=T,decreasing=T)$ix]
  
  
  new_hub_df[,'Hub'] <- factor(new_hub_df[,'Hub'],levels=sorted_hubs)
  new_hub_df[,'Condition'] <- factor(new_hub_df[,'Condition'],levels=sorted_conditions)
  new_hub_df <- dplyr::filter(new_hub_df, Hub %in% names(which(hub_count_list > 0)))
  new_hub_df[,'value'] <- as.factor(new_hub_df[,'value'])
  #new_hub_df <- as.data.frame(t(new_hub_df))
  
  
  
  cols <- color_function(10)
  myplot <- ggplot2::ggplot(data = new_hub_df, ggplot2::aes(y=Hub, x=Condition, fill=value)) + 
    ggplot2::geom_tile(ggplot2::aes(fill = value),color=border_colour,size=border_size) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_manual(
      labels = legend_labels,
      values=c('-1'=cols[2],'0'=nonsig_colour,'1'=cols[9]),
      guide = ggplot2::guide_legend(title='')) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(
      #face="bold",
      size=15,
      vjust=0),
      axis.title.y = ggplot2::element_text(size=20
                                           #face="bold"
      ),
      axis.title.x = ggplot2::element_text(size=50
                                           #face="bold"
      ),
      
      axis.text.x  = ggplot2::element_text(angle=90,
                                           hjust=1,
                                           vjust=0.5,
                                           size=10),
      axis.text.y = ggplot2::element_text(size=10),
      
      legend.text = ggplot2::element_text(size=ggplot2::rel(1)),
      legend.position=legend_position,
      legend.justification=legend_justification,
      legend.key.width=unit(2, "lines"),
      legend.key.height=unit(2, "lines"),
      plot.margin = unit(c(0.5, 7, 0, 0), "lines"))
  
  plot(myplot)
}

bias_over_conditions <- function(hub_df,
                                 p_cutoff=0.05,
                                 title="Proportion of concerted hubs\nover different conditions"){
  conditions <- levels(hub_df$Condition)
  hub_proportion <- sapply(conditions,function(condition){
    subset_df <- dplyr::filter(hub_df,Condition==condition)
    all_hubs <- unique(subset_df$Hub)
    cond_hubs <- unique(filter(subset_df,q.value.BH. < p_cutoff)$Hub)
    return(length(cond_hubs)/length(all_hubs))
  })
  
  hub_proportion <- as.data.frame(hub_proportion)
  hub_proportion <- cbind(rownames(hub_proportion),hub_proportion)
  colnames(hub_proportion) <- c('Condition','Proportion of Concerted Hubs')
  
  order <- sort(hub_proportion$`Proportion of Concerted Hubs`,
                index.return=T,
                decreasing=T)$ix
  #hub_proportion <- hub_proportion[order,]
    
  hub_proportion$Condition <- factor(hub_proportion$Condition,
                                     levels=hub_proportion$Condition[order])
  #levels(hub_proportion$Condition) <- hub_proportion$Condition[order]
  
  
  ggplot2::ggplot(data=hub_proportion,ggplot2::aes(x=Condition,
                                 y=`Proportion of Concerted Hubs`,width=.7)) + 
    ggplot2::ggtitle(title) +
    ggplot2::ylab('Proportion of Concerted Hubs\n') + 
    ggplot2::xlab('\nConditions') +
    ggplot2::geom_bar(stat="identity",fill="snow4",colour="black") +
    ggplot2::scale_y_continuous(expand = c(0,0),
                                limits=c(0,max(hub_proportion$`Proportion of Concerted Hubs`)*1.05)) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white"),
      legend.position=c(0.8,0.8),
      #legend.text = ggplot2::element_text(size=15),
      legend.title = ggplot2::element_text(size=20,hjust=0),
      text = ggplot2::element_text(size=20),
      axis.text.x = ggplot2::element_text(angle=90, hjust=1,vjust=0.5),
      axis.line.x = ggplot2::element_line(size = 1, linetype = "solid", colour = "black"),
      axis.line.y = ggplot2::element_line(size = 1, linetype = "solid", colour = "black"))
  
}
