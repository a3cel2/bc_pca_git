require(data.table)


keep_only_large_connected_components <- function(graph, size){
  components <- igraph::clusters(graph)
  largest_component <- which(components$csize >= size)
  vertices_in_other_components <- names(which(components$membership != largest_component))
  new_graph <- igraph::delete.vertices(graph,vertices_in_other_components)
  return(new_graph)
}

#Gets edges matching a given node
get_connected_edges <- function(universe_edgelist,added_node){
  query <- universe_edgelist[,1] %in% added_node | universe_edgelist[,2] %in% added_node
  return(universe_edgelist[query, , drop = F])
}


#A hack to remove rows
row_subtract  <- function(a1,a2)
{
  a1.vec <- apply(a1, 1, paste, collapse = "")
  a2.vec <- apply(a2, 1, paste, collapse = "")
  a1.without.a2.rows <- a1[!a1.vec %in% a2.vec,]
  return(a1.without.a2.rows)
}

#How many extra nodes will be added to the component?
get_extra_nodes <- function(universe_edgelist,
                            used_nodes,
                            connected=T){
  if(connected == T){
    query <- universe_edgelist[,1] %in% used_nodes | universe_edgelist[,2] %in% used_nodes
  }
  else{
    query <- 1:nrow(universe_edgelist)
  }
  #Which are valid nodes?
  #If starting, all nodes are valid
  if(length(used_nodes) > 0){
    remaining_edgelist <- universe_edgelist[query,]
    remaining_nodelist <- unique(as.vector(remaining_edgelist))
    remaining_nodelist <- setdiff(remaining_nodelist,used_nodes)
  }
  #If nodes already used, have to get the connected nodes
  else{
    remaining_nodelist <- unique(as.vector(universe_edgelist))
  }
  
  #How much will adding a given node increase our connected component size?
  extra_genes <- sapply(remaining_nodelist,function(node){
    query <- universe_edgelist[,1] %in% node | universe_edgelist[,2] %in% node
    extra_genes <- unique(as.vector(universe_edgelist[query,]))
    return(length(setdiff(extra_genes,used_nodes)))
  })
  return(extra_genes)
}


#Creates a component of n nodes from an edgelist by randomly sampling connected nodes and extracting all of their components
#If connected is set to true, components are all connected
simulate_component_nodewise <- function(universe_edgelist,
                                        excluded_nodes,
                                        component_size,
                                        connected=T
                                        ){

  #universe_edgelist <- universe_edgelist[!(universe_edgelist[,1] %in% excluded_nodes | universe_edgelist[,2] %in% excluded_nodes),]
  used_edgelist <- c()
  used_nodelist <- c()
  extra_components <- get_extra_nodes(universe_edgelist,used_nodelist)
  remaining_components <- component_size
  
  #Step for first node is different
  initial_candidate_nodes <- names(extra_components[extra_components <= remaining_components])
  added_node <- sample(initial_candidate_nodes,1)
  added_edges <- get_connected_edges(universe_edgelist,added_node)
  added_nodes <- unique(as.vector(added_edges))
  
  used_edgelist <- rbind(used_edgelist,added_edges)
  used_nodelist <- c(used_nodelist,added_nodes)
  
  remaining_components <- remaining_components - length(added_nodes)
  while(remaining_components > 0){
    extra_components <- get_extra_nodes(universe_edgelist,
                                       used_nodelist,
                                       connected=connected)
    #Cannot use node twice
    candidate_nodes <- names(extra_components[extra_components <= remaining_components])
    candidate_nodes <- setdiff(candidate_nodes,added_node)
    
    #Restart if in dead end
    if(length(candidate_nodes) == 0){
      print('Fucked up')
      used_edgelist <- c()
      used_nodelist <- c()
      added_node <- sample(initial_candidate_nodes,1)
      added_edges <- get_connected_edges(universe_edgelist,added_node)
      added_nodes <- unique(as.vector(added_edges))
      
      used_edgelist <- rbind(used_edgelist,added_edges)
      used_nodelist <- unique(c(used_nodelist,added_nodes))
    }
    else{
      added_node <- sample(candidate_nodes,1)
      added_edges <- get_connected_edges(universe_edgelist,added_node)
      added_nodes <- unique(as.vector(added_edges))
      
      used_edgelist <- rbind(used_edgelist,added_edges)
      used_nodelist <- unique(c(used_nodelist,added_nodes))
    }
    remaining_components <- component_size - length(used_nodelist)
  }
  return(used_edgelist)
}

connected_component_simulation <- function(graph,
                                           universe_graph,
                                           mode='nodewise'){
  connected_component_list <- sort(igraph::clusters(graph)$csize,decreasing=T)
  used_nodes <- c()
  final_graph <- c()
  universe_edgelist <- get.edgelist(universe_graph)
  for(i in connected_component_list){
    #Don't connect to previous components
    universe_edgelist <- universe_edgelist[!(universe_edgelist[,1] %in% used_nodes | universe_edgelist[,2] %in% used_nodes),]
    
    component <- simulate_component_nodewise(universe_edgelist,used_nodes,i)
    used_nodes <- c(used_nodes,unique(as.vector(component)))
    final_graph <- rbind(final_graph,component)
  }
  return(final_graph)
}

#A fast nodewise simulation if connected component preservation does not matter
fast_nodewise_simulation <- function(indexed_universe_edgelist,
                                     n_nodes){
  remaining_nodes <- n_nodes
  remaining_edges <- n_edges
  used_nodes <- c()
  edgelist <- c()
  
  #Not filtering works better
  #to_be_added <- sapply(indexed_universe_edgelist,function(x){length(unique(x))})
  maybe_valid <- names(indexed_universe_edgelist)
  while(remaining_nodes > 0){
    test_edgelist <- rbind(edgelist,indexed_universe_edgelist[[sample(maybe_valid,1)]])
    test_remaining_nodes <- n_nodes - length(unique(as.vector(test_edgelist)))
    while(test_remaining_nodes < 0){
      test_edgelist <- rbind(edgelist,indexed_universe_edgelist[[sample(maybe_valid,1)]])
      test_remaining_nodes <- n_nodes - length(unique(as.vector(test_edgelist)))
    }
    edgelist <- test_edgelist
    remaining_nodes <- n_nodes - length(unique(as.vector(edgelist)))
  }
  return(edgelist)
}

edge_index_universe <- function(universe_edgelist){
  node_index <- list()
  for(i in 1:nrow(universe_edgelist)){
    node_index[[universe_edgelist[i,1]]] <- rbind(node_index[[universe_edgelist[i,1]]],universe_edgelist[i,])
  }
  return(node_index)
}
  
  
edge_uni <- get.edgelist(g_universe)
uni_index <- index_universe(edge_uni)
n_nodes <- length(V(g1))
node_simulation <- sapply(1:10,function(x){
  max(clusters(graph.edgelist(fast_nodewise_simulation(uni_index,85)))$csize)
  })