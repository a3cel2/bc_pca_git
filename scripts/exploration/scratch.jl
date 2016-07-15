function testfun(x,y)
  x*y
end

test_df = readdlm("/Users/Albi/Dropbox/Roth Lab/projects/bc_pca_git/data/scratch/ethanol_universe_test.tsv",'\t',AbstractString)

#Creates an index of a gene universe for quick lookup
function index_universe(universe_array)
  df_index = Dict()
  for(i in 1:size(universe_array)[1])
    pair = test_df[i,:]
    #println(pair)
    for(gene in pair)
      if(get(df_index,gene,0) == 0)
        df_index[gene] = Array(AbstractString,0,2)
      end
      df_index[gene] = [df_index[gene];pair]
    end
  end
  return(df_index)
end


function nodewise_simulation(indexed_universe_edgelist,
  universe_edgelist,
  n_edges)

  remaining_edges = n_edges

  used_nodes = AbstractString[]
  edgelist = Array(AbstractString,0,2)

  potentially_valid_nodes = keys(indexed_universe_edgelist)
  while(remaining_edges > 0)
    
  end
end

indexed_universe_edgelist = index_universe(test_df)
universe_edgelist = test_df
n_edges = 2

println(nodewise_simulation(indexed_universe_edgelist,
universe_edgelist,
n_edges))
# nodewise_simulation <- function(indexed_universe_edgelist,
#                                 potential_edgelist,
#                                      n_edges,
#                                      node_sensitivity=1,
#                                      prob_node=1){
#   remaining_edges <- n_edges
#   used_nodes <- c()
#   edgelist <- c()
#
#   #Not filtering works better
#   to_be_added <- sapply(indexed_universe_edgelist,function(x){length(unique(x))})
#   maybe_valid <- names(indexed_universe_edgelist)
#   while(remaining_edges > 0){
#     node_trial <- runif(1)
#     new_edges <- sample_edges(node_trial,indexed_universe_edgelist,potential_edgelist,node_sensitivity,prob_node,maybe_valid)
#     test_edgelist <- unique(rbind(edgelist,new_edges))
#     test_remaining_edges <- n_edges - nrow(test_edgelist)
#     while(test_remaining_edges < 0 ){
#       new_edges <- sample_edges(node_trial,indexed_universe_edgelist,potential_edgelist,node_sensitivity,prob_node,maybe_valid)
#       test_edgelist <- unique(rbind(edgelist,new_edges))
#       test_remaining_edges <- n_edges - nrow(test_edgelist)
#     }
#     edgelist <- test_edgelist
#
#     #Don't remove node if you did an edgewise trial
#     if(prob_node == 1){
#       maybe_valid <- maybe_valid[!(maybe_valid %in% unique(as.vector(edgelist)))]
#     }
#     else if (prob_node < 1){
#       if(node_trial <= prob_node){
#         maybe_valid <- maybe_valid[!(maybe_valid %in% unique(as.vector(edgelist)))]
#       }
#     }
#     remaining_edges <- n_edges - nrow(test_edgelist)
#   }
#   return(edgelist)
# }
