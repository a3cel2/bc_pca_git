devtools::use_package('reshape2')
devtools::use_package('ggplot2')
devtools::use_package('dplyr')
devtools::use_package('Hmisc')

my_color_list <- c(
  rgb(1,0.45,0.25),
  rgb(0.8,0.25,0.25),
  rgb(0,0,0),
  rgb(0.25,0.45,0.8),
  rgb(0.25,0.75,1)
)


keep_only_largest_connected_component <- function(graph){
  components <- igraph::clusters(graph)
  largest_component <- which.max(components$csize)
  vertices_in_other_components <- names(which(components$membership != largest_component))
  new_graph <- igraph::delete.vertices(graph,vertices_in_other_components)
  return(new_graph)
}



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
                                    text_size_constant=30
){
  igraph::V(graph)$color <- node_colour
  igraph::V(graph)$label.cex <- ((node_size/text_size_constant)*4)/sapply(igraph::V(graph)$name,nchar)
  igraph::V(graph)$label.family <- "Arial Black"
  igraph::V(graph)$label.color <- node_text_colour
  igraph::E(graph)$color <- set_colours(edge_attribute,edge_colour_scale,edge_ncolours,edge_colour_min,edge_colour_max)
  igraph::E(graph)$width <- edge_width
  return(graph) 
}




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
                                       seed=102,
                                       edit=F,
                                       load_saved=T,
                                       layout_save_dir='',
                                       layout_file='default.layout',
                                       my_title='Doxorubicin directional connectivity'){
  condition_enhanced <- dplyr::filter(pca_enhanced,Condition==condition)
  condition_depleted <- dplyr::filter(pca_depleted,Condition==condition)
  
  
  #Filter by largest connected component separately
  g1 <- igraph::graph_from_edgelist(as.matrix(condition_enhanced[,c('Gene.1','Gene.2')]),directed=T)
  g2 <- igraph::graph_from_edgelist(as.matrix(condition_depleted[,c('Gene.1','Gene.2')]),directed=T)
  g1_n <- keep_only_largest_connected_component(g1)
  g2_n <- keep_only_largest_connected_component(g2)
  g <- igraph::graph.union(g1_n,g2_n)
  
  #Get changes for new graph
  edges <- igraph::get.edgelist(g)
  colnames(edges) <- c('Gene.1','Gene.2')
  combined_vals <- rbind(condition_enhanced,condition_depleted)
  combined_vals <- merge(edges,combined_vals,by=c('Gene.1','Gene.2'),sort=F)
  
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
  #layout(t(matrix(c(1,2))),widths=c(0.45,4))
  par(oma=c(0,0,0,0),mar=c(0,3,2,0))
  set.seed(seed)
  layout_file_path <- paste(c(layout_save_dir,layout_file),collapse='/')
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
  plot(g,layout=l,vertex.size=node_size,edge.width=edge_width)
  hub_legend_draw(min_fc,max_fc,my_color_list,main_font_size = 0.7,side_font_size = 0.7 ,width=0.2,new_plot=F,x_adjust=-1,y_adjust=-1)
  
  #par(mar=c(0,0,2,0))
  #par(oma=c(0,0,0,0))
  
  par(cex=1.5)
  title(my_title)
}


get_largest_component_from_edgelist <- function(edgelist){
  graph <- igraph::graph_from_edgelist(edgelist)
  return(max(igraph::clusters(graph)$csize))
}

get_largest_component_size <- function(pca_file,
                                       condition){
  edgelist <- as.matrix(dplyr::filter(pca_file,Condition==condition) 
                        %>% dplyr::select(ORF.1,ORF.2))
  return(get_largest_component_from_edgelist(edgelist))
}

make_connectivity_iterations <- function(pca_universe,pca_file,condition,n_iters){
  real_edgelist <- as.matrix(dplyr::filter(pca_file,Condition==condition) 
                             %>% dplyr::select(ORF.1,ORF.2))
  n_edges <- nrow(real_edgelist)
  potential_edgelist <- as.matrix(dplyr::filter(pca_universe,Condition==condition) 
                                  %>% dplyr::select(ORF.1,ORF.2))
  potential_nrow <- nrow(potential_edgelist)
  
  return(sapply(1:n_iters,function(n){
    edgelist <- potential_edgelist[sample(1:potential_nrow,n_edges,replace=F),,drop=F]
    get_largest_component_from_edgelist(edgelist)
  }))
}

network_connectivity_significance <- function(pca_universe,
                                              pca_enhanced,
                                              pca_depleted,
                                              color_list,
                                              excluded_condition_grep = 'CRISPR',
                                              iterations = 1000,
                                              seed=123,
                                              sigval=0.05){
  set.seed(seed)


  conditions <- unique(c(levels(pca_enhanced$Condition),levels(pca_depleted$Condition)))
  conditions <- grep(excluded_condition_grep,conditions,value=T,invert=T)
  
  record_list <- list()
  iteration_list <- list()
  return_matrix <- matrix(nrow=2,ncol=length(conditions),data=1)
  rownames(return_matrix) <- c('enhanced','depleted')
  colnames(return_matrix) <- conditions
  
  for(direction in c('enhanced','depleted')){
    record_list[[direction]] <- list()
    if(direction == 'enhanced'){
      pca_file <- pca_enhanced
    }
    else{
      pca_file <- pca_depleted
    }
    for(condition in conditions){
      largest_component <- get_largest_component_size(pca_file,condition)
      iteration_list[[direction]][[condition]] <- make_connectivity_iterations(pca_universe,pca_file,condition,n_iters=iterations)
      record_list[[direction]][[condition]] <- largest_component
      return_matrix[direction,condition] <- sum(iteration_list[[direction]][[condition]] >= record_list[[direction]][[condition]])/iterations
    }
  }
  
  
  return_matrix['depleted',][which(return_matrix['depleted',] <= sigval)] <- -1.01
  return_matrix['enhanced',][which(return_matrix['enhanced',] <= sigval)] <- 1.01
  return_matrix[which(!(return_matrix %in% c(-1.01,1.01)))] <- 0
 
  return_matrix[return_matrix == -1.01] <- -1
  return_matrix[return_matrix == 1.01] <- 1
  return(return_matrix)

}




connectivity_graph <- function(my_matr,
                               my_color_list,
                               border_colour="grey30",
                               border_size=0.8){
  vals <- sort(apply(my_matr,2,function(x){sum(x == -1)+2*sum(x == 1)}),decreasing=T)
  level_order <- names(vals)
  cols <- grDevices::colorRampPalette(my_color_list)
  cols <- cols(3)
  my_matr.m <- reshape2::melt(my_matr)
  colnames(my_matr.m) <- c('Direction','Condition','value')
  my_matr.m$Condition <- factor(my_matr.m$Condition,
                                levels=level_order)
  my_matr.m$Direction <- factor(my_matr.m$Direction,
                                levels=c('depleted','enhanced'))
  
  myplot <- ggplot2::ggplot(data = my_matr.m, ggplot2::aes(x=Condition, y=Direction, fill=factor(value))) + 
    ggplot2::geom_tile(color=border_colour,size=border_size) +
    ggplot2::scale_fill_manual(values=cols,
                      labels = c("Increased connectivity\n(depleted complexes)",
                                 "Expected connectivity",
                                 "Increased connectivity\n(enhanced complexes)"),
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

connectivity_histogram <- function(pca_universe,
                                   pca_enhanced,
                                   pca_depleted,
                                   condition,
                                   iterations=10000,
                                   seed=123){
  set.seed(seed)
  par(mfrow=c(1,2),
      mar=c(5,4.5,3,1),
      oma=c(0,0,0,0))
  for(direction in c('enhanced','depleted')){
    if(direction == 'enhanced'){
      pca_file <- pca_enhanced
    }
    else{
      pca_file <- pca_depleted
    }
    largest_component <- get_largest_component_size(pca_file,condition)
    shuffled_iters <- make_connectivity_iterations(pca_universe,pca_file,condition,n_iters=iterations)
    #plot(density(shuffled_iters,from=0,to=largest_component+largest_component*0.1),xlim=c(min(shuffled_iters),largest_component+largest_component*0.05))
    my_hist <- hist(shuffled_iters,
         breaks=seq(0,max(shuffled_iters),by=0.5),
         xlim=c(0,25),
         xlab='Connectivity (shuffled)',
         main='',
         ylab='Frequency',
         col='gray20')
    abline(v=largest_component,lwd=2,lty=4,col='red')
    text(largest_component,mean(c(0,max(my_hist$counts))),'Observed',srt=90,adj=c(0.5,1.5),col='gray60')
    mtext(paste(c(Hmisc::capitalize(condition),direction,'\ncomplexes'),collapse=' '),side=3,cex=1.5)
  }
}
