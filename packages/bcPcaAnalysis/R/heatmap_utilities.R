devtools::use_package('dplyr')
devtools::use_package('tidyr')
devtools::use_package('grDevices')
devtools::use_package('gplots')
devtools::use_package('cba')


uncentered_correlation <- function(x,y){
  sd_x_0 <- sqrt((1/length(x))*sum(x^2))
  sd_y_0 <- sqrt((1/length(y))*sum(y^2))
  
  z_x <- x/sd_x_0
  z_y <- y/sd_y_0
  
  return((1/length(x))*sum(z_x*z_y))
}

uncentered_correlation_dist <- function(x){
  x <- t(x)
  as.dist(apply(x,2,function(x1){
    apply(x,2,function(x2){
      return(1 - uncentered_correlation(x1,x2))
    })
  }))
}


process_matrix_for_heatmap <- function(pca_file,
                                       pca_enhanced_calls,
                                       pca_depleted_calls,
                                       centerings,
                                       centering_type){
  
  sig_pca <- rbind(pca_enhanced_calls,pca_depleted_calls)
  
  
  
  
  simplified_pca_file <- dplyr::select(pca_file,Tag,Condition,FC.avg,FC.rep1,FC.rep2,q.val)
  
  simplified_pca_file <- simplified_pca_file %>% 
    dplyr::group_by(Tag,Condition) %>% 
    dplyr::summarize(FC.avg=mean(FC.avg),FC.rep1=mean(FC.rep1),FC.rep2=mean(FC.rep2))
  
  simplified_pca_file <- simplified_pca_file %>% 
    dplyr::filter(Tag %in% sig_pca$Tag)  %>% 
    dplyr::select(Tag,Condition,FC.rep1,FC.rep2)
  
  rep1 <- simplified_pca_file %>% dplyr::select(-FC.rep2) %>% tidyr::spread(Condition,FC.rep1)
  rep2 <- simplified_pca_file %>% dplyr::select(-FC.rep1) %>% tidyr::spread(Condition,FC.rep2)
  
  
  response_matr <- merge(rep1,rep2,by='Tag')
  rownames(response_matr) <- response_matr$Tag
  response_matr <- as.matrix(response_matr[,-1])
  
  if('gene' %in% centerings){
    response_matr <- response_matr - apply(response_matr,1,centering_type)
  }
  if('experiment' %in% centerings){
    response_matr <- response_matr - apply(response_matr,2,centering_type)
  }
  return(response_matr)
}

convert_pca_file_to_heatmap_format <- function(pca_file,
                                               pca_enhanced_calls,
                                               pca_depleted_calls,
                                               color_function=blue_black_orange,
                                               min_heat_val=-1,
                                               max_heat_val=1,
                                               n_breaks=160,
                                               gene_dendrogram_width=0.25,
                                               condition_dendrogram_height=0.15,
                                               label_size=1,
                                               line_width=1,
                                               label_angle=90,
                                               centerings=c('gene','experiment'),
                                               centering_type=median,
                                               draw=T,
                                               output_path='dummy',
                                               filename='dummy',
                                               png_width=1500,
                                               png_height=2000){
  
  tag_to_pair <- list()
  for(i in 1:(nrow(pca_file)/length(unique(pca_file$Condition)))){
    tag_to_pair[as.vector(pca_file[i,]$Tag)] <- as.vector(pca_file[i,]$PPI.short)
  }
  
  blue_black_orange <- grDevices::colorRampPalette(c(
    rgb(1,0.45,0.25),
    rgb(0.8,0.25,0.25),
    rgb(0,0,0),
    rgb(0.25,0.45,0.8),
    rgb(0.25,0.75,1)
  ))
  
  response_matr <- process_matrix_for_heatmap(pca_file,
                                              pca_enh,
                                              pca_depl,
                                              centerings,
                                              centering_type)
  
  col_dist <- uncentered_correlation_dist(t(response_matr))
  col_clust <- hclust(col_dist,method='average')
  
  row_dist <- uncentered_correlation_dist(response_matr)
  row_clust <- hclust(row_dist,method='average')
  col_order <- cba::order.optimal(col_dist,col_clust$merge)$order
  row_order <- cba::order.optimal(row_dist,row_clust$merge)$order
  
  x_labels <- sapply(colnames(response_matr),function(name){strsplit(name,split='\\.')[[1]][1]})
  
  
  if(draw == F){
    grDevices::png(file=paste(c(output_path,filename),collapse='/'),width=png_width,height=png_height)
  }
  par(lwd=line_width)
  heatmap.2(response_matr,
            Rowv=row_order,
            Colv=col_order,
            scale='none',
            labCol=x_labels,
            labRow=tag_to_pair[rownames(response_matr)],
            breaks=seq(min_heat_val,max_heat_val,by=(max_heat_val-min_heat_val)/n_breaks),
            col=color_function(n_breaks),
            distfun=uncentered_correlation_dist,
            trace='none',
            lmat = rbind(c(4,0),c(0,3),c(2,1),c(0,0)),
            lwid = c(gene_dendrogram_width,1),
            lhei = c(0.1,condition_dendrogram_height,1,0.2),
            cexCol = label_size,
            srtCol = label_angle,
            key=T,
            key.par=list(cex.lab=label_size*0.8,cex.axis=label_size*0.8,mar=c(0,2,2,2),mgp=c(4,2,0)),
            key.xlab='Complex Log2(R)',
            key.title='',
            density.info = 'none',
            key.ylab='')
  if(draw == F){
    dev.off()
  }
  #col_median <- apply(response_matr,2,median)
  #row_median <- apply(response_matr,1,median)
  #response_matr <- re
  
}

condition_summary_barplot <- function(pca_enhanced,
                                      pca_depleted,
                                      color_function=gplots::redgreen,
                                      excluded_conditions='CRISPR',
                                      draw=T,
                                      output_path='dummy',
                                      filename='dummy',
                                      my_width=10,
                                      my_height=7){
  pca_enh <- pca_enhanced
  pca_enh$Condition <- as.vector(pca_enh$Condition)
  for(condition in excluded_conditions){
    pca_enh <- filter(pca_enh,!(grepl(condition,Condition)))
  }
  
  pca_depl <- pca_depleted
  pca_depl$Condition <- as.vector(pca_depl$Condition)
  for(condition in excluded_conditions){
    pca_depl <- filter(pca_depl,!(grepl(condition,Condition)))
  }
  
  enh_table <- table(pca_enh$Condition)
  depl_table <- table(pca_depl$Condition)
  
  summary_table <- merge(enh_table,depl_table,by='Var1',all=T)
  colnames(summary_table) <- c('Condition','Enhanced','Depleted')
  
  #Sort
  bar_order <- sort(apply(summary_table[,2:3],1,sum),index.return=T,decreasing=T)$ix
  summary_table$Condition <- factor(summary_table$Condition,
                                    levels=summary_table$Condition[bar_order])
  
  #Make wide
  summary_table <- tidyr::gather(summary_table,Condition)
  colnames(summary_table) <- c('Condition','Direction','Frequency')
  
  #Sort again, enhanced_first
  summary_table$Direction<- factor(summary_table$Direction,
                                    levels=c('Enhanced','Depleted'))
  
  #summary_table[,2] <- as.factor(summary_table[,2])
  #ggplot(data=summary_table, aes(x=Condition, y=Frequency, fill=Direction)) +
  #  geom_bar(stat="identity", position=position_dodge()) + scale_fill_manual(values=rev(blue_black_orange(2)))
  
  my_plot <- ggplot2::ggplot(data=summary_table, aes(x=Condition, y=Frequency, fill=Direction, width=.8)) +
    ggplot2::geom_bar(stat="identity", position=position_dodge(), colour='black') + 
    ggplot2::scale_y_continuous(expand = c(0,0), limits=c(0,max(summary_table$Frequency)*1.05)) +
    ggplot2::scale_fill_manual(values=rev(color_function(10)[c(2,9)])) +
    ggplot2::ylab('Dynamic Complexes') +
    ggplot2::theme(
          panel.background = element_rect(fill = "white"),
          legend.position=c(0.8,0.8),
          legend.text = element_text(size=15),
          legend.title = element_text(size=20,hjust=0),
          text = element_text(size=25),
          axis.text.x = element_text(angle=90, hjust=1),
          axis.line.x = element_line(size = 1, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 1, linetype = "solid", colour = "black"))
  if(draw==T){
    my_plot
  }
  else{
    filename <- paste(c(output_path,filename),collapse='/')
    ggplot2::ggsave(plot=my_plot,filename=filename,width=my_width,height=my_height,units='in')
  }
}
