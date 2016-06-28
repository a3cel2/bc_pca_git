devtools::use_package('dplyr')
devtools::use_package('grDevices')
devtools::use_package('gplots')
devtools::use_package('cba')


convert_pca_file_to_heatmap_format <- function(pca_file,
                                               pca_enhanced_calls,
                                               pca_depleted_calls,
                                               color_function=blue_black_orange,
                                               min_heat_val=-1,
                                               max_heat_val=1,
                                               n_breaks=160,
                                               gene_dendrogram_width=0.25,
                                               condition_dendrogram_height=0.15,
                                               label_size=5,
                                               line_width=5,
                                               label_angle=90,
                                               centerings=c('gene','experiment'),
                                               centering_type=median,
                                               draw=T,
                                               output_path='dummy',
                                               filename='dummy'){
  
  
  blue_black_orange <- grDevices::colorRampPalette(c(
    rgb(1,0.45,0.25),
    rgb(0.8,0.25,0.25),
    rgb(0,0,0),
    rgb(0.25,0.45,0.8),
    rgb(0.25,0.75,1)
  ))
  
  sig_pca <- rbind(pca_enh,pca_depl)
  #sig_pca <- dplyr::filter(sig_pca, Gene.1 == 'MUP1' | Gene.2 == 'MUP1')
  
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
  
  col_dist <- uncentered_correlation_dist(t(response_matr))
  col_clust <- hclust(col_dist,method='average')
  
  row_dist <- uncentered_correlation_dist(response_matr)
  row_clust <- hclust(row_dist,method='average')
  col_order <- cba::order.optimal(col_dist,col_clust$merge)$order
  row_order <- cba::order.optimal(row_dist,row_clust$merge)$order
  
  x_labels <- sapply(colnames(response_matr),function(name){strsplit(name,split='\\.')[[1]][1]})
  
  #layout(mat = rbind(c(0,3),c(2,1),c(0,4)), widths = c(0.5,4), heights = c(0.5,4,0.1))
  
  if(draw == F){
    grDevices::png(file=paste(c(output_path,filename),collapse='/'),width=1500,height=2000)
  }
  par(lwd=line_width)
  heatmap.2(response_matr,
            Rowv=row_order,
            Colv=col_order,
            scale='none',
            labCol=x_labels,
            labRow='',
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
            key.xlab='Log2(R)',
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
