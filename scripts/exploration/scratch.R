library(dplyr)
library(xlsx)

hub_enrichment_file = '/Users/Albi/Dropbox/barcoded-PCA/2015-09-02/Data for Figure 3D.xlsx'
hub_df <- xlsx::read.xlsx2(hub_enrichment_file,sheetName = "Sheet1", colClasses=c("character","character",rep("numeric",5)))
# new_hub_df <- select(hub_df,Hub,Condition,q.value.BH.)
# new_hub_df[,3] <- as.numeric(new_hub_df[,3] < 0.05)
# new_hub_df[,3] <- new_hub_df[,3]*sign(hub_df$Delta)
# new_hub_df[,3][is.na(hub_df[,3])] <- 0
# colnames(new_hub_df)[3] <- 'value'
# 
# 
# #Number of nonzero conditions
# sorted_new_hub_df <- as.data.frame(new_hub_df %>% dplyr::group_by(Hub) %>% dplyr::summarize(nzero=sum(abs(value))))
# hub_count_list <- unlist(apply(sorted_new_hub_df,1,function(x){
#   retval <- list()
#   retval[x[1]] <- as.numeric(x[2])
#   #names(retval) <- x[1]
#   return(retval)
# }))
# 
# sorted_hubs <- sorted_new_hub_df[,'Hub'][sort(sorted_new_hub_df[,'nzero'],index.return=T,decreasing=T)$ix]
# 
# sorted_new_hub_df <- as.data.frame(new_hub_df %>% dplyr::group_by(Condition) %>% dplyr::summarize(nzero=sum(abs(value))))
# sorted_conditions <- sorted_new_hub_df[,'Condition'][sort(sorted_new_hub_df[,'nzero'],index.return=T,decreasing=F)$ix]
# 
# 
# new_hub_df[,'Hub'] <- factor(new_hub_df[,'Hub'],levels=sorted_hubs)
# new_hub_df[,'Condition'] <- factor(new_hub_df[,'Condition'],levels=sorted_conditions)
# new_hub_df <- dplyr::filter(new_hub_df, Hub %in% names(which(hub_count_list > 0)))
# 
# legend_labels=c('-1'="Depletion-Biased\nInteraction Changes",
#                 '0'="Unbiased\nInteraction Changes",
#                 '1'="Accumulation-Biased\nInteraction Changes")
# 
 my_color_list <- c(
   rgb(1,0.45,0.25),
   rgb(0.8,0.25,0.25),
   rgb(0,0,0),
   rgb(0.25,0.45,0.8),
   rgb(0.25,0.75,1)
 )
 blue_black_orange <- grDevices::colorRampPalette(my_color_list)
# 
# cols <- blue_black_orange(3)
# border_colour <- 'black'
# border_size <- 0.25
# new_hub_df[,'value'] <- as.factor(new_hub_df[,'value'])
# #new_hub_df[,'value'][new_hub_df[,'value'] == 0] <- NA
# myplot <- ggplot2::ggplot(data = new_hub_df, ggplot2::aes(x=Hub, y=Condition, fill=value)) + 
#   ggplot2::geom_tile(ggplot2::aes(fill = value),color=border_colour,size=border_size) +
#   ggplot2::coord_equal() +
#   ggplot2::scale_fill_manual(
#     labels = legend_labels,
#     values=c('-1'=cols[1],'0'='grey10','1'=cols[3]),
#     guide = ggplot2::guide_legend(title='')) +
#   ggplot2::theme(axis.title.x = ggplot2::element_text(
#                                     #face="bold",
#                                     size=15,
#                                     vjust=0),
#         axis.title.y = ggplot2::element_text(size=15
#                                     #face="bold"
#                                     ),
#         
#         axis.text.x  = ggplot2::element_text(angle=90,
#                                  hjust=1,
#                                  vjust=0.5),
#         legend.text = ggplot2::element_text(size=ggplot2::rel(1.1)),
#         legend.position='bottom')
# 
#   #scale_colour_manual(values=cols)
#   #scale_fill_gradientn(colors=cols, guide="colorbar",na.value = "grey10")
# #ggplot2::scale_colour_gradient2()#colours=cols)
# myplot

hub_bias_heatmap <- function(hub_df,
                             color_function,
                             legend_labels=c('-1'="Depletion-Biased\nInteraction Changes",
                                             '0'="Unbiased\nInteraction Changes",
                                             '1'="Accumulation-Biased\nInteraction Changes"),
                             border_colour = 'black',
                             border_size = 0.25,
                             legend_position='bottom',
                             nonsig_colour='grey10'
  
){
  
  #Format hub dataframe for plotting
  new_hub_df <- select(hub_df,Hub,Condition,q.value.BH.)
  new_hub_df[,3] <- as.numeric(new_hub_df[,3] < 0.05)
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
  
  sorted_hubs <- sorted_new_hub_df[,'Hub'][sort(sorted_new_hub_df[,'nzero'],index.return=T,decreasing=T)$ix]
  
  sorted_new_hub_df <- as.data.frame(new_hub_df %>% dplyr::group_by(Condition) %>% dplyr::summarize(nzero=sum(abs(value))))
  sorted_conditions <- sorted_new_hub_df[,'Condition'][sort(sorted_new_hub_df[,'nzero'],index.return=T,decreasing=F)$ix]
  
  
  new_hub_df[,'Hub'] <- factor(new_hub_df[,'Hub'],levels=sorted_hubs)
  new_hub_df[,'Condition'] <- factor(new_hub_df[,'Condition'],levels=sorted_conditions)
  new_hub_df <- dplyr::filter(new_hub_df, Hub %in% names(which(hub_count_list > 0)))
  new_hub_df[,'value'] <- as.factor(new_hub_df[,'value'])
  
  cols <- color_function(3)
  myplot <- ggplot2::ggplot(data = new_hub_df, ggplot2::aes(x=Hub, y=Condition, fill=value)) + 
    ggplot2::geom_tile(ggplot2::aes(fill = value),color=border_colour,size=border_size) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_manual(
      labels = legend_labels,
      values=c('-1'=cols[1],'0'=nonsig_colour,'1'=cols[3]),
      guide = ggplot2::guide_legend(title='')) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(
      #face="bold",
      size=15,
      vjust=0),
      axis.title.y = ggplot2::element_text(size=15
                                           #face="bold"
      ),
      
      axis.text.x  = ggplot2::element_text(angle=90,
                                           hjust=1,
                                           vjust=0.5),
      legend.text = ggplot2::element_text(size=ggplot2::rel(1.1)),
      legend.position=legend_position)
  
  plot(myplot)
}