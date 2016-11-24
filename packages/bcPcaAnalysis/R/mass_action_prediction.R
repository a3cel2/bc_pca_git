devtools::use_package('dplyr')
devtools::use_package('metap')
devtools::use_package('grDevices')
devtools::use_package('Cairo')
devtools::use_package('igraph')
devtools::use_package('gplots')
devtools::use_package('cba')

#protein_abundance_file = "/Users/Albi/Dropbox/Roth Lab/projects/bc_pca_git/data/paxdb_abundance.tsv"
#expression_file = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.14.txt'
#pca_universe = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.6.txt'




#' A simple mass action prediction for protein complex level changes
#'
#' @param c1 concentration of the first protein in the PPI pair, in arbitrary units
#' @param c2 concentration of the second protein in the PPI pair, in arbitrary units
#' @param r1 ratio of concentration change of the first protein
#' @param r2 ratio of concentration change of the second protein
#' @param k affinity estimate (in same units as c1,c2), set to min(c1,c2)/concentration factor by default
#' @param concentration_factor used for unknown Kd assignment, defaults to 20
#' @param mode either analytic for an exact solution (default) or iterative
#'
#' @return a vector representing the predicted complex level ratio
mass_action_predictor <- function(c1,
                                  c2,
                                  fc1,
                                  fc2,
                                  k=NA,
                                  concentration_factor=20,
                                  mode='analytic'){

  if(is.na(k)){
    k = min(c(c1,c2))/concentration_factor
  }
  if(mode == 'analytic'){
    old_pair_analytic = 0.5*(((-1)*sqrt((c1**2)-(2*c1*(c2-k))+(c2+k)**2))+c1+c2+k)
    c1 <- c1*(fc1)
    c2 <- c2*(fc2)
    new_pair_analytic = 0.5*(((-1)*sqrt((c1**2)-(2*c1*(c2-k))+(c2+k)**2))+c1+c2+k)
    return(new_pair_analytic/old_pair_analytic)
  }

  else if(mode == 'iterative'){
    f1 = 0
    f2 = 0
    f3 = 0
    f4 = 0

    for(i in 1:100){
      f1 <- c1/(1 + (f2/k) + (f4/k))
      f2 <- c2/(1 + (f1/k) + (f3/k))
      f3 <- c1/(1 + (f2/k) + (f4/k))
      f4 <- c2/(1 + (f1/k) + (f3/k))
    }
    old_pair <- (1/k)*f1*f2

    c1 <- c1*(fc1)
    c2 <- c2*(fc2)
    c3 <- c1*(fc1)
    c4 <- c2*(fc2)


    f1 = 0
    f2 = 0
    f3 = 0
    f4 = 0

    for(i in 1:100){
      f1 <- c1/(1 + (f2/k) + (f4/k))
      f2 <- c2/(1 + (f1/k) + (f3/k))
      f3 <- c1/(1 + (f2/k) + (f4/k))
      f4 <- c2/(1 + (f1/k) + (f3/k))
    }
    new_pair <- (1/k)*f1*f2

    return((new_pair/old_pair))
  }
}


#A simple function to look up abundance of both orfs in a pair,
#assigns median abundance in the abundance file if not present
get_orf_pair_abundance <- function(pair,abundance_table){
  abundance1 <- abundance_table[pair[1],]
  abundance2 <- abundance_table[pair[2],]
  if(is.na(abundance1)){
    abundance1 <- median(abundance_table[,1])
  }
  if(is.na(abundance2)){
    abundance2 <- median(abundance_table[,1])
  }
  return(c(abundance1,abundance2))
}


#Groups UPtag and DNtag measurements from a PCA file,
#Combines multiple measurements for the same PPI measured twice
#returns a new data frame
merge_pca_file <- function(pca_calls,condition){
  merged_pca_calls_up <- dplyr::filter(pca_calls, Condition == condition, UP.DN=='uptag')
  merged_pca_calls_dn <- dplyr::filter(pca_calls, Condition == condition, UP.DN=='downtag')

  if(!identical(rownames(merged_pca_calls_up),rownames(merged_pca_calls_dn))){
    stop('PCA file is not evenly divided into uptag and uptag, please provide unprocessed data')
  }

  new_fc_avg <- sapply(1:nrow(merged_pca_calls_up),function(i){
    mean(c(merged_pca_calls_up$FC.avg[i],merged_pca_calls_dn$FC.avg[i]))
  })

  new_q_val <- sapply(1:nrow(merged_pca_calls_up),function(i){
    unlist(metap::sumz(c(merged_pca_calls_up$q.val[i],merged_pca_calls_dn$q.val[i])))[2]
  })

  dmso_index <- grep('DMSO',colnames(merged_pca_calls_up))
  dmso_abundance_up <- sapply(1:nrow(merged_pca_calls_up),function(i){
    return(mean(as.numeric(merged_pca_calls_up[i,dmso_index])))
  })

  dmso_abundance_dn <- sapply(1:nrow(merged_pca_calls_up),function(i){
    return(mean(as.numeric(merged_pca_calls_dn[i,dmso_index])))
  })


  merged_pca_calls <- data.frame(ORF.1=merged_pca_calls_up$ORF.1,
                                   ORF.2=merged_pca_calls_up$ORF.2,
                                   AUC=merged_pca_calls_up$AUC,
                                   DMSO.UP=dmso_abundance_up,
                                   DMSO.DN=dmso_abundance_dn,
                                   FC.UP=merged_pca_calls_up$FC.avg,
                                   FC.DN=merged_pca_calls_dn$FC.avg,
                                   FC.avg=new_fc_avg,
                                   q.val=new_q_val)

  return(merged_pca_calls)
}


#Given an mRNA expression file, obtains expression changes for ORFs in a vector or matrix (output arranged same format as input)
get_orf_mrna_changes <- function(pairs,
                                      expression_file,
                                      expression_control_regexp,
                                      expression_condition_regexp,
                                      paired_index_data=F){
  ctrl_index <- grep(expression_control_regexp,
                     colnames(expression_file))
  cond_index <- grep(expression_condition_regexp,
                     colnames(expression_file))
  if(is.null(nrow(pairs))){
    return(sapply(1:length(pairs),function(i){
    mean(as.matrix(expression_file[pairs[i],cond_index]))/
      mean(as.matrix(expression_file[pairs[i],ctrl_index]))
    }))
  } else {
    return(sapply(1:ncol(pairs),function(i){apply(expression_file[pairs[,i],cond_index],1,mean)/apply(expression_file[pairs[,i],ctrl_index],1,mean)}))
  }
}


#Simplifies expression file into a numeric matrix with one measurement per ORF (mean aggregation)
#and ORF names used as indeces
simplify_expression_file <- function(expression_file){
  expression_file <- expression_file %>% dplyr::select(-ID_REF)
  expression_file <- aggregate(expression_file[,-1],list(expression_file$ORF),mean,na.rm=T)
  rownames(expression_file) <- expression_file[,1]
  expression_file <- expression_file[,-1]
  return(expression_file)
}


#' Provides mass-action based predictions for bcPCA output in a specific format
#'
#' @param pca_file filename corresponding to output from the PCA experiment,
#' unfiltered, with both UPtag and DNtag measurements available
#' @param abundance_file filename corresponding to the abundance of each ORF in the PCA
#' file, first column of file should be the ORF ID, the second column should be a measure of
#' abundance (in arbitrary units)
#' @param expression_file mRNA expression data from which the predictions are going to be made,
#' must have an ORF column corresponding to the ORF, and the experiments must be organized in columns
#' @param condition the condition in the pca_file which is to be correlated
#' @param expression_control_regexp regular expression used to match the control condition in the mRNA file
#' @param expression_condition_regexp regular expression used to match the condition of interesti n the mRNA file
#'
#' @return a data frame with the original data combined with the predictions
pca_ma_prediction <- function(
  pca_file,
  abundance_file,
  expression_file,
  condition,
  expression_control_regexp='Ethanol.0h',
  expression_condition_regexp='Ethanol.4h'){

  #Read appropriate fies
  pca_file <- read.csv(pca_file,head = T,stringsAsFactors = F, sep = '\t')
  expression_file <- read.table(expression_file,head=T,stringsAsFactors = F)
  abundance_file <- read.table(abundance_file,head=F,stringsAsFactors = F, row.names=1)

  #Process files accordingly
  merged_pca_calls <- merge_pca_file(pca_file,condition=condition)
  expression_file <- simplify_expression_file(expression_file)


  orf_pairs <- merged_pca_calls[,c('ORF.1','ORF.2')]
  output_df <- data.frame()
  for(i in 1:nrow(orf_pairs)){
    pair <- orf_pairs[i,]
    pair <- as.vector(as.matrix(pair))
    #print(pair)
    abundances <- get_orf_pair_abundance(pair,abundance_file)
    mRNA_changes <- get_orf_mrna_changes(pair,expression_file,expression_control_regexp,expression_condition_regexp)
    prediction <- log2(mass_action_predictor(abundances[1],abundances[2],mRNA_changes[1],mRNA_changes[2]))

    output <- data.frame(ORF1=pair[1],
                         ORF2=pair[2],
                         AUC=merged_pca_calls[i,'AUC'],
                         bcPCA_FC.UP=merged_pca_calls[i,'FC.UP'],
                         bcPCA_FC.DN=merged_pca_calls[i,'FC.DN'],
                         bcPCA_FC.AVG=merged_pca_calls[i,'FC.avg'],
                         bcPCA_qVal=merged_pca_calls[i,'q.val'],
                         bcPCA_DMSO.UP=merged_pca_calls[i,'DMSO.UP'],
                         bcPCA_DMSO.DN=merged_pca_calls[i,'DMSO.DN'],
                         PaxDB_Abundance_ORF1=abundances[1],
                         PaxDB_Abundance_ORF2=abundances[2],
                         mRNAFC_ORF1=mRNA_changes[1],
                         mRNAFC_ORF2=mRNA_changes[2],
                         Log2_MA_prediction=prediction)
    output_df <- rbind(output_df,output)
  }
  return(output_df)
}

#' Plots mRNA-based bcPCA predictions compared to quantitative bc-PCA measurements
#'
#' @param my_predictions prediction data frame, output from pca_ma_prediction
#' @param output_path where the plot is to be written
#' @param filename the file name of the plot'
pca_ma_prediction_plot <- function(my_predictions,
                                   output_path,
                                   draw=F,
                                   filename = 'bcPCA_mRNA_predictions.pdf',
                                   prediction_colname = 'Log2_MA_prediction',
                                   measurement_colname = 'bcPCA_FC.AVG',
                                   point_colours = rgb(0.15,0.2,0.3,0.3),
                                   outline_colours = rgb(0,0,0,0),
                                   point_size = 0.8,
                                   label_size = 1.2,
                                   xlimits = c(-5,5),
                                   ylimits = c(-4,2),
                                   relative_text_position_x = 0.2,
                                   relative_text_position_y = 0.8,
                                   x_axis_label = 'mRNA Predicted Log2(R)',
                                   y_axis_label = 'Measured Log2(R)',
                                   long_tick_length = 0.04,
                                   short_tick_length = 0.02
                                   ){
  if(draw==F){
    CairoPDF(file=paste(c(output_path,filename),collapse='/'),width=6,height=6)
  }
  x <- my_predictions[,prediction_colname]
  y <- my_predictions[,measurement_colname]
  plot(x,
       y,
       xlab=x_axis_label,
       ylab=y_axis_label,
       pch=16,
       col=point_colours,
       xlim=xlimits,
       ylim=ylimits,
       axes=F,
       cex=point_size,
       cex.lab=label_size)

  points(x,
         y,
         cex=point_size,
         col=outline_colours)

  axis(side=1,
       tck=-long_tick_length,
       at=seq(min(xlimits),max(xlimits),2),
       labels=T)
  axis(side=1,
       tck=-short_tick_length,
       at=seq(min(xlimits),max(xlimits),1),
       labels=F)


  axis(side=2,
       las=1,
       tck=-long_tick_length,
       at=seq(min(ylimits),max(ylimits),2),
       labels=T)
  axis(side=2,
       las=1,
       tck=-short_tick_length,
       at=seq(min(ylimits),max(ylimits),1),
       labels=F)

  text_x_pos <- quantile(xlimits,probs=c(relative_text_position_x))
  text_y_pos <- quantile(ylimits,probs=c(relative_text_position_y))

  text(paste(c('r = ',format(cor(x,y,use='pair'),digits=2)),collapse=''),
       x=text_x_pos,
       y=text_y_pos,cex=label_size)

  abline(lm(y~x),col='red')

  if(draw==F){
    dev.off()
  }
}


#' Plots precisions for mass-action based complex predictions as a function of cutoff
#'
#' @param my_predictions output table from pca_ma_prediction
#' @param output_path where the plot will be saved if draw is False
#' @param filename the filename of the plot if draw is False
#' @param draw if True, plots the file, if False, saves a CairoPDF object
#' @param prediction_cutoffs the resolution of the plot, defaults to full
#' @param p_cutoff the adjusted p value cutoff for calling significant interactions
#' @param effect_size_cutoff the adjusted effect size cutoff for caling significant interactions
#' @param bottom_predition_limit lowest effect size plotted
#' @param top_predition_limit highest effect size plotted
#' @param line_width
#' @param bootstrap_iters how many bootstraps to draw the 95% polygon
#' @param enhanced_colour
#' @param depleted_colour
#' @param polygon_transparency
#' @param label_size
#' @param legend_x_position
#' @param legend_y_position
#' @param xlabel
#' @param ylabel
#' @param legend_labels
pca_ma_precision_plot <- function(my_predictions,
                                        output_path,
                                        filename,
                                        draw=F,
                                        prediction_cutoffs = 'AUTO',
                                        p_cutoff = 0.05,
                                        effect_size_cutoff = 0.25,
                                        bottom_predition_limit = -2.5,
                                        top_predition_limit = 2.5,
                                        line_width=3,
                                        bootstrap_iters=1000,
                                        enhanced_colour=rgb(0.75,0.1,0.1),
                                        depleted_colour=rgb(0.1,0.1,0.75),
                                        polygon_transparency=0.3,
                                        label_size=1.3,
                                        legend_x_position=-1.2,
                                        legend_y_position=0.85,
                                        xlabel='mRNA Predicted Log2(R) Cutoff',
                                        ylabel='Precision',
                                        legend_labels=c('Accumulated complexes',
                                                        'Depleted complexes')
                                        ){

  percent_correct_predictions <- function(values,
                                          labels,
                                          prediction_cutoffs,
                                          mode){

    retlist <- sapply(prediction_cutoffs,function(cutoff){
      if(mode=='enhanced'){
        crit <- values >= cutoff
      }
      if(mode=='depleted'){
        crit <- values <= cutoff
      }
      sum(labels[crit],na.rm=T)/sum(crit,na.rm=T)
    })
    return(retlist)
  }

  create_bootstrap_polygon <- function(predictions,
                                       values,
                                       prediction_cutoffs,
                                       mode,
                                       boostrap_iters,
                                       quantile_upper=0.95,
                                       quantile_lower=0.05){

    bootstrapped_vals <- sapply(1:bootstrap_iters,function(iter){
      strap <- sample(1:length(predictions),replace=T)
      percent_correct_predictions(predictions[strap],
                                  values[strap],prediction_cutoffs,mode=mode)
    })

    lower_y <- apply(bootstrapped_vals,1,function(x){quantile(x,probs=quantile_lower,na.rm=T)})
    upper_y <- apply(bootstrapped_vals,1,function(x){quantile(x,probs=quantile_upper,na.rm=T)})

    return(cbind(c(prediction_cutoffs,rev(prediction_cutoffs)),
                 c(lower_y,rev(upper_y))))
  }

  mRNA_predictions <- my_predictions$Log2_MA_prediction
  if(prediction_cutoffs == 'AUTO'){
    prediction_cutoffs <- sort(unique(mRNA_predictions))
    prediction_cutoffs <- prediction_cutoffs[prediction_cutoffs >= bottom_predition_limit &
                                               prediction_cutoffs <= top_predition_limit]
  }
  enhanced <-  my_predictions$bcPCA_qVal <= p_cutoff & my_predictions$bcPCA_FC.AVG >= effect_size_cutoff
  depleted <-  my_predictions$bcPCA_qVal <= p_cutoff & my_predictions$bcPCA_FC.AVG <= -(effect_size_cutoff)

  #Plot real values
  depleted_cutoffs <- prediction_cutoffs[prediction_cutoffs <= 0]
  enhanced_cutoffs <- prediction_cutoffs[prediction_cutoffs >= 0]

  depleted_bpc_precision <- percent_correct_predictions(mRNA_predictions,depleted,depleted_cutoffs,mode='depleted')
  enhanced_bpc_precision <- percent_correct_predictions(mRNA_predictions,enhanced,enhanced_cutoffs,mode='enhanced')

  if(draw == F){
    CairoPDF(file=paste(c(output_path,filename),collapse='/'),width=6,height=6)
  }
  plot(depleted_cutoffs,
       depleted_bpc_precision,
       type='l',
       xlab=xlabel,
       ylab=ylabel,
       lwd=line_width,
       col=depleted_colour,
       xlim=c(bottom_predition_limit+0.5,top_predition_limit-0.5),
       cex.lab=label_size)
  lines(enhanced_cutoffs,enhanced_bpc_precision,
        lwd=line_width,
        col=enhanced_colour)

  #Add error bars
  bootstrap_error_depleted <- create_bootstrap_polygon(mRNA_predictions,
                                                        depleted,
                                                        depleted_cutoffs,
                                                        mode='depleted',
                                                        boostrap_iters)

  bootstrap_error_enhanced <- create_bootstrap_polygon(mRNA_predictions,
                                                       enhanced,
                                                       enhanced_cutoffs,
                                                       mode='enhanced',
                                                       boostrap_iters)

  polygon(bootstrap_error_enhanced,border=NA,col=adjustcolor(enhanced_colour, alpha.f = polygon_transparency))
  polygon(bootstrap_error_depleted,border=NA,col=adjustcolor(depleted_colour, alpha.f = polygon_transparency))

  legend(x=legend_x_position,
         y=legend_y_position,
         legend=legend_labels,
         fill=c(enhanced_colour,depleted_colour)
         )

  if(draw == F){
    dev.off()
  }
}



#' Sets colours based on an attribute
#'
#' @param attribute the attribute to be mapped, a vector of values
#' @param color_list a list of colours to be made into a colour ramp
#' @param ncolors number of colours in the ramp
#' @param minimum minimum value in the colour ramp
#' @param maximum maximum value in the colour ramp
#'
#' @return a list of colours, linearly mapped to the colour ramp on the interval between minimum and maximum
set_colours <- function(attribute,color_list,ncolors,minimum,maximum){
  colour_scheme <- grDevices::colorRampPalette(color_list)(ncolors)
  breaks <- seq(minimum,maximum,(maximum-minimum)/(ncolors-1))
  colours <- colour_scheme[sapply(attribute,function(val){which.min(abs(val-breaks))[1]})]
  return(colours)
}


map_gene_names <- function(gene_names){
  sapply(gene_names,function(name){
    proposed_name <- org.Sc.sgd.db::org.Sc.sgdGENENAME[[name]][1]
    if(is.na(proposed_name)){
      return(name)
    }
    return(proposed_name)
  })
}

reverse_map_gene_names <- function(orf_ids){
  sapply(orf_ids,function(orf){
    proposed_orf <- org.Sc.sgd.db::org.Sc.sgdCOMMON2ORF[[orf]][1]
    if(is.null(proposed_orf)){
      return(orf)
    }
    return(proposed_orf)
  })
}

hub_legend_draw <- function(min_val,
                            max_val,
                            color_list,
                            ncolors=1000,
                            main_font_size=2,
                            side_font_size=2,
                            new_plot=T,
                            width=1,
                            height=1,
                            x_adjust=0,
                            y_adjust=0){
  if(new_plot == T){
    plot.new()
  }
  colors <- rev(grDevices::colorRampPalette(color_list)(ncolors))

  xleft <- 0.1 + x_adjust
  middle <- 0.5 + y_adjust
  for(i in 1:(length(colors))){
    #rect(
    #  xleft,
    #  middle-(((i)/ncolors)*height)/2+0.25*height,
    #  xleft+0.6*width,
    #  middle-(((i-1.1)/ncolors)*height)/2+0.25*height,
    #  col=colors[i],
    #  border=NA)
    lines(c(xleft,xleft+0.6*width),
          c(middle-(((i)/ncolors)*height)/2+0.25*height,
            middle-(((i)/ncolors)*height)/2+0.25*height),
          col=colors[i],
          lwd=1/height,
          pch=15)

    #lines(c(0.1+0.02*width,0.1-0.02*width+0.6*width),
    #      c(0.5-((i/ncolors)*height)/2+0.25*height,
    #        0.5-((i/ncolors)*height)/2+0.25*height),
    #      lwd=3,
    #      col=colors[i])
  }

   #rect(xleft,middle-0.25*height,xleft+0.6*width,middle+0.25*height,lwd = 1)
   text(xleft+0.74*width+0.01,middle+0.25*height,max_val,xpd=T,cex=side_font_size*height)
   text(xleft+0.74*width+0.01,middle-0.25*height,min_val,xpd=T,cex=side_font_size*height)
   text(xleft+0.74*width+0.01,middle,mean(c(max_val,min_val)),xpd=T,cex=side_font_size*height)
   text(xleft-0.1*width,middle+0.05+0.25*height,'Log2(R)',xpd=T,cex=main_font_size*height,adj=0)
}

hub_comparison_graph <- function(my_predictions,
                                 hub_name,
                                 hub_name_mode = 'common',
                                 output_path='test',
                                 filename='test.pdf',
                                 draw=F,
                                 color_list = c('red','black','green'),
                                 titles=c("mRNA\nPredictions",'BC-PCA\nMeasurements'),
                                 title_size=3,
                                 title_offset=1.35,
                                 ncolors=1000,
                                 node_expr_color_limits=c(-1,1),
                                 edge_expr_color_limits=c(-1,1),
                                 pca_color_limits=c(-1,1),
                                 q_val_cutoff = 0.05,
                                 effect_size_cutoff = 0.25,
                                 default_node_color=rgb(0.3,0.3,0.3),
                                 edge_width=7,
                                 node_size=20,
                                 graph_seed=1234,
                                 layout_algorithm=igraph::layout.kamada.kawai){
  set.seed(graph_seed)
  if(typeof(hub_name) == 'character'){
    hub_name <- strsplit(hub_name,split=',')[[1]]
  }
  if(hub_name_mode == 'common'){
    hub_name <- reverse_map_gene_names(hub_name)
  }
  hub_predictions <- filter(my_predictions,ORF1 %in% hub_name | ORF2 %in% hub_name,
                            bcPCA_qVal <= q_val_cutoff,
                            abs(bcPCA_FC.AVG) >= effect_size_cutoff)

  draw_comparison_network(hub_predictions,
               color_list,
               ncolors,
               titles,
               title_size,
               title_offset,
               node_expr_color_limits,
               edge_expr_color_limits,
               pca_color_limits,
               default_node_color,
               edge_width,
               node_size,
               graph_seed,
               draw,
               output_path,
               filename,
               layout_algorithm)
}




draw_comparison_network <- function(my_predictions,
                         color_list,
                         ncolors=100,
                         titles,
                         title_size=3,
                         title_offset=1.35,
                         node_expr_color_limits=c(-1,1),
                         edge_expr_color_limits=c(-1,1),
                         pca_color_limits=c(-1,1),
                         default_node_color=rgb(0.3,0.3,0.3),
                         edge_width=7,
                         node_size=20,
                         graph_seed=1234,
                         draw=T,
                         output_path,
                         filename,
                         layout_algorithm=igraph::layout.kamada.kawai){

  create_orf_expr_list <- function(my_predictions){
    orf1_expr <- dplyr::select(my_predictions,ORF1,mRNAFC_ORF1)
    orf2_expr <- dplyr::select(my_predictions,ORF2,mRNAFC_ORF2)
    orf_expr_list <- list()
    for(i in 1:nrow(orf1_expr)){
      orf_expr_list[as.vector(orf1_expr[i,1])] <- orf1_expr[i,2]
      orf_expr_list[as.vector(orf2_expr[i,1])] <- orf2_expr[i,2]
    }
    return(orf_expr_list)
  }
  #Change names and predict expression
  my_predictions$ORF1 <- map_gene_names(as.vector(my_predictions$ORF1))
  my_predictions$ORF2 <- map_gene_names(as.vector(my_predictions$ORF2))
  orf_expr_list <- create_orf_expr_list(my_predictions)


  #Initialize graph
  orf_graph <- igraph::graph_from_edgelist(as.matrix(my_predictions[,c('ORF1','ORF2')]),directed=F)
  igraph::V(orf_graph)$color <- default_node_color
  igraph::V(orf_graph)$label.cex <- ((node_size/17)*4)/sapply(igraph::V(orf_graph)$name,nchar)
  igraph::V(orf_graph)$label.family="Arial Black"
  igraph::V(orf_graph)$label.font=2
  #Colour by node expression
  vertex_expressions <-
    log2(unlist(orf_expr_list[as.vector(igraph::V(orf_graph))]))
  igraph::V(orf_graph)$log2_node_expr <- vertex_expressions
  igraph::V(orf_graph)$expr_colour <-
    set_colours(
      vertex_expressions,color_list,ncolors,node_expr_color_limits[1],node_expr_color_limits[2]
    )

  l <- layout_algorithm(orf_graph)

  igraph::E(orf_graph)$log2_ma_prediction <-
    my_predictions$Log2_MA_prediction
  igraph::E(orf_graph)$predicted_color <-
    set_colours(
      my_predictions$Log2_MA_prediction,color_list,ncolors,edge_expr_color_limits[1],edge_expr_color_limits[2]
    )
  igraph::E(orf_graph)$observed_color <-
    set_colours(
      my_predictions$bcPCA_FC.AVG,color_list,ncolors,pca_color_limits[1],pca_color_limits[2]
    )



  label_by_colour <- function(node_colours,brightness_threshold=80){
    sapply(node_colours,function(colour){
      colour <- col2rgb(colour)
      brightness <- sum(colour*c(0.2126,0.7152,0.0722))
      if(brightness >= brightness_threshold){
        return('black')
      }
      return('white')
    })
  }

  if(draw == F){
    Cairo::CairoPDF(file=paste(c(output_path,filename),collapse='/'),width=12,height=6)
  }
  par(oma=c(0,0,0,0),mar=c(0,0,6,0))
  layout(t(matrix(c(1,2,3))),widths=c(0.7,4,4))
  hub_legend_draw(node_expr_color_limits[1],
                  node_expr_color_limits[2],
                  color_list,
                  ncolors)
  plot(
    orf_graph,edge.width = edge_width,
    vertex.size = node_size,vertex.color=igraph::V(orf_graph)$expr_colour,
    edge.color = igraph::E(orf_graph)$predicted_color,
    vertex.label.color= label_by_colour(igraph::V(orf_graph)$expr_colour),
    layout = l
  )
  #Not sure how else to add a large title
  text(0,title_offset,titles[1],cex=title_size,xpd=T)

  plot(
    orf_graph,edge.width = edge_width,
    vertex.size = node_size,
    edge.color = igraph::E(orf_graph)$observed_color,vertex.color=igraph::V(orf_graph)$expr_colour,
    vertex.label.color= label_by_colour(igraph::V(orf_graph)$expr_colour),
    layout = l
  )
  text(0,title_offset,titles[2],cex=title_size,xpd=T)

  if(draw==F){
    dev.off()
  }
}

#Compares reproducibility between two mRNA measurements
mRNA_comparison <- function(pca_file,
                            expression_file,
                            condition,
                            condition_regexp_vec,
                            control_regexp_vec){
  pca_file <- read.csv(pca_file,head = T,stringsAsFactors = F, sep = '\t')
  expression_file <- read.table(expression_file,head=T,stringsAsFactors = F)
  expression_file <- simplify_expression_file(expression_file)

  orf_pairs <- pca_file%>% dplyr::filter(Condition==condition) %>% dplyr::select(ORF.1,ORF.2)
  orfs <- unique(unlist(orf_pairs))

  expr_matr <- c()
  for(i in 1:length(control_regexp_vec)){
    expr_matr <- cbind(expr_matr,get_orf_mrna_changes(orfs,expression_file,control_regexp_vec[i],condition_regexp_vec[i]))
  }
  rownames(expr_matr) <- orfs
  return(expr_matr)
}

#Compares the reproducibility between two mRNA measurements, either for raw measurements
#Or for derived predictions
mRNA_comparison_graph <- function(pca_file,
                                  expression_file,
                                  condition='ethanol',
                                  condition_regexp_vec=c('Ethanol.4h','Ethanol.12h'),
                                  control_regexp_vec=c('Ethanol.0h','Ethanol.0h'),
                                  xlabel='4h mRNA Expression log10(R)',
                                  ylabel='12h mRNA Expression log10(R)',
                                  point_colours = rgb(0.15,0.2,0.3,0.3),
                                  point_size = 0.8,
                                  label_size = 1.2,
                                  text_position_x = -0.7,
                                  text_position_y = 1.2,
                                  xlimits=c(-1.5,1.5),
                                  ylimits=c(-1.5,1.5),
                                  output_path='dummy',
                                  filename='dummy',
                                  draw=T
                                ){
  expr_matr <- mRNA_comparison(pca_file,expression_file,condition,condition_regexp_vec,control_regexp_vec)
  if(draw == F){
    Cairo::CairoPDF(file=paste(c(output_path,filename),collapse='/'),width=6,height=6)
  }
  par(bg=NA)
  plot(
    log10(expr_matr[,1]),log10(expr_matr[,2]),col = point_colours,pch = 16,xlab =
      xlabel,ylab = ylabel,cex.lab = label_size,cex = point_size, xlim = xlimits, ylim = ylimits
  )
  text(paste(c('r = ',format(cor(log10(expr_matr[,1]),log10(expr_matr[,2]),use='pair'),digits=2)),collapse=''),
       x=text_position_x,
       y=text_position_y,cex=label_size)
  if(draw == F){
    dev.off()
  }
}

#' Performs an error analysis of the mRNA based BCA predictions, based on the reproducibility
#' of UP and DN tag measurements in the PCA experiment and (by default)
#' the reproducibility of the 4h and 12h mRNA based measurements
#' @param pca_file the location of a PCA file
#' @param expression_file the location of an expression file
#' @param abundance_file the location of a file encoding protein abundance
#' @param condition condition tested, so far expression only works on ethanol
#' @param q_val_cutoff q value cutoff for PCA experiment
#' @param fc_cutoff log2(fold change) cutoff for PCA
#' @param expression_condition_regexp_vec regular expression vector that matches the conditions in the expression experiment to be compared
#' @param expression_control_regexp_vec regular expression vector that matches the reference conditions in the expression experiment by which the conditions are to be compared to
#'
#' @return a matrix of noisy mRNA-based PCA measurements
pca_mRNA_error_analysis <- function(pca_file,
                                    expression_file,
                                    abundance_file,
                                    condition='ethanol',
                                    q_val_cutoff=1,
                                    fc_cutoff=0,
                                    up_dn_diff_cutoff=Inf,
                                    expression_condition_regexp_vec=c('Ethanol.4h','Ethanol.1h'),
                                    expression_control_regexp_vec=c('Ethanol.0h','Ethanol.0h'),
                                    iters=100,
                                    my_seed=1234
){
  set.seed(my_seed)
  abundance_file <- read.table(abundance_file,head=F,stringsAsFactors = F, row.names=1)


  #mRNA reproducibility, in log2 space
  expr_matr <- mRNA_comparison(pca_file,expression_file,condition,expression_condition_regexp_vec,expression_control_regexp_vec)


  #PCA reproducibility
  pca_file <- read.csv(pca_file,head = T,stringsAsFactors = F, sep = '\t')
  #
  merged_pca_calls <- merge_pca_file(pca_file,condition=condition)
  merged_pca_calls <- dplyr::filter(merged_pca_calls,q.val <= q_val_cutoff,
                                    abs(FC.avg) >= fc_cutoff,
                                    abs(FC.UP - FC.DN) <= up_dn_diff_cutoff)

  merged_ppi_index <- apply(merged_pca_calls[,c('ORF.1','ORF.2')],1,function(x){paste(x,collapse='_')})
  filtered_pca_file <- dplyr::filter(pca_file,Condition==condition)
  filtered_pca_ppi_index <- apply(filtered_pca_file[,c('ORF.1','ORF.2')],1,function(x){paste(x,collapse='_')})

  filtered_pca_file <- filtered_pca_file[filtered_pca_ppi_index %in% merged_ppi_index,]
  filtered_pca_file_up <- dplyr::filter(filtered_pca_file,UP.DN=='uptag')
  filtered_pca_file_dn <- dplyr::filter(filtered_pca_file,UP.DN=='downtag')

  filtered_upvals <- c(filtered_pca_file_up$FC.rep1,filtered_pca_file_up$FC.rep2)
  filtered_dnvals <- c(filtered_pca_file_dn$FC.rep1,filtered_pca_file_dn$FC.rep2)

  #pca_error model
  up_dn_cor <- cor(filtered_upvals,filtered_dnvals)
  up_dn_error_sd <- sqrt((1/up_dn_cor)-1)
  rep_cor <- cor(filtered_pca_file$FC.rep1,filtered_pca_file$FC.rep2)
  rep_error_sd <- sqrt((1/rep_cor)-1)

  #Process mRNA expression file
  used_orfs <- unique(c(as.vector(merged_pca_calls$ORF.1),as.vector(merged_pca_calls$ORF.2)))
  expr_matr <- expr_matr[used_orfs,]
  expr_cor <- cor(log2(expr_matr[,1]),log2(expr_matr[,2]),use='pair')
  expr_error_sd <- sqrt((1/expr_cor)-1)

  expression_file <- read.table(expression_file,head=T,stringsAsFactors = F)
  expression_file <- simplify_expression_file(expression_file)



  orf_pairs <- merged_pca_calls %>% dplyr::select(ORF.1,ORF.2)
  master_abundance_vec <- c()
  master_mRNA_change_vec <- c()
  for(i in 1:nrow(orf_pairs)){
    pair <- orf_pairs[i,]
    pair <- as.vector(as.matrix(pair))
    #print(pair)
    abundances <- get_orf_pair_abundance(pair,abundance_file)

    mRNA_changes <-
        get_orf_mrna_changes(pair,
                             expression_file,
                             expression_control_regexp_vec[1],expression_condition_regexp_vec[1])
    master_abundance_vec <- rbind(master_abundance_vec,abundances)
    master_mRNA_change_vec <- rbind(master_mRNA_change_vec,mRNA_changes)
  }

  make_error_prone_predictions <- function(master_abundance_vec,
                               master_mRNA_change_vec,
                               expr_error_sd,
                               up_dn_error_sd,
                               rep_error_sd,
                               n_iters=iters){

  iter_matr <- c()
  nrows <- nrow(master_mRNA_change_vec)
  mrna_sd <- sd(log2(master_mRNA_change_vec),na.rm=T)
   for(i in 1:n_iters){
     predictions <- c()
     #Add mRNA noise
     #mRNA_change_vec <- 2^(log2(master_mRNA_change_vec)+rnorm(length(master_mRNA_change_vec),sd=sd(master_abundance_vec)*expr_error_sd))
     mRNA_change_vec <- 2^(log2(master_mRNA_change_vec) + sapply(1:2,function(x){rnorm(nrows,sd=mrna_sd*expr_error_sd)}))
     #print(mRNA_change_vec)
     #print(nrows)
     #print(mrna_sd)
     #stop()
     for(j in 1:nrow(mRNA_change_vec)){
       prediction <- mass_action_predictor(master_abundance_vec[j,1],master_abundance_vec[j,2],mRNA_change_vec[j,1],mRNA_change_vec[j,2])
       #Add final outcome noise
       prediction <- log2(prediction)
       predictions <- c(predictions,prediction)
     }
     noisy_predictions <- predictions + rnorm(length(predictions),sd=sd(predictions,na.rm=T)*up_dn_error_sd)
     noisy_predictions <- noisy_predictions + rnorm(length(predictions),sd=sd(predictions,na.rm=T)*rep_error_sd)
     iter_matr <- cbind(iter_matr,noisy_predictions)#*mean(abs(predictions),na.rm=T)/pca_mean)
   }
   return(iter_matr)
  }

  return(make_error_prone_predictions(master_abundance_vec,master_mRNA_change_vec,expr_error_sd,up_dn_error_sd,rep_error_sd))
}


pca_mRNA_error_analysis_graph <- function(predictions,
                                    error_analysis,
                                    hist_colour='grey',
                                    border_colour='grey39',
                                    line_colour='red',
                                    seed=32,
                                    line_width=2,
                                    text_adjustment_factor=0,
                                    text_size=1.5,
                                    title='',
                                    x_label='Estimated Variance Explained by Model',
                                    draw=T,
                                    output_path=dummy,
                                    filename=dummy
){
  set.seed(seed)
  error_sim_cor_matr <- cor(error_analysis,use='pair')
  error_sim_cors <- error_sim_cor_matr[upper.tri(error_sim_cor_matr)]
  iters <- length(error_sim_cors)
  cor_matr <- sapply(1:iters,function(iter){
    sample_vec <- sample(1:nrow(predictions),replace=T)
    model_cor <- cor(predictions$Log2_MA_prediction[sample_vec],
        predictions$bcPCA_FC.AVG[sample_vec],use='pair')
    return(1/(1/(model_cor^2) - 1/(error_sim_cors^2) + 1))
  })


  beginning <- quantile(cor_matr,probs=c(0.05))
  middle <- median(cor_matr)
  end <- quantile(cor_matr,probs=c(0.95))

  if(draw == F){
    Cairo::CairoPDF(file=paste(c(output_path,filename),collapse='/'),width=6,height=6)
  }
  my_hist <-
    hist(cor_matr,
         breaks = sqrt(length(cor_matr)),
         col = hist_colour,
         border = border_colour,
         main=title,
         xlab=x_label,
         cex.lab=text_size,
         xlim=c(0,1)
         #xlim=c(min(cor_matr),min(max(cor_matr),1)),
         )
  lines(c(beginning,beginning),c(0,max(my_hist$counts)),col=line_colour,lwd=line_width/2,lty=2,xpd=F)
  lines(c(middle,middle),c(0,max(my_hist$counts)),col=line_colour,lwd=line_width,xpd=F)
  lines(c(end,end),c(0,max(my_hist$counts)),col=line_colour,lwd=line_width/2,lty=2,xpd=F)

  textval <- paste(c('Median variance explained = ',format(middle*100,digits=2),'%\n'),collapse='')
  text(middle-abs(middle-min(my_hist$breaks))*text_adjustment_factor,max(my_hist$counts),textval,cex=text_size,xpd=T)
  if(draw == F){
    dev.off()
  }
}


monochromatic_prediction_accuracy_graph <- function(my_predictions,
                                                    hub_df,
                                                    condition='ethanol',
                                                    q_value_cutoff=0.05,
                                                    effect_size_cutoff=0.25,
                                                    sig_cap_width=0.04,
                                                    height_buffer=0.015,
                                                    line_width=2,
                                                    output_path='dummy',
                                                    filename='dummy',
                                                    pdf_width=5,
                                                    pdf_height=9,
                                                    draw=F,
                                                    text_size=1.5){


  hub_df <- dplyr::filter(hub_df,Condition==condition)
  nonsig_hubs <- as.vector(as.matrix(hub_df %>% dplyr::filter(q.value.BH. > q_value_cutoff) %>% dplyr::select(Hub)))
  sig_hubs <- as.vector(as.matrix(hub_df %>% dplyr::filter(q.value.BH. <= q_value_cutoff) %>% dplyr::select(Hub)))

  nonsig_hubs <- reverse_map_gene_names(nonsig_hubs)
  sig_hubs <- reverse_map_gene_names(sig_hubs)

  sig_hub_predictions <- dplyr::filter(my_predictions,ORF1 %in% sig_hubs | ORF2 %in% sig_hubs, bcPCA_qVal <= q_value_cutoff, abs(bcPCA_FC.AVG) >= effect_size_cutoff, !is.na(Log2_MA_prediction))
  nonsig_hub_predictions <- dplyr::filter(my_predictions,ORF1 %in% nonsig_hubs | ORF2 %in% nonsig_hubs, bcPCA_qVal <= q_value_cutoff, abs(bcPCA_FC.AVG) >= effect_size_cutoff, !is.na(Log2_MA_prediction))

  sig_successes <- sum(sign(sig_hub_predictions$Log2_MA_prediction) == sign(sig_hub_predictions$bcPCA_FC.AVG))
  sig_failures <- sum(sign(sig_hub_predictions$Log2_MA_prediction) != sign(sig_hub_predictions$bcPCA_FC.AVG))
  nonsig_successes <- sum(sign(nonsig_hub_predictions$Log2_MA_prediction) == sign(nonsig_hub_predictions$bcPCA_FC.AVG))
  nonsig_failures <- sum(sign(nonsig_hub_predictions$Log2_MA_prediction) != sign(nonsig_hub_predictions$bcPCA_FC.AVG))

  #print(sig_successes)
  #print(sig_failures)

  test_statistic <- fisher.test(rbind(c(sig_successes,sig_failures),c(nonsig_successes,nonsig_failures)))

  sig_conf_interval <- binom.test(sig_successes,sig_successes+sig_failures)$conf.int
  nonsig_conf_interval <- binom.test(nonsig_successes,nonsig_successes+nonsig_failures)$conf.int

  top <- max(c(sig_conf_interval,nonsig_conf_interval))
  top <- top + top*height_buffer

  if(draw == F){
    Cairo::CairoPDF(file=paste(c(output_path,filename),collapse='/'),width=pdf_width,height=pdf_height)
  }
  par(lwd = line_width, cex=text_size)
  barCenters <- barplot(height=c(sig_successes/(sig_successes+sig_failures),nonsig_successes/(nonsig_successes+nonsig_failures)),
          ylim=c(0.5,1),
          xpd=F,
          space=0.5,
          col = RColorBrewer::brewer.pal(12,'Set3')[c(3,5)],
          ylab='mRNA Prediction Accuracy')

  #Custom make labels
  text(x=barCenters[1],y=0.46,c('Concerted\nhubs'),xpd=T)
  text(x=barCenters[2],y=0.46,c('Other\nhubs'),xpd=T)

  abline(h=0.5,lwd=line_width)
  #Add error bars
  segments(barCenters,c(sig_conf_interval[1],nonsig_conf_interval[1]),barCenters,c(sig_conf_interval[2],nonsig_conf_interval[2]),lwd=line_width)
  #Cap error bars
  #for(i in 1:length(barCenters)){
  for(func in c(max,min)){
    lines(c(barCenters[1]+sig_cap_width,barCenters[1]-sig_cap_width),c(rep(func(sig_conf_interval),2)),lwd=line_width)
  }
  for(func in c(max,min)){
    lines(c(barCenters[2]+sig_cap_width,barCenters[2]-sig_cap_width),c(rep(func(nonsig_conf_interval),2)),lwd=line_width)
  }

#
  #Add comparison line
  lines(barCenters,c(top,top),lwd=line_width/2)
  lines(x=rep(barCenters[1],2),y=c(top-height_buffer/2,top),lwd=line_width/2)
  lines(x=rep(barCenters[2],2),y=c(max(nonsig_conf_interval)+height_buffer/2,top),lwd=line_width/2)

  text(mean(barCenters),top+top*height_buffer,paste(c('p = ',format(test_statistic$p.value,digits=2,scientific=T)),collapse=''),xpd=T,cex=1/1.3)
  if(draw == F){
    dev.off()
  }
}



#hub_enrichment_file = '/Users/Albi/Dropbox/barcoded-PCA/2015-09-02/Data for Figure 3D.xlsx'
#hub_df <- xlsx::read.xlsx2(hub_enrichment_file,sheetName = "Sheet1", colClasses=c("character","character",rep("numeric",5)))

#my_predictions <- pca_ma_prediction(pca_universe,protein_abundance_file,expression_file,'ethanol',expression_condition_regexp='Ethanol.4h')
#my_predictions <- filter(my_predictions,abs(bcPCA_FC.UP - bcPCA_FC.DN) <= 1)

#plot_prediction_precision(my_predictions)
#plot(percent_correct_predictions(my_predictions,diminished_bpcs,,'less_than'),type='l')
#lines(percent_correct_predictions(my_predictions,enhanced_bpcs,i,'greater_than')type='l')#
