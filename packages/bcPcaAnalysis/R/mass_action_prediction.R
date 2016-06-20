devtools::use_package('dplyr')
devtools::use_package('metap')
devtools::use_package('grDevices')
devtools::use_package('Cairo')
#require(dplyr)


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
  merged_pca_calls_up <- dplyr::filter(pca_calls, Condition == condition,UP.DN=='uptag')
  merged_pca_calls_dn <- dplyr::filter(pca_calls, Condition == condition,UP.DN=='downtag')
  
  if(!identical(rownames(merged_pca_calls_up),rownames(merged_pca_calls_dn))){
    stop('PCA file is not evenly divided into uptag and uptag, please provide unprocessed data')
  }
  
  new_fc_avg <- sapply(1:nrow(merged_pca_calls_up),function(i){
    mean(c(merged_pca_calls_up$FC.avg[i],merged_pca_calls_dn$FC.avg[i]))
  })
  
  new_q_val <- sapply(1:nrow(merged_pca_calls_up),function(i){
   #Stouffer.log(c(merged_pca_calls_up$q.val[i],merged_pca_calls_dn$q.val[i]))
    unlist(metap::sumz(c(merged_pca_calls_up$q.val[i],merged_pca_calls_dn$q.val[i])))[2]
    
  })
  
  
  merged_pca_calls <- data.frame(ORF.1=merged_pca_calls_up$ORF.1,
                                   ORF.2=merged_pca_calls_up$ORF.2,
                                   FC.UP=merged_pca_calls_up$FC.avg,
                                   FC.DN=merged_pca_calls_dn$FC.avg,
                                   FC.avg=new_fc_avg,
                                   q.val=new_q_val)
  
#    aggregate_merge_p1 <- aggregate(.~ORF.1+ORF.2,
#                                    data=merged_pca_calls[,c('ORF.1',
#                                                             'ORF.2',
#                                                             'FC.UP',
#                                                             'FC.DN',
#                                                             'FC.avg')]
#                                    ,mean)
#    
#    aggregate_merge_p2 <- aggregate(.~ORF.1+ORF.2,
#                                    data=merged_pca_calls[,c('ORF.1',
#                                                             'ORF.2',
#                                                             'q.val')]
#                                    ,function(x){
#                                      if(length(x)==1){
#                                        return(x)
#                                      }
#                                      unlist(metap::sumz(x))[2]
#                                      })
#    
#    merged_pca_calls <- cbind(aggregate_merge_p1,
#                              aggregate_merge_p2)
  
  return(merged_pca_calls)
}

#Given an mRNA expression file, obtains expression changes for ORFs in a pair
#Paired index data parameter asks if control indexes correspond to condition indexes
#e.g. first column of control is same experiment as first column of experiment, etc
get_orf_pair_mRNA_changes <- function(pair,
                                      expression_file,
                                      expression_control_regexp,
                                      expression_condition_regexp,
                                      paired_index_data=F){
  ctrl_index <- grep(expression_control_regexp,
                     colnames(expression_file))
  cond_index <- grep(expression_condition_regexp,
                     colnames(expression_file))
  
  pair_exp <- c()
  for(orf in pair){
    orf_expression <- dplyr::filter(expression_file,ORF==orf)
    ctrl_orf_expression <- as.matrix(orf_expression[,ctrl_index,drop=F])
    cond_orf_expression <- as.matrix(orf_expression[,cond_index,drop=F])
    fc <- c()
    if(nrow(orf_expression) >= 1){
      for(i in 1:nrow(orf_expression)){
        if(paired_index_data){
          for(j in 1:length(ctrl_index)){
            measurement_fc <- cond_orf_expression[i,j]/ctrl_orf_expression[i,j]
            fc <- c(fc, measurement_fc)  
          }
        }
        else{
          measurement_fc <- mean(cond_orf_expression[i,])/mean(ctrl_orf_expression[i,])  
          fc <- c(fc, measurement_fc)
        }
      }
      pair_exp <- c(pair_exp,mean(fc))
    }
    else{
      pair_exp <- c(pair_exp,NA)
    }
  }
  return(pair_exp)
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
  
  merged_pca_calls <- merge_pca_file(pca_file,condition=condition)
  
  
  orf_pairs <- merged_pca_calls[,c('ORF.1','ORF.2')]
  output_df <- data.frame()
  for(i in 1:nrow(orf_pairs)){
    pair <- orf_pairs[i,]
    pair <- as.vector(as.matrix(pair))
    #print(pair)
    abundances <- get_orf_pair_abundance(pair,abundance_file)
    mRNA_changes <- get_orf_pair_mRNA_changes(pair,expression_file,expression_control_regexp,expression_condition_regexp)
    prediction <- log2(mass_action_predictor(abundances[1],abundances[2],mRNA_changes[1],mRNA_changes[2]))
    
    output <- data.frame(ORF1=pair[1],
                         ORF2=pair[2],
                         bcPCA_FC.UP=merged_pca_calls[i,'FC.UP'],
                         bcPCA_FC.DN=merged_pca_calls[i,'FC.DN'],
                         bcPCA_FC.AVG=merged_pca_calls[i,'FC.avg'],
                         bcPCA_qVal=merged_pca_calls[i,'q.val'],
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
                                   point_colours = rgb(0.08,0.08,0.35,0.3),
                                   outline_colours = rgb(0,0,0,0.1),
                                   point_size = 0.5,
                                   label_size = 1.3,
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
  if(draw==F){
    dev.off()
  }
}


pca_ma_precision_plot <- function(my_predictions,
                                        output_path,
                                        filename,
                                        prediction_cutoffs = 'AUTO',
                                        p_cutoff = 0.05,
                                        effect_size_cutoff = 0.25,
                                        bottom_predition_limit = -2.5,
                                        top_predition_limit = 2.5,
                                        line_width=3,
                                        bootstrap_iters=1000,
                                        enhanced_colour=rgb(1,0,0),
                                        depleted_colour=rgb(0,0,1),
                                        polygon_transparency=0.3,
                                        label_size=1.3,
                                        xlabel='mRNA Predicted Log2(R) Cutoff',
                                        ylabel='Precision'
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
  
  depleted_bpc_precision <- percent_correct_predictions(mRNA_predictions,depleted,prediction_cutoffs,mode='depleted')
  enhanced_bpc_precision <- percent_correct_predictions(mRNA_predictions,enhanced,prediction_cutoffs,mode='enhanced')
  
  CairoPDF(file=paste(c(output_path,filename),collapse='/'),width=6,height=6)
  plot(prediction_cutoffs,
       depleted_bpc_precision,
       type='l',
       xlab=xlabel,
       ylab=ylabel,
       lwd=line_width,
       col=depleted_colour,
       xlim=c(bottom_predition_limit+0.5,top_predition_limit-0.5),
       cex.lab=label_size)
  lines(prediction_cutoffs,enhanced_bpc_precision,
        lwd=line_width,
        col=enhanced_colour)
  
  #Add error bars
  bootstrap_error_depleted <- create_bootstrap_polygon(mRNA_predictions,
                                                        depleted,
                                                        prediction_cutoffs,
                                                        mode='depleted',
                                                        boostrap_iters)
  
  bootstrap_error_enhanced <- create_bootstrap_polygon(mRNA_predictions,
                                                       enhanced,
                                                       prediction_cutoffs,
                                                       mode='enhanced',
                                                       boostrap_iters)
  
  polygon(bootstrap_error_enhanced,border=NA,col=adjustcolor(enhanced_colour, alpha.f = polygon_transparency))
  polygon(bootstrap_error_depleted,border=NA,col=adjustcolor(depleted_colour, alpha.f = polygon_transparency))
  dev.off()
}

#my_predictions <- pca_ma_prediction(pca_universe,protein_abundance_file,expression_file,'ethanol',expression_condition_regexp='Ethanol.4h')
#my_predictions <- filter(my_predictions,abs(bcPCA_FC.UP - bcPCA_FC.DN) <= 1)

#plot_prediction_precision(my_predictions)
#plot(percent_correct_predictions(my_predictions,diminished_bpcs,,'less_than'),type='l')
#lines(percent_correct_predictions(my_predictions,enhanced_bpcs,i,'greater_than')type='l')#