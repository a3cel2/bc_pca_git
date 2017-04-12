library(devtools)
library(Cairo)
library(dplyr)
library(grDevices)
library(xlsx)
library(rmarkdown)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

#Parts of analysis to run
#Options:
#Atorvastatin enrichment
#GO enrichment
#Frequency Perturbation
#Connectivity
to_run <- c()#'Expression PCA')#'Make figures')#'RBD2','Make figures')#c('Frequency Perturbation')#c('Atorvastatin enrichment')

#Markdown directory
figure_path <- "../results/rmarkdown_figures"

#Package containing necessary scripts
devtools::load_all('../packages/bcPcaAnalysis')
devtools::document('../packages/bcPcaAnalysis')

#Global parameters
pca_universe = '../data/external/Additional.file.6.txt'
pca_enhanced_calls = '../data/external/Additional.file.10.txt'
pca_depleted_calls = '../data/external/Additional.file.11.txt'
go_association_file = '../data/funcassociate_go_associations.txt'

protein_abundance_file = '../data/paxdb_abundance.tsv'
protein_abundance_file_msonly = '../data/paxdb_abundance_msms_godoy.txt'
expression_file = '../data/external/Additional.file.14.txt'

hub_enrichment_file = '../data/external/Data for Figure 3D.xlsx'
rbd2_mtx_data_file <- '../data/external/Additional.file.15.txt'
rbd2_no_mtx_data_file <- '../data/external/Additional.file.16.txt'


output_path <- '../results/master_output'
saved_parameter_path <- '../data/script_input'
heatmap_cluster_file <- 'cluster1_heatmap.tsv'
figure_path <- '../results/rmarkdown_figures'

#Color scale for a lot of stuff
my_color_list <- c(
  rgb(1,0.45,0.25),
  rgb(0.8,0.25,0.25),
  rgb(0,0,0),
  rgb(0.25,0.45,0.8),
  rgb(0.25,0.75,1)
)
blue_black_orange <- grDevices::colorRampPalette(my_color_list)



##Maslov simulations, shuffle Kds
expression_pca_output_path <- paste(c(output_path,'expression_pca'),collapse='/')
dir.create(expression_pca_output_path, showWarnings = FALSE)

for(time_point in c('30min','1h','4h','12h')){
  condition_regexp <- paste(c('Ethanol',time_point),collapse='.')
  
  #We'll try everything with shuffled Kds
  shuffle_status <- c(T,F)
  for(shuffle in shuffle_status){
    
    if(shuffle){
      filename <- paste(c('bcPCA_mRNA_predictions_maslov_shuffled',time_point,'.pdf'),collapse='')
    }else{
      filename <- paste(c('bcPCA_mRNA_predictions_maslov_',time_point,'.pdf'),collapse='')
    }
    predictions <- pca_ma_prediction(pca_universe,
                                     protein_abundance_file,
                                     expression_file,
                                     'ethanol',
                                     expression_condition_regexp=condition_regexp,
                                     prediction_method = "maslov",
                                     shuffled_kds = shuffle)
    pca_ma_prediction_plot(predictions,expression_pca_output_path,filename=filename)
    
    #Also do independend shuffled ones
    if(shuffle){
      filename <- paste(c('bcPCA_mRNA_predictions_shuffled',time_point,'.pdf'),collapse='')
      predictions <- pca_ma_prediction(pca_universe,
                                       protein_abundance_file,
                                       expression_file,
                                       'ethanol',
                                       expression_condition_regexp=condition_regexp,
                                       shuffled_kds = shuffle)
      pca_ma_prediction_plot(predictions,expression_pca_output_path,filename=filename)
    }
  }
}


#Baseline growth plot
expression_pca_baseline_output_path <- paste(c(output_path,'expression_pca_baseline'),collapse='/')
dir.create(expression_pca_baseline_output_path, showWarnings = FALSE)
#This doesn't matter, it's just so I can use the same function as above

shuffle_status <- c(T,F)
prediction_methods <- c('maslov','independent')
datasets <- list('whole'=protein_abundance_file,
                 'ms'=protein_abundance_file_msonly)
for(prediction_method in prediction_methods){
  for(dataset_name in names(datasets)){
    for(shuffle in shuffle_status){
    predictions <- pca_ma_prediction(pca_universe,
                                     abundance_file = datasets[[dataset_name]],
                                     expression_file,
                                     condition = 'ethanol',
                                     shuffled_kds = shuffle,
                                     prediction_method = prediction_method,
                                     return_only_baseline = T)
    if(shuffle){
      filename <- paste(c('bcPCA_baseline_predictions_shuffled',dataset_name,'_',prediction_method,'.pdf'),collapse='')
    }else{
      filename <- paste(c('bcPCA_baseline_predictions_',dataset_name,'_',prediction_method,'.pdf'),collapse='')
    }
    outfile <- paste(c(expression_pca_baseline_output_path,filename),collapse='/')
    CairoPDF(file=outfile,width=6,height=6)
    
    #Populate the title
    main_text_vector <- c('Predictions - ')
    if(prediction_method == 'maslov'){
      main_text_vector <- paste(c(main_text_vector,'Competitive binding model'),collapse='')
    }else if(prediction_method =='independent'){
      main_text_vector <- paste(c(main_text_vector,'Independent binding model'),collapse='')
    }
    if(shuffle){
      main_text_vector <- c(main_text_vector,'Shuffled estimated Kds')
    }else{
      main_text_vector <- c(main_text_vector,'Estimated Kds')
    }
    
    if(dataset_name == 'ms'){
      main_text_vector <- c(main_text_vector,'MS data')
    } else if(dataset_name == 'whole'){
      main_text_vector <- c(main_text_vector,'All data')
    }
    
    plot(predictions$Log2_MA_prediction,
         predictions$AUC,
         main=main_text_vector,
         pch=16,
         col=rgb(0,0,0,0.3),
         cex=0.7,
         xlab = 'Predicted abundance [log10 ppm]',
         ylab = 'Isogenic culture [log10 AUC]')
    
    min_pred <- min(predictions$Log2_MA_prediction,na.rm=T)
    max_pred <- max(predictions$Log2_MA_prediction,na.rm=T)
    
    min_auc <- min(predictions$AUC,na.rm=T)
    max_auc <- max(predictions$AUC,na.rm=T)
    text(
      x = min_pred + (max_pred - min_pred)/10,
      y = max_auc - (max_auc - min_auc)/10,
      labels=paste(c('r=',format(cor(predictions$Log2_MA_prediction,predictions$AUC,use='pair'),digits=2)),collapse=''),
      cex=2
    )
    abline(lm(predictions$AUC~predictions$Log2_MA_prediction),col='red',lwd=2,lty=5)
    dev.off()
    }
  }
}

#Spuriously significant results
pca_univ <- read.table(pca_universe,head=T,sep='\t')
spurious_sig_output_path <- paste(c(output_path,'spurious_significant'),collapse='/')
dir.create(spurious_sig_output_path, showWarnings = FALSE)
nsigs <- sapply(1:100,function(i){
  test <- randomize_pca_file(pca_univ)
  nsigs <- get_n_sig_genes(test)
  write(paste(c(nsigs,'spurious'),collapse=' '), stderr())
  return(nsigs)
})
outfile <- paste(c(spurious_sig_output_path,'spurious_hist.pdf'),collapse='/')
CairoPDF(file=outfile,width=6,height=6)
hist(nsigs,breaks=100,col='black',xlab='Spuriously Significant Results',main='')
dev.off()