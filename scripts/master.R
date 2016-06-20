library(devtools)
library(Cairo)
library(dplyr)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

#Parts of analysis to run
#Options:
#Atorvastatin enrichment
#GO enrichment

to_run <- c('Expression PCA')#c('Atorvastatin enrichment')

#Package containing necessary scripts
devtools::load_all('../packages/bcPcaAnalysis')
devtools::document('../packages/bcPcaAnalysis')

#Global parameters
pca_universe = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.6.txt'
pca_enhanced_calls = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.11.txt'
pca_depleted_calls = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.10.txt'
go_association_file = '/Users/Albi/Dropbox/Roth Lab/projects/bc_pca_git/data/funcassociate_go_associations.txt'

protein_abundance_file = "/Users/Albi/Dropbox/Roth Lab/projects/bc_pca_git/data/paxdb_abundance.tsv"
expression_file = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.14.txt'

output_path <- '../results/master_output'


##Check for enrichment of isoprenylation motif under atorvastatin
#Fails, did not bother giving output
if('Atorvastatin enrichment' %in% to_run){
  print('Checking for CAAX motif enrichment amongst atorvastatin decreased drugs')
  print(given_motif_enrichment(
    universe_file = pca_universe,
    #Enhanced interactions
    pca_call_file = pca_enhanced_calls,
    motif_regexp = 'c(g|a|v|l|i){2,2}(m|s|q|a|c|l|e)\\*',
    condition = 'atorvastatin'
  ))
}

#Check for GO enrichment in each condition, enhanced and depleted, nodewise and edgewise
if('GO enrichment' %in% to_run){
  all_conditions <- read.csv(pca_universe,sep='\t',stringsAsFactors=F)
  all_conditions <- unique(all_conditions$Condition)
  funcassociate_output <- bcpca_funcassociate_analysis(pca_universe,
                                                       pca_enhanced_calls,
                                                       pca_depleted_calls,
                                                       all_conditions)
}

if('Expression PCA' %in% to_run){
  expression_pca_output_path <- paste(c(output_path,'expression_pca'),collapse='/')
  my_predictions <- pca_ma_prediction(pca_universe,protein_abundance_file,expression_file,'ethanol',expression_condition_regexp='Ethanol.4h')
  
  dir.create(expression_pca_output_path, showWarnings = FALSE)
  
  
  ##For quantitative
  #All predictions
  pca_ma_prediction_plot(my_predictions,expression_pca_output_path)
  
  #Remove overly noisy measurements
  filtered_predictions <- filter(my_predictions,
                                 abs(bcPCA_FC.UP - bcPCA_FC.DN) <= 1)
  pca_ma_prediction_plot(filtered_predictions,
                         expression_pca_output_path,
                         filename='bcPCA_predictions_filtered.pdf')
  
  #Only significant
  filtered_predictions <- filter(my_predictions,
                                 bcPCA_qVal < 0.05,
                                 abs(bcPCA_FC.AVG) > 0.25)
  pca_ma_prediction_plot(filtered_predictions,
                         expression_pca_output_path,
                         filename='bcPCA_predictions_significant.pdf')
  
  #Significant and filtered
  filtered_predictions <- filter(my_predictions,
                                 abs(bcPCA_FC.UP - bcPCA_FC.DN) <= 1,
                                 bcPCA_qVal < 0.05,
                                 abs(bcPCA_FC.AVG) > 0.25)
  pca_ma_prediction_plot(filtered_predictions,
                         expression_pca_output_path,
                         filename='bcPCA_predictions_significant_filtered.pdf')
  
  ##For categorical
  pca_ma_precision_plot(my_predictions,expression_pca_output_path,filename='bcPCA_precision')

}