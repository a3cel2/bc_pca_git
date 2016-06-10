library(devtools)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


to_run <- c()#c('Atorvastatin enrichment')

#Package containing necessary scripts
devtools::load_all('../packages/bcPcaAnalysis')
devtools::document('../packages/bcPcaAnalysis')

#Global parameters
pca_universe = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.6.txt'
pca_enhanced_calls = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.11.txt'
pca_depleted_calls = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.10.txt'
go_assoc
##Check for enrichment of isoprenylation motif under atorvastatin

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

