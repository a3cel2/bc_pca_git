this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

source('genome_analysis_utilities.R')


##Check for enrichment of atorvastatin motif
print('Checking for CAAX motif amongst atorvastatin decreased drugs')
print(given_motif_enrichment(
  universe_file = '/Users/Albi/Dropbox/barcoded-PCA/2015-07-29/Additional.file.3.txt',
  pca_call_file = '/Users/Albi/Dropbox/barcoded-PCA/2015-07-29/Additional.file.8.txt',
  motif_regexp = 'c(g|a|v|l|i){2,2}(m|s|q|a|c|l|e)\\*',
  drug = 'atorvastatin',
  direction = 'decreased'
))