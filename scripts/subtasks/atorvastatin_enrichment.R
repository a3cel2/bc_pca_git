library(seqinr)

print('Checking for CAAX motif enrichment amongst atorvastatin decreased drugs')
print(given_motif_enrichment(
  universe_file = pca_universe,
  #Enhanced interactions
  pca_call_file = pca_enhanced_calls,
  motif_regexp = 'c(g|a|v|l|i){2,2}(m|s|q|a|c|l|e)\\*',
  condition = 'atorvastatin'
))