source('genome_analysis_utilities.R')

prenylation_motif_regexp <- '(g|a|v|l|i){2,2}(m|s|q|a|c|l|e)\\*'
#prenylation_motif_regexp <- 'g'
pca_call_file <- '/Users/Albi/Dropbox/barcoded-PCA/2015-07-29/Additional.file.8.txt'
pca_calls <- read.csv(pca_call_file,sep='\t')

#Orf Universe
all_orfs <- unique(unlist(bpc_notation_to_yeast_orf(as.vector(pca_calls$BPC))))

#Found in condition
filtered_pca_calls <- pca_calls[pca_calls$Condition =='atorvastatin' & 
            pca_calls$Effect.on.growth =='decreased',]

atorvastatin_enhanced_orfs <- unique(as.vector(bpc_notation_to_yeast_orf(as.vector(filtered_pca_calls$BPC))))

all_prenylated_orfs <- get_genes_matching_regex_motif(motif_regexp=prenylation_motif_regexp)
all_prenylated_orfs <- intersect(all_orfs,all_prenylated_orfs )

atorvastatin_prenylated_orfs <- intersect(atorvastatin_enhanced_orfs,all_prenylated_orfs)

stat_table <- rbind(c(length(all_prenylated_orfs),length(all_orfs)),
                  c(length(atorvastatin_prenylated_orfs),length(atorvastatin_enhanced_orfs)))

print(fisher.test(stat_table))