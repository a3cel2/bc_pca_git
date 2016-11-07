library(knitr)

freq_perturb_output_path <- paste(c(output_path,'freq_peturb'),collapse='/')
dir.create(freq_perturb_output_path, showWarnings = FALSE)

pca_univ <- read.table(pca_universe,head=T,sep='\t')
pca_enh <- read.table(pca_enhanced_calls,head=T,sep='\t')
pca_depl <- read.table(pca_depleted_calls,head=T,sep='\t')

pca_perturbed <- rbind(pca_enh,pca_depl)

frequent_ppis <- get_frequent_ppis(pca_perturbed,4)
frequent_genes <- get_genes_from_ppis(pca_perturbed,frequent_ppis)

all_ppis <- get_frequent_ppis(pca_univ,0)
all_genes <- get_genes_from_ppis(pca_univ,all_ppis)

go_enrichment <- funcassociate(frequent_genes,all_genes,order_mode='unordered')
outfile <- paste(c(freq_perturb_output_path,'go_table_frequent_perturbations.html'),collapse='/')
write(kable(format_funcassociate(go_enrichment),format='html',caption='Go Enrichment amongst proteins participating in frequently dynamic complexes'),outfile)