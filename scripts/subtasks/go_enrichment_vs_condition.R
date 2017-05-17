go_output_path <- paste(c(output_path,'go_enrich'),collapse='/')
dir.create(go_output_path, showWarnings = FALSE)

convert_costanzo_matlab_data('../data/go_maps/go_bp_140819.mat','../data/go_maps/costanzo_go.txt')
convert_sgd_slim_go_data('../data/go_maps/go_slim_sgd_genewise.tab.txt','../data/go_maps/go_slim_sgd_termwise.txt')
#ppi_freq_plot(pca_univ,pca_enh,pca_depl)
all_conditions <- read.csv(pca_universe,sep='\t',stringsAsFactors=F)
all_conditions <- unique(all_conditions$Condition)
#stop()
funcassociate_output <- bcpca_funcassociate_analysis(pca_universe,
                                                     pca_enhanced_calls,
                                                     pca_depleted_calls,
                                                     all_conditions,
                                                     reps=10000,
                                                     network_modes='nodewise')

funcassociate_table <- format_bc_pca_funcassociate(funcassociate_output,conditions = all_conditions)
outfile <- paste(c(go_output_path,'go_table_all_experiments.html'),collapse='/')
write(knitr::kable(funcassociate_table,format='html',caption='Go Enrichment amongst all_experiments'),outfile)