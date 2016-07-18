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

to_run <- c()#('Connectivity')#c('Atorvastatin enrichment')

#Markdown directory


#Package containing necessary scripts
devtools::load_all('../packages/bcPcaAnalysis')
devtools::document('../packages/bcPcaAnalysis')

#Global parameters
pca_universe = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.6.txt'
pca_enhanced_calls = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.10.txt'
pca_depleted_calls = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.11.txt'
go_association_file = '/Users/Albi/Dropbox/Roth Lab/projects/bc_pca_git/data/funcassociate_go_associations.txt'

protein_abundance_file = "/Users/Albi/Dropbox/Roth Lab/projects/bc_pca_git/data/paxdb_abundance.tsv"
expression_file = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.14.txt'

hub_enrichment_file = '/Users/Albi/Dropbox/barcoded-PCA/2015-09-02/Data for Figure 3D.xlsx'

output_path <- '../results/master_output'
saved_parameter_path <- '../data/script_input'

#Color scale for a lot of stuff
my_color_list <- c(
  rgb(1,0.45,0.25),
  rgb(0.8,0.25,0.25),
  rgb(0,0,0),
  rgb(0.25,0.45,0.8),
  rgb(0.25,0.75,1)
)
blue_black_orange <- grDevices::colorRampPalette(my_color_list)

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

if('Heatmap' %in% to_run){
  heatmap_output_path <- paste(c(output_path,'heatmaps'),collapse='/')
  dir.create(heatmap_output_path, showWarnings = FALSE)
  
  pca_file <- read.table(pca_universe,head=T,sep='\t')
  pca_enh <- read.table(pca_enhanced_calls,head=T,sep='\t')
  pca_depl <- read.table(pca_depleted_calls,head=T,sep='\t')
  
  convert_pca_file_to_heatmap_format(pca_file=pca_file,
                                     pca_enhanced_calls = pca_enh,
                                     pca_depleted_calls = pca_depl,
                                     output_path=heatmap_output_path,
                                     filename='all_sig_bpcs.png',
                                     draw=F,
                                     label_size=5,
                                     line_width=3
                                     )
  
  condition_summary_barplot(pca_enh,
                            pca_depl,
                            color_function=blue_black_orange,
                            draw=F,
                            output_path=heatmap_output_path,
                            filename='number_dynamic_bpcs.pdf')
}


if('Connectivity' %in% to_run){
  connectivity_output_path <- paste(c(output_path,'connectivity'),collapse='/')
  dir.create(connectivity_output_path, showWarnings = FALSE)
  
  pca_universe <- read.table(pca_universe,head=T,sep='\t')
  pca_enhanced <- read.table(pca_enhanced_calls,head=T,sep='\t')
  pca_depleted <- read.table(pca_depleted_calls,head=T,sep='\t')
  hub_df <- xlsx::read.xlsx2(hub_enrichment_file,sheetName = "Sheet1", colClasses=c("character","character",rep("numeric",5)))
  
  stop()
  Cairo::CairoPDF(file=paste(c(connectivity_output_path,'doxorubicin_connnectivity.pdf'),collapse='/'),width=5,height=6)
  network_connectivity_graph(pca_universe,
                             pca_enhanced,
                             pca_depleted,
                             'doxorubicin',
                             my_color_list,
                             #Set this to True to edit the network layout
                             edit=F,
                             node_size=10,
                             edge_width=4,
                             layout_save_dir='../data/script_input/',
                             layout_file='doxo_connectivity.layout')
  dev.off()
  
  Cairo::CairoPDF(file=paste(c(connectivity_output_path,'ethanol_connnectivity.pdf'),collapse='/'),width=6,height=6)
  network_connectivity_graph(pca_universe,
                             pca_enhanced,
                             pca_depleted,
                             'ethanol',
                             my_color_list,
                             #Set this to True to edit the network layout
                             edit=F,
                             load_saved=T,
                             layout_save_dir='../data/script_input/',
                             layout_file='ethanol_connectivity.layout',
                             my_title='Ethanol directional connectivity',
                             node_size=5,
                             edge_width=2)
  dev.off()
  
  connectivity_iterations <- 1000
  
  Cairo::CairoPDF(file=paste(c(connectivity_output_path,'component_size_significance.pdf'),collapse='/'),width=12,height=3)
  component_sig_matrix <- network_simulation_significance(pca_universe,
                                               pca_enhanced,
                                               pca_depleted,
                                               iterations=connectivity_iterations)
  connectivity_graph(component_sig_matrix,my_color_list)
  dev.off()
  
  Cairo::CairoPDF(file=paste(c(connectivity_output_path,'component_size_significance_nodewise.pdf'),collapse='/'),width=12,height=3)
  component_sig_matrix_nodewise <- network_simulation_significance(pca_universe,
                                                               pca_enhanced,
                                                               pca_depleted,
                                                               iterations=connectivity_iterations,
                                                               mode = 'nodewise')
  connectivity_graph(component_sig_matrix_nodewise,my_color_list)
  dev.off()
  
  Cairo::CairoPDF(file=paste(c(connectivity_output_path,'density_significance.pdf'),collapse='/'),width=12,height=3)
  density_sig_matrix <- network_simulation_significance(pca_universe,
                                                             pca_enhanced,
                                                             pca_depleted,
                                                             iterations=connectivity_iterations,
                                                             metric='density')
  connectivity_graph(density_sig_matrix,my_color_list,legend_labels=c('-1'="Decreased density",
                                                                           '0'="Expected density",
                                                                           '1'="Increased density"))
  dev.off()
  
  
  Cairo::CairoPDF(file=paste(c(connectivity_output_path,'density_significance_nodewise.pdf'),collapse='/'),width=12,height=3)
  connectivity_sig_matrix <- network_simulation_significance(pca_universe,
                                                               pca_enhanced,
                                                               pca_depleted,
                                                               mode = 'nodewise',
                                                               iterations=connectivity_iterations,
                                                               metric='density')
  connectivity_graph(connectivity_sig_matrix,my_color_list,legend_labels=c('-1'="Decreased density",
                                                                                         '0'="Expected density",
                                                                                         '1'="Increased density"))
  dev.off()
  
  Cairo::CairoPDF(file=paste(c(connectivity_output_path,'component_size_significance_search.pdf'),collapse='/'),width=9,height=6)
  component_size_sig_search_matrix <- network_simulation_significance_node_edge_search_matrix(pca_universe = pca_universe,
                                                                                            pca_enhanced = pca_enhanced,
                                                                                            pca_depleted = pca_depleted,
                                                                                            node_probs=c(0:10)/10,
                                                                                            iterations=connectivity_iterations,
                                                                                            load_saved=T,
                                                                                            save_output=T,
                                                                                            save_directory=saved_parameter_path,
                                                                                            save_filename='parameter_search_component_size.tsv')
  node_edge_search_heatmap(component_size_sig_search_matrix,blue_black_orange)
  dev.off()
  
  Cairo::CairoPDF(file=paste(c(connectivity_output_path,'density_significance_search.pdf'),collapse='/'),width=9,height=6)
  density_sig_search_matrix <- network_simulation_significance_node_edge_search_matrix(pca_universe = pca_universe,
                                                                                              pca_enhanced = pca_enhanced,
                                                                                              pca_depleted = pca_depleted,
                                                                                              node_probs=c(0:10)/10,
                                                                                              iterations=connectivity_iterations,
                                                                                              load_saved=T,
                                                                                              save_output=T,
                                                                                              metric='density',
                                                                                              save_directory=saved_parameter_path,
                                                                                              save_filename='parameter_search_density.tsv')
  node_edge_search_heatmap(density_sig_search_matrix,
                           blue_black_orange,
                           legend_labels=c('-1'="Decreased\ndensity",
                                           '0'="Expected\ndensity",
                                           '1'="Increased\ndensity"))
  dev.off()
  
  
  Cairo::CairoPDF(file=paste(c(connectivity_output_path,'hub_bias_heatmap.pdf'),collapse='/'),width=10,height=5)
  hub_bias_heatmap(hub_df,blue_black_orange,nonsig_colour = 'grey90')
  dev.off()
  
  Cairo::CairoPDF(file=paste(c(connectivity_output_path,'connnectivity_hist_doxo.pdf'),collapse='/'),width=10,height=5)
  connectivity_histogram(pca_universe,
                         pca_enhanced,
                         pca_depleted,
                         'doxorubicin')
  dev.off()
  
  
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
  #Compare mRNA expression reproducibility
  mRNA_comparison_graph(pca_universe,expression_file,output_path=expression_pca_output_path,filename='mRNA_4h_12h_comparison',draw=F)
  #Compare tag reproducibility
  pca_ma_prediction_plot(my_predictions,
                         prediction_colname='bcPCA_FC.UP',
                         measurement_colname='bcPCA_FC.DN',
                         x_axis_label = "bcPCA Log2(R) UPTAG",
                         y_axis_label = "bcPCA Log2(R) DNTAG",
                         output_path=expression_pca_output_path,
                         filename='UPtag_DNtag_comparison',
                         xlimits=c(-4,2),
                         ylimits=c(-4,2))
  
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
  
  #Compare error_prone replicates
  error_analysis <- pca_mRNA_error_analysis(pca_universe,expression_file,protein_abundance_file,iters=1000)
  pca_mRNA_error_analysis_graph(error_analysis,output_path=expression_pca_output_path,filename='replicate_comparison',draw=F)
  
  ##For categorical
  pca_ma_precision_plot(my_predictions,expression_pca_output_path,filename='bcPCA_precision')
  
  #Compare different hubs
  hub_comparison_graph(my_predictions,'HXT1,HSP30',output_path=expression_pca_output_path,filename='HXT1_HSP30_comparison',color_list=my_color_list)
  hub_comparison_graph(my_predictions,'LSP1',
                       node_size=35,
                       edge_width=10,
                       output_path=expression_pca_output_path,
                       filename='LSP1_comparison',
                       title_offset=1.41,
                       color_list=my_color_list,
                       layout_algorithm=igraph::layout.fruchterman.reingold)
  
  #Compare mass action accuracy for monochromatic hubs
  hub_df <- xlsx::read.xlsx2(hub_enrichment_file,sheetName = "Sheet1", colClasses=c("character","character",rep("numeric",5)))
  monochromatic_prediction_accuracy_graph(my_predictions,
                                          hub_df,
                                          output_path=expression_pca_output_path,
                                          filename='hub_comparison')
  
  
  
  
}

#Make Figures
rmarkdown::render("manuscript_latex/figures/figures.Rmd", "pdf_document")
