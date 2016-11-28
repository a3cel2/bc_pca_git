connectivity_output_path <- paste(c(output_path,'connectivity'),collapse='/')
dir.create(connectivity_output_path, showWarnings = FALSE)

#Number of simulations for estimating connectivity
connectivity_iterations <- 1000

pca_universe <- read.table(pca_universe,head=T,sep='\t')
pca_enhanced <- read.table(pca_enhanced_calls,head=T,sep='\t')
pca_depleted <- read.table(pca_depleted_calls,head=T,sep='\t')
hub_df <- xlsx::read.xlsx2(hub_enrichment_file,sheetName = "Sheet1", colClasses=c("character","character",rep("numeric",5)))

#Show connectivity patterns in doxorubicin
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
#stop()

#Show connectivity patterns in ethanol
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


#Estimate the probability of getting the given component size distributions with randomly sampling edges
Cairo::CairoPDF(file=paste(c(connectivity_output_path,'component_size_significance.pdf'),collapse='/'),width=12,height=3)
component_sig_matrix <- network_simulation_significance(pca_universe,
                                                        pca_enhanced,
                                                        pca_depleted,
                                                        iterations=connectivity_iterations)
connectivity_graph(component_sig_matrix,my_color_list)
dev.off()

#Estimate the probability of getting the given component size distributions with randomly sampling nodes
Cairo::CairoPDF(file=paste(c(connectivity_output_path,'component_size_significance_nodewise.pdf'),collapse='/'),width=12,height=3)
component_sig_matrix_nodewise <- network_simulation_significance(pca_universe,
                                                                 pca_enhanced,
                                                                 pca_depleted,
                                                                 iterations=connectivity_iterations,
                                                                 mode = 'nodewise')
connectivity_graph(component_sig_matrix_nodewise,my_color_list)
dev.off()

#Estimate the probability of getting the given graph density distributions with randomly sampling edges
# Cairo::CairoPDF(file=paste(c(connectivity_output_path,'density_significance.pdf'),collapse='/'),width=12,height=3)
# density_sig_matrix <- network_simulation_significance(pca_universe,
#                                                       pca_enhanced,
#                                                       pca_depleted,
#                                                       iterations=connectivity_iterations,
#                                                       metric='density')
# connectivity_graph(density_sig_matrix,my_color_list,legend_labels=c('-1'="Decreased density",
#                                                                     '0'="Expected density",
#                                                                     '1'="Increased density"))
# dev.off()
# 
# #Estimate the probability of getting the given graph density distributions with randomly sampling nodes
# Cairo::CairoPDF(file=paste(c(connectivity_output_path,'density_significance_nodewise.pdf'),collapse='/'),width=12,height=3)
# connectivity_sig_matrix <- network_simulation_significance(pca_universe,
#                                                            pca_enhanced,
#                                                            pca_depleted,
#                                                            mode = 'nodewise',
#                                                            iterations=connectivity_iterations,
#                                                            metric='density')
# connectivity_graph(connectivity_sig_matrix,my_color_list,legend_labels=c('-1'="Decreased density",
#                                                                          '0'="Expected density",
#                                                                          '1'="Increased density"))
# dev.off()


#Estimate the probability of getting the given graph component size distributions with randomly sampling differing proportions
#of nodes and edges

#We only want to plot things with more than n_changes
n_changes <- 50
enh_counts <- table(pca_enhanced$Condition)
depl_counts <- table(pca_depleted$Condition)
conditions <- unique(c(names(enh_counts),names(depl_counts)))
change_counts <- c()
for(condition in conditions){
  change_counts[condition] <- sum(c(enh_counts[condition],depl_counts[condition]),na.rm=T)
}
conditions_plotted <- names(which(change_counts > n_changes))

Cairo::CairoPDF(file=paste(c(connectivity_output_path,'component_size_significance_search.pdf'),collapse='/'),width=9,height=5.33)
component_size_sig_search_matrix <- network_simulation_significance_node_edge_search_matrix(pca_universe = pca_universe,
                                                                                            pca_enhanced = pca_enhanced,
                                                                                            pca_depleted = pca_depleted,
                                                                                            node_probs=c(0:10)/10,
                                                                                            iterations=100,#connectivity_iterations,
                                                                                            load_saved=F,
                                                                                            save_output=T,
                                                                                            save_directory=saved_parameter_path,
                                                                                            save_filename='parameter_search_component_size.tsv')

processed_component_size_sig_search_matrix <- process_node_edge_sig_matrix(component_size_sig_search_matrix,conditions_plotted)
name_order <- names(sort(apply(processed_component_size_sig_search_matrix,1,sum),decreasing = T))
processed_component_size_sig_search_matrix <- processed_component_size_sig_search_matrix[name_order,]
node_edge_search_heatmap(processed_component_size_sig_search_matrix,
                         conditions_plotted=conditions_plotted,
                         my_color_function=blue_black_orange)
dev.off()

#Estimate the probability of getting the given graph density distributions with randomly sampling differing proportions
#of nodes and edges
Cairo::CairoPDF(file=paste(c(connectivity_output_path,'density_significance_search.pdf'),collapse='/'),width=9,height=5.33)
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

processed_density_sig_search_matrix <- process_node_edge_sig_matrix(density_sig_search_matrix,conditions_plotted)
processed_density_sig_search_matrix <- processed_density_sig_search_matrix[name_order,]
node_edge_search_heatmap(processed_density_sig_search_matrix,
                         conditions_plotted=conditions_plotted,
                         my_color_function=blue_black_orange,
                         legend_labels=c('-1'="Increased density (p < 0.05)",
                                         '0'="Expected density",
                                         '1'="Decreased density (p < 0.05)",
                                         '2'="Most consistent simulation"))
dev.off()


#Shows which hubs have a bias of accumulating or depleting interactions
Cairo::CairoPDF(file=paste(c(connectivity_output_path,'hub_bias_heatmap.pdf'),collapse='/'),width=6,height=8.5)
hub_bias_heatmap(hub_df,blue_black_orange,nonsig_colour = 'grey90')
dev.off()

#Shows bias over different conditions
Cairo::CairoPDF(file=paste(c(connectivity_output_path,'hub_bias_over_conds.pdf'),collapse='/'),width=8,height=8)
bias_over_conditions(hub_df)
dev.off()

#For doxorubicin, shows the component size distribution of randomly sampling edges
Cairo::CairoPDF(file=paste(c(connectivity_output_path,'connnectivity_hist_AA_edge.pdf'),collapse='/'),width=10,height=5)
connectivity_histogram(pca_universe,
                       pca_enhanced,
                       pca_depleted,
                       'AA-mixture',
                       xlim=c(0,90))
dev.off()

Cairo::CairoPDF(file=paste(c(connectivity_output_path,'connnectivity_hist_AA_node.pdf'),collapse='/'),width=10,height=5)
connectivity_histogram(pca_universe,
                       pca_enhanced,
                       pca_depleted,
                       'AA-mixture',xlim=c(0,90),sampling_mode='nodewise',prob_node=1)
dev.off()
