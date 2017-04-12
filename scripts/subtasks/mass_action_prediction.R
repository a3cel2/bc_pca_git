expression_pca_output_path <- paste(c(output_path,'expression_pca'),collapse='/')
dir.create(expression_pca_output_path, showWarnings = FALSE)

my_predictions <- pca_ma_prediction(pca_universe,
                                    protein_abundance_file,
                                    expression_file,
                                    'ethanol',
                                    expression_condition_regexp='Ethanol.4h')


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

#For different times
for(time_point in c('30min','1h','4h','12h')){
  condition_regexp <- paste(c('Ethanol',time_point),collapse='.')
  filename <- paste(c('bcPCA_mRNA_predictions_',time_point,'.pdf'),collapse='')
  predictions <- pca_ma_prediction(pca_universe,
                                   protein_abundance_file,
                                   expression_file,
                                   'ethanol',
                                   expression_condition_regexp=condition_regexp)
  predictions <- filter(predictions,
                                abs(bcPCA_FC.UP - bcPCA_FC.DN) <= measurement_diff_cutoff)
                                #bcPCA_qVal < q_val_cutoff,
                                #abs(bcPCA_FC.AVG) >= effect_size_cutoff)
  pca_ma_prediction_plot(predictions,expression_pca_output_path,filename=filename)
}

#Remove overly noisy
filtered_predictions <- filter(my_predictions,
                               abs(bcPCA_FC.UP - bcPCA_FC.DN) <= measurement_diff_cutoff)#,
                               #bcPCA_qVal < q_val_cutoff,
                               #abs(bcPCA_FC.AVG) > effect_size_cutoff)
pca_ma_prediction_plot(filtered_predictions,
                       expression_pca_output_path,
                       filename='bcPCA_predictions_filtered.pdf')


#Compare error_prone replicates
error_analysis <- pca_mRNA_error_analysis(pca_universe,
                                          expression_file,
                                          protein_abundance_file,
                                          iters=100,
                                          q_val_cutoff=1,
                                          fc_cutoff=0,
                                          up_dn_diff_cutoff=measurement_diff_cutoff)
pca_mRNA_error_analysis_graph(filtered_predictions,
                              error_analysis,
                              output_path=expression_pca_output_path,
                              filename='replicate_comparison_filtered',draw=F)



filtered_sig_predictions <- filter(my_predictions,
                                   abs(bcPCA_FC.UP - bcPCA_FC.DN) <= measurement_diff_cutoff,
                                   bcPCA_qVal < q_val_cutoff,
                                   abs(bcPCA_FC.AVG) > effect_size_cutoff)
error_analysis_sig <- pca_mRNA_error_analysis(pca_universe,
                                          expression_file,
                                          protein_abundance_file,
                                          iters=100,
                                          q_val_cutoff = q_val_cutoff,
                                          fc_cutoff = effect_size_cutoff,
                                          up_dn_diff_cutoff = measurement_diff_cutoff)
pca_mRNA_error_analysis_graph(filtered_sig_predictions,
                             error_analysis_sig,
                             output_path=expression_pca_output_path,
                             filename='replicate_comparison_filtered_sig',draw=F)

##For categorical
cols <- blue_black_orange(10)
pca_ma_precision_plot(my_predictions,
                      expression_pca_output_path,
                      filename='bcPCA_precision',
                      enhanced_colour=cols[8],
                      depleted_colour=cols[3]
)

#Compare different hubs
hub_comparison_graph(my_predictions,
                     'HXT1,HSP30',
                     output_path=expression_pca_output_path,
                     filename='HXT1_HSP30_comparison',
                     color_list=my_color_list,
                     node_size=25)

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



