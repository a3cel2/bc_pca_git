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
                                   line_width=3,
                                   row_names=F
)

inter_filename <- paste(c(saved_parameter_path,heatmap_cluster_file),collapse='/')
inters <- as.vector(read.table(inter_filename)[,1])
convert_pca_file_to_heatmap_format(pca_file=pca_file,
                                   pca_enhanced_calls = pca_enh,
                                   pca_depleted_calls = pca_depl,
                                   output_path=heatmap_output_path,
                                   filename='clusters_sig_bpcs.png',
                                   draw=F,
                                   label_size=4,
                                   png_height=1500,
                                   gene_dendrogram_width=0.12,
                                   #png_width=1700,
                                   row_label_size=3.5,
                                   line_width=3,
                                   interaction_subset=inters,
                                   wide_margins=T,
                                   legend=F
)

condition_summary_barplot(pca_enh,
                          pca_depl,
                          color_function=blue_black_orange,
                          draw=F,
                          output_path=heatmap_output_path,
                          filename='number_dynamic_bpcs.pdf')