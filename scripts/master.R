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
#Frequency Perturbation
#Connectivity
to_run <- c('Make figures')#'Reviewer update')#'Expression PCA')#'Make figures')#'RBD2','Make figures')#c('Frequency Perturbation')#c('Atorvastatin enrichment')

#Markdown directory
figure_path <- "../results/rmarkdown_figures"

#Package containing necessary scripts
devtools::load_all('../packages/bcPcaAnalysis')
devtools::document('../packages/bcPcaAnalysis')

#Global parameters
pca_universe = '../data/external/Additional.file.6.txt'
pca_enhanced_calls = '../data/external/Additional.file.10.txt'
pca_depleted_calls = '../data/external/Additional.file.11.txt'
go_association_file = '../data/funcassociate_go_associations.txt'

protein_abundance_file = '../data/paxdb_abundance.tsv'
protein_abundance_file_msonly = '../data/paxdb_abundance_msms_godoy.txt'
expression_file = '../data/external/Additional.file.14.txt'

hub_enrichment_file = '../data/external/Data for Figure 3D.xlsx'
rbd2_mtx_data_file <- '../data/external/Additional.file.15.txt'
rbd2_no_mtx_data_file <- '../data/external/Additional.file.16.txt'


output_path <- '../results/master_output'
saved_parameter_path <- '../data/script_input'
heatmap_cluster_file <- 'cluster1_heatmap.tsv'
figure_path <- '../results/rmarkdown_figures'

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
#Uninteresting result, did not bother giving output, left here for record
if('Atorvastatin enrichment' %in% to_run){
  source('subtasks/atorvastatin_enrichment.R')
}

#Creates a summary heatmap of the data, emulating Cluster 3.0 output
if('Heatmap' %in% to_run){
  source('subtasks/summary_heatmap.R')
}

#Simulation of connected component sizes and subnetwork density of node-based versus edge-based sampling
if('Connectivity' %in% to_run){
  source('subtasks/connectivity_simulations.R')
}

#Functional enrichment compared to frequency of perturbation
if('Frequency Perturbation' %in% to_run){
  source('subtasks/go_enrichment_vs_frequency.R')
}

#Check for GO enrichment in each condition, enhanced and depleted, nodewise and edgewise
if('GO enrichment' %in% to_run){
  source('subtasks/go_enrichment_vs_condition.R')
}

#Mass-action based modelling of protein complex level changes
if('Expression PCA' %in% to_run){
  #For expression mass action analysis
  q_val_cutoff <- 0.05
  measurement_diff_cutoff <- 1
  effect_size_cutoff <- 0.25
  source('subtasks/mass_action_prediction.R')
}

#Makes a graph of the RBD2 CRISPRi data
if('RBD2' %in% to_run){
  dir.create(paste(c(output_path,'rbd2'),collapse='/'),showWarnings=F)
  
  rbd2_mtx_data <- read.csv(rbd2_mtx_data_file,head=T,sep='\t')
  #rbd2_no_mtx_data <- read.csv(rbd2_no_mtx_data_file,head=T,sep='\t')
  rbd2_data <- cbind(rbd2_mtx_data[,1:9],rbd2_mtx_data[,'FC.avg'])# - rbd2_no_mtx_data[,'FC.avg'])
  colnames(rbd2_data)[ncol(rbd2_data)] <- 'FC.avg'
  rbd2_mtx_data_up <- filter(rbd2_data,UP.DN == 'uptag')
  rbd2_data_dn <- filter(rbd2_data,UP.DN == 'downtag')
  rbd2_data_comb <- rbd2_data#cbind(rbd2_data_up[,c('Gene.1','Gene.2')],(rbd2_data_up$FC.avg+rbd2_data_dn$FC.avg)/2)
  colnames(rbd2_data_comb)[ncol(rbd2_data_comb)] <- 'FC.avg'
  first_neighbours <- filter(rbd2_data_comb,Gene.1 == 'RBD2' | Gene.2 == 'RBD2')[,'FC.avg']
  first_neihbour_genes <- unique(as.vector(unlist(filter(rbd2_data_comb,Gene.1 == 'RBD2' | Gene.2 == 'RBD2')[,c('Gene.1','Gene.2')])))
  second_neighbours <- filter(rbd2_data_comb,(Gene.1 %in% first_neihbour_genes | Gene.2 %in% first_neihbour_genes) & Gene.2 != 'RBD2' & Gene.1 != 'RBD2')[,'FC.avg']
  all_others <- filter(rbd2_data_comb,!(Gene.1 %in% first_neihbour_genes | Gene.2 %in% first_neihbour_genes))[,'FC.avg']
  
  #Make histogram
  nbreaks <- 50
  col1 <- col2rgb(RColorBrewer::brewer.pal(12,'Set3')[3])/255
  col2 <- col2rgb(RColorBrewer::brewer.pal(12,'Set3')[5])/255
  darkening_factor <- 1
  opacity_factor <- 0.7
  col_list <- c(
    rgb(col1[1]*darkening_factor,col1[2]*darkening_factor,col1[3]*darkening_factor,opacity_factor),
    rgb(col2[1]*darkening_factor,col2[2]*darkening_factor,col2[3]*darkening_factor,opacity_factor),
    rgb(0.3,0.3,0.3,opacity_factor)
  )
  
  Cairo::CairoPDF(file=paste(c(output_path,'rbd2','rbd2_knockdown.pdf'),collapse='/'),width=6,height=6)
  real_breaks <- hist(rbd2_data[,ncol(rbd2_data)],plot=F,breaks=nbreaks)$breaks
  all_other_hist <- hist(all_others,plot=F,breaks=real_breaks)
  all_other_hist$counts <- all_other_hist$density
  secondary_hist <- hist(second_neighbours,plot=F,breaks=real_breaks)
  secondary_hist$counts <- secondary_hist$density
  primary_hist <- hist(first_neighbours,plot=F,breaks=real_breaks)
  primary_hist$counts <- primary_hist$density
  par(mar=c(5,5,0,0))
  plot(all_other_hist,col=col_list[3],xaxt='n',xlab=c('Log2(R) [+ATC/-ATC]'),ylab='Density',main='',xlim=c(-4,4),ylim=c(0,1.5),cex.lab=1.3)
  axis(1,at=c(-4:4))
  plot(secondary_hist,col=col_list[2],add=T)
  plot(primary_hist,add=T,col=col_list[1])
  legend('left',legend=c('1° Rbd2 interactions','2° Rbd2 interactions','Other interactions'),
         cex=1,bty='n',fill=col_list)#,pch=22,pt.cex=2.4)
  dev.off()
}

if('Reviewer update' %in% to_run){
  source('subtasks/reviewer_update.R')
  
}

#Creates figures using Rmarkdown
if('Make figures' %in% to_run){
  rmarkdown::render("manuscript_latex/figures/BC-PCA_Figure.1.Rmd", "pdf_document",output_dir=figure_path)
  rmarkdown::render("manuscript_latex/figures/BC-PCA_Figure.2.Rmd", "pdf_document",output_dir=figure_path)
  rmarkdown::render("manuscript_latex/figures/BC-PCA_Figure.3.Rmd", "pdf_document",output_dir=figure_path)
  rmarkdown::render("manuscript_latex/figures/BC-PCA_Figure.4.Rmd", "pdf_document",output_dir=figure_path)
  rmarkdown::render("manuscript_latex/figures/BC-PCA-Figure.S3.Rmd", "pdf_document",output_dir=figure_path)
  rmarkdown::render("manuscript_latex/figures/BC-PCA-Figure.S4.Rmd", "pdf_document",output_dir=figure_path)
}

