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
to_run <- c('Make figures')#c('Frequency Perturbation')#c('Atorvastatin enrichment')

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

protein_abundance_file = "../data/paxdb_abundance.tsv"
expression_file = '../data/external/Additional.file.14.txt'

hub_enrichment_file = '../data/external/Data for Figure 3D.xlsx'

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
  q_val_cutoff <- 1#0.05
  measurement_diff_cutoff <- 1
  effect_size_cutoff <- 0#0.25
  source('subtasks/mass_action_prediction.R')
}

#Creates figures using Rmarkdown
if('Make figures' %in% to_run){
  rmarkdown::render("manuscript_latex/figures/BC-PCA_main_figures.Rmd", "pdf_document",output_dir=figure_path)
  rmarkdown::render("manuscript_latex/figures/BC-PCA-Figure.S4.Rmd", "pdf_document",output_dir=figure_path)
  rmarkdown::render("manuscript_latex/figures/BC-PCA-Figure.S5.Rmd", "pdf_document",output_dir=figure_path)
  rmarkdown::render("manuscript_latex/figures/BC-PCA-Figure.S6.Rmd", "pdf_document",output_dir=figure_path)
  rmarkdown::render("manuscript_latex/figures/BC-PCA-Figure.S7.Rmd", "pdf_document",output_dir=figure_path)
}

