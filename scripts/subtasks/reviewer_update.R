require(RCurl)
require(XML)


##Maslov simulations, shuffled and assigned Kds
expression_pca_output_path <- paste(c(output_path,'expression_pca'),collapse='/')
dir.create(expression_pca_output_path, showWarnings = FALSE)

for(time_point in c('30min','1h','4h','12h')){
  condition_regexp <- paste(c('Ethanol',time_point),collapse='.')
  
  #We'll try everything with shuffled assigned Kds
  shuffle_status <- c(T,F)
  for(shuffle in shuffle_status){
    
    if(shuffle){
      filename <- paste(c('bcPCA_mRNA_predictions_maslov_shuffled',time_point,'.pdf'),collapse='')
    }else{
      filename <- paste(c('bcPCA_mRNA_predictions_maslov_',time_point,'.pdf'),collapse='')
    }
    predictions <- pca_ma_prediction(pca_universe,
                                     protein_abundance_file,
                                     expression_file,
                                     'ethanol',
                                     expression_condition_regexp=condition_regexp,
                                     prediction_method = "maslov",
                                     shuffled_kds = shuffle)
    pca_ma_prediction_plot(predictions,expression_pca_output_path,filename=filename)
    
    #Also do independend shuffled ones
    if(shuffle){
      filename <- paste(c('bcPCA_mRNA_predictions_shuffled',time_point,'.pdf'),collapse='')
      predictions <- pca_ma_prediction(pca_universe,
                                       protein_abundance_file,
                                       expression_file,
                                       'ethanol',
                                       expression_condition_regexp=condition_regexp,
                                       shuffled_kds = shuffle)
      pca_ma_prediction_plot(predictions,expression_pca_output_path,filename=filename)
    }
  }
}


#Baseline growth plot with full equation
expression_pca_baseline_output_path <- paste(c(output_path,'expression_pca_baseline'),collapse='/')
dir.create(expression_pca_baseline_output_path, showWarnings = FALSE)

shuffle_status <- c(T,F)
prediction_methods <- c('maslov','independent')
datasets <- list('whole'=protein_abundance_file,
                 'ms'=protein_abundance_file_msonly)
for(prediction_method in prediction_methods){
  for(dataset_name in names(datasets)){
    for(shuffle in shuffle_status){
    predictions <- pca_ma_prediction(pca_universe,
                                     abundance_file = datasets[[dataset_name]],
                                     expression_file,
                                     condition = 'ethanol',
                                     shuffled_kds = shuffle,
                                     prediction_method = prediction_method,
                                     return_only_baseline = T)
    if(shuffle){
      filename <- paste(c('bcPCA_baseline_predictions_shuffled',dataset_name,'_',prediction_method,'.pdf'),collapse='')
    }else{
      filename <- paste(c('bcPCA_baseline_predictions_',dataset_name,'_',prediction_method,'.pdf'),collapse='')
    }
    outfile <- paste(c(expression_pca_baseline_output_path,filename),collapse='/')
    CairoPDF(file=outfile,width=6,height=6)
    
    #Populate the title
    main_text_vector <- c('Predictions - ')
    if(prediction_method == 'maslov'){
      main_text_vector <- paste(c(main_text_vector,'Competitive binding model'),collapse='')
    }else if(prediction_method =='independent'){
      main_text_vector <- paste(c(main_text_vector,'Independent binding model'),collapse='')
    }
    if(shuffle){
      main_text_vector <- c(main_text_vector,'Shuffled estimated Kds')
    }else{
      main_text_vector <- c(main_text_vector,'Estimated Kds')
    }
    
    if(dataset_name == 'ms'){
      main_text_vector <- c(main_text_vector,'MS data')
    } else if(dataset_name == 'whole'){
      main_text_vector <- c(main_text_vector,'All data')
    }
    
    predictions$Log2_MA_prediction <- log10(2^predictions$Log2_MA_prediction)
    predictions$AUC <- log10(predictions$AUC)
    
    predictions <- predictions[predictions$AUC > log10(8),]
    
    
    plot(predictions$Log2_MA_prediction,
         predictions$AUC,
         main=main_text_vector,
         pch=16,
         col=rgb(0,0,0,0.3),
         cex=0.7,
         xlab = 'Predicted abundance [log10 ppm]',
         ylab = 'Isogenic culture [log10 AUC]')
    
    min_pred <- min(predictions$Log2_MA_prediction,na.rm=T)
    max_pred <- max(predictions$Log2_MA_prediction,na.rm=T)
    
    min_auc <- min(predictions$AUC,na.rm=T)
    max_auc <- max(predictions$AUC,na.rm=T)
    text(
      x = min_pred + (max_pred - min_pred)/10,
      y = max_auc - (max_auc - min_auc)/10,
      labels=paste(c('r=',format(cor(predictions$Log2_MA_prediction,predictions$AUC,use='pair'),digits=2)),collapse=''),
      cex=2
    )
    abline(lm(predictions$AUC~predictions$Log2_MA_prediction),col='red',lwd=2,lty=5)
    dev.off()
    }
  }
}

#Spuriously significant results
pca_univ <- read.table(pca_universe,head=T,sep='\t')

spurious_sig_output_path <- paste(c(output_path,'spurious_significant'),collapse='/')
dir.create(spurious_sig_output_path, showWarnings = FALSE)
nsigs <- sapply(1:100,function(i){
  test <- randomize_pca_file(pca_univ)
  nsigs <- get_n_sig_genes(test)
  write(paste(c(nsigs,'spurious'),collapse=' '), stderr())
  return(nsigs)
})
outfile <- paste(c(spurious_sig_output_path,'spurious_hist.pdf'),collapse='/')
CairoPDF(file=outfile,width=6,height=6)
real_sigs <- log10(dim(read.table(pca_depleted_calls,sep='\t'))[1] + dim(read.table(pca_enhanced_calls,sep='\t'))[1])

all_freq <- hist(log10(nsigs),plot=F,breaks=30)

plot(all_freq,col='black',xlab='log10(Number of Result found Significant)',main='',xlim=c(0,real_sigs))
abline(v=real_sigs,col='red',lwd=2,lty=2)
text(real_sigs-0.02*real_sigs,max(all_freq$counts)/10,'Observed in real data',srt=90,adj=0)
legend('topleft',c('Shuffled data'),fill='black')

dev.off()

#Kinase substrate overlap, no output, have to write overlap_table to something
pca_univ <- read.table(pca_universe,head=T,sep='\t')
pca_ppis <- unique(pca_univ[,c('ORF.1','ORF.2')])
pca_ppi_vec <- as.vector(apply(pca_ppis,1,function(x){paste(sort(x),collapse='_')}))

kin_sub <- read.table('../data/kinase_substrate_2011.tsv',sep='\t',head=T,stringsAsFactors = T)
kin_sub_vec <- as.vector(apply(kin_sub,1,function(x){paste(sort(x),collapse='_')}))

overlap <- intersect(pca_ppi_vec,kin_sub_vec)

overlap_table <- t(sapply(overlap,function(x){
  pair <- strsplit(x,split='_')[[1]]
  return(map_gene_names(pair))
}))


#Effect size distributions compared to Levy 2004
levy_filename <- '../data/levy_2014_ab.tsv'
levy_data <- read.table(levy_filename,sep='\t',head=T)
localization_dist <- log2(levy_data$AB.MEM.AGENT/levy_data$AB.CYTO.AGENT)
localization_dist <- localization_dist[!is.na(localization_dist)]

expression_dist <- -log2(levy_data$AB.GFP.YMD/levy_data$AB.GFP.YPD)
expression_dist <- expression_dist[!is.na(expression_dist)]

predictions <- pca_ma_prediction(pca_universe,
                                 protein_abundance_file,
                                 expression_file,
                                 'ethanol',
                                 prediction_method = "independent",
                                 shuffled_kds = F)

density_output_path <- paste(c(output_path,'levy_comparison'),collapse='/')
dir.create(density_output_path, showWarnings = FALSE)
outfile <- paste(c(density_output_path,'effect_size_dist.pdf'),collapse='/')

CairoPDF(file=outfile,width=8.5,height=6)
d1 <- density(predictions$bcPCA_FC.AVG,adj=1.5,from=-6,to=6)
plot(d1,type='n',main='',xlab='Effect Size [Log2(R)]',cex.lab=1.5)
polygon(d1,col=rgb(0,0,0.5,0.3),border=NA)#
lines(d1,lwd=4,col=rgb(0,0,0.5,0.8))

d2 <- density(localization_dist,adj=1.5,from=-6,to=6)
polygon(d2,col=rgb(0.5,0,0,0.3),border=NA)
lines(d2,lwd=4,col=rgb(0.5,0,0,0.8))

d3 <- density(expression_dist,adj=1.5,from=-6,to=6)
polygon(d3,col=rgb(0,0.5,0,0.3),border=NA)
lines(d3,lwd=4,col=rgb(0,0.5,0,0.8))

d4 <- density(predictions$Log2_MA_prediction,adj=1.5,from=-6,to=6)
polygon(d4,col=rgb(0.5,0.7,1,0.3),border=NA)
lines(d4,lwd=4,col=rgb(0.5,0.7,1,0.8))


legend('topleft',legend=c('Membrane vs Cytosol','YPD vs minimal','Ethanol BC-PCA (Observed)','Ethanol BC-PCA (mRNA Predicted)'),fill=c('darkred','darkgreen','darkblue',rgb(0.5,0.7,1)))
dev.off()


#mRNA Levy Comparison
levy_filename <- '../data/levy_2014_ab.tsv'
mRNA_filename <- expression_file
levy_output_path <- paste(c(output_path,'levy_comparison'),collapse='/')
dir.create(levy_output_path, showWarnings = FALSE)
outfile <- paste(c(levy_output_path,'mRNA_0h_vs_levy_abundance.pdf'),collapse='/')

CairoPDF(file=outfile,width=9,height=5)
mRNA_levy_comparison(levy_filename,mRNA_filename)
dev.off()


#Extract available Kds
interpret_kd <- function(kd_string){
  unit_map <- c('nM' = 10^(-9),
                'uM' = 10^(-6))
  split_string <- strsplit(kd_string,split='Kd.')[[1]]
  kd_part <- split_string[length(split_string)]              
  
  kd_num <- substr(kd_part,1,nchar(kd_part)-2)
  kd_units <- substr(kd_part,nchar(kd_part)-1,nchar(kd_part))
  
  return(as.numeric(kd_num)*unit_map[[kd_units]])
  #print(kd_units)
}


pdb_index <- '../data/PDB_kd_list'

pdb_index <- read.csv(pdb_index,sep='\t',head=F,stringsAsFactors = F,row.names=1)

ids <- rownames(pdb_index)
yeast_ids <- c()
for(id in ids){
  webpage <- getURL(sprintf("http://www.rcsb.org/pdb/rest/describeMol?structureId=%s",id))
  webpage <- readLines(tc <- textConnection(webpage))
  close(tc)
  if(length(grep('Saccharomyces',webpage)) > 1){
    yeast_ids <- c(yeast_ids,id)
  }
}


sgd_entries <- c()
for(yeast_id in yeast_ids){
  webpage <- getURL(sprintf("http://www.rcsb.org/pdb/rest/describeMol?structureId=%s",yeast_id))
  webpage <- readLines(tc <- textConnection(webpage)); close(tc)
  xml_page <- xmlParse(webpage,asText = T)
  page_list <- xmlToList(xml_page)
  polymers <- grep('polymer',names(page_list$structureId))
  if(length(polymers)  == 2){
    sgd_ids <- c()
    for(polymer in polymers){
      uniprot_id <- page_list$structureId[[polymer]]$macroMolecule$accession
      webpage <- getURL(sprintf("http://www.uniprot.org/uniprot/%s.txt",uniprot_id))
      webpage <- readLines(tc <- textConnection(webpage)); close(tc)
      sgd_id_line <- webpage[grep('DR   SGD',webpage)]
      split_line <- strsplit(sgd_id_line,split='; ')[[1]]
      sgd_id <- strsplit(split_line[length(split_line)],split='[.]')[[1]]
      sgd_id <- reverse_map_gene_names(sgd_id)
      sgd_ids <- c(sgd_ids,sgd_id)
      #print(sgd_id)
    }
    kd <- pdb_index[yeast_id,]
    #interpreted_kd <- interpret_kd(kd)
    sgd_entries <- rbind(sgd_entries,c(yeast_id,sgd_ids,kd))
    
    #stop()
  }
}

pca_univ <- read.table(pca_universe,head=T,stringsAsFactors = F)
ppis <- unique(pca_univ[,c('ORF.1','ORF.2')])
ppis <- apply(ppis,1,function(ppi){
  paste(sort(ppi),collapse='_')
})

ppis_with_kd <- apply(sgd_entries,1,function(entry){
  paste(sort(entry[2:3]),collapse='_')
})
