#Maintenance
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

devtools::load_all('../../packages/bcPcaAnalysis')
devtools::document('../../packages/bcPcaAnalysis')

pca_universe = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.6.txt'
pca_enhanced_calls = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.10.txt'
pca_depleted_calls = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.11.txt'
abundance_file = "/Users/Albi/Dropbox/Roth Lab/projects/bc_pca_git/data/paxdb_abundance.tsv"

pca_univ <- read.csv(pca_universe,sep='\t',stringsAsFactors=F)
pca_enh <- read.csv(pca_enhanced_calls,sep='\t',stringsAsFactors=F)
pca_depl <- read.csv(pca_depleted_calls,sep='\t',stringsAsFactors=F)
merged_pca_calls <- rbind(pca_enh,pca_depl)
abundance_file <- read.table(abundance_file,head=F,stringsAsFactors = F, row.names=1)

landry <- read.csv('../../data/scratch/landry_data.tsv',sep='\t')
output_df <- data.frame()
for(i in 1:nrow(landry)){
  pair <- landry[i,2:3]
  pair <- as.vector(as.matrix(pair))
  #print(pair)
  abundances <- get_orf_pair_abundance(pair,abundance_file)
  #stop()
  #mRNA_changes <- get_orf_mrna_changes(pair,expression_file,expression_control_regexp,expression_condition_regexp)
  
  #mRNA_changes <- 2^c(mean(as.numeric(esr[pair[1],grep(condition_grep_val,colnames(esr))])),
  #                  mean(as.numeric(esr[pair[2],grep(condition_grep_val,colnames(esr))])))
  
  #mRNA_changes <- 2^(c(as.numeric(diauxic_vals[pair[1],]),as.numeric(diauxic_vals[pair[2],])))
  mRNA_changes <- 2^as.vector(as.matrix(landry[i,6:7]))
  prediction <- log2(mass_action_predictor(abundances[1],abundances[2],mRNA_changes[1],mRNA_changes[2]))
  
  output <- data.frame(ORF1=pair[1],
                       ORF2=pair[2],
                       #bcPCA_FC.UP=merged_pca_calls[i,'FC.UP'],
                       #bcPCA_FC.DN=merged_pca_calls[i,'FC.DN'],
                       #bcPCA_FC.avg=merged_pca_calls[i,'FC.avg'],
                       bcPCA_FC.AVG=landry[i,'I.score'],
                       #bcPCA_qVal=merged_pca_calls[i,'q.val'],
                       PaxDB_Abundance_ORF1=abundances[1],
                       PaxDB_Abundance_ORF2=abundances[2],
                       mRNAFC_ORF1=mRNA_changes[1],
                       mRNAFC_ORF2=mRNA_changes[2],
                       Log2_MA_prediction=prediction)
  output_df <- rbind(output_df,output)
}

stop()


esr <- read.csv('../../data/scratch/complete_dataset.txt',sep='\t')
rownames(esr) <- esr$UID
esr[esr == "#VALUE!"] <- NA
esr[,3:ncol(esr)] <- apply(esr[,3:ncol(esr)],2,function(x){as.numeric(as.vector(x))})

diauxic <- read.csv('../../data/scratch/2010.diauxic.pcl.txt',sep='\t')
#rownames(diauxic) <- diauxic$YORF
diauxic_vals <- diauxic[2:nrow(diauxic),c(1,4:ncol(diauxic))]
diauxic_vals <- summarize(group_by(diauxic_vals,YORF),expr=mean(X20.5.hr))
rownames(diauxic_vals) <- diauxic_vals$YORF
diauxic_vals <- as.data.frame(diauxic_vals)
diauxic_vals <- diauxic_vals[,2,drop=F]

#Set values here
condition <- 'H2O2'
condition_grep_val <- 'DBY7286.*H2O2'

#merged_pca_calls <- filter(merged_pca_calls,Condition==condition)
merged_pca_calls <- merge_pca_file(pca_univ,condition=condition)
retlist <- c()
output_df <- data.frame()
orf_pairs <- merged_pca_calls[,c('ORF.1','ORF.2')]

for(i in 1:nrow(orf_pairs)){
  pair <- orf_pairs[i,]
  pair <- as.vector(as.matrix(pair))
  #print(pair)
  abundances <- get_orf_pair_abundance(pair,abundance_file)
  #stop()
  #mRNA_changes <- get_orf_mrna_changes(pair,expression_file,expression_control_regexp,expression_condition_regexp)
  
  mRNA_changes <- 2^c(mean(as.numeric(esr[pair[1],grep(condition_grep_val,colnames(esr))])),
                    mean(as.numeric(esr[pair[2],grep(condition_grep_val,colnames(esr))])))
  
  #mRNA_changes <- 2^(c(as.numeric(diauxic_vals[pair[1],]),as.numeric(diauxic_vals[pair[2],])))
  prediction <- log2(mass_action_predictor(abundances[1],abundances[2],mRNA_changes[1],mRNA_changes[2]))
  
  output <- data.frame(ORF1=pair[1],
                       ORF2=pair[2],
                       bcPCA_FC.UP=merged_pca_calls[i,'FC.UP'],
                       bcPCA_FC.DN=merged_pca_calls[i,'FC.DN'],
                       bcPCA_FC.avg=merged_pca_calls[i,'FC.avg'],
                       #bcPCA_FC.AVG=mean(c(merged_pca_calls[i,'FC..UP.tag.'],merged_pca_calls[i,'FC..DN.tag.'])),
                       #bcPCA_qVal=merged_pca_calls[i,'q.val'],
                       PaxDB_Abundance_ORF1=abundances[1],
                       PaxDB_Abundance_ORF2=abundances[2],
                       mRNAFC_ORF1=mRNA_changes[1],
                       mRNAFC_ORF2=mRNA_changes[2],
                       Log2_MA_prediction=prediction)
  output_df <- rbind(output_df,output)
}
stop()

ppi_freq_list <- sort(table(c(
  #as.vector(pca_enh$PPI.short)
  as.vector(pca_depl$PPI.short)
)
  ),decreasing=T)


esr <- read.csv('../../data/scratch/complete_dataset.txt',sep='\t')
rownames(esr) <- esr$UID
esr[esr == "#VALUE!"] <- NA
esr[,3:ncol(esr)] <- apply(esr[,3:ncol(esr)],2,function(x){as.numeric(as.vector(x))})
mean_expression_vals <- apply(esr[,3:ncol(esr)],1,function(x){mean(x,na.rm=T)})

esr_genes <- sapply(as.vector(esr$NAME),function(x){strsplit(x,split=' ')[[1]][1]})
esr_positive_genes <- esr_genes[mean_expression_vals > 0]
esr_negative_genes <- esr_genes[mean_expression_vals < 0]

expression_file = '/Users/Albi/Dropbox/barcoded-PCA/2015-08-30/Additional.file.14.txt'
expression_file <- read.table(expression_file,head=T,stringsAsFactors = F)
expression_file <- simplify_expression_file(expression_file)
expression_control_regexp='Ethanol.0h'
expression_condition_regexp='Ethanol.4h'

avg_vals <- apply(expression_file[,grep(expression_condition_regexp,colnames(expression_file)),drop=F],1,mean)/
  apply(expression_file[,grep(expression_control_regexp,colnames(expression_file)),drop=F],1,mean)

plot(log2(avg_vals),apply(esr[names(avg_vals),grep('ethanol',colnames(esr)),drop=F],1,mean),
     xlab='Values Measured by Uli (Log2)',
     ylab='Values Measured by Audrey Gasch')
#mRNA_changes <- get_orf_mrna_changes(pair,expression_file,expression_control_regexp,expression_condition_regexp)

stop()

#esr_ppis <- filter(pca_univ,ORF.1 %in% esr_positive_genes & ORF.2 %in% esr_positive_genes)

frequently_perturbed <- names(which(table(pca_differential$PPI.short) >= 4))

ORFs_frequently_perturbed <- unique(select(filter(pca_univ, PPI.short %in% frequently_perturbed),ORF.1,ORF.2))
ORFs_all_ppis <- unique(select(pca_univ,ORF.1,ORF.2))

ORFs_frequently_perturbed_in_esr <- filter(ORFs_frequently_perturbed,ORF.1 %in% esr_genes & ORF.2 %in% esr_genes)
ORFs_all_ppis_in_esr <- filter(ORFs_all_ppis,ORF.1 %in% esr_genes & ORF.2 %in% esr_genes)

rbind(c(dim(ORFs_frequently_perturbed_in_esr)[1],
        dim(ORFs_frequently_perturbed)[1]),
      c(dim(ORFs_all_ppis_in_esr)[1],
        dim(ORFs_all_ppis)[1])
      )

stop()



ppi_freq_plot <- function(pca_univ,
                          pca_enh,
                          pca_depl,
                          perturbation_categories=c('Static Complexes',
                                                    'Specifically variable\ncomplexes',
                                                    'Frequently variable\ncomplexes')){
  ppi_freq_list <- list()
  all_ppis <- unique(c(pca_univ$PPI.short))
  for(ppi in all_ppis){
    ppi_freq_list[[ppi]] <- 0
  }
  
  for(condition in unique(pca_enh$Condition)){
    condition_ppis <- unique(unlist(rbind(dplyr::filter(pca_enh,Condition==condition),
                                          dplyr::filter(pca_depl,Condition==condition)) %>% dplyr::select(PPI.short)))
    for(ppi in condition_ppis){
      ppi_freq_list[[ppi]] <- ppi_freq_list[[ppi]] + 1
    }
  }
  
  ppi_freq_list <- unlist(ppi_freq_list)
  
  plot_table <- as.data.frame(t(table(ppi_freq_list)))
  colnames(plot_table) <- c('Classification','ChangedConds','Freq')
  plot_table$Classification <- factor(c(perturbation_categories[1],
                                        rep(perturbation_categories[2],3),
                                        rep(perturbation_categories[3],6)),
                                      levels=perturbation_categories)
  
  ggplot(plot_table,ggplot2::aes(x=as.factor(ChangedConds),y=Freq,fill=Classification)) +
    ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge()) +
    ggplot2::scale_y_continuous(expand = c(0,0), limits=c(0,max(plot_table$Freq)*1.05)) +
    ggplot2::scale_fill_manual(values = c('grey','skyblue','darkblue')) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
                   legend.position=c(0.6,0.8),
                   legend.text = ggplot2::element_text(size=7),
                   legend.title = ggplot2::element_text(size=0.01,hjust=0)) +
    ggplot2::ylab('Number of Complexes') + 
    ggplot2::xlab('Frequency of Variability')
}

#frequent_ppi_table <- t(sapply(names(which(ppi_freq_list > 3)),function(x){strsplit(x,split=':')[[1]]}))
# frequent_ppi_table <- t(sapply(names(which(ppi_freq_list > 0 & ppi_freq_list < 4)),function(x){strsplit(x,split=':')[[1]]}))
# frequent_ppis <- igraph::graph_from_edgelist(frequent_ppi_table,directed=F)
# frequent_ppis <- keep_only_largest_connected_component(frequent_ppis)
# 
# frequent_complexes <- reverse_map_gene_names(unique(as.vector(frequent_ppi_table)))
# frequent_complexes_largest <- reverse_map_gene_names(unique(as.vector(igraph::get.edgelist(frequent_ppis))))
# all_complexes <- reverse_map_gene_names(unique(as.vector(c(pca_univ$Gene.1,pca_univ$Gene.2))))
