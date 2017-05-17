randomize_pca_file <- function(pca_file,control_regexp='DMSO',condition_regexp='Treatment'){
  dmso_cols <- grep(control_regexp,colnames(pca_file))
  treatment_cols <- grep(condition_regexp,colnames(pca_file))
  
  both_cols <- c(dmso_cols,treatment_cols)
  
  
  #seed_val <- round(runif(1,max=1000))
  
  pca_file_up <- dplyr::filter(pca_file, UP.DN == 'uptag')
  pca_file_dn <- dplyr::filter(pca_file, UP.DN == 'downtag')
  
  sampled_cols <- t(sapply(1:nrow(pca_file_up),function(x){
    return(sample(both_cols))
  }))
  
  test_up <- t(sapply(1:nrow(sampled_cols),function(i){
    return(as.numeric(pca_file_up[i,sampled_cols[i,]]))
  }))
#  
  test_dn <- t(sapply(1:nrow(sampled_cols),function(i){
    return(as.numeric(pca_file_dn[i,sampled_cols[i,]]))
  }))
  
  
  new_vals <- rbind(test_up,test_dn)
  colnames(new_vals) <- colnames(pca_file)[both_cols]
  #pca_file[,both_cols] <- new_vals
  #pca_file <- new_vals
  pca_file <- cbind(pca_file[,-both_cols],new_vals)
  return(pca_file)
  
}

#randomize_pca_file <- function(pca_file){
#  
#}



get_n_sig_genes <- function(pca_file,p_cutoff = 0.05,control_regexp='DMSO',condition_regexp='Treatment'){
  condsum <- 0
  for(condition in unique(pca_file$Condition)){
    pca_cond <- dplyr::filter(pca_file, Condition == condition)
    pca_cond_up <- dplyr::filter(pca_cond, UP.DN == 'uptag')
    pca_cond_dn <- dplyr::filter(pca_cond, UP.DN == 'downtag')
    dmso_cols <- grep('DMSO',colnames(pca_cond))
    treatment_cols <- grep('Treatment',colnames(pca_cond))
    
    design_matrix <- matrix(nrow=length(dmso_cols)+length(treatment_cols),ncol=2)
    rownames(design_matrix) <- c(colnames(pca_cond)[dmso_cols],colnames(pca_cond)[treatment_cols])
    design_matrix[,1] <- c(rep(1,length(dmso_cols)),rep(0,length(treatment_cols)))
    design_matrix[,2] <- 1 - design_matrix[,1]
    colnames(design_matrix) <- c('DMSO','cond')
    cont.matrix <- makeContrasts(condvsDMSO=cond-DMSO, levels=design_matrix)
    #
    eset_up <- pca_cond_up[,rownames(design_matrix)]
    fit_up <- lmFit(eset_up, design_matrix)
    fit2_up <- contrasts.fit(fit_up, cont.matrix)
    fit2_up <- p.adjust(eBayes(fit2_up)$p.val,method='BH')
    
    eset_dn <- pca_cond_dn[,rownames(design_matrix)]
    fit_dn <- lmFit(eset_dn, design_matrix)
    fit2_dn <- contrasts.fit(fit_dn, cont.matrix)
    fit2_dn <- p.adjust(eBayes(fit2_dn)$p.val,method='BH')
    
    #print(condition)
    condsum <- condsum + sum(fit2_dn < p_cutoff & fit2_up < p_cutoff)
    #print(sum(fit2_dn < p_cutoff & fit2_up < p_cutoff))
    #}
  }
  return(condsum)
}


filter_mtx_defects <- function(pca_file,excluded_strains){
  ret <- c()
  for(condition in unique(pca_file$Condition)){
    pca_cond <- dplyr::filter(pca_file, Condition == condition)
    excl_cond <- dplyr::filter(excluded_strains, Condition == condition)
    ret <- rbind(ret,filter(pca_cond, !(Tag %in% excl_cond$Tag)))
    
  }
  return(ret)
}

