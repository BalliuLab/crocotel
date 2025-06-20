
crossval_helper = function(Ys, X, lengths_y, rownames_y, contexts_vec, out_dir, gene_name, num_folds, alpha = 0.5){
    
  # define some commonly-used variables
  q<-length(Ys)
  m<-ncol(X)
  N<-nrow(X)
  Yhats_tiss<-vector("list", q)
  
  for(i in 1:q){
    Yhats_tiss[[i]]<-matrix(NA, ncol=1, nrow=lengths_y[i], 
                           dimnames = list(rownames_y[[i]], "pred"))
  }
  
  hom_expr_mat<-matrix(NA, nrow = nrow(X), ncol=q)
  rownames(hom_expr_mat)<-rownames(X)
  
  for(i in 1:q){
    hom_expr_mat[rownames(Ys[[i]]),i]<-Ys[[i]]
  }
  
  all_missing<-names(rowMeans(hom_expr_mat, na.rm = T)[which(is.nan(rowMeans(hom_expr_mat, na.rm = T)))])
  remove_inds<-which(rownames(hom_expr_mat) %in% all_missing)
  # These people are not in any tissues, we can just remove them from the data
  if(length(remove_inds)>0){
    #message("These individuals have NA in all Ys: ")
    #message(paste0(all_missing, collapse = ","))
    hom_expr_mat<-hom_expr_mat[-remove_inds,,drop=F]
    X<-X[-remove_inds,]
    message(paste0("Proceeding with ", nrow(X), " individuals after preliminary filtering."))
  }
  
  # prepare for 10-fold cross-validation
  test_inds_idx<-caret::createFolds(y=rowMeans(hom_expr_mat,na.rm=T), k = num_folds)
  if(length(test_inds_idx)<num_folds){
    while(length(test_inds_idx)<num_folds){
      test_inds_idx<-caret::createFolds(y=rowMeans(hom_expr_mat,na.rm=T), k = num_folds)
    }
  }
  test_inds_ids<-lapply(test_inds_idx, function(x) {rownames(X)[x]})
  test_inds<-list()
  for(fold in 1:num_folds){
    for(i in 1:q){
      if(fold ==1){
        test_inds[[i]]<-list()
      }
      test_inds[[i]][[fold]]<-which(rownames(Ys[[i]]) %in% test_inds_ids[[fold]])
    }
  }
  message("Starting cross-validation")
  message("CONTENT temporary file is ", paste0(out_dir,gene_name, "_content_tmp.bk"))
  if(file.exists(paste0(out_dir,gene_name, "_content_tmp.bk"))){
    system(paste0("rm ", paste0(out_dir,gene_name, "_content_tmp.bk")))
  }
  explanatory=as_FBM(X, backingfile=paste0(out_dir,gene_name, "_content_tmp"))
  

  # start cross-validation
  for(cur_fold in 1:num_folds){
    
    train_inds_id<-setdiff(rownames(X), test_inds_ids[[cur_fold]])
    train_inds_id<-intersect(train_inds_id, rownames(hom_expr_mat))
    fold_hom_expr_mat<-hom_expr_mat[train_inds_id,,drop=F]
    
    ## Fit the heterogeneous components
    tiss_betas<-vector("list", q)
    tiss_ints<-vector("list", q)
    
    for(j in 1:q){
      
      tiss_inds=rownames(fold_hom_expr_mat)[which(!is.na(fold_hom_expr_mat[,j]))]
      tiss_response=fold_hom_expr_mat[tiss_inds,j]
      
      tiss_fit<-big_spLinReg(X = explanatory, ind.train = match(tiss_inds, rownames(X)),
                             y.train = tiss_response,K=10, alphas = c(alpha),warn=F)
      
      tiss_beta_vals<-unlist(summary(tiss_fit)$beta[
        which.min(summary(tiss_fit)$validation_loss)])[1:m]
      tiss_tiss_int<-unlist(summary(tiss_fit)$intercept[
        which.min(summary(tiss_fit)$validation_loss)])
      idx_remove<-c()
      if(any(is.na(tiss_beta_vals))){
        idx_remove<-which(is.na(tiss_beta_vals))
        tiss_beta_vals[idx_remove]<-0
      }
      if(length(tiss_beta_vals)<dim(explanatory)[2]){
        new_betas<-rep(0, dim(explanatory)[2])
        new_betas[attr(tiss_fit, "ind.col")]<-tiss_beta_vals
        tiss_beta_vals<-new_betas
      }
      tiss_betas[[j]]<-tiss_beta_vals
      tiss_ints[[j]]<-tiss_tiss_int
      
      
      safe_test_inds<-intersect(rownames(Ys[[j]])[test_inds[[j]][[cur_fold]]], rownames(Ys[[j]]))
      if(length(safe_test_inds) < 1){
        next
      }
      Yhats_tiss[[j]][rownames(Ys[[j]])[test_inds[[j]][[cur_fold]]],]<-X[rownames(Ys[[j]])[test_inds[[j]][[cur_fold]]], ] %*% tiss_betas[[j]] + tiss_ints[[j]]
    }
    

    if(cur_fold < num_folds){
      message("Finished fold: ", cur_fold, " of ", num_folds)
    }
  }
  
  names(Yhats_tiss)<-names(Ys)
  message("removing temporary files")
  system(paste0("rm ", paste0(out_dir,gene_name, "_content_tmp.bk")))
  
  return(list(Yhats_tiss = Yhats_tiss, hom_expr_mat = hom_expr_mat))
  
}

evaluation_helper = function(Ys, hom_expr_mat, Yhats_tiss, contexts_vec, is_GBAT, Yhats_full, out_dir, gene_name){
  message("Calculating r2 performance metrics")
  ## Score each method
  het_cv_pvals<-vector("list", length(contexts_vec)); het_cv_r2s<-vector("list", length(contexts_vec))
  hom_cv_pvals<-vector("list", length(contexts_vec)); hom_cv_r2s<-vector("list", length(contexts_vec))
  het_cv_pvals.herit<-vector("list", length(contexts_vec)); het_cv_r2s.herit<-vector("list", length(contexts_vec))
  full_cv_pvals<-vector("list", length(contexts_vec)); full_cv_r2s<-vector("list", length(contexts_vec))
  tiss_cv_pvals<-vector("list", length(contexts_vec)); tiss_cv_r2s<-vector("list", length(contexts_vec))
  full_weights<-vector("list", length(contexts_vec)); het_scales=vector("list", length(contexts_vec))
  
  for(context in contexts_vec){
    if(context != "AverageContext"){
      if(!is_GBAT){
        index_exp = which(contexts_vec == context)
        index_avg_exp = which(contexts_vec == "AverageContext")
        hom.tmp=Yhats_tiss[[index_exp]]
        # baseline model
        m1<-lm((hom_expr_mat[rownames(Ys[[index_exp]]),index_exp]+hom_expr_mat[rownames(Ys[[index_exp]]),index_avg_exp])~1) ### tests how much an intercept explains total expresison
        # homogeneous model
        m2<-lm((hom_expr_mat[,index_exp]+hom_expr_mat[,index_avg_exp]) ~ Yhats_tiss[[index_avg_exp]]) ### tests how much predicted shared expression of this context explains total expression
        # heterogeneous model
        ## if we learned one for this tissue/context
        ## is the het term heritable:
        hetresponse=hom_expr_mat[rownames(Ys[[index_exp]]),index_exp]
        hetbaseline<-lm(hetresponse ~ 1)
        hetfit=lm(hetresponse ~ Yhats_tiss[[index_exp]]) 
        
        het_test_stat.herit<--2*(logLik(hetbaseline))+2*logLik(hetfit)
        het_cv_pvals.herit[[index_exp]]<-pchisq(het_test_stat.herit,1,lower.tail=F) ## pvalue of the observed specific component of expression to the predicted specific expression
        het_cv_r2s.herit[[index_exp]]<-summary(hetfit)$adj.r.squared
        
        
        m3<-lm((hom_expr_mat[rownames(Ys[[index_exp]]),index_exp]+hom_expr_mat[rownames(Ys[[index_exp]]),index_avg_exp]) ~ Yhats_tiss[[index_exp]][rownames(Ys[[index_exp]]),])
        #het_scales[[index_exp]]=coef(m3)[2]
        # full model
        m4=lm((hom_expr_mat[rownames(Ys[[index_exp]]),index_exp]+hom_expr_mat[rownames(Ys[[index_exp]]),index_avg_exp]) ~ hom_expr_mat[rownames(Ys[[index_exp]]),index_avg_exp] + 
                Yhats_tiss[[index_exp]][rownames(Ys[[index_exp]]),])
        full_values <- predict(m4)
        full_values = data.frame(full_values, check.names = F)
        names(full_values) = "pred"
        Yhats_full[[index_exp]][rownames(full_values), ]<- full_values$pred
        
        # test for signif of full model
        full_test_stat<--2*(logLik(m1))+2*logLik(m4)
        full_cv_pvals[[index_exp]]<-pchisq(full_test_stat,2,lower.tail=F) ## pvalue of the total expression for this context with shared and specific expression
        full_cv_r2s[[index_exp]]<-summary(m4)$adj.r.squared
        full_weights[[index_exp]]<-coef(m4)[2:3]
        # test for signif of het | hom
        het_test_stat<--2*(logLik(m2))+2*logLik(m4) ## pvalue of total against just shared with total against full
        het_cv_pvals[[index_exp]]<-pchisq(het_test_stat,1,lower.tail=F)
        het_cv_r2s[[index_exp]]<-summary(m3)$adj.r.squared
        # test for signif of hom | het
        hom_test_stat<--2*(logLik(m3))+2*logLik(m4) ## pvalue of total against specific and total against shared and specific
        hom_cv_pvals[[index_exp]]=pchisq(hom_test_stat,1,lower.tail=F)
        hom_cv_r2s[[index_exp]]=summary(m2)$adj.r.squared
      }
      # tissue by tissue approach
      if(is_GBAT){
        index_exp = which(contexts_vec == context)
        t1<-lm(hom_expr_mat[rownames(Ys[[index_exp]]),index_exp] ~ Yhats_tiss[[index_exp]][rownames(Ys[[index_exp]]),])
        m1<-lm((hom_expr_mat[rownames(Ys[[index_exp]]),index_exp])~1) ### tests how much an intercept explains total expresison
        tiss_test_stat<--2*(logLik(m1))+2*logLik(t1)
        tiss_cv_pvals[[index_exp]]<-pchisq(tiss_test_stat,1,lower.tail=F)
        tiss_cv_r2s[[index_exp]]<-summary(t1)$adj.r.squared
      }
    }else{
      if(!is_GBAT){
        ## first is the hom term heritable:
        index_exp = which(contexts_vec == context)
        index_avg_exp = which(contexts_vec == "AverageContext")
        hom.tmp=Yhats_tiss[[context]]
        baseline=lm(hom_expr_mat[,index_avg_exp] ~ 1)
        homfit=lm(hom_expr_mat[,index_avg_exp] ~ hom.tmp)
        hom_test_stat.herit<--2*(logLik(baseline))+2*logLik(homfit)
        het_cv_pvals.herit[[index_exp]]=pchisq(hom_test_stat.herit,1,lower.tail=F)
        het_cv_r2s.herit[[index_exp]]=summary(homfit)$adj.r.squared
      }
    }
  }
  if(is_GBAT){
    pvaldf=cbind(tiss_cv_pvals)
    rownames(pvaldf)=names(Ys)
    pvaldf = data.frame(context = rownames(pvaldf), pvaldf, check.names = F)
    r2df=cbind(tiss_cv_r2s)
    rownames(r2df)=names(Ys)
    r2df = data.frame(context = rownames(r2df), r2df, check.names = F)
    fwrite(pvaldf, file = paste0(out_dir, gene_name, ".GBAT.crossval_pvalues.txt"), sep = "\t")
    fwrite(r2df, file = paste0(out_dir, gene_name, ".GBAT.crossval_r2.txt"), sep = "\t")
  }else{
    pvaldf=cbind(het_cv_pvals, het_cv_pvals.herit, hom_cv_pvals, full_cv_pvals)
    rownames(pvaldf)=names(Ys)
    r2df=cbind(het_cv_r2s, het_cv_r2s.herit, hom_cv_r2s, full_cv_r2s)
    rownames(r2df)=names(Ys)
    pvaldf = data.frame(cbind(context = rownames(pvaldf), pvaldf), check.names = F) %>% mutate(across(everything(), ~ map(.x, ~ if (is.null(.x)) NA else .x)))
    r2df = data.frame(cbind(context = rownames(r2df), r2df), check.names = F) %>% mutate(across(everything(), ~ map(.x, ~ if (is.null(.x)) NA else .x)))
    
    fwrite(pvaldf, file = paste0(out_dir, gene_name, ".crocotel.crossval_pvalues.txt"), sep = "\t")
    fwrite(r2df, file = paste0(out_dir, gene_name, ".crocotel.crossval_r2.txt"), sep = "\t")
    
    Yhat_full_mat<-matrix(NA, nrow = nrow(hom_expr_mat), ncol=length(contexts_vec))
    rownames(Yhat_full_mat)<-rownames(hom_expr_mat)
    colnames(Yhat_full_mat)<-contexts_vec
    for(i in 1:length(contexts_vec)){
      context = contexts_vec[i]
      Yhat_full_mat[rownames(Yhats_full[[i]]),i]<-Yhats_full[[i]]
    }
    all_missing<-names(rowMeans(Yhat_full_mat, na.rm = T)[which(is.nan(rowMeans(Yhat_full_mat, na.rm = T)))])
    remove_inds<-which(rownames(Yhat_full_mat) %in% all_missing)
    Yhat_full_mat = data.frame(cbind(id = rownames(Yhat_full_mat), Yhat_full_mat), check.names = F)
    if(length(remove_inds) != 0){
      print("here")
      Yhat_full_mat = Yhat_full_mat[-remove_inds,]
    }
    fwrite(Yhat_full_mat[,!names(Yhat_full_mat) %in% "AverageContext"], file = paste0(out_dir,gene_name,".crocotel.GReX_predictors.txt"), sep = "\t")
  }
  
  message("Done computing evaluation metrics.")
  
}

# function to format Crocotel summary stat ouput file into a format that can be input into treeQTL
format_treeQTL_old = function(input_file, outdir, top_level){
  df = fread(input_file, sep = "\t", data.table = F)
  
  df %>%
    group_by(context) %>%
    arrange(pvalue) %>%
    group_split() %>%
    purrr::walk(function(sub_df) {
      group_name <- unique(sub_df$context)
      if (top_level == "R"){
        sub_df = sub_df %>%
          rename(
            SNP = target,
            gene = regulator,
            t.stat = se,
            p.value = pvalue
          ) 
      }else if(top_level == "T"){
        sub_df = sub_df %>%
          rename(
            SNP = regulator,
            gene = target,
            t.stat = se,
            p.value = pvalue
          ) 
      }else{
        stop("No valid input specified for target or regulator as top level.")
        }
      
      sub_df%>% mutate(FDR = NA) %>% select(SNP, gene, beta, t.stat, p.value, FDR) %>% filter(p.value <= 0.05) %>%
        fwrite(file = paste0(outdir, "all_gene_pairs.", group_name, ".txt"), sep = "\t", na = NA)
    })
  
  df %>%
    group_by(context) %>%
    group_split() %>%
    purrr::walk(function(sub_df) {
      group_name <- unique(sub_df$context)
      
      if (top_level == "R"){
        sub_df = sub_df %>%
          rename(
            SNP = target,
            gene = regulator
          ) 
      }else if(top_level == "T"){
        sub_df = sub_df %>%
          rename(
            SNP = regulator,
            gene = target
          ) 
      }
      
      sub_df %>% group_by(gene) %>% mutate(fam_p = n()) %>% rename(family = gene) %>%
        select(family, fam_p) %>% distinct() %>%
        fwrite(file = paste0(outdir, "n_tests_per_gene.", group_name, ".txt"), sep = ",")
    })
}

format_treeQTL = function(crocotel_dir, top_level, tmp_dir){
  files = list.files(crocotel_dir, full.names = T)
  for(file in files){
    context = sub("\\..*", "", basename(file))
    sub_df = fread(file, sep = "\t", data.table = F)
    if (top_level == "R"){
      sub_df = sub_df %>%
        rename(
          SNP = target,
          gene = regulator,
        ) 
    }else if(top_level == "T"){
      sub_df = sub_df %>%
        rename(
          SNP = regulator,
          gene = target,
        ) 
    }else{
      stop("No valid input specified for target or regulator as top level.")
    }
    sub_df %>% fwrite(file = paste0(tmp_dir, "all_gene_pairs.", context, ".txt"), sep = "\t", na = NA)
    
    sub_df = fread(file, sep = "\t", data.table = F)
    if (top_level == "R"){
      sub_df = sub_df %>%
        rename(
          SNP = target,
          gene = regulator
        ) 
    }else if(top_level == "T"){
      sub_df = sub_df %>%
        rename(
          SNP = regulator,
          gene = target
        ) 
    }
    
    sub_df %>% group_by(gene) %>% mutate(fam_p = n()) %>% rename(family = gene) %>%
      select(family, fam_p) %>% distinct() %>%
      fwrite(file = paste0(outdir, "n_tests_per_gene.", context, ".txt"), sep = ",")
  }

}

# Modified treeQTL function to get eGenes in a multi-context experiment
get_eGenes_multi_tissue_mod = function(crocotel_dir, exp_suffix, out_dir, top_level = "R", level1 = 0.05, level2 = 0.05, level3 = 0.05) {
  
  print(paste("Step 0.1: Computing summary statistics for each context"))
  
  ### set up summary stats per context and number of tests per context
  treeQTL_dir = paste0(out_dir, "/treeQTL_output/")
  dir.create(treeQTL_dir, showWarnings = F)
  tmp_dir = paste0(treeQTL_dir, "/treeQTL_tmp/")
  dir.create(tmp_dir, showWarnings = F)
  
  format_treeQTL(crocotel_dir, top_level, tmp_dir)
  
  crocotel_outfiles = list.files(tmp_dir, pattern = "all_gene_pairs", full.names = T)
  n_SNPs_per_gene_files = list.files(tmp_dir, pattern = "n_tests_per_gene", full.names = T)
  sprintf("Proceeding with %i summary statistic files", length(crocotel_outfiles))
  sprintf("Proceeding with %i tests per gene files", length(n_SNPs_per_gene_files))
  
  
  n_tissue <- length(contexts_vec)
  for (i in 1:n_tissue) {
    cur_tissue_name <- contexts_vec[i]
    
    print(paste("Computing summary statistics for context ", cur_tissue_name, sep = ""))
    n_SNPs_per_gene_this_tissue <- data.frame(fread(input = n_SNPs_per_gene_files[i], header = T), stringsAsFactors = F,check.names = F)
    colnames(n_SNPs_per_gene_this_tissue)=c("family","n_tests")
    n_SNPs_per_gene_this_tissue <- n_SNPs_per_gene_this_tissue[n_SNPs_per_gene_this_tissue$n_tests > 0, ]
    
    gene_simes_cur_tissue <- get_eGenes(n_tests_per_gene = n_SNPs_per_gene_this_tissue, m_eqtl_out = crocotel_outfiles[i], method = "BH", level1 = 1, level2 = 1, silent = TRUE)
    gene_simes_cur_tissue <- merge(gene_simes_cur_tissue, n_SNPs_per_gene_this_tissue, by = "family", all = TRUE)
    gene_simes_cur_tissue$fam_p[which(is.na(gene_simes_cur_tissue$fam_p))] <- 1
    
    if (i == 1) {
      eGene_pvals <- gene_simes_cur_tissue[, c("family", "fam_p")]
      n_SNPs_per_gene_xT <- n_SNPs_per_gene_this_tissue
    } else {
      eGene_pvals <- merge(eGene_pvals, gene_simes_cur_tissue[, c("family", "fam_p")], by = "family", all = TRUE)
      n_SNPs_per_gene_xT <- merge(n_SNPs_per_gene_xT, n_SNPs_per_gene_this_tissue, by = "family", all = TRUE)
    }
    names(eGene_pvals)[i + 1] <- cur_tissue_name
    names(n_SNPs_per_gene_xT)[i + 1] <- cur_tissue_name
  }
  names(eGene_pvals)[1] <- "gene"
  remove(cur_tissue_name, n_SNPs_per_gene_this_tissue, gene_simes_cur_tissue)
  
  print("Step 0.2: Computing summary statistics across contexts")
  col_ind_pvals <- 2:(n_tissue + 1)
  eGene_pvals$simes_p <- apply(eGene_pvals[, col_ind_pvals], 1, TreeQTL:::get_simes_p)
  
  print("Step 1: Selecting eGenes across contexts")
  eGene_xT_qvals <- qvalue(eGene_pvals$simes_p, lambda = 0)$qvalue
  R_G <- sum(eGene_xT_qvals <= level1)
  print(paste("Number of cross-tissue eGenes = ", R_G))
  
  print("Step 2: Selecting contexts in which eGenes are active")
  if(R_G == 0){
    stop("no significant eGenes. Not writing output.")
  }
  q2_adj <- R_G * level2/nrow(eGene_pvals)
  ind_sel_simes <- which(eGene_xT_qvals <= level1)
  sel_eGenes_simes <- eGene_pvals[ind_sel_simes, ]
  rej_simes <- t(1 * apply(sel_eGenes_simes[, c(col_ind_pvals)], 1, TreeQTL:::qsel_by_fam, q2_adj))
  
  print("Step 3: Selecting regulators or targets associated to each gene in each contexts")
  sel_eGenes_simes$n_sel_tissues <- rowSums(rej_simes)
  sel_eGenes_simes$n_tested_tissues <- rowSums(!is.na(sel_eGenes_simes[, col_ind_pvals]))
  
  for (i in 1:n_tissue) {
    cur_tissue_name <- contexts_vec[i]
    print(paste("Selecting regulators or targets for contexts", cur_tissue_name))
    sel_gene_names_this_tissue <- sel_eGenes_simes$gene[which(rej_simes[, i] == 1)]
    sel_gene_info <- n_SNPs_per_gene_xT[which(n_SNPs_per_gene_xT$family %in% sel_gene_names_this_tissue), c(1, i + 1)]
    names(sel_gene_info)[2] <- "n_tests"
    sel_gene_info <- merge(sel_gene_info, sel_eGenes_simes[, c("gene", "n_sel_tissues", "n_tested_tissues")],
                           by.x = "family", by.y = "gene", all.x = TRUE, all.y = FALSE)
    n_sel_per_gene <- TreeQTL:::get_nsel_SNPs_per_gene_tissue_pair(sel_gene_info, cur_tissue_name, crocotel_outfiles[i], R_G, nrow(eGene_pvals),
                                                                   level3 = level3)
    
    print(paste("Total number of associations for context", cur_tissue_name, "=", sum(n_sel_per_gene$n_sel_snp)))
    out_file_name <- paste0(tmp_dir,"/eAssoc_by_gene.", cur_tissue_name, ".", exp_suffix,".txt")
    print(paste("Writing output file", out_file_name))
    if(nrow(n_sel_per_gene) == 0){
      input_df <- n_sel_per_gene
    }else{
      input_df <- data.frame(family = n_sel_per_gene$family, pval = NA, n_sel = n_sel_per_gene$n_sel_snp, check.names = F)
    }
    get_eAssociations(input_df, NULL,
                      crocotel_outfiles[i], out_file_name, by_snp = FALSE, silent = TRUE)
  }
  eGene_xT_sel <- data.frame(gene = sel_eGenes_simes$gene, check.names = F)
  eGene_xT_sel <- cbind(eGene_xT_sel, rej_simes)
  names(eGene_xT_sel)[2:(n_tissue + 1)] <- contexts_vec
  
  
  ## reformat eAssoc files and save to output directory
  eAssoc_files = list.files(tmp_dir, pattern = "eAssoc_by_gene.", full.names = T)
  for(file in eAssoc_files){
    outfile_name = basename(file)
    outfile = paste0(treeQTL_dir, "/", outfile_name)
    sub_df = fread(file, sep = "\t", header = T, data.table = F)
    if (top_level == "R"){
      sub_df = sub_df %>%
        rename(
          SNP = target,
          gene = regulator
        ) 
    }else if(top_level == "T"){
      sub_df = sub_df %>%
        rename(
          SNP = regulator,
          gene = target
        ) 
    }
    fwrite(sub_df, outfile, sep = "\t")
  }
  
  ### remove tmp directory
  unlink(tmp_dir, recursive = TRUE)
  return(eGene_xT_sel)
}


  
  
  