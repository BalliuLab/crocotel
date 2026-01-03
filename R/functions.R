crossval_helper_parallel = function(total_exp_mat, decomp_exp_mat, X,
                                    GReX_outdir, gene_name, is_cxc, num_folds, alpha = 0.5,
                                    n_cores = 4) {
  
  ### create a working matrix based on the decomposed (regular expression if cxc) to use in this function
  context_exp_mat = decomp_exp_mat
  rownames(context_exp_mat) = decomp_exp_mat$id
  context_exp_mat = context_exp_mat %>% dplyr::select(-id)
  ## set up some needed data frames
  
  rownames(total_exp_mat) = total_exp_mat$id
  
  # Drop samples missing in all tissues
  all_missing <- rownames(context_exp_mat)[rowSums(!is.na(context_exp_mat)) == 0]
  if (length(all_missing) > 0) {
    context_exp_mat <- context_exp_mat[!rownames(context_exp_mat) %in% all_missing, , drop = FALSE]
  }
  X <- X[rownames(X) %in% rownames(context_exp_mat), , drop = FALSE]
  message("Proceeding with ", nrow(X), " individuals after filtering.")
  
  # Cross-validation setup
  set.seed(123)  # for reproducibility
  test_inds_idx <- caret::createFolds(rowMeans(context_exp_mat, na.rm = TRUE), k = num_folds)
  test_inds_ids <- lapply(test_inds_idx, function(x) rownames(X)[x])
  
  ## get person IDs that are not NA for the test set in each context
  test_ids <- setNames(lapply(colnames(context_exp_mat), function(context) {
    lapply(test_inds_ids, function(fold_ids) rownames(context_exp_mat)[
      !is.na(context_exp_mat[, context]) & rownames(context_exp_mat) %in% fold_ids
    ])
  }), colnames(context_exp_mat))
  
  # Prepare parallel
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  Yhats_list <- foreach(cur_fold = 1:num_folds, .packages = c("bigstatsr")) %dopar% {
    
    train_inds_id <- setdiff(rownames(X), test_inds_ids[[cur_fold]])
    
    fold_exp_mat <- context_exp_mat[train_inds_id, , drop = FALSE]
    
    ### create a matrix of mxc to store estimated betas for this fold
    beta_hat_mat = data.frame(matrix(NA, ncol = ncol(context_exp_mat), nrow = ncol(X)))
    colnames(beta_hat_mat) = colnames(context_exp_mat)
    intercept_mat = data.frame(matrix(NA, ncol = ncol(context_exp_mat), nrow = 1))
    colnames(intercept_mat) = colnames(context_exp_mat)
    
    # Per-fold FBM (avoid file collisions)
    fold_bkfile <- paste0(GReX_outdir, gene_name, "_crocotel_tmp_fold", cur_fold)
    if(file.exists(paste0(fold_bkfile, ".bk"))){
      file.remove(paste0(fold_bkfile, ".bk"))
    }
    explanatory <- as_FBM(X, backingfile = fold_bkfile)
    
    fold_Yhats <- data.frame(matrix(NA,nrow = nrow(context_exp_mat), ncol = ncol(context_exp_mat),
                                    dimnames = list(rownames(context_exp_mat), colnames(context_exp_mat))), check.names = F)
    
    ## fit each specific component per context
    for (context in colnames(context_exp_mat)) {
      context_ids <- rownames(fold_exp_mat)[!is.na(fold_exp_mat[,context])]
      context_fit <- big_spLinReg(X = explanatory,
                                  ind.train = match(context_ids, rownames(X)),
                                  y.train = fold_exp_mat[context_ids, context],
                                  K = 10, alphas = alpha, warn = FALSE)
      
      best_idx <- which.min(summary(context_fit)$validation_loss)
      beta_vec <- unlist(summary(context_fit)$beta[best_idx])
      intercept <- summary(context_fit)$intercept[[best_idx]]
      
      # Handle NA beta 
      beta_vec[is.na(beta_vec)] <- 0
      if (length(beta_vec) < ncol(X)) {
        full_beta <- rep(0, ncol(X))
        full_beta[attr(tiss_fit, "ind.col")] <- beta_vec
        beta_vec <- full_beta
      }
      
      
      beta_hat_mat[,context] = beta_vec
      intercept_mat[1,context] = intercept
    }
    
    ### predict on the train for each context
    if(is_cxc){
      contexts_prediction = colnames(context_exp_mat)
    }else{
      contexts_prediction = setdiff(colnames(context_exp_mat), "shared")
    }
    for(context in contexts_prediction){
      ## get yhat specific 
      specific_preds = X[train_inds_id,]%*% beta_hat_mat[,context] + intercept_mat[1,context]
      
      ### if we are running crocotel, get yhat shared for this context
      if(!is_cxc){
        shared_preds = X[train_inds_id, ] %*% beta_hat_mat[,"shared"] + intercept_mat[1,"shared"]
      }
      
      ## get test indices that have data for this context
      safe_test_ids = intersect(test_ids[[context]][[cur_fold]], rownames(X))
      
      ### get Yhat specific on the test set
      specific_preds_test = X[safe_test_ids, ] %*% beta_hat_mat[,context] + intercept_mat[1,context]
      fold_Yhats[safe_test_ids, context] <- specific_preds_test
      
      ### if running crocotel, get Yhat shared on the test set and Yhat full
      if(!is_cxc){
        shared_preds_test = X[safe_test_ids, ] %*% beta_hat_mat[,"shared"] + intercept_mat[1,"shared"]
        fold_Yhats[safe_test_ids, "shared"] = shared_preds_test
      }
    }
    
    file.remove(paste0(fold_bkfile, ".bk"))  # Clean up
    return_list = list(fold_Yhats)
    return(return_list)
  }
  stopCluster(cl)
  
  # Combine predictions across parallelized folds
  context_predictions_list = lapply(Yhats_list, `[[`, 1) 
  combined_context_preds = Reduce(function(a, b) dplyr::coalesce(a, b), context_predictions_list)
  combined_context_preds = cbind(id = rownames(combined_context_preds), combined_context_preds)
  
  if(is_cxc){
    method = "cxc"
    fwrite(combined_context_preds, file = paste0(GReX_outdir, gene_name,".", method, ".GReX_predictors.txt"), sep = "\t")
    full = NA
    ## get r2
    evaluation_helper(total_exp_mat, combined_context_preds, GReX_outdir, gene_name, method = method)
  }else{
    shared_col = combined_context_preds[,"shared"]
    ## get crocotel_added
    added = as.data.frame(
      lapply(combined_context_preds[,(2:(ncol(combined_context_preds)-1))], function(col) col + shared_col), 
      check.names = F)
    added = cbind(id = rownames(combined_context_preds), added)
    ## get added r2
    method = "crocotel_added"
    evaluation_helper(total_exp_mat, added, GReX_outdir, gene_name, method = method, combined_context_preds)
    fwrite(added, file = paste0(GReX_outdir, gene_name,".", method, ".GReX_predictors.txt"), sep = "\t")
    ## get crocotel_full
    full = Reduce(
      function(x, y) merge(x, y, by = "rownames", all = TRUE),
        lapply(colnames(combined_context_preds[,(2:(ncol(combined_context_preds)-1))]), function(col){
        model = lm(total_exp_mat[rownames(combined_context_preds),col]~combined_context_preds[,col] + shared_col)
        preds = predict(model)
        df = data.frame(
          rownames = combined_context_preds[!is.na(combined_context_preds[,col]), "id"],
          col = preds,
          check.names = FALSE
        )
        names(df) = c("rownames", col)
        df
      }))
    names(full) = colnames(combined_context_preds[,(2:(ncol(combined_context_preds)-1))])
    full = cbind(id = rownames(combined_context_preds), full)
    ## get full r2
    method = "crocotel"
    evaluation_helper(total_exp_mat, full, GReX_outdir, gene_name, method = method, combined_context_preds)
    fwrite(full, file = paste0(GReX_outdir, gene_name,".", method, ".GReX_predictors.txt"), sep = "\t")
    fwrite(combined_context_preds, file = paste0(GReX_outdir, gene_name,".", method, "_context.GReX_predictors.txt"), sep = "\t")
  }
  
  return(list(full = full, context = combined_context_preds))
}


evaluation_helper = function(total_exp_mat, combined_full_preds, out_dir, gene_name, method, combined_context_preds=NULL){
  message("Calculating r2 performance metrics")
  
  contexts_vec = setdiff(colnames(total_exp_mat), "id")
  full_cv_pvals<-vector("list", length(contexts_vec)); full_cv_r2s<-vector("list", length(contexts_vec))
  names(full_cv_pvals) = contexts_vec; names(full_cv_r2s) = contexts_vec
  shared_cv_pvals<-vector("list", length(contexts_vec)); shared_cv_r2s<-vector("list", length(contexts_vec))
  names(shared_cv_pvals) = contexts_vec; names(shared_cv_r2s) = contexts_vec
  specific_cv_pvals<-vector("list", length(contexts_vec)); specific_cv_r2s<-vector("list", length(contexts_vec))
  names(specific_cv_pvals) = contexts_vec; names(specific_cv_r2s) = contexts_vec
  
  for(context in contexts_vec){
    m1<-lm(total_exp_mat[,context]~1) ### tests how much an intercept explains total expresison
    
    # full model
    m4=lm(total_exp_mat[combined_full_preds$id,context] ~ combined_full_preds[,context])
    
    if(method != "cxc"){
      #shared model
      m2 = lm(total_exp_mat[combined_context_preds$id, context] ~ combined_context_preds[,"shared"])
      #specific model
      m3 = lm(total_exp_mat[combined_context_preds$id, context] ~ combined_context_preds[, context])

      # test for signif of both models
      shared_test_stat<--2*(logLik(m1))+2*logLik(m2)
      shared_cv_pvals[[context]]<-pchisq(shared_test_stat,2,lower.tail=F) ## pvalue of the total expression for this context with shared  expression
      shared_cv_r2s[[context]]<-summary(m2)$adj.r.squared
      
      specific_test_stat<--2*(logLik(m1))+2*logLik(m3)
      specific_cv_pvals[[context]]<-pchisq(specific_test_stat,2,lower.tail=F) ## pvalue of the total expression for this context with specific expression
      specific_cv_r2s[[context]]<-summary(m3)$adj.r.squared
    }
    
    # test for signif of full model
    full_test_stat<--2*(logLik(m1))+2*logLik(m4)
    full_cv_pvals[[context]]<-pchisq(full_test_stat,2,lower.tail=F) ## pvalue of the total expression for this context with shared and specific expression
    full_cv_r2s[[context]]<-summary(m4)$adj.r.squared
    
  }
  
  full_cv_pvals     <- unlist(full_cv_pvals)
  specific_cv_pvals <- unlist(specific_cv_pvals)
  shared_cv_pvals   <- unlist(shared_cv_pvals)
  
  full_cv_r2s       <- unlist(full_cv_r2s)
  specific_cv_r2s   <- unlist(specific_cv_r2s)
  shared_cv_r2s     <- unlist(shared_cv_r2s)
  
  pvaldf=cbind(full_cv_pvals, specific_cv_pvals, shared_cv_pvals)
  r2df=cbind(full_cv_r2s, specific_cv_r2s, shared_cv_r2s)
  
  
  if(method != "cxc"){
    pvaldf <- data.frame(
      context = contexts_vec,
      full = full_cv_pvals,
      specific = specific_cv_pvals,
      shared = shared_cv_pvals,
      check.names = FALSE
    )
    
    r2df <- data.frame(
      context = contexts_vec,
      full = full_cv_r2s,
      specific = specific_cv_r2s,
      shared = shared_cv_r2s,
      check.names = FALSE
    )
  }else{
    pvaldf <- data.frame(
      context = contexts_vec,
      full = full_cv_pvals,
      specific = NA,
      shared = NA,
      check.names = FALSE
    )
    
    r2df <- data.frame(
      context = contexts_vec,
      full = full_cv_r2s,
      specific = NA,
      shared = NA,
      check.names = FALSE
    )
  }

  fwrite(pvaldf, file = paste0(out_dir, gene_name, ".", method, ".crossval_pvalues.txt"), sep = "\t")
  fwrite(r2df, file = paste0(out_dir, gene_name, ".", method, ".crossval_r2.txt"), sep = "\t")
  
  message("Done computing evaluation metrics.")
  
}

# Modified treeQTL function to get eGenes in a multi-context experiment
get_eGenes_multi_tissue_mod = function(crocotel_dir, exp_suffix, out_dir, top_level = "R", level1 = 0.05, level2 = 0.05, level3 = 0.05) {
  
  print(paste("Step 0.1: Computing summary statistics for each context"))
  
  ### set up summary stats per context and number of tests per context
  #tmp_dir = paste0(out_dir, "/treeQTL_tmp/")
  tmp_dir = paste0(tempfile(tmpdir = out_dir), "/")
  dir.create(tmp_dir, showWarnings = F)
  
  format_treeQTL(crocotel_dir, top_level, tmp_dir)
  
  crocotel_outfiles = list.files(tmp_dir, pattern = "all_gene_pairs", full.names = T)
  n_SNPs_per_gene_files = list.files(tmp_dir, pattern = "n_tests_per_gene", full.names = T)
  sprintf("Proceeding with %i summary statistic files", length(crocotel_outfiles))
  sprintf("Proceeding with %i tests per gene files", length(n_SNPs_per_gene_files))
  
  contexts_vec = sub("\\..*", "", basename(list.files(crocotel_dir)))
  print(paste("inferred contexts:", paste(contexts_vec, collapse = ",")))
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
    print("no significant eGenes. writing empty eGene file.")
    empty_df <- as.data.frame(matrix("", nrow = 1, ncol = length(contexts_vec)))
    colnames(empty_df) <- contexts_vec
    empty_df$gene <- ""
    empty_df <- empty_df[, c("gene", contexts_vec)] 
    return(empty_df)
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
    outfile = paste0(out_dir, "/", outfile_name)
    sub_df = fread(file, sep = " ", header = T, data.table = F)
    if (top_level == "R"){
      sub_df = sub_df %>%
        rename(
          target = SNP,
          regulator = gene
        ) 
    }else if(top_level == "T"){
      sub_df = sub_df %>%
        rename(
          regulator = SNP,
          target = gene
        ) 
    }
    fwrite(sub_df, outfile, sep = "\t")
  }
  
  ### remove tmp directory
  unlink(tmp_dir, recursive = TRUE)
  return(eGene_xT_sel)
}



###### fastgxc functions

decompose=function(X,design){
  X = as.matrix(X)
  rep.measures = factor(design)
  if (any(summary(as.factor(rep.measures)) == 1)) 
    stop("A multilevel analysis can not be performed when at least one some sample is not repeated.")
  indiv.names = rownames(X)
  rownames(X) = as.character(rep.measures)
  X.mean.indiv = matrix(apply(X, 2, tapply, rep.measures, mean, na.rm = TRUE), 
                        nrow = length(unique(rep.measures)), 
                        ncol = dim(X)[2], 
                        dimnames = list(levels(as.factor(rep.measures)), colnames(X)))
  Xb = X.mean.indiv[as.character(rep.measures), ]
  Xw = X - Xb
  dimnames(Xw) = list(indiv.names, colnames(X))
  return(list(Xw=Xw,Xb=X.mean.indiv))
}

#' Decomposition Step
#'
#' Function to decompose expression into one shared component and specific components per context to run crocotel direct
#'
#' @param  exp_files - list of expression files for each context
#' @return returns the shared component of expression per individual and each specific expression component for each of the C contexts
#'
decomposition_step = function(exp_files){

  context_names = c()
  
  ### merge the expression files
  exp_all=data.frame(fread(input = exp_files[1], header = T), check.names = F,stringsAsFactors = F)
  cur_context = gsub("_expression.txt", "", basename(exp_files[1]))
  names(exp_all)[-1] = paste(names(exp_all)[-1],cur_context, sep = " - ")
  context_names = c(context_names, cur_context)
  print(paste("Finished merging context",1))
  
  for(i in 2:length(exp_files)){
    
    # Read expression matrix for tissue t
    cur_context = gsub("_expression.txt", "", basename(exp_files[i]))
    exp_t=data.frame(fread(input = exp_files[i], header = T), check.names = F,stringsAsFactors = F)
    colnames(exp_t)[-1] = paste(colnames(exp_t)[-1],cur_context, sep = " - ")
    context_names = c(context_names, cur_context)
    
    # Merge with other tissues
    exp_all = merge(x = exp_all, y = exp_t, by="id", all = TRUE)
    
    print(paste("Finished merging context",i))
  }
  
  # Transpose merged expression matrix to have genes in the columns 
  exp_mat = t(as.matrix(exp_all[, -1]))
  colnames(exp_mat) = exp_all[,1]
  exp_mat = cbind(id = rownames(exp_mat), data.frame(exp_mat, check.names = F))
  print("Finished transposing merged file")
  
  #%%%%%%%%%%%%%%% Sample and context names
  design=sapply(1:nrow(exp_mat), function(i) unlist(strsplit(exp_mat[,1][i], split = " - "))[1])
  context_names=sapply(1:nrow(exp_mat), function(i) unlist(strsplit(exp_mat[,1][i], split = " - "))[2])
  contexts=unique(context_names)
  
  #%%%%%%%%%%%%%%% Decompose expression into homogeneous and heterogeneous context expression
  print("Decomposing data")
  rownames(exp_mat) = exp_mat[,1]
  exp_mat = exp_mat[,-1]
  
  #%%%%%%%%%%%%%%% Print number of genes and samples
  string1 = sprintf("There are %s samples and %s genes. The max number of missing samples for a gene is  %s. The max number of missing genes for a sample is  %s. \n", nrow(exp_mat), ncol(exp_mat),max(colSums(is.na(exp_mat))),max(rowSums(is.na(exp_mat))))
  cat(string1)
  
  dec_exp_all=decompose(X = exp_mat, design = design)
  bexp_all=dec_exp_all$Xb
  wexp_all=dec_exp_all$Xw
  bexp_all[is.nan(bexp_all)]=NA
  wexp_all[is.nan(wexp_all)]=NA
  
  string2 = sprintf("Between individual matrix: There are %s individuals and %s genes. The max number of missing samples for a gene is  %s. The max number of missing genes for a sample is  %s. \n", nrow(bexp_all), ncol(bexp_all),max(colSums(is.na(bexp_all))),max(rowSums(is.na(bexp_all))))
  cat(string2)
  
  string3 = sprintf("Within individual matrix: There are %s samples and %s genes. The max number of missing samples for a gene is  %s. The max number of missing genes for a sample is  %s. \n", nrow(wexp_all), ncol(wexp_all),max(colSums(is.na(wexp_all))),max(rowSums(is.na(wexp_all))))
  cat(string3)
  
  #%%%%%%%%%%%%%%% Save decomposed expression files 
  print("Finished decomposition")
  
  exp_list = list()
  shared_exp = data.table::data.table(t(bexp_all),keep.rownames = T) %>% {setnames(., old = "rn", new = "geneID")[]}
  exp_list[["shared"]] = shared_exp

  for(i in 1:length(contexts)){
    print(contexts[i])
    wexp_t = wexp_all[grep(pattern = paste0(contexts[i],"$"), rownames(wexp_all)),]
    rownames(wexp_t)=gsub(pattern = paste0(" - ",contexts[i]), replacement = "", x = rownames(wexp_t))
    cur_sp_exp = data.table::data.table(t(wexp_t),keep.rownames = T) %>% {setnames(., old = "rn", new = "geneID")[]}
    exp_list[[contexts[i]]] = cur_sp_exp
  }
  
  return(exp_list)
  
}

