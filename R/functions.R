
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
  
  return(Yhats_tiss)
  
}