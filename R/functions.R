
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
        full_values = data.frame(full_values)
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
    pvaldf = data.frame(context = rownames(pvaldf), pvaldf)
    r2df=cbind(tiss_cv_r2s)
    rownames(r2df)=names(Ys)
    r2df = data.frame(context = rownames(r2df), r2df)
    fwrite(pvaldf, file = paste0(out_dir, gene_name, "_GBAT_crossval_pvalues.txt"), sep = "\t")
    fwrite(r2df, file = paste0(out_dir, gene_name, "_GBAT_crossval_r2.txt"), sep = "\t")
  }else{
    pvaldf=cbind(het_cv_pvals, het_cv_pvals.herit, hom_cv_pvals, full_cv_pvals)
    rownames(pvaldf)=names(Ys)
    r2df=cbind(het_cv_r2s, het_cv_r2s.herit, hom_cv_r2s, full_cv_r2s)
    rownames(r2df)=names(Ys)
    pvaldf = data.frame(cbind(context = rownames(pvaldf), pvaldf)) %>% mutate(across(everything(), ~ map(.x, ~ if (is.null(.x)) NA else .x)))
    r2df = data.frame(cbind(context = rownames(r2df), r2df)) %>% mutate(across(everything(), ~ map(.x, ~ if (is.null(.x)) NA else .x)))
    
    fwrite(pvaldf, file = paste0(out_dir, gene_name, "_cstem_crossval_pvalues.txt"), sep = "\t")
    fwrite(r2df, file = paste0(out_dir, gene_name, "_cstem_crossval_r2.txt"), sep = "\t")
    
    Yhat_full_mat<-matrix(NA, nrow = nrow(hom_expr_mat), ncol=length(contexts_vec))
    rownames(Yhat_full_mat)<-rownames(hom_expr_mat)
    colnames(Yhat_full_mat)<-contexts_vec
    for(i in 1:length(contexts_vec)){
      context = contexts_vec[i]
      Yhat_full_mat[rownames(Yhats_full[[i]]),i]<-Yhats_full[[i]]
    }
    all_missing<-names(rowMeans(Yhat_full_mat, na.rm = T)[which(is.nan(rowMeans(Yhat_full_mat, na.rm = T)))])
    remove_inds<-which(rownames(Yhat_full_mat) %in% all_missing)
    Yhat_full_mat = data.frame(cbind(id = rownames(Yhat_full_mat), Yhat_full_mat))
    Yhat_full_mat = Yhat_full_mat[-remove_inds,]
    fwrite(Yhat_full_mat, file = paste0(out_dir,gene_name,"_cstem_full_predictors.txt"), sep = "\t")
  }
  
  message("Done computing evaluation metrics.")
  
}
  
  
  
  
  
  
  
  