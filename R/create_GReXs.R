
### Test parameters
#X_file = "/Users/lkrockenberger/C-STEM/example_data/genos"
#exp_files = list.files("/Users/lkrockenberger/C-STEM/example_data/expression/")
#contexts = exp_files
#exp_files = paste0("/Users/lkrockenberger/C-STEM/example_data/expression/", exp_files)
#out_dir = "/Users/lkrockenberger/C-STEM/example_data/GReXs/"
#snps = "/Users/lkrockenberger/C-STEM/example_data/snps.txt"
#decomposition_dir = "/Users/lkrockenberger/C-STEM/example_data/decomposed_expression/"
#gene_name = "gene1"
#context_thresh = 3
#alpha = 0.5
#num_folds = 10
#run_GBAT = TRUE

#library(data.table)
#library(bigstatsr)
#library(dplyr)
#library(caret)

#source("R/decompose_expression.R")
#source("R/functions.R")


#' Creates GReXs (Genetically regulated predictors of expression) for one gene across contexts
#'
#' @param X_file - a genotype file with the set of individuals in all Y_files and their corresponding cis-SNPs. Individuals are row names, no column names
#' @param exp_files - vector with list of all expression files per context. The individuals (row names) must be a subset of genotype (X_file). Contains an unnamed
#' @param contexts - vector of strings which have context names in them.
#' @param out_dir - output directory for GReXs
#' @param gene_name - identifier for current gene being run. Add this prefix to all the saved results files... necessary to distinguish results where more than one gene-analysis is run.
#' @param decomposition_dir - directory to store decomposed expression files
#' @param context_thresh - minimum number of contexts to run C-STEM on. 
#' @param alpha - The regularization constant. Default is .5 (eNet). Minimum value of 1e-4.
#' @param num_folds - Number of folds for cross-validation
#' @param run_GBAT - Takes values of TRUE or FALSE. Default is FALSE, if TRUE then will run the GBAT* method too.
#' @return writes out a file of predicted expression across individuals and contexts 
#' @export
create_GReXs = function(X_file, exp_files, contexts, out_dir, gene_name, decomposition_dir, context_thresh = 3, alpha = 0.5, num_folds = 10, run_GBAT = FALSE){
  seed = 9000
  set.seed(seed)
  message("Saving cross-validated predictors and performance metrics in ", out_dir)

  message("Reading in files...") 
  suppressWarnings(expr = {X<-fread(file = X_file, sep='\t', data.table=F)})
  X<-as.matrix(data.frame(X, row.names=1, check.names = F))
  
  ###### read in expression and decompose - files are written out to decomposed exp directory
  decompose_expression(exp_files, gene_name, contexts, context_thresh, decomposition_dir)
  ## this assumes that the file name of decomposition dir is saved as "gene.context.etc" and shared is called "Average Context" (output of decompose function)
  contexts_vec = sapply(strsplit(list.files(decomposition_dir), "\\."), "[[", 2)
  
  ### read in expression
  Ys<-vector("list", length = length(list.files(decomposition_dir)))
  names(Ys)<-contexts_vec
  lengths_y = c()
  rownames_y = list()
  Yhats_full<-vector("list", length(Ys))
  for(i in 1:length(Ys)){
    suppressWarnings(expr={ Ys[[i]]<-fread(file = 
                                             paste0(decomposition_dir,list.files(decomposition_dir)[i]), sep='\t', data.table=F)})
    Ys[[i]]<-as.matrix(data.frame(Ys[[i]], row.names=1, check.names = F))
    lengths_y = c(lengths_y, nrow(Ys[[i]]))
    rownames_y[[i]] = rownames(Ys[[i]])
    if(any(is.na(Ys[[i]])) | any(is.nan(Ys[[i]]))){
      remove=unique(c( which(is.na(Ys[[i]])), which(is.nan(Ys[[i]])) ))
      Ys[[i]]=Ys[[i]][-remove,,drop=F]
    }
    
    Yhats_full[[i]]<-matrix(NA, ncol=1, nrow=lengths_y[i], 
                            dimnames = list(rownames_y[[i]], "pred"))
  }
  
  
  crossval_output = crossval_helper(Ys, X, lengths_y, rownames_y, contexts_vec, out_dir, gene_name, num_folds, alpha)
  Yhats_tiss = crossval_output[["Yhats_tiss"]]
  hom_expr_mat = crossval_output[["hom_expr_mat"]]
  
  ## combine Yhats_tiss into large dataframe with individuals as rows and contexts as columns:
  Yhat_tiss_mat<-matrix(NA, nrow = nrow(X), ncol=length(contexts_vec))
  rownames(Yhat_tiss_mat)<-rownames(X)
  colnames(Yhat_tiss_mat)<-contexts_vec
  for(context in contexts_vec){
    Yhat_tiss_mat[rownames(Yhats_tiss[[context]]),context]<-Yhats_tiss[[context]]
  }
  all_missing<-names(rowMeans(Yhat_tiss_mat, na.rm = T)[which(is.nan(rowMeans(Yhat_tiss_mat, na.rm = T)))])
  remove_inds<-which(rownames(Yhat_tiss_mat) %in% all_missing)
  Yhat_tiss_mat = data.frame(Yhat_tiss_mat[-remove_inds,])
  Yhat_tiss_mat = cbind(id = rownames(Yhat_tiss_mat), Yhat_tiss_mat)
  fwrite(Yhat_tiss_mat, file = paste0(out_dir,gene_name,"_cstem_predictors.txt"), sep = "\t")
  evaluation_helper(Ys, hom_expr_mat, Yhats_tiss, contexts_vec, FALSE, Yhats_full, out_dir, gene_name)
  
  ### read in expression for gbat
  if(run_GBAT){
    Ys_gbat<-vector("list", length = length(exp_files))
    names(Ys_gbat)<-sapply(strsplit(exp_files, "/"), tail, 1)
    lengths_y_gbat = c()
    rownames_y_gbat = list()
    for(i in 1:length(Ys_gbat)){
      suppressWarnings(expr={ Ys_gbat[[i]]<-fread(file = 
                                                   exp_files[i], sep='\t', data.table=F)})
      Ys_gbat[[i]]<-as.matrix(data.frame(Ys_gbat[[i]], row.names=1, check.names = F))
      lengths_y_gbat = c(lengths_y_gbat, nrow(Ys_gbat[[i]]))
      rownames_y_gbat[[i]] = rownames(Ys_gbat[[i]])
      if(any(is.na(Ys_gbat[[i]])) | any(is.nan(Ys_gbat[[i]]))){
        remove=unique(c( which(is.na(Ys_gbat[[i]])), which(is.nan(Ys_gbat[[i]])) ))
        Ys_gbat[[i]]=Ys_gbat[[i]][-remove,,drop=F]
      }
    }
    gbat_contexts_vec = contexts_vec[!grepl("AverageContext", contexts_vec)]
    output_gbat = crossval_helper(Ys_gbat, X, lengths_y_gbat, rownames_y_gbat, gbat_contexts_vec, out_dir, gene_name, num_folds, alpha)
    Yhats_gbat = output_gbat[["Yhats_tiss"]]
    hom_expr_mat_gbat = output_gbat[["hom_expr_mat"]]
    
    Yhat_gbat_mat<-matrix(NA, nrow = nrow(X), ncol=length(gbat_contexts_vec))
    rownames(Yhat_gbat_mat)<-rownames(X)
    colnames(Yhat_gbat_mat)<-gbat_contexts_vec
    for(context in gbat_contexts_vec){
      Yhat_gbat_mat[rownames(Yhats_gbat[[context]]),context]<-Yhats_gbat[[context]]
    }
    all_missing<-names(rowMeans(Yhat_gbat_mat, na.rm = T)[which(is.nan(rowMeans(Yhat_gbat_mat, na.rm = T)))])
    remove_inds<-which(rownames(Yhat_gbat_mat) %in% all_missing)
    Yhat_gbat_mat = data.frame(Yhat_gbat_mat[-remove_inds,])
    Yhat_gbat_mat = cbind(id = rownames(Yhat_gbat_mat), Yhat_gbat_mat)
    fwrite(Yhat_gbat_mat, file = paste0(out_dir,gene_name,"_gbat_predictors.txt"), sep = "\t")
    
    evaluation_helper(Ys_gbat, hom_expr_mat_gbat, Yhats_gbat, gbat_contexts_vec, TRUE, NULL, out_dir, gene_name)
  }
  message("Finished computing GReXs and evalutation metrics")
  
}
