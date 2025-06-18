

regress_target_GReX = function(exp_files, Yhats_tiss, outdir){
  dir.create(paste0(outdir, "exp_residualized_GReX/"))
  all_residuals = NULL
  for(file in exp_files){
    context = gsub(".*/", "", file)
    context = sub("\\..*", "", context)
    original_exp = data.frame(fread(file, header = F), check.names = F,stringsAsFactors = F)
    ids <- original_exp[, 1]
    expr <- original_exp[, 2]
    
    cur_grex = Yhats_tiss[[context]]
    
    # Make sure IDs match — reorder or subset if needed
    matching_ids <- intersect(ids, rownames(cur_grex))
    expr_sub <- expr[ids %in% matching_ids]
    grex_sub <- cur_grex[matching_ids,]
    ids_sub <- ids[ids %in% matching_ids]
    
    ## get GReX expression of context
    
    model <- lm(expr_sub ~ grex_sub)
    residuals_vec <- residuals(model)
    
    sub_df = data.frame(id = ids_sub, residual = residuals_vec)
    names(sub_df)[2] = context
    
    # Merge into final residuals data frame
    if (is.null(all_residuals)) {
      all_residuals = sub_df
    } else {
      all_residuals = merge(all_residuals, sub_df, by = "id", all = TRUE)
    }
  }
  outfile = paste0(outdir, "exp_residualized_GReX/", context, ".crocotel.GReX_residuals.txt")
  fwrite(all_residuals, outfile, sep = "\t", quote = F)
}


#' Creates GReXs (Genetically regulated predictors of expression) for one gene across contexts
#'
#' @param genotype_file - a genotype file with the set of individuals in all Y_files and their corresponding cis-SNPs. Individuals are row names, no column names
#' @param exp_files - vector with list of all expression files per context. The individuals (row names) must be a subset of genotype (genotype_file). Contains an unnamed
#' @param out_dir - output directory for GReXs
#' @param gene_name - identifier for current gene being run. Add this prefix to all the saved results files... necessary to distinguish results where more than one gene-analysis is run.
#' @param context_thresh - minimum number of contexts to run crocotel on. 
#' @param alpha - The regularization constant. Default is .5 (eNet). Minimum value of 1e-4.
#' @param num_folds - Number of folds for cross-validation
#' @param method - Takes values of "crocotel" or "cxc". Default is "crocotel", if "cxc" then will not run decomposition and will buil GReXs in a context by context manner.
#' @return writes out a file of predicted expression across individuals and contexts 
#' @export
create_GReXs = function(genotype_file, exp_files, out_dir, gene_name, context_thresh = 3, alpha = 0.5, num_folds = 10, method = "crocotel"){
  seed = 9000
  set.seed(seed)
  dir.create(paste0(outdir, "GReXs"), showWarnings = F)
  decomposition_dir = paste0(out_dir, gene_name, "_decomposed/")
  dir.create(decomposition_dir, showWarnings = F)
  message("Saving decomposed expression in ",  decomposition_dir)
  message("Saving cross-validated predictors and performance metrics in ", out_dir)
  
  message("Reading in files...") 
  suppressWarnings(expr = {X<-fread(file = genotype_file, sep='\t', data.table=F, check.names = F)})
  X<-as.matrix(data.frame(X, row.names=1, check.names = F))
  
  ###### read in expression and decompose - files are written out to decomposed exp directory
  decompose_expression(exp_files, gene_name, context_thresh, decomposition_dir)
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
                                             paste0(decomposition_dir,list.files(decomposition_dir)[i]), sep='\t', data.table=F, check.names = F)})
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
  
  unlink(decomposition_dir, recursive = TRUE)
  crossval_output = crossval_helper(Ys, X, lengths_y, rownames_y, contexts_vec, out_dir, gene_name, num_folds, alpha)
  Yhats_tiss = crossval_output[["Yhats_tiss"]]
  hom_expr_mat = crossval_output[["hom_expr_mat"]]
  ### residualise GReXs from target expression
  regress_target_GReX(exp_files, Yhats_tiss, outdir)
  
  ## combine Yhats_tiss into large dataframe with individuals as rows and contexts as columns:
  Yhat_tiss_mat<-matrix(NA, nrow = nrow(X), ncol=length(contexts_vec))
  rownames(Yhat_tiss_mat)<-rownames(X)
  colnames(Yhat_tiss_mat)<-contexts_vec
  for(context in contexts_vec){
    Yhat_tiss_mat[rownames(Yhats_tiss[[context]]),context]<-Yhats_tiss[[context]]
  }
  all_missing<-names(rowMeans(Yhat_tiss_mat, na.rm = T)[which(is.nan(rowMeans(Yhat_tiss_mat, na.rm = T)))])
  remove_inds<-which(rownames(Yhat_tiss_mat) %in% all_missing)
  if(length(remove_inds) != 0){
    Yhat_tiss_mat = data.frame(Yhat_tiss_mat[-remove_inds,], check.names = F)
  }
  Yhat_tiss_mat = data.frame(cbind(id = rownames(Yhat_tiss_mat), Yhat_tiss_mat))
  #fwrite(Yhat_tiss_mat, file = paste0(out_dir,gene_name,".crocotel_predictors.txt"), sep = "\t")
  evaluation_helper(Ys, hom_expr_mat, Yhats_tiss, contexts_vec, FALSE, Yhats_full, out_dir, gene_name)
  
  ### read in expression for gbat
  if(method == "cxc"){
    Ys_gbat<-vector("list", length = length(exp_files))
    names(Ys_gbat)<-sapply(strsplit(exp_files, "/"), tail, 1)
    lengths_y_gbat = c()
    rownames_y_gbat = list()
    for(i in 1:length(Ys_gbat)){
      suppressWarnings(expr={ Ys_gbat[[i]]<-fread(file = 
                                                   exp_files[i], sep='\t', data.table=F, check.names = F)})
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
    Yhat_gbat_mat = data.frame(Yhat_gbat_mat[-remove_inds,], check.names = F)
    Yhat_gbat_mat = cbind(id = rownames(Yhat_gbat_mat), Yhat_gbat_mat)
    fwrite(Yhat_gbat_mat, file = paste0(out_dir,gene_name,".cxc.predictors.txt"), sep = "\t")
    
    evaluation_helper(Ys_gbat, hom_expr_mat_gbat, Yhats_gbat, gbat_contexts_vec, TRUE, NULL, out_dir, gene_name)
  }
  message("Finished computing GReXs and evalutation metrics")
  
}
