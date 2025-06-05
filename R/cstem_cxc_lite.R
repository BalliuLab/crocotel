

### have to run this separately for gbat and Crocotel but with different parameter settings
#' @export
crocotel_lite = function(regulator_pred_exp_file, target_exp_files, contexts_vec, regulator_gene_name, target_gene_name, outdir, method = "Crocotel", target_cis_pred = F, target_pred_exp_file = NULL, r2_thresh = NULL, regulator_r2_file = NULL){
  dir.create(outdir, showWarnings = F)
  ## get target expression across all contexts
  regulator_exp_mat = fread(regulator_pred_exp_file, sep = "\t", data.table = F, check.names = F)
  
  if(!is.null(r2_thresh)){
    regulator_r2 = fread(regulator_r2_file, sep = "\t", data.table = F, check.names = F)
    ### checks that at least one contexts has r2 > threshold
    r2 = max(regulator_r2$full_cv_r2s, na.rm = T) < r2_thresh
    if(r2){
      print("Regulator GReX did not pass specified R2 threshold in any context. Not running Crocotel for this target regulator pair.")
      return(NULL)
    }
  }
  
  regulator_contexts = names(regulator_exp_mat)[names(regulator_exp_mat) %in% contexts_vec]
  target_contexts = c()
  for(cur_context in contexts_vec){
    cur_file = target_exp_files[grepl(paste0("/",cur_context), target_exp_files)]
    if(length(cur_file) !=0){
      target_contexts = c(target_contexts, cur_context)
    }
  }
  
  ### only run on contexts where target and regulator have information for
  intersected_contexts = intersect(regulator_contexts, target_contexts)
  
  this_gene = bind_rows(lapply(intersected_contexts, function(context_name){
    if(method != "GBAT"){
      print(paste("Running Crocotel lite for gene pair ", regulator_gene_name, " and ", target_gene_name, " in context ", context_name))
    }else{
      print(paste("Running GBAT for gene pair ", regulator_gene_name, " and ", target_gene_name, " in context ", context_name))
    }
    target_exp_vec = fread(target_exp_files[grepl(paste0("/",context_name), target_exp_files)], sep = "\t", data.table = F)
    names(target_exp_vec) = c("id", "target_exp")
    
    regulator_exp_vec = regulator_exp_mat[,c("id", context_name)]
    names(regulator_exp_vec) = c("id", "regulator_pred")
    
    if (target_cis_pred) {
      # if get target cis predicted expression across all contexts
      if(!is.null(target_pred_exp_file)){
        target_cis_pred_mat = fread(target_pred_exp_file, sep = "\t", data.table = F, check.names = F)
      }else{
        stop("Target GReX flag set, but no target GReX file provided. Exiting.")
      }
      target_cis_pred_vec = target_cis_pred_mat[,c("id", context_name)]
      names(target_cis_pred_vec) = c("id", "target_pred")
      lm_df = target_exp_vec %>%
        full_join(regulator_exp_vec, by = "id") %>%
        full_join(target_cis_pred_vec, by = "id")
      
      trans_model <- lm(lm_df[,"target_exp"] ~ lm_df[,"regulator_pred"] + lm_df[,"target_pred"])
      summary_model <- summary(trans_model)
      regulator_pvalue <- summary_model$coefficients[2, 4]
      regulator_beta = summary_model$coefficients[2, 1]
      regulator_se = summary_model$coefficients[2,2]
      regulator_tstat = summary_model$coefficients[2,3]
      #target_pvalue <- summary_model$coefficients[3, 4]
      #df = data.frame(SNP = target_gene_name, gene = regulator_gene_name, beta = regulator_beta, 'se' = regulator_se, 'pvalue' = regulator_pvalue, FDR = NA, context = context_name)#, target_pvalue))
      df = data.frame(target = target_gene_name, regulator = regulator_gene_name, beta = regulator_beta, 'se' = regulator_se, 'pvalue' = regulator_pvalue, context = context_name)#, target_pvalue))
    } else {
      lm_df = target_exp_vec %>%
        full_join(regulator_exp_vec, by = "id") 
      trans_model <- lm(lm_df[,"target_exp"] ~ lm_df[,"regulator_pred"])
      summary_model <- summary(trans_model)
      regulator_pvalue <- summary_model$coefficients[2, 4]
      regulator_beta = summary_model$coefficients[2, 1]
      regulator_se = summary_model$coefficients[2,2]
      #df = data.frame(SNP = target_gene_name, gene = regulator_gene_name, beta = regulator_beta, 'se' = regulator_se, 'pvalue' = regulator_pvalue, FDR = NA, context = context_name)
      df = data.frame(target = target_gene_name, regulator = regulator_gene_name, beta = regulator_beta, 'se' = regulator_se, 'pvalue' = regulator_pvalue, context = context_name)#, target_pvalue))
    }
  }))
  
  if(target_cis_pred){
    file_prefix = ".cis_crocotel_lite.txt"
    if(method == "GBAT"){
      file_prefix = ".cis_gbat.txt"
    }
  }else{
    file_prefix = ".crocotel_lite.txt"
    if(method == "GBAT"){
      file_prefix = ".gbat.txt"
    }
  }
  fwrite(this_gene, file = paste0(outdir, regulator_gene_name, "_", target_gene_name, file_prefix),  sep = "\t")
  return(this_gene)
}




