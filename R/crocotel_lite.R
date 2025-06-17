

### have to run this separately for gbat and Crocotel but with different parameter settings
crocotel_lite_old = function(regulator_pred_exp_file, target_exp_files, contexts_vec, regulator_gene_name, target_gene_name, outdir, method = "Crocotel", target_cis_pred = F, target_pred_exp_file = NULL, r2_thresh = NULL, regulator_r2_file = NULL, write_output = F){
  dir.create(outdir, showWarnings = F)
  ## get target expression across all contexts
  regulator_exp_mat = fread(regulator_pred_exp_file, sep = "\t", data.table = F, check.names = F, header = T)
  
  if(!is.null(r2_thresh)){
    regulator_r2 = fread(regulator_r2_file, sep = "\t", data.table = F, check.names = F, header = T)
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
        target_cis_pred_mat = fread(target_pred_exp_file, sep = "\t", data.table = F, check.names = F, header = T)
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
  if(write_output){
    fwrite(this_gene, file = paste0(outdir, regulator_gene_name, "_", target_gene_name, file_prefix),  sep = "\t") 
  }
  return(this_gene)
}

get_trans_genes = function(gene, geneloc_file, trans_threshold){
  genelocs = fread(geneloc_file, sep = "\t", data.table = F)
  cur_gene_chr = genelocs %>% filter(geneid == gene) %>% select(chr) %>% unlist()
  cur_gene_start = genelocs %>% filter(geneid == gene) %>% select(s1) %>% unlist()
  upstream_pos = cur_gene_start - trans_threshold
  downstream_pos = cur_gene_start + trans_threshold
  same_chr_trans = genelocs %>% filter(chr == chr, geneid != gene) %>% filter(s1 <= upstream_pos & s2 >= downstream_pos) %>% select(geneid) %>% unlist() %>% unname()
  trans_genes = genelocs %>% filter(chr != cur_gene_chr) %>% select(geneid) %>% unlist() %>% unname()
  trans_genes = c(trans_genes, same_chr_trans)
  return(trans_genes)
}

trans_genes = get_trans_genes(gene_name, geneloc_file, trans_threshold)
regulator_pred_exp_file = paste0(workdir, gene_name, ".cstem.full_predictors.txt")
regulator_gene_name = gene_name
output = bind_rows(lapply(trans_genes, function(trans_gene){
  if(trans_gene != gene_name){
    target_pred_exp_file = paste0(workdir, trans_gene, ".cstem.full_predictors.txt")
    if(file.exists(target_pred_exp_file)){
      target_exp_files = list.files(paste0(target_exp_dir, trans_gene, "/"), full.names = T)
      target_gene_name = trans_gene
      
      #df = crocotel_lmm(regulator_pred_exp_file, target_exp_files, contexts_vec, regulator_gene_name, target_gene_name, outdir, target_cis_pred, target_pred_exp_file)
      df = crocotel_lite(regulator_pred_exp_file, target_exp_files, contexts_vec, regulator_gene_name, target_gene_name, outdir, method, target_cis_pred, target_pred_exp_file)
      
    }
  }
}))
fwrite(output, file = paste0(outdir, gene_name, ".", method, ".txt"), sep = "\t")


#' @export
crocotel_lite = function(regulator_pred_exp_file, target_exp_dir, contexts_vec, outdir, r2_thresh = NULL, regulator_r2_file = NULL, write_output = F){
  dir.create(outdir, showWarnings = F)
  ## get target expression across all contexts
  regulator_gene_name = sapply(strsplit(regulator_pred_exp_file, "\\."), "[[", 1)
  
  regulator_exp_mat = fread(regulator_pred_exp_file, sep = "\t", data.table = F, check.names = F, header = T)
  
  if(!is.null(r2_thresh)){
    regulator_r2 = fread(regulator_r2_file, sep = "\t", data.table = F, check.names = F, header = T)
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
      print(paste("Running Crocotel lite for gene ", regulator_gene_name, " in context ", context_name))

    target_exp_vec = fread(target_exp_files[grepl(paste0("/",context_name), target_exp_files)], sep = "\t", data.table = F)
    names(target_exp_vec) = c("id", "target_exp")
    
    regulator_exp_vec = regulator_exp_mat[,c("id", context_name)]
    names(regulator_exp_vec) = c("id", "regulator_pred")
    
    if (target_cis_pred) {
      # if get target cis predicted expression across all contexts
      if(!is.null(target_pred_exp_file)){
        target_cis_pred_mat = fread(target_pred_exp_file, sep = "\t", data.table = F, check.names = F, header = T)
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
  if(write_output){
    fwrite(this_gene, file = paste0(outdir, regulator_gene_name, "_", target_gene_name, file_prefix),  sep = "\t") 
  }
  return(this_gene)
}




