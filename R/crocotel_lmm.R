
#' @export
crocotel_lmm = function(regulator_gene_name, target_gene_name, out_dir, target_exp_file = NULL, GReX_dir = NULL, regress_target_GReX = T, pval_thresh = 1, r2_thresh = NULL, context_dependence = F){
  out_dir_crocotel_lmm = paste0(out_dir, "/crocotel_lmm_output/")
  dir.create(out_dir_crocotel_lmm, showWarnings = F)
  
  if(is.null(GReX_dir)){
    message("inferring GReX directory and regulator GReX file...")
    GReX_dir = paste0(out_dir, "/GReXs/")
    regulator_pred_exp_file = list.files(GReX_dir, pattern = paste0(regulator_gene_name, ".crocotel.GReX_predictors.txt"), full.names = T)
  }
  if(is.null(target_exp_file)){
    if(regress_target_GReX == T){
      message("inferring expression files with residualized target GReX and formatting.")
      target_exp_file = list.files(paste0(out_dir, "/exp_residualized_GReX/"), pattern = paste0(target_gene_name, ".crocotel.GReX_residuals.txt"), full.names = T)
    }else{
      stop("You are running crocotel lmm without regressing out GReX for each target. Cannot infer expression of the target. Please specify an expression file for the target_exp_file parameter.")
    }
  }
  
  ## get target expression across all contexts
  target_output = fread(target_exp_file, sep = "\t", data.table = F, header = T)
  target_exp_mat = pivot_longer(target_output, cols = -id, names_to = "context", values_to = "target_exp")
  
  regulator_exp_mat = fread(regulator_pred_exp_file, sep = "\t", data.table = F, check.names = F, header = T)
  ###### get contexts present for regulator and for target
  regulator_contexts = names(regulator_exp_mat)
  target_contexts = unique(target_exp_mat$context)
  intersected_contexts = intersect(regulator_contexts, target_contexts)
  
  ###### subset regulator and target expression matrices to only include intersected contexts
  regulator_exp_mat = regulator_exp_mat %>% pivot_longer(cols = intersected_contexts,
                                                         names_to = "context",
                                                         values_to = "regulator_pred" )
  target_exp_mat = target_exp_mat %>% filter(context %in% intersected_contexts)
  
  if(!is.null(r2_thresh)){
    regulator_r2_file = list.files(paste0(out_dir, "/GReXs/"), pattern = paste0(regulator_gene_name, ".crocotel.crossval_r2.txt"), full.names = T)
    regulator_r2 = fread(regulator_r2_file, sep = "\t", data.table = F, check.names = F, header = T)
    ### checks that at least one contexts has r2 > threshold
    r2 = max(regulator_r2$full_cv_r2s, na.rm = T) < r2_thresh
    if(r2){
      print("Regulator GReX did not pass specified R2 threshold in any context. Not running Crocotel for this regulator-target pair.")
      return(NULL)
    }
  }
  
  if(regress_target_GReX){
    file_prefix = ".crocotel_lmm_regress.txt"
  }else{
    file_prefix = ".crocotel_lmm.txt"
  }
  
  if(length(intersected_contexts) < 3){
    print("Target and Regulator do not include enough overlapping contexts to run Crocotel lmm, running Crocotel lite for this pair. Overall results will still be saved in lmm file.")
    output_df = bind_rows(lapply(intersected_contexts, function(context_name){
      print(paste("Running Crocotel lite for gene pair ", regulator_gene_name, " and ", target_gene_name, " in context ", context_name))
      target_exp_vec = target_exp_mat %>% filter(context == context_name) %>% select(id, target_exp)
      names(target_exp_vec) = c("id", "target_exp")
      
      #### build dataframe to run trans model - makes sure that all IDs are in the same order across vectors.
      regulator_exp_vec = regulator_exp_mat %>% filter(context == context_name) %>% select(id, regulator_pred)
      names(regulator_exp_vec) = c("id", "regulator_pred")
      
      lm_df = target_exp_vec %>%
        full_join(regulator_exp_vec, by = "id")
      trans_model <- lm(lm_df$target_exp ~ lm_df$regulator_pred)
      summary_model <- summary(trans_model)
      regulator_pvalue <- summary_model$coefficients[2, 4]
      regulator_beta = summary_model$coefficients[2, 1]
      regulator_se = summary_model$coefficients[2,2]
      #df = data.frame(SNP = target_gene_name, gene = regulator_gene_name, beta = regulator_beta, 'se' = regulator_se, 'pvalue' = regulator_pvalue, FDR = NA, context = context_name)
      df = data.frame(target = target_gene_name, regulator = regulator_gene_name, beta = regulator_beta, 'se' = regulator_se, 'pvalue' = regulator_pvalue, context = context_name)
    }))
    
    for(cur_context in intersected_contexts){
      context_df = output_df %>% filter(context == cur_context) %>% select(-context)
      fwrite(context_df, file = paste0(out_dir_crocotel_lmm, cur_context, ".", regulator_gene_name, "_", target_gene_name, file_prefix),
             sep = "\t", quote = F)
    }
  }else{
    print(paste("Running Crocotel lmm for gene pair ", regulator_gene_name, " and ", target_gene_name, " in contexts ", paste0(intersected_contexts, collapse = ",")))
    #### build dataframe to run trans model - makes sure that all IDs are in the same order across vectors.
    trans_model_df = target_exp_mat %>%
      full_join(regulator_exp_mat, by = c("id", "context"))
    trans_model_df$id = factor(trans_model_df$id)
    trans_model_df$context = factor(trans_model_df$context)
    ref_context = intersected_contexts[1]
    trans_model_df$context <- relevel(factor(trans_model_df$context), ref = ref_context)

    trans_model <- suppressMessages(suppressWarnings(lmer(target_exp ~ regulator_pred + context + 
                          regulator_pred:context + (1 | id), data = trans_model_df)))
    if(context_dependence){
      print(paste("Assessing r2 of regulator GReX by context interaction for gene pair ", regulator_gene_name, " and ", target_gene_name))
      null_model = suppressMessages(suppressWarnings(lmer(target_exp ~ regulator_pred + context + (1 | id), data = trans_model_df)))
      context_dependence_pval = anova(trans_model, null_model, test = "LRT")[2,"Pr(>Chisq)"]
      r2_null = r.squaredGLMM(null_model)[,"R2m"]
      r2_full = r.squaredGLMM(trans_model)[,"R2m"]
      
      # Calculate how much variance the interaction explains:
      delta_r2 <- r2_full - r2_null
      context_dependence_df = data.frame(r2 = delta_r2, pvalue = context_dependence_pval)
      fwrite(context_dependence_df, file = paste0(out_dir_crocotel_lmm, regulator_gene_name, "_", target_gene_name, ".reg_GReX_by_context_interaction_r2.", file_prefix),  sep = "\t")
    }
    # extract marginal trends for each predicted exp
    #reg_marginal_trends = suppressMessages(suppressWarnings(emtrends(trans_model, ~ context, var = "regulator_pred", data = trans_model_df))) %>% tidy() %>% mutate(FDR = NA) %>% select(regulator_pred.trend, std.error, p.value, FDR, context)
    #names(reg_marginal_trends) = c("beta", "se", "pvalue", "FDR", "context")
    reg_marginal_trends = suppressMessages(suppressWarnings(emtrends(trans_model, ~ context, var = "regulator_pred", data = trans_model_df))) %>% tidy() %>% select(regulator_pred.trend, std.error, p.value, context)
    names(reg_marginal_trends) = c("beta", "se", "pvalue", "context")
    output_df = cbind(regulator = regulator_gene_name, target = target_gene_name, reg_marginal_trends)
    contexts = unique(output_df$context)
    for(cur_context in contexts){
      context_df = output_df %>% filter(context == cur_context) %>% select(-context)
      fwrite(context_df, file = paste0(out_dir_crocotel_lmm, cur_context, ".", regulator_gene_name, "_", target_gene_name, file_prefix),
             sep = "\t", quote = F)
    }
  }
  print("Finished running crocotel lmm for this pair.")
}

#' @export
format_original_expression_crocotel_lmm = function(workdir, out_dir){
  gene_dirs = list.dirs(path = workdir, full.names = T, recursive = F)
  
  for (gene_dir in gene_dirs) {
    gene_id = basename(gene_dir)
    contexts = list.files(path = gene_dir, full.names = T)
    context_dfs = list()
    
    for (context_file in contexts) {
      context_id = basename(context_file)
      context_id = sub("\\..*", "", context_id)
      df = fread(context_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, data.table = F, check.names = F)
      colnames(df) = c("id", context_id)
      context_dfs[[context_id]] = df
    }
    
    # Merge all context data.frames by sample ID
    merged_df <- Reduce(function(x, y) merge(x, y, by = "id", all = TRUE), context_dfs)
    
    out_file = paste0(out_dir, paste0(gene_id, "_expression.txt"))
    fwrite(merged_df, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}



#regulator_pred_exp_file = "/u/scratch/l/lkrocken/crocotile_example/GReXs/gene1.cstem.full_predictors.txt"
#regulator_r2_file = "/u/scratch/l/lkrocken/crocotile_example/GReXs/gene1.cstem.crossval_r2.txt"
#target_pred_exp_file = "/u/scratch/l/lkrocken/crocotile_example/GReXs/gene2.cstem.full_predictors.txt"
#target_exp_files = "/u/scratch/l/lkrocken/crocotile_example/input_data/gene2"
#target_exp_files = "/u/scratch/l/lkrocken/crocotile_example/input_data/gene2/"
#target_exp_files = list.files(target_exp_files, full.names = T)
#target_exp_files
#contexts_vec = paste0(seq(0,9))
#contexts_vec
#regulator_gene_name = "gene1"
#target_gene_name = "gene2"
#outdir = "/u/scratch/l/lkrocken/"
#target_cis_pred = T
#r2_trhesh = 0.01
#r2_thresh = NULL







