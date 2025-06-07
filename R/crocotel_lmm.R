
get_target_exp = function(target_exp_files, contexts_vec){
  contexts = c()
  for(cur_context in contexts_vec){
    cur_file = target_exp_files[grepl(paste0("/",cur_context), target_exp_files)]
    if(length(cur_file) !=0){
      contexts = c(contexts, cur_context)
    }
  }
  targ_exp = bind_rows(lapply(contexts_vec, function(cur_context){
    cur_file = target_exp_files[grepl(paste0("/",cur_context), target_exp_files)]
    if(length(cur_file) !=0){
      df = fread(cur_file, sep = "\t", data.table = F, check.names = F)
      names(df) = c("id", "target_exp")
      df = df %>% mutate(context = cur_context) %>% select(id, context, target_exp)
    }
  }))
  return(list(targ_exp, contexts))
}

#' @export
crocotel_lmm = function(regulator_pred_exp_file, target_exp_files, contexts_vec, regulator_gene_name, target_gene_name, outdir, target_cis_pred = F, target_pred_exp_file, r2_thresh = NULL, regulator_r2_file = NULL, context_dependence = F){
  dir.create(outdir, showWarnings = F)
  ## get target expression across all contexts
  target_output = get_target_exp(target_exp_files, contexts_vec)
  target_exp_mat = target_output[[1]]
  regulator_exp_mat = fread(regulator_pred_exp_file, sep = "\t", data.table = F, check.names = F, header = T)
  ###### get contexts present for regulator and for target
  regulator_contexts = names(regulator_exp_mat)
  target_contexts = target_output[[2]]
  intersected_contexts = intersect(regulator_contexts, target_contexts)
  
  ###### subset regulator and target expression matrices to only include intersected contexts
  regulator_exp_mat = regulator_exp_mat %>% pivot_longer(cols = intersected_contexts,
                                                         names_to = "context",
                                                         values_to = "regulator_pred" )
  
  ################################################################################
  ##### TAKE THIS OUT AFTER FIXING AVERAGE CONTEXT IN PREDICTOR FILE #############
  ################################################################################
  regulator_exp_mat = regulator_exp_mat %>% select(id, context, regulator_pred)
  target_exp_mat = target_exp_mat %>% filter(context %in% intersected_contexts)
  
  if(!is.null(r2_thresh)){
    regulator_r2 = fread(regulator_r2_file, sep = "\t", data.table = F, check.names = F, header = T)
    ### checks that at least one contexts has r2 > threshold
    r2 = max(regulator_r2$full_cv_r2s, na.rm = T) < r2_thresh
    if(r2){
      print("Regulator GReX did not pass specified R2 threshold in any context. Not running Crocotel for this target regulator pair.")
      return(NULL)
    }
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
      
      if (target_cis_pred) {
        target_cis_pred_vec = target_cis_pred_mat %>% filter(context == context_name) %>% select(id, target_cis_pred)
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
        df = data.frame(target = target_gene_name, regulator = regulator_gene_name, beta = regulator_beta, 'se' = regulator_se, 'pvalue' = regulator_pvalue, context = context_name)
      }
    }))
  }else{
    print(paste("Running Crocotel lmm for gene pair ", regulator_gene_name, " and ", target_gene_name))
    # get target cis predicted expression across all contexts
    target_cis_pred_mat = fread(target_pred_exp_file, sep = "\t", data.table = F, check.names = F, header = T)
    target_cis_pred_mat = target_cis_pred_mat %>% pivot_longer(cols = -id,
                                                               names_to = "context",
                                                               values_to = "target_cis_pred")
    
    #### build dataframe to run trans model - makes sure that all IDs are in the same order across vectors.
    trans_model_df = target_exp_mat %>%
      full_join(regulator_exp_mat, by = c("id", "context")) %>%
      full_join(target_cis_pred_mat, by = c("id", "context"))
    trans_model_df$id = factor(trans_model_df$id)
    trans_model_df$context = factor(trans_model_df$context)
    
    # setup pvalue matricies for target and regulator cis-predicted expression
    # pvalue matrix for cis-genetic predicted target associations with simulated trans expression
    target_assoc_pvalues = data.frame(matrix(0, nrow = 1, ncol = length(intersected_contexts)+1), check.names = F)
    names(target_assoc_pvalues) = c("target_gene", intersected_contexts)
    # pvalue matrix for cis-genetic predicted regulator associations with simulated trans expression
    regulator_assoc_pvalues = data.frame(matrix(0, nrow = 1, ncol = length(intersected_contexts)+2), check.names = F)
    names(regulator_assoc_pvalues) = c("target_gene", "regulator_gene", intersected_contexts)
    
    ref_context = intersected_contexts[1]
    trans_model_df$context <- relevel(factor(trans_model_df$context), ref = ref_context)
    if(target_cis_pred){
      trans_model <- suppressMessages(suppressWarnings(lmer(target_exp ~ regulator_pred + target_cis_pred + context + 
                            regulator_pred:context + target_cis_pred:context + (1 | id), data = trans_model_df)))
      if(context_dependence){
        print(paste("Assessing r2 of regulator GReX by context interaction for gene pair ", regulator_gene_name, " and ", target_gene_name, " in context ", context_name))
        null_model = suppressMessages(suppressWarnings(lmer(target_exp ~ regulator_pred + target_cis_pred + context + target_cis_pred:context + (1 | id), data = trans_model_df)))
        context_dependence_pval = anova(trans_model, null_model, test = "LRT")[2,"Pr(>Chisq)"]
        r2_null = r.squaredGLMM(null_model)[,"R2m"]
        r2_full = r.squaredGLMM(trans_model)[,"R2m"]
        
        # Calculate how much variance the interaction explains:
        delta_r2 <- r2_full - r2_null
        context_dependence_df = data.frame(r2 = delta_r2, pvalue = context_dependence_pval)
        fwrite(context_dependence_df, file = paste0(outdir, regulator_gene_name, "_", target_gene_name, ".reg_GReX_by_context_interaction_r2.txt"),  sep = "\t")
      }
      # extract marginal trends for each predicted exp
      #reg_marginal_trends = suppressMessages(suppressWarnings(emtrends(trans_model, ~ context, var = "regulator_pred", data = trans_model_df))) %>% tidy() %>% mutate(FDR = NA) %>% select(regulator_pred.trend, std.error, p.value, FDR, context)
      #names(reg_marginal_trends) = c("beta", "se", "pvalue", "FDR", "context")
      reg_marginal_trends = suppressMessages(suppressWarnings(emtrends(trans_model, ~ context, var = "regulator_pred", data = trans_model_df))) %>% tidy() %>% select(regulator_pred.trend, std.error, p.value, context)
      names(reg_marginal_trends) = c("beta", "se", "pvalue", "context")
      target_marginal_trends <- suppressMessages(suppressWarnings(emtrends(trans_model, ~ context, var = "target_cis_pred"))) %>% tidy()
      output_df = cbind(target = target_gene_name, regulator = regulator_gene_name, reg_marginal_trends)
    } else {
      trans_model <- suppressMessages(suppressWarnings(lmer(target_exp ~ regulator_pred + context + 
                            regulator_pred:context + (1 | id), data = trans_model_df)))
      # extract marginal trends for each predicted exp
      reg_marginal_trends = suppressMessages(suppressWarnings(emtrends(trans_model, ~ context, var = "regulator_pred", data = trans_model_df))) %>% tidy() %>% select(regulator_pred.trend, std.error, p.value, context)
      names(reg_marginal_trends) = c("beta", "se", "pvalue", "context")
      output_df = cbind(target = target_gene_name, regulator = regulator_gene_name, reg_marginal_trends)
    }
  }
  
  if(target_cis_pred){
    file_prefix = ".cis_crocotel_lmm.txt"
  }else{
    file_prefix = ".crocotel_lmm.txt"
  }
  
  fwrite(output_df, file = paste0(outdir, regulator_gene_name, "_", target_gene_name, file_prefix),  sep = "\t")
  return(output_df)
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







