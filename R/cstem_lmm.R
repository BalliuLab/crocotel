
#library(data.table)
#library(dplyr)
#library(foreach)
#library(tidyr)
#library(lme4)
#library(emmeans)
#library(broom)

#setwd("project-bballiu/C-STEM/r_package")
#regulator_pred_exp_file = "example_data/GReXs/gene1_cstem_full_predictors.txt"
#target_pred_exp_file = "example_data/GReXs/gene1_cstem_predictors.txt"
#target_exp_files = list.files("example_data/expression/")
#contexts_vec = target_exp_files
#target_exp_files = paste0("example_data/expression/", target_exp_files)
#regulator_gene_name = "gene1"
#target_gene_name = "gene2"
#target_cis_pred = TRUE
#outdir = "example_data/trans_output/"


get_target_exp = function(target_exp_files, contexts_vec){
  targ_exp = bind_rows(lapply(contexts_vec, function(cur_context){
    cur_file = target_exp_files[grepl(paste0("/",cur_context), target_exp_files)]
    df = fread(cur_file, sep = "\t", data.table = F, check.names = F)
    names(df) = c("id", "target_exp")
    df = df %>% mutate(context = cur_context) %>% select(id, context, target_exp)
  }))
  return(targ_exp)
}

#' @export
cstem_lmm = function(regulator_pred_exp_file, target_pred_exp_file, target_exp_files, contexts_vec, regulator_gene_name, target_gene_name, outdir, target_cis_pred = T){
  dir.create(outdir, showWarnings = F)
  ## get target expression across all contexts
  regulator_exp_mat = fread(regulator_pred_exp_file, sep = "\t", data.table = F, check.names = F)
  regulator_exp_mat = regulator_exp_mat %>% pivot_longer(cols = -id,
                                                         names_to = "context",
                                                         values_to = "regulator_pred" )
  
  # get target cis predicted expression across all contexts
  target_cis_pred_mat = fread(target_pred_exp_file, sep = "\t", data.table = F, check.names = F)
  target_cis_pred_mat = target_cis_pred_mat %>% pivot_longer(cols = -id,
                                                             names_to = "context",
                                                             values_to = "target_cis_pred")
  
  target_exp_mat = get_target_exp(target_exp_files, contexts_vec)
  trans_model_df = target_exp_mat %>%
    full_join(regulator_exp_mat, by = c("id", "context")) %>%
    full_join(target_cis_pred_mat, by = c("id", "context"))
  trans_model_df$id = factor(trans_model_df$id)
  trans_model_df$context = factor(trans_model_df$context)
  
  # setup pvalue matricies for target and regulator cis-predicted expression
  # pvalue matrix for cis-genetic predicted target associations with simulated trans expression
  target_assoc_pvalues = data.frame(matrix(0, nrow = 1, ncol = length(contexts_vec)+1), check.names = F)
  names(target_assoc_pvalues) = c("target_gene", contexts_vec)
  # pvalue matrix for cis-genetic predicted regulator associations with simulated trans expression
  regulator_assoc_pvalues = data.frame(matrix(0, nrow = 1, ncol = length(contexts_vec)+2), check.names = F)
  names(regulator_assoc_pvalues) = c("target_gene", "regulator_gene", contexts_vec)
  
  ref_context = contexts_vec[1]
  trans_model_df$context <- relevel(factor(trans_model_df$context), ref = ref_context)
  if(target_cis_pred){
    trans_model <- lmer(target_exp ~ regulator_pred + target_cis_pred + context + 
                          regulator_pred:context + target_cis_pred:context + (1 | id), data = trans_model_df)
    # extract marginal trends for each predicted exp
    reg_marginal_trends = emtrends(trans_model, ~ context, var = "regulator_pred", data = trans_model_df) %>% tidy() %>% mutate(FDR = NA) %>% select(regulator_pred.trend, std.error, p.value, FDR, context)
    names(reg_marginal_trends) = c("beta", "se", "pvalue", "FDR", "context")
    target_marginal_trends <- emtrends(trans_model, ~ context, var = "target_cis_pred") %>% tidy()
    output_df = cbind(SNP = target_gene_name, gene = regulator_gene_name, reg_marginal_trends)
  } else {
    trans_model <- lmer(target_exp ~ regulator_pred + context + 
                          regulator_pred:context + (1 | id), data = trans_model_df)
    # extract marginal trends for each predicted exp
    reg_marginal_trends = emtrends(trans_model, ~ context, var = "regulator_pred", data = trans_model_df) %>% tidy() %>% mutate(FDR = NA) %>% select(regulator_pred.trend, std.error, p.value, FDR, context)
    names(reg_marginal_trends) = c("beta", "se", "pvalue", "FDR", "context")
    output_df = cbind(SNP = target_gene_name, gene = regulator_gene_name, reg_marginal_trends)
  }
  if(target_cis_pred){
    file_prefix = ".cis_cstemlmm.txt"
  }else{
    file_prefix = ".cstemlmm.txt"
  }
  fwrite(output_df, file = paste0(outdir, regulator_gene_name, "_", target_gene_name, file_prefix),  sep = "\t")
  return(output_df)
}











