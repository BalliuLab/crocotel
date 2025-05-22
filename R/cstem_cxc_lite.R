#library(data.table)
#library(dplyr)
#library(foreach)

#regulator_pred_exp_file = "/Users/lkrockenberger/C-STEM/example_data/GReXs/gene1_cstem_full_predictors.txt"
#regulator_gbat_pred_exp_file = "/Users/lkrockenberger/C-STEM/example_data/GReXs/gene1_gbat_predictors.txt"
#target_pred_exp_file = "/Users/lkrockenberger/C-STEM/example_data/GReXs/gene1_cstem_predictors.txt"
#target_gbat_pred_exp_file = "/Users/lkrockenberger/C-STEM/example_data/GReXs/gene1_gbat_predictors.txt"
#target_exp_files = list.files("/Users/lkrockenberger/C-STEM/example_data/expression/")
#method = "CSTEM"
#contexts_vec = target_exp_files
#target_exp_files = paste0("/Users/lkrockenberger/C-STEM/example_data/expression/", target_exp_files)
#regulator_gene_name = "gene1"
#target_gene_name = "gene2"
#target_cis_pred = TRUE
#outdir = "/Users/lkrockenberger/C-STEM/example_data/trans_output/"

####################
### have to run this separately for gbat and C-STEM but with different parameter settings
#' @export
cstem_gbat_lite = function(regulator_pred_exp_file, target_pred_exp_file, target_exp_files, contexts_vec, method, regulator_gene_name, target_gene_name, outdir, target_cis_pred = T){
  dir.create(outdir, showWarnings = F)
  ## get target expression across all contexts
  regulator_exp_mat = fread(regulator_pred_exp_file, sep = "\t", data.table = F, check.names = F)
  
  # get target cis predicted expression across all contexts
  target_cis_pred_mat = fread(target_pred_exp_file, sep = "\t", data.table = F, check.names = F)
  
  this_gene = bind_rows(lapply(contexts_vec, function(context_name){
    print(paste("Running c-stem lite for gene pair ", regulator_gene_name, " and ", target_gene_name, " in context ", context_name))
    target_exp_vec = fread(target_exp_files[grepl(paste0("/",context_name), target_exp_files)], sep = "\t", data.table = F)
    names(target_exp_vec) = c("id", "target_exp")
    
    regulator_exp_vec = regulator_exp_mat[,c("id", context_name)]
    names(regulator_exp_vec) = c("id", "regulator_pred")
    target_cis_pred_vec = target_cis_pred_mat[,c("id", context_name)]
    names(target_cis_pred_vec) = c("id", "target_pred")
    lm_df = target_exp_vec %>%
      full_join(regulator_exp_vec, by = "id") %>%
      full_join(target_cis_pred_vec, by = "id")
    
    if (target_cis_pred) {
      trans_model <- lm(lm_df[,"target_exp"] ~ lm_df[,"regulator_pred"] + lm_df[,"target_pred"])
      summary_model <- summary(trans_model)
      regulator_pvalue <- summary_model$coefficients[2, 4]
      regulator_beta = summary_model$coefficients[2, 1]
      regulator_se = summary_model$coefficients[2,2]
      regulator_tstat = summary_model$coefficients[2,3]
      #target_pvalue <- summary_model$coefficients[3, 4]
      #df = data.frame(regulator_gene_name, target_gene_name, context_name, regulator_beta, regulator_se, regulator_pvalue)#, target_pvalue))
      df = data.frame(SNP = target_gene_name, gene = regulator_gene_name, beta = regulator_beta, 'se' = regulator_se, 'pvalue' = regulator_pvalue, FDR = NA, context = context_name)#, target_pvalue))
    } else {
      #trans_model <- lm.fit(cbind(1, regulator_exp_vec), trans_exp_df_vec)
      trans_model2 <- lm(lm_df[,"target_exp"] ~ lm_df[,"regulator_pred"])
      summary_model <- summary(trans_model2)
      regulator_pvalue <- summary_model$coefficients[2, 4]
      regulator_beta = summary_model$coefficients[2, 1]
      regulator_se = summary_model$coefficients[2,2]
      df = data.frame(SNP = target_gene_name, gene = regulator_gene_name, beta = regulator_beta, 'se' = regulator_se, 'pvalue' = regulator_pvalue, FDR = NA, context = context_name)
    }
  }))
  
  if(target_cis_pred){
    file_prefix = ".cis_cstemlite.txt"
    if(method == "GBAT"){
      file_prefix = ".cis_gbat.txt"
    }
  }else{
    file_prefix = ".cstemlite.txt"
    if(method == "GBAT"){
      file_prefix = ".gbat.txt"
    }
  }
  fwrite(this_gene, file = paste0(outdir, regulator_gene_name, "_", target_gene_name, file_prefix),  sep = "\t")
  return(this_gene)
}

#cstem_gbat_lite(regulator_pred_exp_file, target_pred_exp_file, target_exp_files, contexts_vec, FALSE, regulator_gene_name, target_gene_name, outdir, target_cis_pred)
#cstem_gbat_lite(regulator_gbat_pred_exp_file, target_gbat_pred_exp_file, target_exp_files, contexts_vec, TRUE, regulator_gene_name, target_gene_name, outdir, target_cis_pred)






