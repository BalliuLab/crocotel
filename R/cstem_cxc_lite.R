library(data.table)
library(dplyr)
library(foreach)

regulator_pred_exp_file = "/Users/lkrockenberger/C-STEM/example_data/GReXs/gene1_cstem_full_predictors.txt"
regulator_cxc_pred_exp_file = "/Users/lkrockenberger/C-STEM/example_data/GReXs/gene1_CxC_predictors.txt"
target_pred_exp_file = "/Users/lkrockenberger/C-STEM/example_data/GReXs/gene1_cstem_predictors.txt"
target_cxc_pred_exp_file = "/Users/lkrockenberger/C-STEM/example_data/GReXs/gene1_CxC_predictors.txt"
target_exp_files = list.files("/Users/lkrockenberger/C-STEM/example_data/expression/")
run_CxC = T
contexts_vec = target_exp_files
target_exp_files = paste0("/Users/lkrockenberger/C-STEM/example_data/expression/", target_exp_files)
regulator_gene_name = "gene1"
target_gene_name = "gene2"
target_cis_pred = TRUE
outdir = "/Users/lkrockenberger/C-STEM/example_data/trans_output/"

####################
### have to run this separately for CxC and C-STEM but with different parameter settings
cstem_cxc_lite = function(regulator_pred_exp_file, target_pred_exp_file, target_exp_files, contexts_vec, run_CxC, regulator_gene_name, target_gene_name, outdir, target_cis_pred = T){
  ## get target expression across all contexts
  regulator_exp_mat = fread(regulator_pred_exp_file, sep = "\t", data.table = F)
  
  # get target cis predicted expression across all contexts
  target_cis_pred_mat = fread(target_pred_exp_file, sep = "\t", data.table = F)
  
  this_gene = foreach(context = 1:length(contexts_vec), .combine = 'rbind') %dopar% {
    context_name = contexts_vec[context]
    print(paste("Running c-stem lite for gene pair ", regulator_gene_name, " and ", target_gene_name, " in context ", context_name))
    target_exp_vec = fread(target_exp_files[grepl(context_name, target_exp_files)], sep = "\t", data.table = F)
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
      #target_pvalue <- summary_model$coefficients[3, 4]
      return(data.frame(regulator_gene_name, target_gene_name, context_name, regulator_beta, regulator_se, regulator_pvalue))#, target_pvalue))
    } else {
      #trans_model <- lm.fit(cbind(1, regulator_exp_vec), trans_exp_df_vec)
      trans_model2 <- lm(lm_df[,"target_exp"] ~ lm_df[,"regulator_pred"])
      summary_model <- summary(trans_model2)
      regulator_pvalue <- summary_model$coefficients[2, 4]
      regulator_beta = summary_model$coefficients[2, 1]
      regulator_se = summary_model$coefficients[2,2]
      return(data.frame(regulator_gene_name, target_gene_name, context_name, regulator_beta, regulator_se, regulator_pvalue))
    }
  }
  if(target_cis_pred){
    file_prefix = "_cis_cstemlite.txt"
    if(run_CxC){
      file_prefix = "_cis_cxc.txt"
    }
  }else{
    file_prefix = "_cstemlite.txt"
    if(run_CxC){
      file_prefix = "_cxc.txt"
    }
  }
  fwrite(this_gene, file = paste0(outdir, regulator_gene_name, "_", target_gene_name, file_prefix),  sep = "\t")
}

cstem_cxc_lite(regulator_pred_exp_file, target_pred_exp_file, target_exp_files, contexts_vec, FALSE, regulator_gene_name, target_gene_name, outdir, target_cis_pred)
cstem_cxc_lite(regulator_cxc_pred_exp_file, target_cxc_pred_exp_file, target_exp_files, contexts_vec, TRUE, regulator_gene_name, target_gene_name, outdir, target_cis_pred)






