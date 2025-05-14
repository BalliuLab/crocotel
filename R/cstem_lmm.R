library(data.table)
library(dplyr)
library(tidyr)




# Script that implements CSTEM-lmm
# lmm based association testing between simulated trans and CONTENT predictions
# also includes optional target-cis predicted regression covariate

args=commandArgs(trailingOnly = T)

scenario <- args[1]
num_contexts <- args[2]
num_sim_genes <- args[3]
method <- args[4]
out_dir <- args[5]
content_dir <- args[6]
work_dir <- args[7]
simulated_trans_exp <- args[8]
target_cis_pred <- args[9]

#############
# test inputs
work_dir <- "/u/scratch/n/nmohamm/cstem_new_sims/test"
scenario <- 1
num_contexts <- 3
num_sim_genes <- 10
method <- "GBAT"
content_dir <- "CONTENT"
simulated_trans_exp <- "trans"
out_dir <- "/u/scratch/n/nmohamm/cstem_new_sims/test_cstem_outputs"
target_cis_pred <- TRUE
#############

format_lmm_data <- function(simulated_trans_data_df, content_reg_preds, content_target_preds, num_contexts) {
	# take in the simulated trans data and content predictions and
	# pivot longer to melding all contexts into 1 column
	# ultimately so data can easily be passed into lmm
	
	# merge trans expression across all contexts into 1 column and create
	# column to designate which context
	long_sim_data <- simulated_trans_data_df %>%
		pivot_longer(cols = everything(), names_to = "context", values_to = "trans_exp")
	# merge all contexts into one vector
	reg_pred_vec <- unlist(content_reg_preds, use.names = FALSE)
	target_pred_vec <- unlist(content_target_preds, use.names = FALSE)
	long_sim_data$reg_predicted_exp <- reg_pred_vec
	long_sim_data$target_predicted_exp <- target_pred_vec
	# add indiv names
	long_sim_data$indivs <- rep(seq_len(nrow(simulated_trans_data_df)), num_contexts)
	long_sim_data$context <- factor(long_sim_data$context)

	return(long_sim_data)
}

# load in trans exp files
trans_exp_files <- list.files(file.path(work_dir, scenario, simulated_trans_exp), full.names = TRUE, pattern = "^[0-9]")
trans_exp_data <- lapply(trans_exp_files, fread, sep = "\t", header = TRUE, data.table = FALSE)

# setup pvalue matricies for target and regulator cis-predicted expression
# pvalue matrix for cis-genetic predicted target associations with simulated trans expression
target_assoc_pvalues <- data.frame(matrix(0, nrow = num_sim_genes, ncol = num_contexts))
# pvalue matrix for cis-genetic predicted regulator associations with simulated trans expression
regulator_assoc_pvalues <- array(0, dim = c(num_sim_genes, num_sim_genes, num_contexts))


# loop through all pairs of regulator and target genes for association testing
for (gene1 in 1:num_sim_genes) {
	for (gene2 in 1:num_sim_genes) {
	        # get CONTENT preds
	        content_preds <- get_content_preds(gene1, gene2, scenario, work_dir, content_dir)
	        regulator_yhats <- content_preds$gene1
	        target_yhats <- content_preds$gene2
	
	        # set Y hats as either full or tissue specific for either CSTEM or GBAT mode
	        if (method == "CSTEM") {
	                regulator_exp <- regulator_yhats$Yhats_full
	                target_exp <- target_yhats$Yhats_full
	        } else {
	                regulator_exp <- regulator_yhats$Yhats_tiss
	                target_exp <- target_yhats$Yhats_tiss
	        }
	
	        # load in simulated trans expression data for current gene
		trans_exp_df <- trans_exp_data[[gene_num]]
	
		# create lmm compatible input data
		lmm_inputs <- format_lmm_data(trans_exp_df, regulator_exp, target_exp, num_contexts)
	
		# fit model depending on if we are including target predictions
		if (target_cis_pred) {
			# model both the regulator and target predicted effect as a 3 way interaction with context
			#trans_model <- lmer(trans_exp ~ reg_predicted_exp * target_predicted_exp * relevel(context, ref = 1) + (1 | indivs), data = lmm_inputs)
			trans_model <- lmer(trans_exp ~ reg_predicted_exp + target_predicted_exp + relevel(context, ref = 1) + 
			                    reg_predicted_exp:relevel(context, ref = 1) + 
			                    target_predicted_exp:relevel(context, ref = 1) + (1 | indivs), data = lmm_inputs)
			# extract marginal trends for each predicted exp
			reg_marginal_trends <- emtrends(trans_model, ~ context, var = "reg_predicted_exp") %>% tidy()
			target_marginal_trends <- emtrends(trans_model, ~ context, var = "target_predicted_exp") %>% tidy()
			# add p values to respective matrix
			regulator_assoc_pvalues[gene1, gene2, ] <- reg_marginal_trends %>% pull(p.value)
			target_assoc_pvalues[gene1, gene2, ] <- target_marginal_trends %>% pull(p.value)
		} else {
			trans_model <- lmer(trans_exp ~ reg_predicted_exp * relevel(context, ref = 1) + (1 | indivs), data = lmm_inputs)
			# calculate average signal of reg_predicted_exp in each context
			marginal_trends <- emtrends(trans_model, ~ context, var = "reg_predicted_exp") %>% tidy()
			# add p values from each context into matrix
			regulator_assoc_pvalues[gene1, gene2, ] <- marginal_trends %>% pull(p.value)
		}
			
	}
}

# assign names to assoc matrices
dimnames(regulator_assoc_pvalues) <- list(
  gene1 = 1:num_sim_genes, 
  gene2 = 1:num_sim_genes, 
  context = 1:num_contexts
)

dimnames(target_assoc_pvalues) <- list(
  gene1 = 1:num_sim_genes, 
  gene2 = 1:num_sim_genes, 
  context = 1:num_contexts
)

# convert to tibble for export
regulator_pvalues_df <- as_tibble(as.table(regulator_assoc_pvalues), .name_repair = "minimal")
target_pvalues_df <- as_tibble(as.table(target_assoc_pvalues), .name_repair = "minimal")

# Rename columns for clarity
colnames(regulator_pvalues_df) <- c("gene1", "gene2", "context", "p_value")
colnames(target_pvalues_df) <- c("gene1", "gene2", "context", "p_value")

if (target_cis_pred) {
        fwrite(target_pvalues_df, file.path(out_dir, paste0(method, "_lmm_target_scenario_", scenario, ".txt")), sep = "\t", row.names = T, col.names = T, quote = F)
}
fwrite(regulator_pvalues_df, file.path(out_dir, paste0(method, "_lmm_regulator_scenario_", scenario,".txt")), sep = "\t", row.names = T, col.names = T, quote = F)

########################################

library(data.table)
library(dplyr)
library(foreach)
library(tidyr)
library(lme4)
library(emmeans)
library(broom)

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


get_target_exp = function(target_exp_files, contexts_vec){
  targ_exp = foreach(context = 1:length(contexts_vec), .combine = 'rbind') %dopar% {
    cur_context = contexts_vec[context]
    cur_file = target_exp_files[grepl(cur_context, target_exp_files)]
    df = fread(cur_file, sep = "\t", data.table = F)
    names(df) = c("id", "target_exp")
    df = df %>% mutate(context = cur_context) %>% select(id, context, target_exp)
    return(df)
  }
  return(targ_exp)
}

####################
### have to run this separately for CxC and C-STEM but with different parameter settings
cstem_cxc_lmm = function(regulator_pred_exp_file, target_pred_exp_file, target_exp_files, contexts_vec, run_CxC, regulator_gene_name, target_gene_name, outdir, target_cis_pred = T){
  ## get target expression across all contexts
  regulator_exp_mat = fread(regulator_pred_exp_file, sep = "\t", data.table = F)
  regulator_exp_mat = regulator_exp_mat %>% pivot_longer(cols = -id,
                                                         names_to = "context",
                                                         values_to = "regulator_pred" )
  
  # get target cis predicted expression across all contexts
  target_cis_pred_mat = fread(target_pred_exp_file, sep = "\t", data.table = F)
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
  target_assoc_pvalues = data.frame(matrix(0, nrow = 1, ncol = length(contexts_vec)+1))
  names(target_assoc_pvalues) = c("target_gene", contexts_vec)
  # pvalue matrix for cis-genetic predicted regulator associations with simulated trans expression
  regulator_assoc_pvalues = data.frame(matrix(0, nrow = 1, ncol = length(contexts_vec)+2))
  names(regulator_assoc_pvalues) = c("target_gene", "regulator_gene", contexts_vec)
  
  if(target_cis_pred){
    ref_context = contexts_vec[1]
    trans_model <- lmer(target_exp ~ regulator_pred + target_cis_pred + relevel(context, ref = ref_context) + 
                          regulator_pred:relevel(context, ref = ref_context) + 
                          target_cis_pred:relevel(context, ref = ref_context) + (1 | id), data = trans_model_df)
    # extract marginal trends for each predicted exp
    reg_marginal_trends <- emtrends(trans_model, ~ context, var = "regulator_pred") %>% tidy()
    target_marginal_trends <- emtrends(trans_model, ~ context, var = "target_cis_pred") %>% tidy()
    # add p values to respective matrix
    regulator_assoc_pvalues[gene1, gene2, ] <- reg_marginal_trends %>% pull(p.value)
    target_assoc_pvalues[gene1, gene2, ] <- target_marginal_trends %>% pull(p.value)
  } else {
    trans_model <- lmer(target_exp ~ regulator_pred + relevel(context, ref = ref_context) + 
      regulator_pred:relevel(context, ref = ref_context) + (1 | id), data = trans_model_df)
    # calculate main effect of predictor
    marginal_trends <- emtrends(trans_model, ~ context, var = "regulator_pred") %>% tidy()
    # add p values from each context into matrix
    regulator_assoc_pvalues[gene1, gene2, ] <- marginal_trends %>% pull(p.value)
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











