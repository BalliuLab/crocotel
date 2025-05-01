library(data.table)
library(foreach)
library(doParallel)
library(tidyr)
library(readr)
library(dplyr)

# Script that implements CSTEM-lite
# lm based association testing between simulated trans and CONTENT predictions
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
target_cis_pred <- FALSE
#############

load_rdata <- function(filepath) {
        # load in CONTENT R objects into independent envs so more
        # than one gene can be parsed at once
        env <- new.env()  # Create an isolated environment
        loaded_names <- load(filepath, envir = env)  # Load into the new environment
        obj_list <- mget(loaded_names, envir = env)  # Extract objects as a list
        return(obj_list)
}

get_content_preds <- function(gene1, gene2, scenario, work_dir, content_dir) {
        # TODO: remove hardcoding for formatted_expresison, really should just take in path
        # gene1 will always be regulator preds
        gene1_content_path <- file.path(work_dir, scenario, content_dir, "regulator",  paste0(gene1, "_crossval_predictors"))
        # gene2 will always be target preds
        gene2_content_path <- file.path(work_dir, scenario, content_dir, "target",  paste0(gene2, "_crossval_predictors"))

        gene1_content_preds <- load_rdata(gene1_content_path)
        gene2_content_preds <- load_rdata(gene2_content_path)

        return(list(gene1 = gene1_content_preds, gene2 = gene2_content_preds))
}

# setup pvalue matricies for target and regulator cis-predicted expression
# pvalue matrix for cis-genetic predicted target associations with simulated trans expression
regulator_assoc_pvalues <- array(0, dim = c(num_sim_genes, num_sim_genes, num_contexts))
# pvalue matrix for cis-genetic predicted regulator associations with simulated trans expression
target_assoc_pvalues <- array(0, dim = c(num_sim_genes, num_sim_genes, num_contexts))

trans_exp_data <- list()
for (gene1 in 1:num_sim_genes) {
  trans_exp_data[[gene1]] <- list()
  for (context in 0:(num_contexts - 1)) {
    context_filename <- paste0(context, "_cis.txt")
    trans_exp_data[[gene1]][[context_filename]] <- fread(file.path(work_dir, scenario, gene1, simulated_trans_exp, context_filename), header=F, data.table=F)[,-1]
  }
}

# Run parallel processing for gene pairs
results <- foreach(gene1 = 1:num_sim_genes, .combine = 'rbind', .packages = c("data.table")) %:%
            foreach(gene2 = 1:num_sim_genes, .combine = 'rbind') %dopar% {
                
                # Get predicted expression values
		print(paste("GENE 1: ", gene1))
		print(paste("GENE 2: ", gene2))
                content_preds <- get_content_preds(gene1, gene2, scenario, work_dir, content_dir)
                regulator_yhats <- content_preds$gene1
                target_yhats <- content_preds$gene2

                regulator_exp <- if (method == "CSTEM") regulator_yhats$Yhats_full else regulator_yhats$Yhats_tiss
                target_exp <- if (method == "CSTEM") target_yhats$Yhats_full else target_yhats$Yhats_tiss

                # Process each context in parallel
                foreach(context = 0:(num_contexts - 1), .combine = 'rbind') %dopar% {
		    print(paste("CONTEXT: ", context))
                    context_filename <- paste0(context, "_cis.txt")
                    regulator_exp_vec <- regulator_exp[[context_filename]]
                    target_exp_vec <- target_exp[[context_filename]]
                    
                    # Load pre-stored trans expression data
                    trans_exp_df_vec <- trans_exp_data[[gene1]][[context_filename]]

                    if (target_cis_pred) {
                        # Use matrix regression for speed
                        #trans_model <- lm.fit(cbind(1, regulator_exp_vec, target_exp_vec), trans_exp_df_vec)
			trans_model <- lm(trans_exp_df_vec ~ regulator_exp_vec + target_exp_vec)
                        summary_model <- summary(trans_model)
                        regulator_pvalue <- summary_model$coefficients[2, 4]
                        target_pvalue <- summary_model$coefficients[3, 4]
                        return(list(gene1, gene2, context + 1, regulator_pvalue, target_pvalue))
                    } else {
                        #trans_model <- lm.fit(cbind(1, regulator_exp_vec), trans_exp_df_vec)
			trans_model <- lm(trans_exp_df_vec ~ regulator_exp_vec)
                        summary_model <- summary(trans_model2)
                        regulator_pvalue <- summary_model$coefficients[2, 4]
                        return(list(gene1, gene2, context + 1, regulator_pvalue, NA))
                    }
                }
            }

print(str(r))
print(r)
# Fill result matrices from parallel results
for (r in results) {
    regulator_assoc_pvalues[r[[1]], r[[2]], r[[3]]] <- r[[4]]
    if (!is.na(r[[5]])) {
        target_assoc_pvalues[r[[1]], r[[2]], r[[3]]] <- r[[5]]
    }
}

# Clean up parallel workers
stopCluster(cl)
