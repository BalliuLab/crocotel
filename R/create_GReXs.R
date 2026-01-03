

regress_target_GReX = function(gene_name, exp_files, Yhats_tiss, outdir){
  dir.create(paste0(outdir, "exp_residualized_GReX/"), showWarnings = F)
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
    
    model <- lm(expr_sub ~ grex_sub, na.action = na.exclude)
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
  outfile = paste0(outdir, "exp_residualized_GReX/", gene_name, ".crocotel.GReX_residuals.txt")
  fwrite(all_residuals, outfile, sep = "\t", quote = F)
}


#' Creates GReXs (Genetically regulated predictors of expression) for one gene across contexts
#'
#' @param gene_name - identifier for current gene being run. Add this prefix to all the saved results files... necessary to distinguish results where more than one gene-analysis is run.
#' @param out_dir - output directory for GReXs
#' @param genotype_file - a genotype file with the set of individuals in all Y_files and their corresponding cis-SNPs. Individuals are row names, no column names
#' @param exp_files - vector with list of all expression files per context. The individuals (row names) must be a subset of genotype (genotype_file). Contains an unnamed
#' @param context_thresh - minimum number of contexts to run crocotel on. 
#' @param alpha - The regularization constant. Default is .5 (eNet). Minimum value of 1e-4.
#' @param num_folds - Number of folds for cross-validation
#' @param method - Takes values of "crocotel" or "cxc". Default is "crocotel", if "cxc" then will not run decomposition and will build GReXs in a context by context manner.
#' @return writes out a file of predicted expression across individuals and contexts 
#' @export
create_GReXs = function(gene_name, out_dir, genotype_file = NULL, exp_files = NULL, context_thresh = 2, alpha = 0.5, num_folds = 10, impute = F, run_cxc = F, save_decomposition = F){
  seed = 9000
  set.seed(seed)
  GReX_outdir = paste0(out_dir, "GReXs/")
  dir.create(GReX_outdir, showWarnings = F)
  message("Saving cross-validated predictors and performance metrics in ", GReX_outdir)
  
  if(is.null(genotype_file)){
    message("inferring genotype file...")
    genotype_file = paste0(out_dir, "/crocotel_formatted_data/",gene_name,"_genotypes.txt")
    if(!file.exists(genotype_file)){
      stop(paste0("genotype file not specified and inferred genotype file: ", genotype_file ," does not exit. exiting."))
    }
  }
  if(is.null(exp_files)){
    message("inferring expression files...")
    expression_directory=paste0(out_dir, "/crocotel_formatted_data/",gene_name,"/")
    exp_files = list.files(expression_directory, full.names = T)
    if(!file.exists(exp_files[1])){
      stop(paste0("expression files not specified and inferred expression file: ", exp_files[1] ," does not exit. exiting."))
    }
  }
  
  message("Reading in genotype and expression files...") 
  suppressWarnings(expr = {X<-fread(file = genotype_file, sep='\t', data.table=F, check.names = F)})
  X<-as.matrix(data.frame(X, row.names=1, check.names = F))
  
  ### read in expression files
  total_exp=bind_rows(lapply(exp_files, function(cur_file){
    context = basename(cur_file)
    context = sub("\\..*", "", context)
    context = sub("*.txt", "", context)
    context = sub("*.tsv", "", context)
    cur_df = data.frame(fread(cur_file, header = F, sep = "\t"), check.names = T, stringsAsFactors = F)
    names(cur_df) = c("id", gene_name)
    final_df = data.frame(context = context, cur_df, check.names = F)
    final_df
  }))
  
  total_exp_mat = total_exp %>% pivot_wider(names_from = "context", values_from = gene_name) %>% as.data.frame()
  
  # check if imputation is needed
  if(impute){
    if(sum(is.na(total_exp_mat)) == 0){
      message("imputation not needed")
      impute = F
    }else{
      run_cxc = T
    }
  }
  
  if(run_cxc){
    message("Running cxc")
    crossval_output_cxc = crossval_helper_parallel(total_exp_mat, total_exp_mat, X, GReX_outdir, gene_name, TRUE, num_folds, alpha)
  }
  
  ### decompose expression and run crocotel
  decomp_exp_mat = decompose_expression(total_exp, gene_name, context_thresh, if(save_decomposition) data_dir = GReX_outdir)
  message("Running crocotel")
  crossval_output = crossval_helper_parallel(total_exp_mat, decomp_exp_mat, X, GReX_outdir, gene_name, FALSE, num_folds, alpha)
  
  message("Finished computing GReXs and evalutation metrics")
  
  #### impute using cxc if flag is set
  if(impute){
    message("imputing crocotel GReX")
    crocotel_grex = crossval_output$full
    cxc_grex = crossval_output_cxc$context
    crocotel_grex[is.na(crocotel_grex)] = cxc_grex[is.na(crocotel_grex)]
    # write out imputed results
    method = "crocotel_imputed"
    fwrite(crocotel_grex,  file = paste0(GReX_outdir, gene_name,".", method, "_context.GReX_predictors.txt"), sep = "\t")
  }
  
}
