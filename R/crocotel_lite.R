

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

format_GReX_for_association = function(GReX_dir, context, r2_genes, tmp_dir){
  plan(multisession, workers = parallel::detectCores() - 1)
  
  ### change this after changing other functions
  files = list.files(GReX_dir, pattern = "GReX")
  
  # Function to process one file
  process_file <- function(f, GReX_dir, context) {
    gene_id <- sub("\\.crocotel.*", "", basename(f))
    dt <- fread(paste0(GReX_dir, f), sep = "\t", header = T)
    
    result <- list()
    ct = context
    if (ct %in% names(dt)) {
      sub_dt <- dt %>% select(id, value = get("ct"))
      sub_dt[, gene := gene_id]
      casted <- dcast(sub_dt, gene ~ id, value.var = "value")
      result[[ct]] <- casted
    }
    result
  }
  
  # Process all files in parallel
  all_results <- future_lapply(files, process_file, GReX_dir, context)
  
  # Reorganize by cell type

  ct = context
  ct_list <- lapply(all_results, function(res) res[[ct]])
  ct_combined <- rbindlist(ct_list, fill = TRUE)
  
  if(!is.null(r2_genes)){
    if(nrow(ct_combined != 0)){
      ct_combined = ct_combined %>% filter(gene %in% r2_genes[[context]])
    }
  }
  
  fwrite(data.frame(ct_combined), paste0(tmp_dir, context, ".txt"), sep = "\t", na = "NA", quote = F)

}

get_genes_passing_r2 = function(GReX_dir, r2_thresh){
  all_r2_files = list.files(GReX_dir, pattern = "crocotel.crossval_r2", full.names = T)
  context_gene_list = list()
  
  for(file in all_r2_files){
    regulator_r2 = fread(file, sep = "\t", data.table = F, check.names = F, header = T)
    gene_id = sub("\\..*", "", basename(file))
    selected_contexts = regulator_r2 %>% filter(full_cv_r2s > r2_thresh)
    for(context in selected_contexts$context){
      context_gene_list[[context]] = c(context_gene_list[[context]], gene_id)
    }
  }
  return(context_gene_list)
}

#' @export
crocotel_lite = function(context, geneloc_file, out_dir, exp_files = NULL, GReX_dir = NULL, regress_target_GReX = T, pval_thresh = 1, r2_thresh = NULL){
  out_dir_crocotel_lite = paste0(out_dir, "/crocotel_lite_output/")
  dir.create(out_dir_crocotel_lite, showWarnings = F)
  ## create temp dir to store input matrixEQTL files
  tmp_dir = paste0(out_dir_crocotel_lite, "/MEQTL_input/")
  dir.create(tmp_dir, showWarnings = F)
  
  if(!is.null(r2_thresh)){
    r2_genes = get_genes_passing_r2(GReX_dir, r2_thresh)
  }else{
    r2_genes = NULL
  }
  
  #### write out formatted GReX files in MatrixEQTL format
  if(is.null(GReX_dir)){
    message("inferring GReX directory...")
    GReX_dir = paste0(out_dir, "/GReXs/")
  }
  if(is.null(exp_files)){
    if(regress_target_GReX == T){
      message("inferring expression files with residualized target GReX and formatting.")
      tmp_dir_regressed = paste0(out_dir_crocotel_lite, "/MEQTL_input/regressed_exp/")
      dir.create(tmp_dir_regressed, showWarnings = F)
      format_GReX_for_association(paste0(out_dir, "/exp_residualized_GReX/"), context, NULL, tmp_dir_regressed)
    }else{
      stop("You are running crocotel lite without regressing out GReX for each target. Please specify a vector of expression files in the `exp_files` parameter.")
    }
  }
  format_GReX_for_association(GReX_dir, context, r2_genes, tmp_dir)
    
  
  ### unchanging parameters across contexts
  geneloc = fread(geneloc_file, sep = "\t", data.table = F)
  snpsloc = geneloc[,c(1:3)]
  names(snpsloc) = c("SNP", "chr", "pos")
  # Only associations significant at this level will be saved
  pvOutputThreshold_tra = pval_thresh;
  # Distance for local gene-SNP pairs
  cisDist = 1e6; ## hard coded cis distance of 1MB
  file_prefix = ".crocotel_lite.txt"
  # Covariates file name
  # Set to character() for no covariates
  covariates_file_name = character()
  # Error covariance matrix
  # Set to numeric() for identity.
  errorCovariance = numeric();
  useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  
  # Genotype file name
  SNP_file_name = paste0(tmp_dir, context, ".txt");
  
  # Gene expression file name
  if(!is.null(exp_files)){
    pattern <- paste0("(^|[^A-Za-z0-9])", context, "([^A-Za-z0-9]|$)")
    # Match against basenames
    matching_file <- exp_files[grepl(pattern, basename(exp_files))]
    # Check for match
    if (length(matching_file) == 0) {
      stop(paste0("No file found matching context exactly: ", context))
    }
    expression_file_name <- matching_file[1]
  }else{
    expression_file_name = paste0(tmp_dir_regressed, context, ".txt");
  }
  
  
  
  ## Load genotype data
  
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);
  
  ## Load gene expression data
  
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }
  
  # Output file name
  output_file_name_cis = tempfile();
  output_file_name_tra = tempfile();
  
  
  ## Run the analysis
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    pvOutputThreshold.cis = 0,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    snpspos = snpsloc,
    genepos = geneloc,
    cisDist = cisDist,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  output = me$all$eqtls
  output = output %>% filter(snps != gene) %>% rename(regulator = "snps", target = "gene")
  outfile = paste0(out_dir_crocotel_lite, context, file_prefix)
  fwrite(output, file = outfile, sep = "\t", quote = F)
  print(paste0("finished analysis association mapping for context ", context))
  
  ## remove uneeded directories
  unlink(tmp_dir, recursive = T)
  unlink(tmp_dir_regressed, recursive = T)
}




