
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
  plan(multisession, workers = 1)
  
  ### change this after changing other functions
  files = list.files(GReX_dir, pattern = "GReX")
  
  # Function to process one file
  process_file <- function(f, GReX_dir, context) {
    gene_id <- sub("\\.crocotel.*", "", basename(f))
    dt <- fread(paste0(GReX_dir, f), sep = "\t", header = T, check.names = F)
    
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
  
  fwrite(data.frame(ct_combined, check.names = F), paste0(tmp_dir, context, ".txt"), sep = "\t", na = "NA", quote = F)

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
crocotel_lite = function(context, geneloc_file, out_dir, exp_files = NULL, GReX_dir = NULL, regress_target_GReX = T, pval_thresh = 1, r2_thresh = NULL, cisDist = 1e6){
  out_dir_crocotel_lite = paste0(out_dir, "/crocotel_lite_output/")
  dir.create(out_dir_crocotel_lite, showWarnings = F)
  ## create temp dir to store input matrixEQTL files
  #tmp_dir = paste0(out_dir_crocotel_lite, "/MEQTL_input/")
  tmp_dir = paste0(tempfile(tmpdir = out_dir_crocotel_lite), "/")
  dir.create(tmp_dir, showWarnings = F)
  
  #### write out formatted GReX files in MatrixEQTL format
  if(is.null(GReX_dir)){
    message("inferring GReX directory...")
    GReX_dir = paste0(out_dir, "/GReXs/")
  }
  if(is.null(exp_files)){
    if(regress_target_GReX == T){
      message("inferring expression files with residualized target GReX and formatting.")
      #tmp_dir_regressed = paste0(out_dir_crocotel_lite, "/MEQTL_input/regressed_exp/")
      tmp_dir_regressed = paste0(tmp_dir, "/regressed_exp/")
      dir.create(tmp_dir_regressed, showWarnings = F)
      format_GReX_for_association(paste0(out_dir, "/exp_residualized_GReX/"), context, NULL, tmp_dir_regressed)
    }else{
      stop("You are running crocotel lite without regressing out GReX for each target. Please specify a vector of expression files in the `exp_files` parameter.")
    }
  }
  
  if(!is.null(r2_thresh)){
    r2_genes = get_genes_passing_r2(GReX_dir, r2_thresh)
  }else{
    r2_genes = NULL
  }
  
  format_GReX_for_association(GReX_dir, context, r2_genes, tmp_dir)
    
  
  ### unchanging parameters across contexts
  geneloc = fread(geneloc_file, sep = "\t", data.table = F)
  snpsloc = geneloc[,c(1:3)]
  names(snpsloc) = c("SNP", "chr", "pos")
  # Only associations significant at this level will be saved
  pvOutputThreshold_tra = pval_thresh;
  # Distance for local gene-SNP pairs
  #cisDist = 1e6; ## hard coded cis distance of 1MB
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
  genos = fread(SNP_file_name, sep = "\t", data.table = F)
  
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
  
  ## Load gene expression data
  
  gene_mat = fread(expression_file_name, sep = "\t", data.table = F)
  gene_mat_formatted = as.matrix(gene_mat[,-1])
  rownames(gene_mat_formatted) = gene_mat$id
  ### remove individuals with all NAs
  gene_mat_formatted = data.frame(gene_mat_formatted, check.names = F) %>% select_if(~ !any(is.na(.)))
  ### remove same individuals from genotype matrix
  genos_formatted = genos[,colnames(gene_mat_formatted)]
  rownames(geos_formatted) = genos$id
  
  gene = SlicedData$new();
  gene$CreateFromMatrix(as.matrix(gene_mat_formatted))
  
  ## Load genotype data
  snps = SlicedData$new();
  snps$CreateFromMatrix(as.matrix(genos_formatted))
  ## Load genotype data
  
  #snps = SlicedData$new();
  #snps$fileDelimiter = "\t";      # the TAB character
  #snps$fileOmitCharacters = "NA"; # denote missing values;
  #snps$fileSkipRows = 1;          # one row of column labels
  #snps$fileSkipColumns = 1;       # one column of row labels
  #snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  #snps$LoadFile(SNP_file_name);
  
  ## Load gene expression data
  
  #gene = SlicedData$new();
  #gene$fileDelimiter = "\t";      # the TAB character
  #gene$fileOmitCharacters = "NA"; # denote missing values;
  #gene$fileSkipRows = 1;          # one row of column labels
  #gene$fileSkipColumns = 1;       # one column of row labels
  #gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  #gene$LoadFile(expression_file_name);
  
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
  
 
  ## make sure individual IDS are the same
  snps$ColumnSubsample(match(gene$columnNames, snps$columnNames))
  
  ## Run the analysis
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = 1,
    snpspos = snpsloc,
    genepos = geneloc,
    cisDist = cisDist,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  unlink(output_file_name_tra);
  unlink(output_file_name_cis);
  
  output = me$trans$eqtls
  output = output %>% filter(snps != gene) %>% rename(regulator = "snps", target = "gene")
  outfile = paste0(out_dir_crocotel_lite, context, file_prefix)
  fwrite(output, file = outfile, sep = "\t", quote = F)
  print(paste0("finished analysis association mapping for context ", context))
  
  ## remove uneeded directories
  unlink(tmp_dir, recursive = T)
  if (dir.exists(tmp_dir_regressed)) {
    unlink(tmp_dir_regressed, recursive = TRUE)
  }
}




