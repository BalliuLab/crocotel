
#' Test one regulator-target pair across contexts with a linear mixed model
#'
#' Fits \code{target ~ regulator_GReX * context + (1 | id)} and reports the
#' per-context marginal trend of the regulator's GReX, writing one file per context.
#' Falls back to per-context \code{lm} when fewer than 3 contexts are shared by the
#' regulator and target. Per-pair output is merged by
#' \code{\link{concat_crocotel_lmm_files}}, which must be called with a matching
#' \code{file_suffix}.
#'
#' @param regulator_gene_name regulator gene identifier
#' @param target_gene_name target gene identifier
#' @param out_dir crocotel output directory; inputs are inferred from its
#'   conventional subdirectories when file paths are not given
#' @param target_exp_file optional wide-format target expression file (\code{id} plus
#'   one column per context, with a header). When supplied it is used as-is and
#'   \code{regress_target_GReX} has no effect on which expression is read.
#' @param GReX_dir optional directory of GReX predictor files; defaults to
#'   \code{out_dir/GReXs/}
#' @param regress_target_GReX controls which target expression is inferred when
#'   \code{target_exp_file} is NULL. \code{TRUE} reads the GReX-residualized
#'   expression from \code{out_dir/exp_residualized_GReX/}, so the target's own cis
#'   genetics are controlled for. \code{FALSE} (default) reads the unresidualized
#'   per-context expression written by \code{format_data} to
#'   \code{out_dir/crocotel_formatted_data/<target_gene_name>/}; note that leaves the
#'   target's cis genetics uncontrolled.
#' @param pval_thresh only associations with p <= this value are written. Default 1
#'   (write everything). Thresholding is safe for downstream treeQTL correction
#'   because \code{record_n_tests_per_gene} supplies design-based family sizes.
#' @param r2_thresh optional minimum cross-validated R2 for the regulator's GReX in at
#'   least one context; the pair is skipped if unmet
#' @param context_dependence if TRUE, also test the regulator-by-context interaction
#'   and write its LRT p-value and delta R2 to
#'   \code{out_dir/crocotel_lmm_output/context_dependence/}
#' @param file_suffix suffix for the per-pair output files. Must match the
#'   \code{file_suffix} passed to \code{concat_crocotel_lmm_files}.
#' @export
crocotel_lmm = function(regulator_gene_name, target_gene_name, out_dir, target_exp_file = NULL, GReX_dir = NULL, regress_target_GReX = F, pval_thresh = 1, r2_thresh = NULL, context_dependence = F, file_suffix = ".crocotel_lmm.txt"){
  out_dir_crocotel_lmm = paste0(out_dir, "/crocotel_lmm_output/")
  dir.create(out_dir_crocotel_lmm, showWarnings = F)
  
  if(is.null(GReX_dir)){
    message("inferring GReX directory and regulator GReX file...")
    GReX_dir = paste0(out_dir, "/GReXs/")
    regulator_pred_exp_file = list.files(GReX_dir, pattern = paste0(regulator_gene_name, ".crocotel.GReX_predictors.txt"), full.names = T)
  }
  ## resolve which target expression to read. An explicit target_exp_file is always
  ## read as wide format (id + one column per context, with a header). Otherwise
  ## regress_target_GReX picks between the GReX-residualized wide file and
  ## format_data's per-context files for this gene.
  target_exp_files = NULL
  if(is.null(target_exp_file)){
    if(regress_target_GReX == T){
      message("inferring expression file with residualized target GReX...")
      target_exp_file = list.files(paste0(out_dir, "/exp_residualized_GReX/"), pattern = paste0(target_gene_name, ".crocotel.GReX_residuals.txt"), full.names = T)
      if(length(target_exp_file) == 0){
        stop("could not infer residualized expression for target ", target_gene_name, " in ",
             out_dir, "/exp_residualized_GReX/. Run create_GReXs(regress_tGReX = T), set ",
             "regress_target_GReX = F, or pass target_exp_file explicitly.")
      }
    }else{
      exp_dir = paste0(out_dir, "/crocotel_formatted_data/", target_gene_name, "/")
      target_exp_files = list.files(exp_dir, full.names = T)
      if(length(target_exp_files) == 0){
        stop("could not infer unresidualized expression for target ", target_gene_name, " in ",
             exp_dir, ". Run format_data(), set regress_target_GReX = T, or pass ",
             "target_exp_file explicitly.")
      }
      message("inferring unresidualized per-context target expression from ", exp_dir)
    }
  }

  ## get target expression across all contexts, in long format
  if(is.null(target_exp_file)){
    ## per-context files from format_data(): headerless, columns are (id, expression),
    ## context taken from the filename as elsewhere in the package
    target_exp_mat = bind_rows(lapply(target_exp_files, function(cur_file){
      context = sub("\\..*", "", basename(cur_file))
      cur_df = fread(cur_file, header = F, sep = "\t", data.table = F)
      if(ncol(cur_df) < 2){
        stop("expected 2 columns (id, expression) in ", cur_file, " but found ", ncol(cur_df))
      }
      data.frame(id = cur_df[[1]], context = context, target_exp = cur_df[[2]],
                 stringsAsFactors = F, check.names = F)
    }))
  }else{
    target_output = fread(target_exp_file, sep = "\t", data.table = F, header = T)
    target_exp_mat = pivot_longer(target_output, cols = -id, names_to = "context", values_to = "target_exp")
  }
  
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
    r2 = max(regulator_r2$full, na.rm = T) < r2_thresh
    if(r2){
      print("Regulator GReX did not pass specified R2 threshold in any context. Not running Crocotel for this regulator-target pair.")
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
      
      lm_df = target_exp_vec %>%
        full_join(regulator_exp_vec, by = "id")
      trans_model <- lm(lm_df$target_exp ~ lm_df$regulator_pred)
      summary_model <- summary(trans_model)
      regulator_pvalue <- summary_model$coefficients[2, 4]
      regulator_beta = summary_model$coefficients[2, 1]
      regulator_se = summary_model$coefficients[2,2]
      regulator_statistic = regulator_beta/regulator_se
      regulator_FDR = NA
      #df = data.frame(SNP = target_gene_name, gene = regulator_gene_name, beta = regulator_beta, 'se' = regulator_se, 'pvalue' = regulator_pvalue, FDR = NA, context = context_name)
      df = data.frame(regulator = regulator_gene_name, target = target_gene_name, beta = regulator_beta, 'se' = regulator_se, 'statistic' = regulator_statistic, 'pvalue' = regulator_pvalue, 'FDR' = regulator_FDR, context = context_name)
    }))
    
    if(pval_thresh < 1) output_df = output_df %>% filter(pvalue <= pval_thresh)
    for(cur_context in unique(output_df$context)){
      context_df = output_df %>% filter(context == cur_context) %>% select(-context)
      fwrite(context_df, file = paste0(out_dir_crocotel_lmm, cur_context, ".", regulator_gene_name, "_", target_gene_name, file_suffix),
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
      ## written to a subdirectory: format_treeQTL scans every non-directory file in
      ## crocotel_lmm_output/ and would read this one as a context named <reg>_<tar>
      cd_dir = paste0(out_dir_crocotel_lmm, "context_dependence/")
      dir.create(cd_dir, showWarnings = F)
      fwrite(context_dependence_df, file = paste0(cd_dir, regulator_gene_name, "_", target_gene_name, ".reg_GReX_by_context_interaction_r2", file_suffix),  sep = "\t")
    }
    # extract marginal trends for each predicted exp
    #reg_marginal_trends = suppressMessages(suppressWarnings(emtrends(trans_model, ~ context, var = "regulator_pred", data = trans_model_df))) %>% tidy() %>% mutate(FDR = NA) %>% select(regulator_pred.trend, std.error, p.value, FDR, context)
    #names(reg_marginal_trends) = c("beta", "se", "pvalue", "FDR", "context")
    reg_marginal_trends = suppressMessages(suppressWarnings(emtrends(trans_model, ~ context, var = "regulator_pred", data = trans_model_df))) %>% tidy() %>% select(regulator_pred.trend, std.error, p.value, context)
    names(reg_marginal_trends) = c("beta", "se", "pvalue", "context")
    reg_marginal_trends$statistic = reg_marginal_trends$beta/reg_marginal_trends$se
    reg_marginal_trends$FDR = NA
    reg_marginal_trends = reg_marginal_trends %>% select(beta, se, statistic, pvalue, FDR, context)
    output_df = cbind(regulator = regulator_gene_name, target = target_gene_name, reg_marginal_trends)
    if(pval_thresh < 1) output_df = output_df %>% filter(pvalue <= pval_thresh)
    contexts = unique(output_df$context)
    for(cur_context in contexts){
      context_df = output_df %>% filter(context == cur_context) %>% select(-context)
      fwrite(context_df, file = paste0(out_dir_crocotel_lmm, cur_context, ".", regulator_gene_name, "_", target_gene_name, file_suffix),
             sep = "\t", quote = F)
    }
  }
  print("Finished running crocotel lmm for this pair.")
}

#' Record the total number of trans tests per gene for crocotel_lmm
#'
#' crocotel_lmm tests one regulator-target pair per call, so no single call sees
#' the whole design. Run this once (e.g. alongside \code{concat_crocotel_lmm_files})
#' to record, per context and per family role, the total number of trans tests each
#' gene participates in over the full regulator x target trans grid. These counts are
#' independent of any p-value threshold applied to stored output, so treeQTL family
#' sizes stay correct. Files are written in the same format treeQTL consumes (CSV;
#' columns \code{family, fam_p}) under \code{out_dir/<method>_output/n_tests_per_gene/},
#' one per context and role (\code{<context>.R.n_tests_per_gene.txt},
#' \code{<context>.T.n_tests_per_gene.txt}).
#'
#' Per-context gene membership for regulators is inferred from which contexts each
#' gene's GReX predictor file contains, read cheaply from file headers (one line per
#' file). For targets it follows \code{regress_target_GReX}, matching whichever
#' expression crocotel_lmm read. If \code{r2_thresh} is given, regulators are filtered
#' exactly as crocotel_lmm does: kept only if their maximum full-model cross-validated
#' R2 across contexts is >= r2_thresh (this reads each gene's small crossval_r2 file,
#' which also supplies its contexts). This is a one-time O(#genes) pass of light reads.
#'
#' @param out_dir crocotel output directory (reads out_dir/GReXs/ and, depending on
#'   regress_target_GReX, out_dir/exp_residualized_GReX/ or out_dir/crocotel_formatted_data/)
#' @param geneloc_file gene location file: geneid, chr, start, end
#' @param dist trans distance (bp); pairs within this window on the same chromosome are cis and not counted. Default 1e6.
#' @param method output subdirectory to write into ("crocotel_lmm" or "crocotel_lite"). Default "crocotel_lmm".
#' @param r2_thresh optional regulator R2 threshold, matching crocotel_lmm's filter (keep regulator if max full R2 across contexts >= r2_thresh). Default NULL (no filtering).
#' @param regress_target_GReX source for the target gene/context set; must match the
#'   value passed to \code{crocotel_lmm}. TRUE reads \code{out_dir/exp_residualized_GReX/}
#'   (contexts from each file's header), FALSE (default) reads the per-gene directories
#'   under \code{out_dir/crocotel_formatted_data/} (contexts from filenames). If this
#'   disagrees with the value used for the association run, the recorded family sizes
#'   describe a different target set than the one actually tested and treeQTL's
#'   correction is silently miscalibrated.
#' @return invisibly, the vector of contexts written
#' @export
record_n_tests_per_gene = function(out_dir, geneloc_file, dist = 1e6, method = "crocotel_lmm", r2_thresh = NULL, regress_target_GReX = F){
  nt_dir = paste0(out_dir, "/", method, "_output/n_tests_per_gene/")
  dir.create(nt_dir, showWarnings = FALSE, recursive = TRUE)

  geneloc = fread(geneloc_file, sep = "\t", data.table = F)
  GReX_dir = paste0(out_dir, "/GReXs/")
  exp_residualized_dir = paste0(out_dir, "/exp_residualized_GReX/")

  reg_files = list.files(GReX_dir, pattern = "\\.crocotel\\.GReX_predictors\\.txt$", full.names = TRUE)

  header_contexts = function(f) setdiff(strsplit(readLines(f, n = 1L), "\t", fixed = TRUE)[[1]], "id")

  ## regulators: contexts present per gene (+ max R2 if filtering).
  reg_by_ctx = list()
  if(is.null(r2_thresh)){
    if(length(reg_files) == 0) stop("No regulator GReX predictor files found in ", GReX_dir)
    for(f in reg_files){
      g = strip_suffix(f, ".crocotel.GReX_predictors.txt")
      for(ctx in header_contexts(f)) reg_by_ctx[[ctx]] = c(reg_by_ctx[[ctx]], g)
    }
  } else {
    ## the crossval_r2 file gives both the gene's contexts and its R2 in one read
    r2_files = list.files(GReX_dir, pattern = "\\.crocotel\\.crossval_r2\\.txt$", full.names = TRUE)
    if(length(r2_files) == 0) stop("No crossval_r2 files found in ", GReX_dir)
    for(f in r2_files){
      g = strip_suffix(f, ".crocotel.crossval_r2.txt")
      d = fread(f, sep = "\t", data.table = F)
      if(!length(d$full) || max(d$full, na.rm = TRUE) < r2_thresh) next   # regulator fails threshold
      for(ctx in as.character(d$context)) reg_by_ctx[[ctx]] = c(reg_by_ctx[[ctx]], g)
    }
  }

  ## targets: contexts present per gene (not R2-filtered). The target set must match
  ## the expression crocotel_lmm actually read, so it follows the same switch.
  tar_by_ctx = list()
  if(regress_target_GReX){
    tar_files = list.files(exp_residualized_dir, pattern = "\\.crocotel\\.GReX_residuals\\.txt$", full.names = TRUE)
    if(length(tar_files) == 0) stop("No residualized target expression files found in ", exp_residualized_dir)
    for(f in tar_files){
      g = strip_suffix(f, ".crocotel.GReX_residuals.txt")
      for(ctx in header_contexts(f)) tar_by_ctx[[ctx]] = c(tar_by_ctx[[ctx]], g)
    }
  }else{
    formatted_dir = paste0(out_dir, "/crocotel_formatted_data/")
    tar_dirs = list.dirs(formatted_dir, full.names = TRUE, recursive = FALSE)
    if(length(tar_dirs) == 0) stop("No per-gene expression directories found in ", formatted_dir)
    for(d in tar_dirs){
      g = basename(d)
      for(ctx in sub("\\..*", "", basename(list.files(d)))) tar_by_ctx[[ctx]] = c(tar_by_ctx[[ctx]], g)
    }
  }

  contexts = intersect(names(reg_by_ctx), names(tar_by_ctx))
  message("Recording n_tests_per_gene for contexts: ", paste(contexts, collapse = ", "))

  for(context in contexts){
    nt = count_trans_tests(reg_by_ctx[[context]], tar_by_ctx[[context]], geneloc, dist)
    write_n_tests_per_gene(nt, nt_dir, context)
  }
  message("Wrote n_tests_per_gene files to ", nt_dir)
  invisible(contexts)
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







