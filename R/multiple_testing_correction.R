 

get_top_pairs = function(crocotel_dir, top_level) {
  file_paths = list.files(crocotel_dir, pattern = "crocotel", full.names = TRUE)
  pair_min_pvals <- list()
  
  if (top_level == "R") {
    for (file in file_paths) {
      df <- fread(file, select = c("regulator", "target", "pvalue"))
      # Take min p-value per pair within this file
      df_min <- df[, .SD[which.min(pvalue)], by = regulator]
      pair_min_pvals[[file]] <- df_min
    }
    # Combine all minimal p-values per pair across contexts
    all_min <- rbindlist(pair_min_pvals)
    # Take the global minimum p-value per pair
    top_pairs = all_min[, .SD[which.min(pvalue)], by = regulator]
    top_pairs = top_pairs %>% mutate(reg_targ = paste0(regulator, ":", target)) %>% select(reg_targ, pvalue)
    return(top_pairs$reg_targ)
  } else if(top_level == "T"){
    for (file in file_paths) {
      df <- fread(file, select = c("regulator", "target", "pvalue"))
      # Take min p-value per pair within this file
      df_min <- df[, .SD[which.min(pvalue)], by = target]
      pair_min_pvals[[file]] <- df_min
    }
    # Combine all minimal p-values per pair across contexts
    all_min <- rbindlist(pair_min_pvals)
    # Take the global minimum p-value per pair
    top_pairs = all_min[, .SD[which.min(pvalue)], by = target]
    top_pairs = top_pairs %>% mutate(reg_targ = paste0(regulator, ":", target)) %>% select(reg_targ, pvalue)
    return(top_pairs$reg_targ)
  }else {
    stop("No valid input specified for target or regulator as top level.")
  }
}

format_crocotel_for_mash = function(crocotel_dir){
  files <- list.files(crocotel_dir, pattern = "crocotel", full.names = TRUE)
  
  # Initialize empty lists for betas and SEs
  beta_list <- list()
  se_list <- list()
  
  for (i in seq_along(files)) {
    df <- fread(files[i], select = c("regulator", "target", "statistic", "beta"))
    # Calculate SE
    df[, se := beta / statistic]
    df[, reg_targ := paste0(regulator, ":", target)]
    
    # Create a context-specific name
    context = sub("\\..*", "", basename(files[i]))
    print(paste0("processing file for current inferred context: ", context))
    
    # Store beta and se separately
    beta_list[[context]] <- df[, .(reg_targ, value = beta)]
    setnames(beta_list[[context]], "value", context)
    
    se_list[[context]] <- df[, .(reg_targ, value = se)]
    setnames(se_list[[context]], "value", context)
  }
  
  merge_all <- function(dt_list) {
    Reduce(function(x, y) merge(x, y, by = c("reg_targ"), all = TRUE), dt_list)
  }
  
  # Create final wide-format beta and se data.tables
  beta_total <- merge_all(beta_list)
  se_total <- merge_all(se_list)
  
  setkey(beta_total, reg_targ)
  setkey(se_total, reg_targ)
  
  # Optional: sanity check
  stopifnot(identical(beta_total$reg_targ, se_total$reg_targ))
  
  return(list(beta = beta_total, se = se_total))
}


#' @export
multiple_testing_correction = function(crocotel_dir, out_dir, fdr_thresh = 0.05, method = "treeQTL", top_level = "R"){
  method_outdir = paste0(out_dir, "/", method, "_output/")
  dir.create(method_outdir, showWarnings = F)
  exp_suffix = gsub("_output", "", basename(crocotel_dir))
  if(top_level == "R"){
    output_prefix = "eRegulators"
  }else if(top_level == "T"){
    output_prefix = "eTargets"
  }else{
    stop("No valid input specified for target or regulator as top level.")
  }
  if(method == "treeQTL"){
    level1 = fdr_thresh
    level2 = fdr_thresh
    level3 = fdr_thresh
    
    eGenes = get_eGenes_multi_tissue_mod(crocotel_dir = crocotel_dir, 
                                         exp_suffix = exp_suffix,
                                         out_dir = method_outdir,
                                         top_level = top_level,
                                         level1 = level1, level2 = level2, level3 = level3)
    fwrite(eGenes, file = paste0(method_outdir, output_prefix, ".", exp_suffix, ".txt"), sep = "\t")
  }
  if(method == "mashr"){
    data = format_crocotel_for_mash(crocotel_dir)
    betas = data.frame(data[["beta"]], check.names = F)
    rownames(betas) = betas$reg_targ
    betas = as.matrix(betas %>% select(-reg_targ))
    ses = data.frame(data[["se"]], check.names = F)
    rownames(ses) = ses$reg_targ
    ses = as.matrix(ses %>% select(-reg_targ))
    
    mash_data = mash_set_data(betas, ses)
    
    ### set up canonical covariance matrix
    U.c = cov_canonical(mash_data)  
    
    ### set up data driven covariance matrix
    pairs = get_top_pairs(crocotel_dir, top_level)
    indices = which(rownames(betas) %in% pairs)
    U.pca = cov_pca(mash_data,5,subset=indices)
    
    #### apply extreme deconvolution
    U.ed = cov_ed(mash_data, U.pca, subset=indices)
    
    #### run mash
    m = mash(mash_data, c(U.c,U.ed))
    sig_results = get_lfsr(m)
    sig_beta = get_pm(m) 
    sig_se = get_psd(m)
    
    ## pivot to long format
    sig_results = as.data.frame(sig_results) %>% mutate(pair = rownames(.)) %>%
      separate(pair, into = c("regulator", "target"), sep = ":") %>%
      pivot_longer(
        cols = -c(regulator, target),
        names_to = "context",
        values_to = "p.value"
      ) 
    beta_long <- as.data.frame(sig_beta) %>% mutate(pair = rownames(.)) %>%
      separate(pair, into = c("regulator", "target"), sep = ":") %>%
      pivot_longer(
        cols = -c(regulator, target),
        names_to = "context",
        values_to = "beta"
      )
    se_long <- as.data.frame(sig_se) %>% mutate(pair = rownames(.)) %>%
      separate(pair, into = c("regulator", "target"), sep = ":") %>%
      pivot_longer(
        cols = -c(regulator, target),
        names_to = "context",
        values_to = "se"
      )

    sig_results = sig_results %>%
      left_join(beta_long, by = c("regulator", "target", "context")) %>%
      left_join(se_long, by = c("regulator", "target", "context")) %>%
      filter(p.value <= 0.05) %>% select(regulator, target, beta, se, p.value) %>%
      arrange(p.value)
    
    fwrite(sig_results, file = paste0(method_outdir, output_prefix, ".", exp_suffix, ".txt"), sep = "\t")
  }
}

#' Merge per-pair crocotel_lmm output files for one context
#'
#' crocotel_lmm writes one file per regulator-target pair per context. This merges
#' them into a single p-value-sorted \code{<context><file_suffix>} in the same
#' directory. Per-pair inputs are deleted only after the merged file is written, and
#' the merged output is never itself consumed as an input.
#'
#' @param directory crocotel out_dir (the parent of the output subdirectory)
#' @param context context to merge
#' @param file_suffix suffix of the per-pair input files; must match the value used
#'   by crocotel_lmm. Default ".crocotel_lmm.txt".
#' @param to_concatenate optional file listing input files, one per line; if NULL,
#'   inputs are discovered in the output subdirectory
#' @param subdir subdirectory of \code{directory} to merge within
#' @param keep_inputs if TRUE, leave the per-pair files in place after merging
#' @param regress_target_GReX deprecated; only for merging output written by older
#'   versions, which named files ".crocotel_lmm_regress.txt". Use file_suffix instead.
#' @return invisibly, the path to the merged file
#' @export
concat_crocotel_lmm_files <- function(directory = ".", context,
                                      file_suffix = ".crocotel_lmm.txt",
                                      to_concatenate = NULL,
                                      subdir = "crocotel_lmm_output",
                                      keep_inputs = FALSE,
                                      regress_target_GReX = NULL) {
  if (!is.null(regress_target_GReX)) {
    warning("`regress_target_GReX` is deprecated; pass `file_suffix` instead.")
    if (isTRUE(regress_target_GReX)) file_suffix = ".crocotel_lmm_regress.txt"
  }
  if (missing(context) || length(context) != 1L || !nzchar(context))
    stop("`context` must be a single non-empty context name.")

  out_subdir = file.path(normalizePath(directory, mustWork = TRUE), subdir)
  if (!dir.exists(out_subdir)) stop("output subdirectory does not exist: ", out_subdir)
  out_file = file.path(out_subdir, paste0(context, file_suffix))

  ## ---- collect inputs -------------------------------------------------------
  if (!is.null(to_concatenate)) {
    if (!file.exists(to_concatenate)) stop("to_concatenate file not found: ", to_concatenate)
    files = readLines(to_concatenate)
    files = files[nzchar(files)]
    files = ifelse(startsWith(files, "/"), files, file.path(out_subdir, basename(files)))
  } else {
    nm = list.files(out_subdir)
    ## per-pair files are "<context>.<reg>_<tar><file_suffix>"; the merged output is
    ## "<context><file_suffix>", with nothing between, so the length test excludes it
    keep = startsWith(nm, paste0(context, ".")) &
           endsWith(nm, file_suffix) &
           nchar(nm) > nchar(context) + 1L + nchar(file_suffix)
    files = file.path(out_subdir, nm[keep])
  }
  files = setdiff(files[file.exists(files)], out_file)

  ## a killed parallel job can leave zero-byte files behind; skip rather than abort
  sizes = file.info(files)$size
  empty = is.na(sizes) | sizes == 0
  if (any(empty)) {
    warning(sprintf("skipping %d empty file(s), e.g. %s",
                    sum(empty), basename(files[which(empty)[1]])))
    files = files[!empty]
  }
  if (length(files) == 0L)
    stop("no per-pair files found for context '", context, "' with suffix '",
         file_suffix, "' in ", out_subdir)
  message(sprintf("merging %d per-pair files for context %s", length(files), context))

  ## ---- read under a column contract ----------------------------------------
  ## use.names = TRUE also aligns the regulator/target column order difference
  ## between crocotel_lmm's LMM path and its <3-context fallback path
  expected = c("regulator", "target", "beta", "se", "statistic", "pvalue", "FDR")
  merged = rbindlist(lapply(files, function(f) {
    d = fread(f, sep = "\t", header = TRUE, data.table = TRUE)
    if (!setequal(names(d), expected))
      stop("unexpected columns in ", f, ": ", paste(names(d), collapse = ", "),
           "\nexpected: ", paste(expected, collapse = ", "))
    d
  }), use.names = TRUE)

  setorder(merged, pvalue, na.last = TRUE)
  fwrite(merged, out_file, sep = "\t", quote = FALSE, na = "NA")

  ## ---- only now is it safe to drop the inputs -------------------------------
  if (!keep_inputs) {
    if (!file.exists(out_file)) stop("merged file was not written; inputs left in place")
    unlink(files)
  }
  message(sprintf("wrote %s (%d rows from %d pairs)", out_file, nrow(merged), length(files)))
  invisible(out_file)
}

# Legacy shell-based merge, preserved verbatim as a fallback for merges too large to
# hold in memory (sort(1) spills to disk where rbindlist cannot). NOT exported, and
# NOT currently usable as-is: the regress_target_GReX flag is never interpolated into
# the script, headers are not stripped before sorting, `uniq` can drop duplicate data
# rows, inputs are deleted before the merge is finalized, and `find -printf` is
# GNU-only. Fix those before relying on it.
concat_crocotel_lmm_files_shell <- function(directory = ".", context, regress_target_GReX = F, to_concatenate = NULL) {
  tmp_dir = paste0(tempfile(tmpdir = paste0(directory, "crocotel_lmm_output/")), "/")
  bash_script <- sprintf('
    cd "%s"
    cd "crocotel_lmm_output/"
    tmp_outdir="%s"
    mkdir -p "$tmp_outdir"
    echo "tmp directory created: $tmp_outdir"
    file_suffix="crocotel_lmm.txt"
    
    ##### have to fix this it does not do anything right now ###
    if [ "$regress_target_GReX" = true ]; then
      file_suffix=".crocotel_lmm_regress.txt"
    fi
    ############################################################
    
    context="%s"
    out_file="${context}.${file_suffix}"
    tmp_merged="${tmp_outdir}${out_file}"
    
    to_concatenate="%s"
    if [ -f "$to_concatenate" ]; then
      echo "valid input concatenation file. proceeding with concatenation"
    else
      echo "finding files to concatenate for context $context"
      find . -maxdepth 1 -type f -name "${context}.*${file_suffix}" -printf "%%p\n" > "${tmp_outdir}${context}_to_concatenate.txt"
      to_concatenate="${tmp_outdir}${context}_to_concatenate.txt"
    fi
    
    #xargs -a "$to_concatenate" awk \'FNR==1 && NR!=1 { next } { print }\' > $tmp_merged
    while IFS= read -r f; do
      if [ -f "$f" ]; then
        cat "$f" >> "$tmp_merged"
        rm -f "$f"                     
      fi  
    done < "$to_concatenate"

    # Sort by 6th column (p-value) ascending, keeping header
    cat "$tmp_merged" | sort -k6,6g | uniq > "${tmp_merged}.sorted"
    cat "${tmp_merged}.sorted" > "${out_file}"

    rm "$tmp_merged" "${tmp_merged}.sorted"
    echo "Wrote $out_file"
    echo "removing tmp directory"
    rmdir "$tmp_outdir"
  ', normalizePath(directory, mustWork = TRUE), tmp_dir, context, to_concatenate)
  
  system(bash_script)
}



