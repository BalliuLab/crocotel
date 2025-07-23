 

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

#' @export
concat_crocotel_lmm_files <- function(directory = ".", regress_target_GReX = TRUE) {
  bash_script <- sprintf('
    cd "%s"
    cd "crocotel_lmm_output/"
    tmp_outdir="tmp_files/"
    mkdir -p "$tmp_outdir"
    file_suffix="crocotel_lmm.txt"
    if [ "$regress_target_GReX" = true ]; then
      file_suffix=".crocotel_lmm_regress.txt"
    fi

    #for prefix in $(ls *${file_suffix} | cut -d. -f1 | sort -u); do
    for prefix in $(find . -maxdepth 1 -type f -name "*${file_suffix}" | sed "s|^\\./||" | cut -d. -f1 | sort -u); do
      out_file="${prefix}.${file_suffix}"
      tmp_merged="${tmp_outdir}${out_file}"
      first=1
      for file in ${prefix}.*${file_suffix}; do
        if [ $first -eq 1 ]; then
          cat "$file" > "$tmp_merged"
          first=0
        else
          tail -n +2 "$file" >> "$tmp_merged"
        fi
        #rm "$file"
      done

      # Sort by 6th column (p-value) ascending, keeping header
      header=$(head -n 1 "$tmp_merged")
      tail -n +2 "$tmp_merged" | sort -k6,6g > "${tmp_merged}.sorted"
      echo "$header" | cat - "${tmp_merged}.sorted" > "${out_file}"

      rm "$tmp_merged" "${tmp_merged}.sorted"
      echo "Wrote $out_file"
    done

    rmdir "$tmp_outdir"
  ', normalizePath(directory, mustWork = TRUE))
  
  system(bash_script)
}





