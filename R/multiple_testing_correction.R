 
#' @export
multiple_testing_correction = function(crocotel_sum_stats, contexts_vec, fdr_thresh, outdir, method = "treeQTL", top_level = "R"){
  dir.create(outdir, showWarnings = F)
  exp_suffix = strsplit(crocotel_sum_stats, "\\.")[[1]][2]
  if(method == "treeQTL"){
    level1 = fdr_thresh
    level2 = fdr_thresh
    level3 = fdr_thresh
    
    eGenes = get_eGenes_multi_tissue_mod(crocotel_sum_stats = crocotel_sum_stats, 
                                         contexts_vec = contexts_vec, 
                                         exp_suffix = exp_suffix,
                                         outdir = outdir,
                                         top_level = top_level,
                                         level1 = level1, level2 = level2, level3 = level3)
    fwrite(eGenes, file = paste0(outdir, "eGenes.", exp_suffix, ".txt"))
  }
  if(method == "mash"){
    print("This hasn't been implemented yet! Come back later")
  }
  
}