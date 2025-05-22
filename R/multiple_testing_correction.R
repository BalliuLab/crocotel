 
#' @export
multiple_testing_correction = function(m_eqtl_outfiles, n_sSNPs_per_gene_files, contexts_vec, fdr_thresh, outdir, method = "treeQTL", top_level = "R", exp_suffix){
  dir.create(outdir, showWarnings = F)
  if(method == "treeQTL"){
    level1 = fdr_thresh
    level2 = fdr_thresh
    level3 = fdr_thresh
    
    eGenes = get_eGenes_multi_tissue_mod(m_eqtl_outfiles = m_eqtl_outfiles, 
                                         n_SNPs_per_gene_files = n_SNPs_per_gene_files, 
                                         contexts_vec = contexts_vec, 
                                         level1 = level1, level2 = level2, level3 = level3, 
                                         exp_suffix = exp_suffix)
    fwrite(eGenes, file = paste0(outdir, "eGenes.", exp_suffix, ".txt"))
  }
  if(method == "mash"){
    print("This hasn't been implemented yet! Come back later")
  }
  
}