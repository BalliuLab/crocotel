 

multiple_testing_correction = function(input_dir, contexts_vec, fdr_thresh, outdir, method = "treeQTL", top_level = "T"){
  if(method == "treeQTL"){
    level1 = fdr_thresh
    level2 = fdr_thresh
    level3 = fdr_thresh
    exp_suffix = ""
    
    eGenes = get_eGenes_multi_tissue_mod(m_eqtl_out_dir = input_dir,
                                         treeQTL_dir=outdir,
                                         tissue_names=contexts_vec,
                                         level1 = level1, level2 = level1, level3 = level1,
                                         exp_suffix=exp_suffix)
    fwrite(eGenes, file = paste0(outdir, "test.txt"))
  }
  if(method == "mash"){
    
  }
  
}