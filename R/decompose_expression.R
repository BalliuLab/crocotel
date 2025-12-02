
## this implementation requires that the order of contexts is the same order as the list of files provided - will try to optimize this such that this does not have to be the case.
#' @export
decompose_expression = function(total_exp, gene_name, context_thresh, data_dir = NULL){
  #%%%%%%%%%%%%%%% Decompose expression into homogeneous and heterogeneous context expression
  print("Decomposing data")
  if(context_thresh < 2){
    print("Context threshold is too low. Filtering for individuals that are present in at least 2 contexts.")
    context_thresh = 2
  }
  
  ### keep only samples that have repeats across number of contexts specified by threshold
  ids_to_keep = names(which(table(total_exp$id) >= context_thresh))
  
  expression = total_exp %>% filter(id %in% ids_to_keep)
  design = factor(expression$id)
  contexts=as.character(unique(expression$context))
  X = scale(x = as.matrix(expression[,-c(1:2)]), center = T, scale = F)
  
  indiv.names = expression$id
  rownames(X) = expression$id
  colnames(X) = gene_name
  
  ## calculate shared
  X.mean.indiv = matrix(apply(X, 2, tapply, design, mean, na.rm = TRUE),
                        nrow = length(unique(design)),
                        ncol = dim(X)[2],
                        dimnames = list(levels(as.factor(design)), colnames(X)))
  Xb = X.mean.indiv[as.character(design), ]
  Xw = X - Xb
  dimnames(Xw) = list(indiv.names, colnames(X))
  
  Xw = data.frame(id=expression$id,context=expression$context, Xw)
  X.mean.indiv_df = cbind(data.frame(X.mean.indiv), id = rownames(X.mean.indiv)) %>% rename("shared" = gene_name)
  
  ### reformat shared and specific into a matrix
  decomp_exp_mat_sp = Xw %>% pivot_wider(names_from = "context", values_from = gene_name) %>% as.data.frame()
  decomp_exp_mat = merge(decomp_exp_mat_sp, X.mean.indiv_df, by = "id")
  
  ## if output directory is provided, write out decoposed files
  if(!is.null(data_dir)){
    decomp_exp_file_name= paste0(data_dir, gene_name, ".decomposed_expression.txt")
    
    fwrite(decomp_exp_mat,
           file = decomp_exp_file_name, quote = F, row.names = F,
           col.names = T, append = F, sep = '\t')
    print("Saved decomposition matrix")
  }
  
  return(decomp_exp_mat)
}



#### test functions:
#exp_files = list.files("/Users/lkrockenberger/C-STEM/example_data/expression/")
#contexts = exp_files
#exp_files = paste0("/Users/lkrockenberger/C-STEM/example_data/expression/", exp_files)
#data_dir = "/Users/lkrockenberger/C-STEM/example_data/decomposed_expression/"
#context_thresh = 3
#gene = "gene1"

#decompose_expression(exp_files, gene, contexts, context_thresh, data_dir)
