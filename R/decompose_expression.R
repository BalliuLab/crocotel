
## this implementation requires that the order of contexts is the same order as the list of files provided - will try to optimize this such that this does not have to be the case.
#' @export
decompose_expression = function(exp_files, gene, contexts, context_thresh, data_dir){
  # Read expression matrix for Context t and merge with other Contexts 
  exp_all=data.frame(fread(input = exp_files[1], header = F), check.names = F,stringsAsFactors = F)
  names(exp_all) = c("id", gene)
  exp_all = exp_all %>% mutate(context = contexts[1]) %>% select("id", "context", all_of(gene))
  print(paste("Finished merging context",1))
  
  for(i in 2:length(exp_files)){
    
    # Read expression matrix for Context t
    exp_t=data.frame(fread(input = exp_files[i], header = F), check.names = F,stringsAsFactors = F)
    #### remove NA rows
    exp_t = na.omit(exp_t)
    names(exp_t) = c("id", gene)
    exp_t = exp_t %>% mutate(context = contexts[i]) %>% select("id", "context", all_of(gene))
    
    # Merge with other Contexts
    exp_all = rbind(exp_all, exp_t)
    
    print(paste("Finished merging context",i))
  }
  
  
  #%%%%%%%%%%%%%%% Decompose expression into homogeneous and heterogeneous context expression
  print("Decomposing data")
  #if(context_thresh < 2){
  #  print("Context threshold is too low. Filtering for individuals that are present in at least 2 contexts.")
  #  context_thresh = 2
  #}
  shared_exp_file_name= paste0(data_dir, gene, ".AverageContext.decomposed_expression.txt")
  spec_exp_file_name= paste0(data_dir, gene, ".", contexts,".decomposed_expression.txt")
  
  ### keep only samples that have repeats across number of contexts specified by threshold
  ids_to_keep = names(which(table(exp_all$id) >= context_thresh))
  
  expression = exp_all %>% filter(id %in% ids_to_keep)
  design = factor(expression$id)
  contexts=as.character(unique(expression$context))
  #X = scale(x = as.matrix(expression[,-c(1:2)]), center = T, scale = F)
  X = as.matrix(expression %>%
                  group_by(context) %>%
                  mutate(
                    scaled = scale(!!sym(gene))
                  ) %>% ungroup() %>% select(scaled))
  
  indiv.names = expression$id
  rownames(X) = expression$id
  colnames(X) = gene
  
  ## calculate shared
  X.mean.indiv = matrix(apply(X, 2, tapply, design, mean, na.rm = TRUE),
                        nrow = length(unique(design)),
                        ncol = dim(X)[2],
                        dimnames = list(levels(as.factor(design)), colnames(X)))
  Xb = X.mean.indiv[as.character(design), ]
  Xw = X - Xb
  dimnames(Xw) = list(indiv.names, colnames(X))
  
  
  fwrite(x = data.table(X.mean.indiv,keep.rownames = T) %>% rename("id" = rn),
         file = shared_exp_file_name, quote = F, row.names = F,
         col.names = T, append = F, sep = '\t')
  print("Saved shared expression matrix")
  
  
  Xw = data.frame(id=expression$id,context=expression$context, Xw)
  
  for(j in 1:length(contexts)){
    
    wexp_t = data.frame(Xw[Xw$context == contexts[j],-2],row.names = 1)
    
    fwrite(x = data.table(wexp_t,keep.rownames = T) %>% rename("id" = rn),
           file = spec_exp_file_name[j],quote = F, row.names = F,
           col.names = T, append = F, sep = '\t')
    
    print(paste0("Saving (specific) expression matrix for context: ",contexts[j]))
    
  }
  
}



#### test functions:
#exp_files = list.files("/Users/lkrockenberger/C-STEM/example_data/expression/")
#contexts = exp_files
#exp_files = paste0("/Users/lkrockenberger/C-STEM/example_data/expression/", exp_files)
#data_dir = "/Users/lkrockenberger/C-STEM/example_data/decomposed_expression/"
#context_thresh = 3
#gene = "gene1"

#decompose_expression(exp_files, gene, contexts, context_thresh, data_dir)
