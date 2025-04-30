

decompose=function(expression, shared_exp_file_name, spec_exp_file_name, genos){

  design = factor(expression$id)
  contexts=as.character(unique(expression$context))

  if (any(summary(as.factor(design)) == 1))
    stop("A multilevel analysis can not be performed when at least one some sample is not repeated.")

  X = scale(x = as.matrix(expression[,-c(1:2)]), center = T, scale = F)

  indiv.names = rownames(X)
  rownames(X) = as.character(design)

  X.mean.indiv = matrix(apply(X, 2, tapply, design, mean, na.rm = TRUE),
                        nrow = length(unique(design)),
                        ncol = dim(X)[2],
                        dimnames = list(levels(as.factor(design)), colnames(X)))
  Xb = X.mean.indiv[as.character(design), ]
  Xw = X - Xb
  dimnames(Xw) = list(indiv.names, colnames(X))


  fwrite(x = data.table(t(X.mean.indiv[colnames(genos),]),keep.rownames = T) %>% {setnames(., old = "rn", new = "geneID")[]},
         file = shared_exp_file_name, quote = F, row.names = F,
         col.names = T, append = F, sep = '\t')
  print("Saved shared expression matrix")


  Xw = data.frame(id=expression$id,context=expression$context, Xw)

  for(j in 1:length(contexts)){

    wexp_t = data.frame(Xw[Xw$context == contexts[j],-2],row.names = 1)

    fwrite(x = data.table(t(wexp_t[colnames(genos),]),keep.rownames = T) %>% {setnames(., old = "rn", new = "geneID")[]},
           file = spec_exp_file_name[j],quote = F, row.names = F,
           col.names = T, append = F, sep = '\t')

    print(paste0("Saving (specific) expression matrix for context: ",contexts[j]))

  }

}

decompose_expression = function(SNP_file_name, exp_files, contexts, outfile_suffix){
    #%%%%%%%%%%%%%%% Read genotype matrix, SNPs in rows, samples in columns
    genos=data.frame(fread(file = SNP_file_name, nrows = 5),row.names = 1, check.names = F)

    # Read expression matrix for Context t and merge with other Contexts 
    exp_all=data.frame(fread(input = exp_files[1], header = T), check.names = F,stringsAsFactors = F)
    colnames(exp_all)[-1] = paste(colnames(exp_all)[-1],contexts[1], sep = " - ")
    print(paste("Finished merging context",1))

    for(i in 2:length(exp_files)){

    # Read expression matrix for Context t
    exp_t=data.frame(fread(input = exp_files[i], header = T), check.names = F,stringsAsFactors = F)
    colnames(exp_t)[-1] = paste(colnames(exp_t)[-1],contexts[i], sep = " - ")

    # Merge with other Contexts
    exp_all = merge(x = exp_all, y = exp_t, by="gene_id", all = TRUE)

    print(paste("Finished merging context",i))
    }

    # Transpose merged expression matrix to have genes in the columns 
    exp_all=t(data.frame(exp_all,row.names = 1,check.names = F))
    print("Finished transposing merged file")

    #%%%%%%%%%%%%%%% Sample and context names
    indv_contexts=matrix(unlist(strsplit(rownames(exp_all), split = " - ")), ncol = 2,byrow = T)
    exp_all = data.frame(id=indv_contexts[,1],context=indv_contexts[,2], exp_all)


    #%%%%%%%%%%%%%%% Decompose expression into homogeneous and heterogeneous context expression
    print("Decomposing data")
    shared_exp_file_name= paste0(data_dir, "AverageContext.", outfile_suffix)
    spec_exp_file_name= paste0(data_dir, contexts,".", outfile_suffix)
    decompose(expression=exp_all, shared_exp_file_name, spec_exp_file_name, genos)

}



