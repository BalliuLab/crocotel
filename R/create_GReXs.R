## THIS IS NOT A FUNCTION YET - it is the code directly used from C-STEM simulation pipeline

### Test parameters
X_file = "/Users/lkrockenberger/C-STEM/example_data/genos"
exp_files = list.files("/Users/lkrockenberger/C-STEM/example_data/expression/")
contexts = exp_files
exp_files = paste0("/Users/lkrockenberger/C-STEM/example_data/expression/", exp_files)
out_dir = "/Users/lkrockenberger/C-STEM/example_data/GReXs/"
snps = "/Users/lkrockenberger/C-STEM/example_data/snps.txt"
decomposition_dir = "/Users/lkrockenberger/C-STEM/example_data/decomposed_expression/"
gene_name = "gene1"
context_thresh = 3
alpha = 0.5
num_folds = 10
run_CxC = TRUE

library(data.table)
library(bigstatsr)
library(dplyr)
library(caret)

source("R/decompose_expression.R")
source("R/functions.R")

## parameters:
# X_file - a genotype file with the set of individuals in all Y_files and their corresponding cis-SNPs. Individuals are row names, no column names
# exp_files - vector with list of all expression files per context. The individuals (row names) must be a subset of genotype (X_file). Contains an unnamed
# contexts - vector of strings which have context names in them.
# column with the expression (or if cov_file_dir is NULL, the residual expression for a given context.")
# out_dir - genotype (X_file). Contains an unnamed column with the expression (or if cov_file_dir is NULL, the residual expression for a given context.")
# snps - A tab-delimited txt file containing information from your .bed or other genotype file. Contains 6 columns and the number of rows corresponds to the number of snps. 
#No column names or row names. Col2 must be rsIDs for TWAS to work. 
#Col1 chromosome
#Col2 rsID
#Col3 location CM (this doesn't really matter for TWAS)
#               Col4 location/locus on chromosome
#               Col5 allele1
#               Col6 allele2
# gene_name - identifier for current gene being run. Add this prefix to all the saved results files... necessary to distinguish results where more than one gene-analysis is run.
# decomposition_dir - directory to store decomposed expression files
# context_thresh - minimum number of contexts to run C-STEM on. (must be >=2 )
# alpha - The regularization constant. Default is .5 (eNet). Minimum value of 1e-4.
# num_folds - Number of folds for cross-validation
# run_CxC - Takes values of TRUE or FALSE. Default is FALSE, if TRUE then will run the CxC method too.
create_GReXs = function(X_file, exp_files, contexts, out_dir, snps, gene_name, decomposition_dir = NULL, context_thresh = 3, alpha = 0.5, num_folds = 10, run_CxC = FALSE){
  seed = 9000
  set.seed(seed)
  message("Saving cross-validated predictors and performance metrics in ", out_dir)

  message("Reading in files...") 
  suppressWarnings(expr = {X<-fread(file = X_file, sep='\t', data.table=F)})
  X<-as.matrix(data.frame(X, row.names=1, check.names = F))
  
  ###### read in expression and decompose - files are written out to decomposed exp directory
  decompose_expression(exp_files, gene_name, contexts, context_thresh, decomposition_dir)
  ## this assumes that the file name of decomposition dir is saved as "gene.context.etc" and shared is called "Average Context" (output of decompose function)
  contexts_vec = sapply(strsplit(list.files(decomposition_dir), "\\."), "[[", 2)
  
  ### read in expression
  Ys<-vector("list", length = length(list.files(decomposition_dir)))
  names(Ys)<-contexts_vec
  lengths_y = c()
  rownames_y = list()
  Yhats_full<-vector("list", length(Ys))
  for(i in 1:length(Ys)){
    suppressWarnings(expr={ Ys[[i]]<-fread(file = 
                                             paste0(decomposition_dir,list.files(decomposition_dir)[i]), sep='\t', data.table=F)})
    Ys[[i]]<-as.matrix(data.frame(Ys[[i]], row.names=1, check.names = F))
    lengths_y = c(lengths_y, nrow(Ys[[i]]))
    rownames_y[[i]] = rownames(Ys[[i]])
    if(any(is.na(Ys[[i]])) | any(is.nan(Ys[[i]]))){
      remove=unique(c( which(is.na(Ys[[i]])), which(is.nan(Ys[[i]])) ))
      Ys[[i]]=Ys[[i]][-remove,,drop=F]
    }
    
    Yhats_full[[i]]<-matrix(NA, ncol=1, nrow=lengths_y[i], 
                            dimnames = list(rownames_y[[i]], "pred"))
  }
  
  
  Yhats_tiss = crossval_helper(Ys, X, lengths_y, rownames_y, contexts_vec, out_dir, gene_name, num_folds, alpha)
  fwrite(Yhats_tiss, file = paste0(out_dir,gene_name,"_cstem_predictors.txt", sep = "\t"))
  evaluation_helper(Yhats_tiss, contexts_vec, run_CxC, Yhats_full, out_dir, gene_name)
  
  ### read in expression for CxC
  if(run_CxC){
    Ys_cxc<-vector("list", length = length(exp_files))
    names(Ys_cxc)<-sapply(strsplit(exp_files, "/"), tail, 1)
    lengths_y_cxc = c()
    rownames_y_cxc = list()
    for(i in 1:length(Ys_cxc)){
      suppressWarnings(expr={ Ys_cxc[[i]]<-fread(file = 
                                                   exp_files[i], sep='\t', data.table=F)})
      Ys_cxc[[i]]<-as.matrix(data.frame(Ys_cxc[[i]], row.names=1, check.names = F))
      lengths_y_cxc = c(lengths_y_cxc, nrow(Ys_cxc[[i]]))
      rownames_y_cxc[[i]] = rownames(Ys_cxc[[i]])
      if(any(is.na(Ys_cxc[[i]])) | any(is.nan(Ys_cxc[[i]]))){
        remove=unique(c( which(is.na(Ys_cxc[[i]])), which(is.nan(Ys_cxc[[i]])) ))
        Ys_cxc[[i]]=Ys_cxc[[i]][-remove,,drop=F]
      }
    }
    Yhats_cxc = crossval_helper(Ys_cxc, X, lengths_y_cxc, rownames_y_cxc, out_dir, gene_name, num_folds, alpha)
    fwrite(Yhats_cxc, file = paste0(out_dir,gene_name,"_CxC_predictors.txt", sep = "\t"))
    evaluation_helper(Yhats_cxc, contexts_vec, run_CxC, NULL, out_dir, gene_name)
  }
  
}

