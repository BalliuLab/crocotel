#' Creates GReXs (Genetically regulated predictors of expression) for one gene across contexts
#'
#' @param exp_files - vector of file names for each context's expression for all genes and individuals 
#' 
#' @return writes out a file of predicted expression across individuals and contexts 
#' @export
format_data = function(exp_files, geneloc_file, snpsloc_file, genotypes_file, out_dir, cis_window = 1e6){
  out_dir = paste0(out_dir, "/crocotel_formatted_data/") 
  dir.create(out_dir, showWarnings = F)
  print("inferring context names from input expression files...")
  
  contexts = sub("\\..*$", "", basename(exp_files))
  print(paste0("inferred contexts are: ", paste(contexts, collapse = ", ")))
  
  gene_loc = fread(geneloc_file, sep = "\t", data.table = F)
  snps_loc = fread(snpsloc_file, sep = "\t", data.table = F)
  genotypes = fread(genotypes_file, sep = "\t", data.table = F)
  ## center and scale genotypes
  SNP_ids = genotypes[,1]
  genotypes = as.matrix(genotypes[, -1])
  row_means = rowMeans2(genotypes, na.rm = TRUE)
  row_sds   = rowSds(genotypes, na.rm = TRUE)
  genotypes = (genotypes - row_means) / row_sds
  genotypes = data.frame(SNP = SNP_ids, genotypes)
  
  get_gene_genotypes = function(chrom, upstream_pos, downstream_pos, snps_loc, genotypes){
    cur_snps = snps_loc %>% filter(chrom == chr) %>% filter(pos >= upstream_pos & pos <= downstream_pos) %>% dplyr::select(snp) %>% unlist() %>% unname()
    cur_genos = genotypes %>% filter(SNP %in% cur_snps)
    cur_genos = data.frame(t(cur_genos[,-1]))
    cur_genos = cbind(id = rownames(cur_genos), cur_genos)
    return(cur_genos)
  }
  
  ###### format each expression file such that each gene has one directory and one context file per directory
  for(i in 1:length(exp_files)){
    file = exp_files[i]
    df = fread(file, sep = "\t", data.table = F)
    genes = df$id
    context = sub("\\..*$", "", basename(file))
    for(gene in genes){
      #### get gene TSS position (start position of gene)
      names(gene_loc) = c("geneid", "chr", "s1", "s2")
      names(snps_loc) = c("snp", "chr", "pos")
      start_pos = gene_loc %>% filter(geneid == gene) %>% dplyr::select(s1) %>% unlist() %>% unname()
      chrom = gene_loc %>% filter(geneid == gene) %>% dplyr::select(chr) %>% unlist() %>% unname()
      upstream_pos = (start_pos - (cis_window))
      downstream_pos = (start_pos + (cis_window))
      
      # get cis-SNP genotypes for this gene
      if (!file.exists(paste0(out_dir, gene, "_genotypes.txt"))) {
        gene_genotypes = get_gene_genotypes(chrom, upstream_pos, downstream_pos, snps_loc, genotypes)
        fwrite(gene_genotypes, file = paste0(out_dir, gene, "_genotypes.txt"), sep = "\t", col.names = F, row.names = F)
      }
      
      dir.create(paste0(out_dir, gene), showWarnings = F)
      gene_df = df %>% filter(id == gene)
      gene_df = gene_df[,-1]
      gene_df = data.frame(t(gene_df))
      gene_df = cbind(id = rownames(gene_df), gene_df)
      ##### match expression ids with genotype ids
      gene_df = gene_df[gene_genotypes$id,]
      
      if(nrow(gene_genotypes) == 0){
        next
      }
      
      fwrite(gene_df, file = paste0(out_dir, gene, "/", context, ".txt"), sep = "\t", col.names = F, row.names = F)
    }
    print(paste0("finished formatting context ", i))
    
  }
}