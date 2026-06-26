#' Format MatrixEQTL-style inputs into per-gene crocotel directories
#'
#' Splits genome-wide expression / genotype inputs into one directory per gene
#' (one expression file per context) plus a per-gene cis-SNP genotype file, the
#' layout consumed by \code{create_GReXs}.
#'
#' @param exp_files - vector of expression files, one per context (rows = genes, columns = individuals)
#' @param geneloc_file - gene location file: geneid, chr, start, end
#' @param snpsloc_file - SNP location file: snp, chr, pos
#' @param genotypes_file - genotype matrix: first column "SNP", remaining columns individuals
#' @param out_dir - output directory; files are written under out_dir/crocotel_formatted_data/
#' @param cis_window - cis distance (bp) around the gene TSS used to select SNPs. Default 1e6.
#' @return writes a per-gene genotype file and per-(gene,context) expression file
#' @export
format_data = function(exp_files, geneloc_file, snpsloc_file, genotypes_file, out_dir, cis_window = 1e6){
  out_dir = paste0(out_dir, "/crocotel_formatted_data/")
  dir.create(out_dir, showWarnings = FALSE)
  message("inferring context names from input expression files...")

  contexts = sub("\\..*$", "", basename(exp_files))
  message("inferred contexts are: ", paste(contexts, collapse = ", "))

  gene_loc  = data.table::fread(geneloc_file,   sep = "\t", data.table = FALSE)
  snps_loc  = data.table::fread(snpsloc_file,   sep = "\t", data.table = FALSE)
  genotypes = data.table::fread(genotypes_file, sep = "\t", data.table = FALSE)
  names(gene_loc) = c("geneid", "chr", "s1", "s2")
  names(snps_loc) = c("snp", "chr", "pos")

  ## Individual order and SNP order are taken from the genotype file and are the
  ## same for every gene; all output files follow these orders (as before).
  geno_ind_order = names(genotypes)[-1]
  geno_snp_order = genotypes$SNP
  rank_in_file = seq_along(geno_snp_order)
  names(rank_in_file) = geno_snp_order

  ## Index gene positions once (was an O(G) dplyr filter per gene per context).
  gene_pos = gene_loc[!duplicated(gene_loc$geneid), c("geneid", "chr", "s1")]
  rownames(gene_pos) = gene_pos$geneid

  ## Keyed SNP locations for binary-search cis range queries, and a SNP-keyed
  ## genotype matrix for fast subsetting (was a full %in% scan per gene).
  snps_dt = data.table::as.data.table(snps_loc)
  data.table::setkey(snps_dt, chr, pos)
  geno_dt = data.table::as.data.table(genotypes)
  data.table::setkey(geno_dt, SNP)
  pos = snp = SNP = NULL  # silence R CMD check NSE notes

  ## cis-genotype block for one gene, in the genotype file's individual/SNP
  ## order. Returns NULL when the gene has no cis SNPs (caller skips it).
  extract_cis_genos = function(gene){
    if (!gene %in% rownames(gene_pos)) return(NULL)
    chrom     = gene_pos[gene, "chr"]
    start_pos = gene_pos[gene, "s1"]
    up   = start_pos - cis_window
    down = start_pos + cis_window
    cur_snps = snps_dt[.(chrom), nomatch = 0L][pos >= up & pos <= down, snp]
    if (length(cur_snps) == 0) return(NULL)
    sub = geno_dt[.(cur_snps), nomatch = 0L]
    sub = sub[order(rank_in_file[sub$SNP])]              # restore genotype-file SNP order
    m = t(as.matrix(sub[, -1, with = FALSE]))            # individuals x SNPs
    data.frame(id = geno_ind_order, m, check.names = FALSE, stringsAsFactors = FALSE)
  }

  skip_genes = character(0)   # genes with no cis SNPs (warned once, skipped everywhere)

  for (i in seq_along(exp_files)) {
    df = data.table::fread(exp_files[i], sep = "\t", data.table = FALSE)
    context = contexts[i]
    genes = df$id
    expr_mat = as.matrix(df[, -1, drop = FALSE])         # genes x individuals, row j == genes[j]
    expr_inds = colnames(expr_mat)
    ind_idx = match(geno_ind_order, expr_inds)           # NA where a genotype individual is absent
    present     = !is.na(ind_idx)                         # genotype individuals present in this context
    present_ids = geno_ind_order[present]                 # only these are written (no NA-id padding)
    present_col = ind_idx[present]

    for (g in seq_along(genes)) {
      gene = genes[g]
      if (gene %in% skip_genes) next
      geno_file = paste0(out_dir, gene, "_genotypes.txt")

      # Extract + write genotypes exactly once per gene; skip genes with no cis SNPs.
      if (!file.exists(geno_file)) {
        cur_genos = extract_cis_genos(gene)
        if (is.null(cur_genos)) {
          warning("Skipping gene '", gene, "': no cis SNPs within ", format(cis_window, scientific = FALSE),
                  " bp of its TSS.", call. = FALSE, immediate. = TRUE)
          skip_genes = c(skip_genes, gene)
          next
        }
        data.table::fwrite(cur_genos, file = geno_file,
                           sep = "\t", col.names = FALSE, row.names = FALSE)
      }

      dir.create(paste0(out_dir, gene), showWarnings = FALSE)

      # Expression row -> one value per individual PRESENT in this context's expression.
      # Do NOT pad absent genotype individuals: padding wrote rows with an empty id,
      # which collide as duplicate (id, context) keys in create_GReXs' pivot_wider
      # whenever the genotype individuals are a superset of a context's expression
      # individuals -- the normal ragged multi-tissue case (e.g. GTEx full panel vs
      # per-tissue expression).
      gene_df = data.frame(
        id    = present_ids,
        value = expr_mat[g, present_col],
        check.names = FALSE, stringsAsFactors = FALSE
      )
      data.table::fwrite(gene_df, file = paste0(out_dir, gene, "/", context, ".txt"),
                         sep = "\t", col.names = FALSE, row.names = FALSE)
    }
    message("finished formatting context ", i)
  }
}
