#' crocotel: Context-Resolved Trans-eQTL Mapping via Genetically Regulated
#' Expression
#'
#' @keywords internal
#' @importFrom stats ave coef cor lm lm.fit quantile rbinom rnorm runif sd var
#' @importFrom utils read.table write.table
"_PACKAGE"

# Quiet R CMD check's "no visible binding for global variable" NOTE on
# data.table non-standard-evaluation symbols (., .N, .SD, :=) and the column
# names referenced inside data.table [i, j, by] expressions. These are not real
# globals; this declaration only suppresses the false-positive NOTE.
utils::globalVariables(c(
  ".", ".N", ".SD", ":=", "..triplet_cols",
  "g1", "g2", "gs", "gene", "gene_id", "chr", "start", "end", "context",
  "n_pairs", "reg", "tgt", "rchr", "tchr", "strength", "snps", "snp", "pos",
  "snp_chr", "partner", "sub", "lo", "hi", "pvalue",
  "accept", "bh_threshold", "is_active", "is_eGene", "k", "m_tests",
  "max_accept_k", "n_sel_tissues", "n_tested", "n_tested_tissues",
  "q1", "q2", "simes_min_term", "simes_p", "simes_p_xC"
))
