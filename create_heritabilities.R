args <- commandArgs(trailingOnly = TRUE)

numContexts <- as.integer(args[1])
filepath_out <- args[2]
cis_shared <- as.numeric(args[3])
total_cis <- as.numeric(args[4])
trans_herit <- as.numeric(args[5])

get_cis_context_heritabilities <- function(numContexts, shared, total_cis_herit) {
  contexts_herit <- rep(0, 3)
  for(item in 1:length(contexts_herit)) {
    specific <- runif(n = 1, min = 0.05, max = total_cis_herit - shared)
    contexts_herit[item, 1] <- shared 
    contexts_herit[item, 2] <- specific
    contexts_herit[item, 3] <- (shared + specific) 
  }
  contexts_herit
}

get_trans_contexts_heritabilities <- function(numContexts, max_trans_herit) {
  contexts_herit <- rep(0, 1)
  for(item in 1:length(contexts_herit)) {
    herit <- runif(n = 1, min = max_trans_herit, max_trans_herit) 
    contexts_herit[item, 0] <- herit 
  }
  contexts_herit
}

get_null_cis_context_heritabilities <- function(numContexts) {
  contexts_herit <- matrix(0, nrow = numContexts, ncol = 3)
  contexts_herit
}

get_null_trans_context_heritabilities(numContexts) {
  contexts_herit <- matrix(0, nrow = numContexts, ncol = 1)
  contexts_herit
}

cis_contexts_herit <- get_cis_context_heritabilities(numContexts, cis_shared, total_cis)
trans_contexts_herit <- get_trans_context_heritabilities(numContexts, trans_herit)
null_cis_herit <- get_null_cis_context_heritabilities(numContexts)
null_trans_herit <- get_null_trans_contexts_heritabilities(numContexts)

cis_herit_file <- file(paste0(filepath_out, numContexts, "-", total_cis + "_cis_heritability.txt"), "w")
trans_herit_file <- file(paste0(filepath_out, numContexts, "_", trans_herit, "_trans_heritability.txt"), "w")
null_cis_herit_file <- file(paste0(filepath_out, numContexts, "_", "null_cis_heritability.txt"), "w")
null_trans_herit_file <- file(paste0(filepath_out, numContexts, "_", "null_trans_heritability.txt"), "w")


write.table(cis_contexts_herit, file = cis_herit_file, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(trans_contexts_herit, file = trans_herit_file, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(null_cis_herit, file = null_cis_herit_file, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(null_trans_herit, file = null_trans_herit_file, sep = "\t", row.names = FALSE, col.names = FALSE)

close(cis_herit_file)
close(trans_herit_file)
close(null_cis_herit_file)
close(null_trans_herit_file)
