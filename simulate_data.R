# The script is used for simulating data for C-STEM projects. 
# It reads in parameters from the command line and uses the functions from `Functions.R`.
# Generate by: Lena Krockenberger, Feiyang Huang, Vivian Yee
# Date [2024-03-14]
# sample command: Rscript simulate_data.R 1000 1000 10 0.01 0.3 0.2 10 10 3 "effect_sizes/" "expression/" "trans_expression/" "cov_output/" 3 0.05 0.15 0.2
# (make sure the snp_mafs.txt file is in the same directory as the simulate.R file)

library(MASS)
library(data.table)
library(fs)
source("simulation_functions.R")


#setwd("~/Documents/BalliuLab/C-STEM/simulation")
# Read parameters from command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) < 17) {
    stop("Insufficient number of arguments!")
}

# Assign parameters
# n: number of individuals, m: length of SNP array, numContexts: number of contexts
# p: probability a SNP is causal, lam: probability a SNP is context specific
# rho: the correlation between contexts, num_iterations: number of iterations
n <- as.integer(args[1])
m <- as.integer(args[2])
numContexts <- as.integer(args[3])
p <- as.numeric(args[4])
lam <- as.numeric(args[5])
rho <- as.numeric(args[6])
num_iterations <- as.integer(args[7])
seed <- as.numeric(args[8])

# t: number of target gene
# ntransT: number of Tissues that have trans effect
t <- as.numeric(args[9])

effects_output = args[10]
cis_expression_dir = args[11]
trans_exp_dir = args[12]
cov_output_dir = args[13]
ntransT = as.numeric(args[14])


# shared cis heritability
cis_shared = as.numeric(args[15])
# total cis heritability
total_cis = as.numeric(args[16])
# total trans heritability
trans_herit = as.numeric(args[17])

# Set working directory and input given parameters
#setwd("~/BalliuLab/CSTEM")

set.seed(seed)

mafs <- generate_mafs(m)
generate_mafs_file <- fwrite(as.data.frame(mafs), "snp_mafs.txt", sep = "\t", col.names = FALSE)


contexts_cis_herit <- get_cis_context_heritabilities(numContexts, cis_shared, total_cis)
trans_herit <- get_trans_context_heritabilities(numContexts, trans_herit)

cis_herit_file <- fwrite(as.data.frame(contexts_cis_herit), "cis_heritability.txt", sep="\t", col.names=FALSE)
trans_herit_file <- fwrite(as.data.frame(trans_herit), "trans_heritability.txt", sep="\t", col.names=FALSE)


covariance = get_trans_covariance(numContexts, trans_herit, rho)
is_trans_effect_vec <- is_trans_effect(ntransT, numContexts, seed)
genotype_matrix <- makeGenoMatrix(n, m, mafs)
genotype_matrix_file <- fwrite(as.data.frame(genotype_matrix), "genotype_matrix.txt", sep = "\t", col.names = FALSE)


# Run simulations and save data
time_start <- Sys.time()
for (i in 1:num_iterations) {
    set.seed(seed)
    cis_covariance <- get_cis_covariance(numContexts, contexts_cis_herit, rho)
    data <- simulate_one_gene(numContexts, n, m, p, lam, mafs, contexts_cis_herit, cis_covariance, genotype_matrix, seed)
    # Share effect and specific effect
    data_share_eff <- as.matrix(data[[1]])
    data_specific_eff <- as.matrix(data[[2]])
    
    # CIS expression
    data_exp <- as.matrix(data[[3]])
    
    # Genotype matrix
    Geno <- as.matrix(data[[4]])
    
    
    # Format data and save
    format_exp(i, n, numContexts, data_exp, "~/expression")
    format_geno(n, m, "~/geno", Geno, i)
    if (!is.null(seed)) {
        seed <- seed + 1
    }
}

simulateTransExpression(cis_expression_dir, n, covariance,contexts_cis_herit, trans_herit, is_trans_effect_vec, effects_output, trans_exp_dir)

time_end <- Sys.time()
print(paste("Time taken:", time_end - time_start))



