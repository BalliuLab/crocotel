# The script is used for simulating data for C-STEM projects. 
# It reads in parameters from the command line and uses the functions from `Functions.R`.
# Generate by: Lena Krockenberger, Feiyang Huang, Vivian Yee
# Date [2024-03-14]
# sample command: Rscript simulate.R 1000 1000 10 0.01 0.3 0.2 snp_mafs.txt
# (make sure the snp_mafs.txt file is in the same directory as the simulate.R file)

library(MASS)
library(data.table)
library(fs)
source("functions.R")


#setwd("~/Documents/BalliuLab/C-STEM/simulation")
# Read parameters from command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) < 8) {
    stop("Insufficient number of arguments! Usage: Rscript simulate.R n m numContexts p lam rho num_sim")
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
maf_file <- args[7]
num_iterations <- as.integer(args[8])
seed <- ifelse(length(args) > 9, as.integer(args[9]), NULL)

# t: number of target gene

#t < as.numeric(args[8])
#context_cis_herit_file = args[9]
#context_herit_file = args[10]
#effects_output = args[11]
#cis_expression_dir = args[12]
#trans_exp_dir = args[13]
#cov_output_dir = args[14]
#ntransT = as.numeric(args[15])

# Set working directory and input given parameters
#setwd("~/BalliuLab/CSTEM")
mafs <- read.table(maf_file, header = FALSE)
mafs <- as.matrix(mafs)
cis_hearitabilities <- read.table("cis_heritability.txt", header = FALSE)


# Run simulations and save data
time_start <- Sys.time()
for (i in 1:num_iterations) {
    cis_covariance <- get_cis_covariance(numContexts, cis_hearitabilities, rho)
    data <- simulate_one_gene(numContexts, n, m, p, lam, mafs, cis_hearitabilities, cis_covariance, seed)
    # Share effect and specific effect
    data_share_eff <- as.matrix(data[[1]])
    data_specific_eff <- as.matrix(data[[2]])
    
    # CIS expression
    data_exp <- as.matrix(data[[3]])
    
    # Genotype matrix
    Geno <- as.matrix(data[[4]])
    
    
    # Format data and save
    #format_exp(iteration, n, numContexts, data_exp, "~/expression")
    format_geno(n, m, "~/geno", Geno, iteration)
    if (!is.null(seed)) {
        seed <- seed + 1
    }
}
time_end <- Sys.time()
print(paste("Time taken:", time_end - time_start))



