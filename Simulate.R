# The script is used for simulating data for C-STEM projects. 
# It reads in parameters from the command line and uses the functions from `Functions.R`.
# Generate by: Lena Krockenberger, Feiyang Huang, Vivian Yee
# Date [2024-03-14]
# sample command: Rscript simulate.R 1000 1000 10 0.01 0.3 0.2 10 snp_mafs.txt
# (make sure the snp_mafs.txt file is in the same directory as the simulate.R file)

library(MASS)
source("functions.R")
set.seed(12)

# Read parameters from command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) < 16) {
    stop("Insufficient number of arguments! Usage: Rscript simulate.R n m numContexts p lam rho num_sim")
}

# Assign parameters
# n: number of individuals, m: length of SNP array, numContexts: number of contexts
# p: probability a SNP is causal, lam: the heritability of the gene
# rho: the correlation between contexts, num_sim: number of simulations
n <- as.integer(args[1])
m <- as.integer(args[2])
numContexts <- as.integer(args[3])
p <- as.numeric(args[4])
lam <- as.numeric(args[5])
rho <- as.numeric(args[6])
num_sim <- as.integer(args[7])
maf_file <- args[8]

# t: number of target gene
t < as.numeric(args[9])
context_cis_herit_file = args[10]
context_herit_file = args[11]
effects_output = args[12]
cis_expression_dir = args[13]
trans_exp_dir = args[14]
cov_output_dir = args[15]
ntransT = as.numeric(args[16])

# Set working directory and input given parameters
#setwd("~/BalliuLab/CSTEM")
mafs <- read.table(maf_file, header = FALSE)
mafs <- as.matrix(mafs)
cis_hearitabilities <- read.table("cis_heritability.txt", header = FALSE)

simulate_and_save <- function(iteration) {
    cis_covariance <- get_cis_covariance(numContexts, cis_hearitabilities, rho)
    data <- simulate_one_gene(numContexts, n, m, p, lam, mafs, cis_hearitabilities, cis_covariance)
    
    # Share effect and specific effect
    data_share_eff <- as.matrix(data[[1]])
    data_specific_eff <- as.matrix(data[[2]])
    
    # CIS expression
    data_exp <- as.matrix(data[[3]])
    
    # Genotype matrix
    Geno <- as.matrix(data[[4]])
    
    # Format data
    
    
    # Write data to files
  
}

# Run simulations and save data
time_start <- Sys.time()
for (i in 1:num_sim) {
    simulate_and_save(i)
    
}
time_end <- Sys.time()
print(paste("Time taken:", time_end - time_start))



