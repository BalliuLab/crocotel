# The script contains the functions used to simulate the data for for C-STEM projects
# Generate by: Lena Krockenberger, Feiyang Huang, Vivian Yee
# Date [2024-03-14]

library(MASS)
set.seed(12)

# contexts_variance: the variance of each context, 1-h^2-h_c^2
# h^2: the context shared heritability 
# h_c^2: the context specific heritability
# context_herit: the heritability of each context
get_cis_covariance <- function(numContexts, contexts_herit, rho) {
    rho_mat <- matrix(1, nrow = numContexts, ncol = numContexts)
    diag(rho_mat) <- 1/rho
    rho_mat <- rho * rho_mat
    contexts_variance <- rep(0, numContexts)
    for (item in 1:nrow(contexts_herit)) {
        total_herit <- contexts_herit[item, 3]
        contexts_variance[item] <- 1 - total_herit
    }
    sigma <- sqrt(contexts_variance)
    covariance <- rho_mat * (sigma %*% t(sigma))
    return(covariance)

}


# Genotype matrix: m x n matrix of genotypes for n individuals and m SNPs
# Input params - maf of SNP, number of individuals to simulate genotype 
# Output - vector of 0, 1, or 2 corresponding to 1 SNP for n number of individuals
# (should all have same MAF)
simGeno <- function(maf, n, numAlleles = 2) {
    geno <- rbinom(n, size = numAlleles, prob = maf)
    return(list(geno, maf))
}


# Input params - number of individuals n, number of SNPs m, min and max MAF
# Outputs SNP array of n x m matrix of genotypes
makeGenoMatrix <- function(n, m, mafs) {
    genotypes <- matrix(0, nrow = n, ncol = m)
    
    for (i in 1:m) {
        genotypes[, i] <- simGeno(mafs[i, 1], n)[[1]]
    }
    
    # center and scale genotypes
    genotypes1 <- scale(genotypes, center = TRUE, scale = TRUE)
    
    return(genotypes1)
}



# Input params - m is the number of SNPs, p is the probability a SNP is causal
# Outputs a m list of 0s and 1s indicating 
# whether the SNP at that index has a causal cis genetic effect
get_causal_snps <- function(m, p) {
    snp_list <- as.numeric(rbinom(m, size = 1, prob = p))
    return(snp_list)
}


# Input params - 
# heritability: the cis_shared heritability, 
# snp_list: the list of causal snps from get_causal_snps
# causal: the number of causal snps causal=M*p,p is the prob that a SNP is causal (0.01)
# Outputs: a list of shared effect sizes, 1 by m list of shared effect sizes
get_shared_effects <- function(heritability, snp_list) {
    causal <- sum(snp_list == 1)
    for (item in 1:length(snp_list)) {
        if (snp_list[item] == 1) {
            effect_size <- rnorm(1, mean = 0, sd = (sqrt(heritability/causal))[1,1])
            snp_list[item] <- effect_size
        }
    }
    return(snp_list)
}

# Input params - 
# lam: probability a SNP is context specific
# heritability: context specific heritability 
# snp_list: a list of causal snps from get_causal_snps
# num_causal: the number of causal snps, num_causal=M*p*lam
# Outputs a list of specific effect sizes, 1 by m list of specific effect sizes
get_specific_effects <- function(lam, heritability, snp_list) {
    for (item in 1:length(snp_list)) {
        if (snp_list[item] == 1) {
            specific <- rbinom(1, size = 1, prob = lam)
            snp_list[item] <- specific
        }
    }
    num_causal <- sum(snp_list == 1)
    for (item in 1:length(snp_list)) {
        if (snp_list[item] == 1) {
            effect_size <- rnorm(1, mean = 0, sd = sqrt(heritability/num_causal))
            snp_list[item] <- effect_size
        } else {
            snp_list[item] <- 0
        }
    }
    return(snp_list)
}

# Input params - 
# genotypes: makeGenoMatrix
# shared effects for the gene: get_shared_effects
# specific effects for the gene: get_specific_effects
# covariance matrix: get_cis_covariance  
# context number: fixed
# context: the context we are going to simulate the expression for
# Outputs: 
# expression vector for the gene and context (simulated cis expression for one gene and one context)
# Can be deleted
get_cis_expression <- function(genotypes, shared_effects, specific_effects, covariance, context) {
    residuals <- mvrnorm(n = nrow(covariance), mu = rep(0, ncol(covariance)), Sigma = covariance)[,context]
    expression <- genotypes %*% shared_effects + genotypes %*% specific_effects + residuals
    return(expression)
}

# Simulates one cis gene in all contexts
simulate_one_gene <- function(numContexts, n, m, p, lam, mafs, cis_heritabilities, cis_covariance) {
    genotype_matrix <- makeGenoMatrix(n, m, mafs)
    causal_snps <- get_causal_snps(m, p)
    shared_effects <- get_shared_effects(cis_heritabilities[1], causal_snps)
    
    context_cis_expression <- vector("list", numContexts)
    context_specific_effects <- vector("list", numContexts)
    context_specific_residuals <- vector("list", numContexts)
    
    residuals <- mvrnorm(n = nrow(cis_covariance), mu = rep(0, ncol(cis_covariance)), Sigma = cis_covariance)
    specific_effects <- matrix(0, nrow = numContexts, ncol = m)
    cis_expression <- matrix(0, nrow = numContexts, ncol = n)
    for (context in 1:numContexts) {
        specific_effects[context,] <- get_specific_effects(lam, cis_heritabilities[context, 2], causal_snps)
        #cis_expression <- get_cis_expression(genotype_matrix, shared_effects, specific_effects, cis_covariance, context)
        cis_expression[context,] <- genotype_matrix %*% shared_effects + genotype_matrix %*% specific_effects[context,] + residuals[context,]
        
    }
    return(list(shared_effects, specific_effects, cis_expression, genotype_matrix, residuals))
}













