# The script contains the functions used to simulate the data for for C-STEM projects
# Generate by: Lena Krockenberger, Feiyang Huang, Vivian Yee
# Date [2024-03-14]

library(data.table)
library(MASS)


generate_mafs <- function(x) {
    mafs <- runif(x, 0.05, 0.5)
    mafs <-as.matrix(mafs)
    return(mafs)
}



get_cis_context_heritabilities <- function(numContexts, shared, total_cis_herit) {
    contexts_herit <- matrix(0, nrow=numContexts, ncol=3)
    for (i in 1:nrow(contexts_herit)) {
        specific <- runif(1, 0.05, (total_cis_herit - shared))
        contexts_herit[i,] <- c(shared, specific, shared + specific)
    }
    return(contexts_herit)
}


get_trans_context_heritabilities <- function(numContexts, trans_herit) {
    contexts_herit <- matrix(trans_herit, nrow = numContexts, ncol = 1)
    return(contexts_herit)
}

get_null_cis_context_heritabilities <- function(numContexts) {
    return(matrix(0, nrow=numContexts, ncol=3))
}


get_null_trans_context_heritabilities <- function(numContexts) {
    return(matrix(0, nrow=numContexts, ncol=1))
}


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
            effect_size <- rnorm(1, mean = 0, sd = (sqrt(heritability/causal)))
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
simulate_one_gene <- function(numContexts, n, m, p, lam, mafs, cis_heritabilities, cis_covariance, genotype_matrix,seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    causal_snps <- get_causal_snps(m, p)
    shared_effects <- get_shared_effects(cis_heritabilities[,1], causal_snps)
    
    context_cis_expression <- vector("list", numContexts)
    context_specific_effects <- vector("list", numContexts)
    context_specific_residuals <- vector("list", numContexts)
    
    residuals <- mvrnorm(n = nrow(cis_covariance), mu = rep(0, ncol(cis_covariance)), Sigma = cis_covariance)
    specific_effects <- matrix(0, nrow = numContexts, ncol = m)
    cis_expression <- matrix(0, nrow = numContexts, ncol = n)
    for (context in 1:numContexts) {
        specific_effects[context, ] <- get_specific_effects(lam, cis_heritabilities[context, 2], causal_snps)
        # cis_expression <- get_cis_expression(genotype_matrix, shared_effects, specific_effects, cis_covariance, context)
        cis_expression[context, ] <- genotype_matrix %*% shared_effects + 
            genotype_matrix %*% specific_effects[context, ] + residuals[context, ]
    }
    return(list(shared_effects, specific_effects, cis_expression, genotype_matrix, residuals))
}



# Formatting the file
format_exp <- function(iteration, num_individuals, numContexts, cis_exp, scratch_dir){
    dir_base <- file.path(scratch_dir, iteration)
    if (!dir.exists(dir_base)) {
        dir.create(dir_base, recursive = TRUE)
    }
    
    dir_create(scratch_dir)
    cis_exp = t(cis_exp)
    individuals = c()
    
    for(i in 0:(num_individuals-1)){
        individuals = c(individuals, paste0("N", i))
    }
   
    for(context in 0:(numContexts-1)){
        
        to_write = data.frame(cis_exp[,(context+1)])
        rownames(to_write) = individuals
        file_name = paste0(dir_base, "/", context, "_cis.txt")
        fwrite(to_write, file_name, row.names = TRUE, col.names = FALSE, 
               quote = FALSE, sep = "\t")
       
    }
    
    
}

format_geno <- function(num_individuals, numGeno, geno_dir, genotypes, iteration){
    dir_base <- file.path(geno_dir, iteration)
    if (!dir.exists(dir_base)) {
        dir.create(dir_base, recursive = TRUE)
    }
    individuals = c()
    ## create directories for each gene
    for(i in 0:(num_individuals-1)){
        individuals = c(individuals, paste0("N", i))
    }
    for(index in 0:(numGeno-1)){
        to_write = data.frame(genotypes)
        rownames(to_write) = individuals
        file_name = file.path(dir_base, "genotypes.txt")  
        fwrite(to_write, file_name, row.names = TRUE, col.names = FALSE, 
               quote = FALSE, sep = "\t")
        
    }
}

### trans_simulation functions ###

# input params - number of contexts, context heritabilities, rho_mat: correlation matrix
# output: context covariance
get_trans_covariance <- function(numContexts, contexts_herit, rho){
    rho_mat = matrix(rho, nrow = numContexts, ncol = numContexts)
    diag(rho_mat) = 1
    contexts_variance = rep(0, numContexts)
    for(item in 1:nrow(contexts_herit)){
        herit = contexts_herit[item,1]
        contexts_variance[item] = 1-herit
    }
    sigma = sqrt(contexts_variance)
    covariance = rho_mat * (matrix(sigma, nrow = numContexts)  %*%  matrix(sigma, nrow = 1))
    return(covariance)
}


# returns a 1 by m list of 0s and 1s indicating whether the target trans gene at that index has a causal trans genetic 
# effect from the regulator gene
# input params - t is the number of target genes, p is the probability a target gene has a trans effect
# outputs: list of 0s and 1s 
get_causal_genes = function(t, p){
    gene_list = rbinom(t, 1, p)
    return(gene_list)
}





# returns simulated trans expression for one gene and one context 
# input params - predicted_exp: predicted regulator gene expression, 
# t: number of target genes
# context covariance matrix
# context number
# alt: indicates whether to has trans effect or not
# outputs expression vector for target genes
get_trans_expression = function(predicted_exp, t, covariance, context, trans_herit, alt, residuals){
    target_trans_expression = matrix(0, nrow = length(predicted_exp), ncol = t)
    for(gene in 1:t){
        residuals_vec = residuals[,context]
        var_residual = var(residuals_vec)
        var_cisExp = var(predicted_exp, na.rm = T)
        if(alt){
            trans_effect = sqrt((var_residual*trans_herit)/(var_cisExp-(var_cisExp*trans_herit)))
        }
        else{
            trans_effect = 0
            residuals_vec = scale(residuals_vec)
        }
        expression = predicted_exp*trans_effect + residuals_vec
        target_trans_expression[,gene] = expression
    }
    
    return(list(target_trans_expression, residuals, trans_effect))
}


# Set the trans effects
is_trans_effect <- function(ntransT, numContexts, seed) {
    set.seed(seed)
    is_trans_effect_vec <- rep(0, numContexts)
    is_trans_effect_vec <- rbinom(numContexts, 1, ntransT/numContexts)
    return(is_trans_effect_vec)
}


# numSamples is the same as number of individuals
simulateTransExpression <- function(cis_expression_dir, numSamples, covariance,
                                    contexts_cis_herit, trans_herit, is_trans_effect_vec, effects_output, trans_exp_dir) {
    cis_exp_genes = list.files(cis_expression_dir)
    for (gene in cis_exp_genes) {
        reg_gene_name <- gene
        contexts <- list.files(paste0(cis_expression_dir, gene, "/"))
        trans_exp_df <- data.frame(matrix(nrow = numSamples, ncol = length(contexts)))
        residuals <- mvrnorm(numSamples, rep(0, nrow(covariance)), covariance) # dimension: numSamples x nrow(covariance)
        colnames(trans_exp_df) <- sapply(contexts, function(c) { strsplit(c, "_")[[1]][1] })
        
        for (c in contexts) {
            context <- as.numeric(strsplit(c, "_")[[1]][1]) + 1 # Assuming context filenames are like '0_filename'
            cis_exp <- fread(paste0(cis_expression_dir, gene, "/", c), sep = "\t", data.table = FALSE)
            rownames(cis_exp) <- cis_exp$V1
            
            cis_herit <- contexts_cis_herit[context, 3] # Column 3 is the total heritability
            is_trans_effect <- is_trans_effect_vec[context]
            trans_expression <- get_trans_expression(cis_exp$V2, length(cis_exp_genes), covariance, context, 
                                                     trans_herit, is_trans_effect, residuals)
            
            trans_exp_df[, context] <- trans_expression[[1]] # target_trans_expression
            residuals_output <- trans_expression[[2]]
            effects_size <- trans_expression[[3]] 
  
            effects_output_file = paste0(effects_output, reg_gene_name, "_", c-1, "_trans_effects.txt")
            residuals_output_file = paste0(effects_output, reg_gene_name, "_", c-1, "_trans_residuals.txt")
            #fwrite(data.frame(trans_expression[[3]]), effects_output_file, sep = "\t", quote = F, row.names = F, col.names = F)
            #fwrite(data.frame(trans_expression[[2]]), residuals_output_file, sep = "\t", quote = F, row.names = F, col.names = F)
            
    
        }
        
        trans_output_file <- paste0(trans_exp_dir, reg_gene_name, "_trans_exp.txt")
        fwrite(trans_exp_df, trans_output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      
    }
}












