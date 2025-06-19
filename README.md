# _crocotel_ R package


### Install _crocotel_ via Github:
##### Note: qvalue and TreeQTL must be installed from source before crocotel is installed 
```

# Install qvalue
if (!requireNamespace("qvalue", quietly = TRUE)) {
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("qvalue")
}

# Install treeQTL
install.packages("http://bioinformatics.org/treeqtl/TreeQTL_2.0.tar.gz", repos = NULL, type = "source")

# Install crocotel
devtools::install_github("BalliuLab/crocotel", dependencies = TRUE)
library(crocotel)
```
### Preliminary step: Format the input data
#### This step takes in input files formatted exactly as in MatrixEQTL (1 expression files for each context, 1 genotype file, 1 snpsloc file, and 1 geneloc file)
```
exp_files=list.files("crocotel_example/input_data/", pattern = ".exp.txt")
geneloc_file="crocotel_example/input_data/geneloc.txt"
snpsloc_file="crocotel_example/input_data/snpsloc.txt"
genotypes_file="crocotel_example/input_data/all_genotypes.txt"
out_dir="crocotel_example/"


```

### Step 1: Build cis Genetically Regulated eXpression componentS (GReXs)
#### This step builds cross-validated cis genetic predictors of expression for a gene across all contexts (e.g. cell types and tissues) using elastic net regularized regression.

#### example code for one gene:
```
gene_name="gene1"
expression_directory=paste0("crocotile_example/input_data/",gene_name,"/")
exp_files = paste0(expression_directory, list.files(expression_directory))
genotype_file = paste0("crocotile_example/input_data/",gene_name,"_genotypes.txt") # maybe change this one

out_dir = "crocotile_example/"

context_thresh = 3 # minimum # of contexts a gene has to have expression on 
alpha = 0.5 # elastic net mixture parameter 
num_folds = 10 # number of folds for cross-validation 
method = "crocotel"

create_GReXs(genotype_file, exp_files, out_dir, gene_name, context_thresh, alpha, num_folds, method)
```

### Step 2: Run regulator-target associations 
##### Option 1 (faster but less powerful): Lite version which tests all regulator-target pairs simultaneously in each context using ultra-fast linear regression 
```
regulator_pred_exp_file = "crocotile_example/GReXs/gene1.cstem.full_predictors.txt"
target_pred_exp_file = "crocotile_example/GReXs/gene2.cstem.full_predictors.txt"
target_exp_files = list.files("crocotile_example/input_data/gene2/", full.names = T)
method = "Crocotel"
contexts_vec = as.character(seq(0,9))
regulator_gene_name = "gene1"
target_gene_name = "gene2"
target_cis_pred = TRUE
outdir = "crocotile_example/trans_output/"

crocotel_lite(regulator_pred_exp_file, target_exp_files, contexts_vec, regulator_gene_name, target_gene_name, outdir, method, target_cis_pred, target_pred_exp_file)
```

##### Option 2 (slower but more powerful): LMM version which tests each regulator-target pairs separately but jointly models all contexts for each pair 
```
regulator_pred_exp_file = "crocotile_example/GReXs/gene1.cstem.full_predictors.txt"
target_pred_exp_file = "crocotile_example/GReXs/gene2.cstem.full_predictors.txt"
target_exp_files = list.files("crocotile_example/input_data/gene2/", full.names = T)
contexts_vec = as.character(seq(0,9))
regulator_gene_name = "gene1"
target_gene_name = "gene2"
target_cis_pred = TRUE
outdir = "crocotile_example/trans_output/"

crocotel_lmm(regulator_pred_exp_file, target_exp_files, contexts_vec, regulator_gene_name, target_gene_name, outdir, target_cis_pred, target_pred_exp_file)
```

### Step 3: Multiple testing correction
#### example using crocotel lmm output
```
crocotel_sum_stats = "crocotile_example/all_gene_pairs.crocotel_lmm.txt"
contexts_vec = as.character(seq(0,9))
fdr_thresh = 0.05
outdir = "crocotile_example/treeQTL_output/"
method = "mashr"

multiple_testing_correction(crocotel_sum_stats, contexts_vec, fdr_thresh, outdir, method)

```

### Important notes:
1. Input file names for total expression should be in the following format with the following filenames "contextName.txt":
    
2. Input file names for expression of each gene must have its own directory with no other files. The names of the files should be "contextName.txt" or "contextName" with no other leading or trailing strings.









