# C-STEM R package


### Necessary packages:
```
library(data.table) 
library(dplyr) 
library(reshape2) 
library(magrittr) 
library(bigstatsr)
library(caret) 
library(devtools) 
library(purrr) 
library(foreach)
library(tidyr)
library(lme4)
library(emmeans)
library(broom)
library(TreeQTL)
library(mashr)
```

### Install Crocotile via github:
##### Note: qvalue and TreeQTL must be installed from source before Crocotile is installed 
```
install.packages("qvalue")
install.packages("", repos = NULL, type = "source")
devtools::install_github("BalliuLab/C-STEM", dependencies = TRUE)
library(cstem)
```

### Step 1: Build GReXs
#### This step decomposes expression and builds cross-validated cis genetic predictors of expression using elastic net regularized regression.

#### example code for one gene:
```
X_file = "crocotile_example/input_data/gene1_genotypes.txt"
exp_files = list.files("crocotile_example/input_data/gene1/")
exp_files = paste0("crocotile_example/input_data/gene1/", exp_files)
contexts = as.character(seq(0,9))
out_dir = "crocotile_example/GReXs/"
gene_name = "gene1"
context_thresh = 3
alpha = 0.5
num_folds = 10
run_GBAT = FALSE

create_GReXs(X_file, exp_files, contexts, out_dir, gene_name, context_thresh, alpha, num_folds, run_GBAT)
```

### Step 2a: Run Crocotile lite
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

cstem_gbat_lite(regulator_pred_exp_file, target_pred_exp_file, target_exp_files, contexts_vec, method, regulator_gene_name, target_gene_name, outdir, target_cis_pred)
```

### Step 2b: Run Crocotile lmm
```
regulator_pred_exp_file = "crocotile_example/GReXs/gene1.cstem.full_predictors.txt"
target_pred_exp_file = "crocotile_example/GReXs/gene2.cstem.full_predictors.txt"
target_exp_files = list.files("crocotile_example/input_data/gene2/", full.names = T)
contexts_vec = as.character(seq(0,9))
regulator_gene_name = "gene1"
target_gene_name = "gene2"
target_cis_pred = TRUE
outdir = "crocotile_example/trans_output/"

cstem_lmm(regulator_pred_exp_file, target_pred_exp_file, target_exp_files, contexts_vec, regulator_gene_name, target_gene_name, outdir, target_cis_pred)
```

### Step 3: Multiple testing correction
#### example using crocotile lmm output
```
m_eqtl_outfiles = list.files("crocotile_example/treeQTL_input/", pattern = "all_gene_pairs.crocotile_lmm", full.names = T)
n_SNPs_per_gene_files = list.files("crocotile_example/treeQTL_input/", pattern = "n_tests_per_gene.crocotile_lite", full.names = T)
contexts_vec = as.character(seq(0,9))
fdr_thresh = 0.05
outdir = "crocotile_example/treeQTL_output/"
method = "treeQTL"
exp_suffix = "crocotile_lmm"

multiple_testing_correction(m_eqtl_outfiles, n_SNPs_per_gene_files, contexts_vec, fdr_thresh, outdir, method, top_level, exp_suffix)

```

### Important notes:
1. Input file names for expression of each gene must have its own directory with no other files. The names of the files should be "contextName.txt" or "contextName" with no other leading or trailing strings.









