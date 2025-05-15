# C-STEM R package


### Necessary packages:
```
library(data.table) 
library(dplyr) 
library(reshape2) 
library(data.table) 
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
```

### Install C-STEM via github:
```
devtools::install("BalliuLab/FastGxC"")
```

### Step 1: Build GReXs
#### This step decomposes expression and builds cross-validated cis genetic predictors of expression using elastic net regularized regression.

#### example code:
## gene 1:
```
X_file = "example_data/gene1.genos"
exp_files = list.files("example_data/expression/gene1/")
contexts = exp_files
exp_files = paste0("example_data/expression/gene1/", exp_files)
out_dir = "example_data/GReXs/"
decomposition_dir = "example_data/decomposed_expression/"
gene_name = "gene1"
context_thresh = 3
alpha = 0.5
num_folds = 10
run_GBAT = FALSE

create_GReXs(X_file, exp_files, contexts, out_dir, gene_name, decomposition_dir, context_thresh, alpha, num_folds, run_GBAT)
```

## gene 2:
```
X_file = "example_data/gene2.genos"
exp_files = list.files("example_data/expression/gene2/")
contexts = exp_files
exp_files = paste0("example_data/expression/gene2/", exp_files)
out_dir = "example_data/GReXs/"
decomposition_dir = "example_data/decomposed_expression/"
gene_name = "gene2"
context_thresh = 3
alpha = 0.5
num_folds = 10
run_GBAT = FALSE

create_GReXs(X_file, exp_files, contexts, out_dir, gene_name, decomposition_dir, context_thresh, alpha, num_folds, run_GBAT)
```

### Step 2: Run C-STEM lite
```
regulator_pred_exp_file = "example_data/GReXs/gene1_cstem_full_predictors.txt"
regulator_gbat_pred_exp_file = "example_data/GReXs/gene1_gbat_predictors.txt"
target_pred_exp_file = "example_data/GReXs/gene1_cstem_predictors.txt"
target_gbat_pred_exp_file = "example_data/GReXs/gene1_gbat_predictors.txt"
target_exp_files = list.files("example_data/expression/")
run_GBAT = F
contexts_vec = target_exp_files
target_exp_files = paste0("example_data/expression/", target_exp_files)
regulator_gene_name = "gene1"
target_gene_name = "gene2"
target_cis_pred = TRUE
outdir = "example_data/trans_output/"

cstem_gbat_lite(regulator_pred_exp_file, target_pred_exp_file, target_exp_files, contexts_vec, run_GBAT, regulator_gene_name, target_gene_name, outdir, target_cis_pred = T)
```

### Step 3: Run C-STEM lmm
```
regulator_pred_exp_file = "example_data/GReXs/gene1_cstem_full_predictors.txt"
target_pred_exp_file = "example_data/GReXs/gene1_cstem_predictors.txt"
target_exp_files = list.files("example_data/expression/")
contexts_vec = target_exp_files
target_exp_files = paste0("example_data/expression/", target_exp_files)
regulator_gene_name = "gene1"
target_gene_name = "gene2"
target_cis_pred = TRUE
outdir = "example_data/trans_output/"

cstem_lmm(regulator_pred_exp_file, target_pred_exp_file, target_exp_files, contexts_vec, regulator_gene_name, target_gene_name, outdir, target_cis_pred = T){
```










