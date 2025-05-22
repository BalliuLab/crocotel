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
devtools::install_github("BalliuLab/C-STEM")
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

### Step 2a: Run C-STEM lite
```
regulator_pred_exp_file = "example_data/GReXs/gene1_cstem_full_predictors.txt"
regulator_gbat_pred_exp_file = "example_data/GReXs/gene1_gbat_predictors.txt"
target_pred_exp_file = "example_data/GReXs/gene2_cstem_full_predictors.txt"
target_gbat_pred_exp_file = "example_data/GReXs/gene2_gbat_predictors.txt"
target_exp_files = list.files("example_data/expression/gene2/")
run_GBAT = F
contexts_vec = target_exp_files
target_exp_files = paste0("example_data/expression/gene2/", target_exp_files)
regulator_gene_name = "gene1"
target_gene_name = "gene2"
target_cis_pred = TRUE
outdir = "example_data/trans_output/"

cstem_gbat_lite(regulator_pred_exp_file, target_pred_exp_file, target_exp_files, contexts_vec, run_GBAT, regulator_gene_name, target_gene_name, outdir, target_cis_pred)
```

### Step 2b: Run C-STEM lmm
```
regulator_pred_exp_file = "example_data/GReXs/gene1_cstem_full_predictors.txt"
target_pred_exp_file = "example_data/GReXs/gene2_cstem_full_predictors.txt"
target_exp_files = list.files("example_data/expression/gene2/")
contexts_vec = target_exp_files
target_exp_files = paste0("example_data/expression/gene2/", target_exp_files)
regulator_gene_name = "gene1"
target_gene_name = "gene2"
target_cis_pred = TRUE
outdir = "example_data/trans_output/"

cstem_lmm(regulator_pred_exp_file, target_pred_exp_file, target_exp_files, contexts_vec, regulator_gene_name, target_gene_name, outdir, target_cis_pred)
```

### Step 3: Multiple testing correction
```

```









