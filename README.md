# Crocotel R package


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
library(MuMIn)
```

### Install Crocotel via github:
##### Note: qvalue and TreeQTL must be installed from source before Crocotel is installed 
```
install.packages("qvalue")
install.packages("", repos = NULL, type = "source")
devtools::install_github("BalliuLab/Crocotel", dependencies = TRUE)
library(Crocotel)
```

### Step 1: Build GReXs
#### This step decomposes expression and builds cross-validated cis genetic predictors of expression using elastic net regularized regression.

#### example code for one gene:
```
X_file = "crocotile_example/input_data/gene1_genotypes.txt"
exp_files = list.files("crocotile_example/input_data/gene1/")
exp_files = paste0("crocotile_example/input_data/gene1/", exp_files)
out_dir = "crocotile_example/GReXs/"
gene_name = "gene1"
context_thresh = 3
alpha = 0.5
num_folds = 10
run_GBAT = FALSE

create_GReXs(X_file, exp_files, out_dir, gene_name, context_thresh, alpha, num_folds, run_GBAT)
```

### Step 2a: Run Crocotel lite
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

### Step 2b: Run Crocotel lmm
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
#### example using Crocotel lmm output
```
crocotel_sum_stats = "crocotile_example/all_gene_pairs.crocotel_lmm.txt"
contexts_vec = as.character(seq(0,9))
fdr_thresh = 0.05
outdir = "crocotile_example/treeQTL_output/"
method = "mashr"

multiple_testing_correction(crocotel_sum_stats, contexts_vec, fdr_thresh, outdir, method)

```

### Important notes:
1. Input file names for expression of each gene must have its own directory with no other files. The names of the files should be "contextName.txt" or "contextName" with no other leading or trailing strings.









