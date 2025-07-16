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
#### This step takes in input files formatted exactly as in MatrixEQTL (1 expression file for each context, 1 genotype file, 1 snpsloc file, and 1 geneloc file)
#### The output creates a directory called "crocotel_formatted_data/" inside the specified output directory with a subdirectory for the expression of each gene in each context. 
```
exp_files=list.files("crocotel_example/input_data/", pattern = ".exp.txt", full.names = T)
geneloc_file="crocotel_example/input_data/geneloc.txt"
snpsloc_file="crocotel_example/input_data/snpsloc.txt"
genotypes_file="crocotel_example/input_data/all_genotypes.txt"

# output directory where all crocotel final and intermediate output will be stored. Only set this directory once.
out_dir="crocotel_example/"

format_data(exp_files, geneloc_file, snpsloc_file, genotypes_file, out_dir)
```

### Step 1: Build cis Genetically Regulated eXpression componentS (GReXs)
#### This step builds cross-validated cis genetic predictors of expression for a gene across all contexts (e.g. cell types and tissues) using elastic net regularized regression.
#### If the preliminary format data step was run, only the gene name and the crocotel "out_dir" set above is a required parameter. 
#### Filenames can also be passed in for all parameters if they were not created by the above function. 

#### example code for one gene:
```
out_dir="crocotel_example/"
gene_name = "gene1"
context_thresh = 3 # minimum # of contexts a gene has to have expression on 
alpha = 0.5 # elastic net mixture parameter 
num_folds = 10 # number of folds for cross-validation 
method = "crocotel" # alternative is "cxc"

create_GReXs(gene_name, out_dir, 
  context_thresh = context_thresh,
  alpha = alpha,
  num_folds = num_folds,
  method = method)
  
```

### Step 2: Run regulator-target associations 
##### Option 1 (faster but less powerful): Lite version which tests all regulator-target pairs simultaneously in each context using ultra-fast linear regression 
##### The p-value output threshold 
```
out_dir="crocotel_example/"
context = "0"
geneloc_file = "crocotel_example/input_data/geneloc.txt"

crocotel_lite(context, geneloc_file, out_dir)
```

##### Option 2 (slower but more powerful): LMM version which tests each regulator-target pairs separately but jointly models all contexts for each pair
```

regulator_gene_name = "gene1"
target_gene_name = "gene2"
out_dir = "crocotel_example/"

crocotel_lmm(regulator_gene_name = regulator_gene_name, target_gene_name = target_gene_name, out_dir = out_dir)
```

##### The implementation of crocotel lmm requires user parallelization to run across all regulator-target pairs. Once association for all pairs is run, use this function to concatenate output files for multiple testing correction

```
out_dir = "crocotel_example/"
regress_target_GReX = T
concat_crocotel_lmm_files(directory = out_dir, regress_target_GReX = regress_target_GReX)
```


### Step 3: Multiple testing correction
#### example using crocotel lite output with treeQTL:


```
method = "treeQTL"
crocotel_dir = "crocotel_example/crocotel_lite_output/"
out_dir = "crocotel_example/"
top_level = "R"

multiple_testing_correction(crocotel_dir = crocotel_dir, out_dir = out_dir, method = method, top_level = top_level)

```

#### example using crocotel lite output with mashr:
```
method = "mashr"
crocotel_dir = "crocotel_example/crocotel_lite_output/"
out_dir = "crocotel_example/"

multiple_testing_correction(crocotel_dir = crocotel_dir, out_dir = out_dir, method = method)

```

### Important notes:
1. Input file names for total expression should be in the following format with the following filenames "contextName.txt":
    
2. Input file names for expression of each gene must have its own directory with no other files. The names of the files should be "contextName.txt" or "contextName" with no other leading or trailing strings.









