# Load R packages

library(data.table)
library(reshape2)
library(data.table)
library(reshape2)
library(magrittr)

# load functions
source('R/decompose_expression.R')
source('R/create_GReXs.R')

# load data
list.files('example_data')
Adrenal_Gland=read.table('example_data/Adrenal_Gland',header = F)

