###### make sample data into correct format

library(data.table)
library(dplyr)

files = list.files("example_data/")
files = files[!grepl(".R", files)]
files = files[!grepl("genos", files)]
files = files[!grepl("snp", files)]

final_df = fread(paste0("example_data/", files[1]), sep = "\t", data.table = F)
names(final_df) = c("ind_ID", files[1])
for(i in 2:length(files)){
    cur_file = fread(paste0("example_data/", files[i]), sep = "\t", data.table = F)
    names(cur_file) = c("ind_ID", files[i])
    final_df = merge(final_df, cur_file, by = "ind_ID", all = TRUE)
}
fwrite(final_df, "example_data/expression.txt", sep = "\t", quote = F)