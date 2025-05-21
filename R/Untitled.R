library(data.table)

expression = list.files("C-STEM/example_data/expression/gene1/")
outdir = "C-STEM/example_data/expression/gene2/"

num_effect_sizes = length(expression)


effect_sizes = rnorm(n = 17, mean = 0, sd = 0.2)
for(i in 1:length(expression)){
  file = expression[i]
  cur_df = fread(paste0("C-STEM/example_data/expression/gene1/",file), sep = "\t", data.table = F)
  new_expression = (cur_df$V2 * effect_sizes[i]) + rnorm(1, mean = 0, sd = 1)
  out_df = data.frame(id = cur_df$V1, new_expression)
  fwrite(out_df, file = paste0(outdir, file), sep = "\t", row.names = F, col.names = F)
}

################################################################################
### simulate data and run mash

library(ashr)
library(mashr)
set.seed(1)
simdata = simple_sims(500,5,1)

data = mash_set_data(simdata$Bhat, simdata$Shat)
U.c = cov_canonical(data)  
print(names(U.c))

m.c = mash(data, U.c)
head(get_significant_results(m.c))



