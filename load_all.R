# load_all.R
# -----------
# Sources all simulation functions. Run this once at the top of any
# analysis script or interactive session before using the simulation pipeline.
#
# Usage:
#   source("load_all.R")

local({
  r_dir <- file.path(dirname(sys.frame(1)$ofile), "R")
  if (!dir.exists(r_dir))
    r_dir <- file.path(getwd(), "R")   # fallback when sourced interactively

  files <- c(
    "utils.R",
    "generate_genotypes.R",
    "simulate_regulator_expression.R",
    "simulate_target_expression.R",
    "simulate_expression.R",
    "write_simulated_genotypes.R",
    "write_simulated_expression.R"
  )

  for (f in files) {
    path <- file.path(r_dir, f)
    if (!file.exists(path))
      stop("Cannot find ", path, "\nMake sure load_all.R is in the project root.")
    source(path)
  }

  message("Loaded: ", paste(files, collapse = ", "))
})
