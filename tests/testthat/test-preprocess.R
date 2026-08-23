# preprocess_expression: per-gene assembly happy path + the two loud-failure
# guards (duplicate gene IDs after version stripping; non-numeric cells
# silently coerced to NA).

.write_ctx_file <- function(path, mat) {
  df <- data.frame(gene_id = rownames(mat), mat, check.names = FALSE)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

test_that("preprocess_expression assembles per-gene matrices (happy path)", {
  td <- file.path(tempdir(), "prep_ok"); unlink(td, recursive = TRUE)
  dir.create(td, recursive = TRUE)
  set.seed(2)
  m1 <- matrix(rnorm(6), 3, 2, dimnames = list(c("gA", "gB", "gC"),
                                               c("i1", "i2")))
  m2 <- matrix(rnorm(6), 3, 2, dimnames = list(c("gA", "gB", "gC"),
                                               c("i2", "i3")))
  .write_ctx_file(file.path(td, "liver.txt"),  m1)
  .write_ctx_file(file.path(td, "spleen.txt"), m2)
  out <- preprocess_expression(
    expr_files = c(file.path(td, "liver.txt"), file.path(td, "spleen.txt")),
    output_dir = file.path(td, "genes"), gene_id_col = "gene_id",
    verbose = FALSE)
  gA <- readRDS(file.path(td, "genes", "gA.rds"))
  expect_setequal(rownames(gA), c("i1", "i2", "i3"))
  expect_identical(colnames(gA), c("liver", "spleen"))
  expect_equal(gA["i1", "liver"], m1["gA", "i1"])
  expect_true(is.na(gA["i1", "spleen"]))    # i1 unmeasured in spleen
})

test_that("duplicate gene IDs after version stripping stop loudly", {
  td <- file.path(tempdir(), "prep_dup"); unlink(td, recursive = TRUE)
  dir.create(td, recursive = TRUE)
  m <- matrix(rnorm(4), 2, 2,
              dimnames = list(c("ENSG000123.4", "ENSG000123.7"),
                              c("i1", "i2")))
  .write_ctx_file(file.path(td, "liver.txt"), m)
  expect_error(
    preprocess_expression(expr_files = file.path(td, "liver.txt"),
                          output_dir = file.path(td, "genes"),
                          gene_id_col = "gene_id",
                          strip_version = TRUE, verbose = FALSE),
    "collapsed multiple entries")
})

test_that("non-numeric expression cells stop loudly instead of becoming NA", {
  td <- file.path(tempdir(), "prep_txt"); unlink(td, recursive = TRUE)
  dir.create(td, recursive = TRUE)
  writeLines(c("gene_id\ti1\ti2",
               "gA\t0.5\thigh",        # a text cell
               "gB\t1.2\t0.3"),
             file.path(td, "liver.txt"))
  expect_error(
    preprocess_expression(expr_files = file.path(td, "liver.txt"),
                          output_dir = file.path(td, "genes"),
                          gene_id_col = "gene_id", verbose = FALSE),
    "non-numeric")
})
