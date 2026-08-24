# crocotel

Context-resolved **trans**-eQTL mapping via genetically regulated expression
(GReX). crocotel builds cross-validated GReX predictors with a CONTENT-style
shared / context-specific decomposition and an elastic net, then tests each
regulator gene's GReX against every distal target gene across contexts
(tissues, cell types) — either per context (`run_trans_eqtl()`) or jointly
across contexts (`run_trans_lmm()`) — with hierarchical treeQTL-style FDR
control, validated exactly against TreeQTL 2.0 in both hierarchy orientations.

## Installation

```r
# install.packages("remotes")
remotes::install_github("BalliuLab/crocotel")
```

Dependencies (`MatrixEQTL`, `bigsnpr`, `glmnet`, `data.table`) install from
CRAN automatically. On a cluster without compilers, `install_deps.sh` sets up
a conda environment with everything pre-built. Requires R >= 4.1.

## Pipeline at a glance

```
expression (per context)          genotypes (PLINK bed / plink2 .traw)
        |                                    |
  preprocess_expression()            prepare_genotypes()
        |                                    |
        +--------------+---------------------+
                       |
              fit_grex_gene() / fit_grex_batch()     one gene per fit /
                       |                             a batch of genes per task
              assemble_grex_matrices()               per-context matrices
                       |                             + target-eligibility file
            +----------+----------+
            |                     |
    run_trans_eqtl()        run_trans_lmm()
    (per-context OLS)      (cross-context LMM)
            |                     |
            +----------+----------+
                       |
              apply_crossmap_post()                  cross-mappability filter
                       |
                   run_fdr()                         3-level hierarchical FDR
```

GReX is fit one gene at a time, and the fits are independent — so this is the
step you parallelize. `fit_grex_gene()` fits a single gene;
`fit_grex_batch()` fits a *batch* of genes inside one job, spreading them over
that job's cores (`mc.cores`). The intended cluster pattern is therefore two
levels: an array task per batch of genes, and `mc.cores` genes at a time
within each task.

## Quick start

```r
library(crocotel)

## 1. One file per gene from per-context expression tables
preprocess_expression(
  expr_files  = c("liver.txt.gz", "muscle.txt.gz", "blood.txt.gz"),
  output_dir  = "expr_by_gene",
  gene_id_col = "gene_id")

## 2. Convert genotypes to a memory-mapped backing (once)
prepare_genotypes("genotypes/cohort")          # cohort.bed/.bim/.fam

## 3. Fit GReX -- one array task per batch of genes, mc.cores genes per task
fit_grex_batch(
  gene_ids       = my_genes,                   # this task's batch
  expr_dir       = "expr_by_gene",
  plink_prefix   = "genotypes/cohort",
  gene_locations = "gene_locations.txt",       # gene_id, chr, start, end
  output_dir     = "grex",
  mc.cores       = as.integer(Sys.getenv("NSLOTS", "1")))

## 4. Assemble per-context matrices (also decides target eligibility)
assemble_grex_matrices(
  grex_dir = "grex", expr_dir = "expr_by_gene", output_dir = "matrices")

## 5. Trans scan -- per-context OLS ...
run_trans_eqtl(matrix_dir = "matrices", gene_locations = "gene_locations.txt",
               output_dir = "trans", method = "crocotel")

##    ... or the joint cross-context test (same inputs, same output layout)
run_trans_lmm(matrix_dir = "matrices", gene_locations = "gene_locations.txt",
              output_dir = "trans_lmm")

## 6. Cross-mappability filter (required on real data), then FDR
apply_crossmap_post("trans", "crocotel", regulator = "gene",
                    cross_map_file = "crossmap_pairs.txt.gz",
                    gene_locations = "gene_locations.txt")
run_fdr(trans_dir = "trans", output_dir = "fdr", method = "crocotel")
```

`run_fdr()` writes the discoveries at three levels: `eTargets_*` (which genes
are trans-regulated at all), `eTarget_context_*` (in which contexts), and
`triplets_*` (by which regulator).

## The two trans tests

| Function | Test | Notes |
|---|---|---|
| `run_trans_eqtl(method="crocotel")` | per-context OLS on GReX | one test per (regulator, target, context); `target_response = "residualized"` (default) removes the target's own cis signal first |
| `run_trans_lmm()` | joint cross-context score test (het-CS noise model) | pools evidence across contexts, so it is better powered when an effect is shared across them; the het-CS noise model is fixed, not a user knob |

Both consume the same `matrix_dir`, obey the same regulator gate and
target-eligibility file, and feed the same `apply_crossmap_post()` /
`run_fdr()` steps — so their discoveries are directly comparable. The LMM
writes `trans_lmm_<ctx>.tsv`, so its downstream calls take `method = "lmm"`.

## Key concepts

- **GReX** — the genetic component of a gene's expression, predicted from cis
  variants by cross-validated elastic net; testing GReX (not observed
  expression) against distal genes protects trans hits from environmental
  confounding. The decomposition splits it into a **shared** (cross-context)
  and **context-specific** part; the fitted effect sizes and SNP support are
  saved on every record.
- **Target eligibility** — decided once at assembly (`expressed_targets.rds`):
  a gene is a target in a context only with >= 30 observed expression values
  there (`min_obs_per_ctx`). Every scanner obeys the same file, so FDR
  families are identical across methods by construction.
- **Regulator gate** — a gene is admitted as a regulator in a context only if
  its GReX is significantly predictive there (`grex_gate_mode = "pval"`,
  default), clears a prediction-R2 floor (`"r2"`, the GBAT criterion), or
  both.
- **Cross-mappability** — on real data, `run_fdr()` refuses to run unless the
  (regulator, target) pairs confounded by sequence similarity (direct +
  LD-halo proximity) have been removed by `apply_crossmap_post()`. Pass
  `crossmap = NA` only when no cross-mappability applies by construction
  (simulations).
- **Hierarchy orientation** — the default tree asks "which genes are
  trans-regulated?" (eTargets). `hierarchy = "regulator"` flips it to "which
  genes are trans-regulators?" (eRegulators); `write_n_tests()` converts an
  existing scan to the other orientation in seconds, no re-scan.

## Simulating data

`simulate_expression()` generates multi-context regulator + target expression
with known trans effects (`n_active_contexts` controls in how many contexts
each signal is active); `write_simulated_genotypes()` /
`write_simulated_expression()` write it in exactly the pipeline's input
format. See the examples in their help pages — they form one coherent flow.

## Comparator methods

The package also ships the published baselines crocotel is benchmarked
against, so a comparison runs on identical eligibility rules, identical FDR
families, and identical inputs rather than on a reimplementation.

| Function | Test | Regulator unit | Notes |
|---|---|---|---|
| `run_trans_eqtl(method="cbc")` | per-context OLS on context-by-context GReX | gene | GBAT-style comparator |
| `run_trans_eqtl_snp(snp_method="lead")` | per-context OLS on each gene's lead cis-SNP | gene | GBAT "top cis-eQTL" baseline |
| `run_trans_eqtl_snp(snp_method="genome_wide")` | per-context OLS on every variant | variant | GTEx/eQTLGen-style baseline |

`run_trans_eqtl_snp()` needs no GReX: it takes its predictors straight from
the genotypes and reads `matrix_dir` only for the raw per-context expression
(`expr_<ctx>.rds`) and the target-eligibility sidecar. To build that
`matrix_dir` without ever fitting GReX, call `assemble_grex_matrices()` in
expression-only mode -- leave `grex_dir` and `grex_list` unset and pass only
`expr_dir` -- which writes the expression matrices and the same
`expressed_targets.rds` sidecar the GReX assembly would.

## License

MIT. Please file issues at
<https://github.com/BalliuLab/crocotel/issues>.
