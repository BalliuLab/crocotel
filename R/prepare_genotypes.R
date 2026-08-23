# prepare_genotypes.R
# -------------------
# One-time conversion of a genotype source to the bigSNP backing files
# (.rds / .bk) required by load_genotypes(). Two input formats:
#   * "bed"    - PLINK bed/bim/fam hardcalls (0/1/2) via bigsnpr::snp_readBed.
#   * "dosage" - a plink2 --export A-transpose (.traw) dosage matrix (imputed
#                data, e.g. GTEx/OneK1K). Dosages are stored in the same
#                bigSNP FBM.code256 structure via bigsnpr::CODE_DOSAGE, so all
#                downstream code (load_genotypes, fit_grex_gene) is unchanged.
#
# The upstream conversion to the ingestable intermediate stays OUTSIDE this
# package (a heavy external-tool step): plink turns a VCF into .bed for the
# hardcall path; plink2 turns an imputed VCF's DS field into a .traw for the
# dosage path (see simulation_study/prep_dosage_genotypes.sh). This function
# only builds the mmap FBM from that intermediate.


#' Prepare genotype files (hardcalls or dosages) for use with crocotel
#'
#' Converts a genotype source to the bigSNP backing files (.rds and .bk)
#' required by \code{load_genotypes()}. One-time per dataset; once the backing
#' files exist, all parallel GReX-fitting jobs read from them concurrently.
#'
#' \code{format = "bed"} converts a PLINK bed/bim/fam fileset (0/1/2 hardcalls)
#' via \code{bigsnpr::snp_readBed}. \code{format = "dosage"} ingests a plink2
#' \code{--export A-transpose} matrix (\code{.traw}: one row per variant, a
#' \code{CHR SNP (C)M POS COUNTED ALT} prefix then one column per sample) of
#' imputed dosages, storing them in the same \code{FBM.code256} structure via
#' \code{bigsnpr::CODE_DOSAGE} (dosages quantized to ~0.01, the standard imputed
#' precision). The \code{.traw} is read in variant chunks so genome-wide inputs
#' never load fully into memory.
#'
#' @param plink_prefix Character. Output path prefix for the bigSNP backing
#'   files (\code{<plink_prefix>.rds} / \code{.bk}). For \code{format = "bed"}
#'   this is also the input PLINK prefix (reads \code{<plink_prefix>.bed}).
#' @param format       Character. \code{"bed"} (default) or \code{"dosage"}.
#' @param dosage_file  Character or \code{NULL}. Path to the plink2
#'   \code{.traw} dosage matrix. Required when \code{format = "dosage"}.
#' @param chunk_size   Integer. Number of variants read per chunk in the dosage
#'   path. Default 20000 (~\code{chunk_size} x n_individuals doubles in memory
#'   at a time). Only used for \code{format = "dosage"}.
#'
#' @param verbose      Logical. Print the "backing already exists" skip
#'   message. Rebuild/broken-backing messages always print. Default
#'   \code{TRUE}.
#'
#' @return Invisibly returns the path to the .rds backing file. An existing
#'   backing is reused only when it is complete (\code{.rds} and
#'   \code{.bk}) and not older than the genotype source; a source
#'   regenerated after the backing was built triggers an automatic rebuild
#'   (so a stale backing can never silently serve outdated genotypes), and
#'   a broken backing is rebuilt when the source is available or reported
#'   with instructions when it is not.
#'
#' @examples
#' \dontrun{
#' # hardcalls
#' prepare_genotypes("/scratch/gtex/geno")
#' # imputed dosages (after: plink2 --vcf imputed.vcf.gz dosage=DS \\
#' #   --export A-transpose --out /scratch/gtex/geno_dosage)
#' prepare_genotypes("/scratch/gtex/geno", format = "dosage",
#'                   dosage_file = "/scratch/gtex/geno_dosage.traw")
#' }
#' @export
prepare_genotypes <- function(plink_prefix,
                              format      = c("bed", "dosage"),
                              dosage_file = NULL,
                              chunk_size  = 20000L,
                              verbose     = TRUE) {

  format <- match.arg(format)

  if (!requireNamespace("bigsnpr", quietly = TRUE))
    stop("Package 'bigsnpr' is required: install.packages('bigsnpr')")

  # ------------------------------------------------------------------
  # Backing validity + currency guard (THE single authority; the loaders
  # call this instead of converting inline). A backing is served only when
  # it is COMPLETE (.rds AND .bk) and NOT older than the genotype source:
  #   * source newer than backing -> the genotypes were regenerated; the
  #     old backing would silently serve OUTDATED genotypes -> rebuild.
  #   * one backing file missing (scratch purge, partial copy) -> broken;
  #     rebuild when a source exists, otherwise stop with instructions.
  # ------------------------------------------------------------------
  rds_file <- paste0(plink_prefix, ".rds")
  bk_file  <- paste0(plink_prefix, ".bk")
  src_file <- if (format == "bed") paste0(plink_prefix, ".bed")
              else dosage_file
  has_src  <- !is.null(src_file) && file.exists(src_file)
  have_rds <- file.exists(rds_file)
  have_bk  <- file.exists(bk_file)

  if (have_rds && have_bk) {
    if (has_src && file.mtime(src_file) > file.mtime(rds_file)) {
      message("Genotype source (", basename(src_file), ") is newer than ",
              "the bigSNP backing; rebuilding the backing (the old one ",
              "would silently serve outdated genotypes).")
      unlink(c(rds_file, bk_file))
    } else {
      if (verbose)
        message("bigSNP backing files already exist, skipping conversion: ",
                rds_file)
      return(invisible(rds_file))
    }
  } else if (have_rds || have_bk) {
    missing_f <- if (have_rds) bk_file else rds_file
    if (has_src) {
      message("Broken bigSNP backing (missing ", basename(missing_f),
              "); rebuilding from ", basename(src_file), ".")
      unlink(c(rds_file, bk_file))
    } else {
      stop("Broken bigSNP backing at ", plink_prefix, ": ",
           basename(missing_f), " is missing (scratch purge or partial ",
           "copy?) and no genotype source (.bed / dosage file) is ",
           "available to rebuild from. Delete the leftover backing file ",
           "and regenerate the genotypes.")
    }
  }

  if (format == "bed") {
    bed_file <- paste0(plink_prefix, ".bed")
    if (!file.exists(bed_file))
      stop("PLINK bed file not found: ", bed_file)
    message("Converting PLINK to bigSNP format (one-time operation)...")
    bigsnpr::snp_readBed(bed_file,
                         backingfile = sub("\\.bed$", "", bed_file))
    message("Done. Backing files written to: ",
            sub("\\.bed$", "", bed_file), ".{rds,bk}")
    return(invisible(rds_file))
  }

  # ------------------------------------------------------------------
  # Dosage path: .traw -> FBM.code256 (CODE_DOSAGE)
  # ------------------------------------------------------------------
  if (is.null(dosage_file) || !file.exists(dosage_file))
    stop("format = 'dosage' requires an existing dosage_file (plink2 .traw).")
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required for the dosage path.")
  code <- bigsnpr::CODE_DOSAGE
  # dosage value -> byte code: nearest representable value; NA -> an NA byte.
  keepv   <- which(!is.na(code)); cv <- code[keepv]; cb <- keepv - 1L
  o       <- order(cv); cv <- cv[o]; cb <- cb[o]
  cuts    <- (cv[-1] + cv[-length(cv)]) / 2
  na_byte <- which(is.na(code))[1] - 1L
  encode  <- function(d) {
    b <- cb[findInterval(d, cuts) + 1L]
    b[is.na(d)] <- na_byte
    as.integer(b)
  }

  # Compressed dosage files are rejected outright: the chunked fread(skip=)
  # ingest below addresses byte offsets in the DEcompressed stream while any
  # line count on the raw file would see compressed bytes -- silently
  # inconsistent. Decompress first.
  if (grepl("\\.(gz|bz2|xz|zst)$", dosage_file, ignore.case = TRUE))
    stop("dosage_file appears to be compressed (", dosage_file, "). ",
         "Decompress it first (e.g. `gunzip`) -- the chunked ingest reads ",
         "the file by line offsets and does not support compressed input.")

  # header: CHR SNP (C)M POS COUNTED ALT <sample1> ...
  hdr      <- strsplit(readLines(dosage_file, n = 1L), "\t")[[1]]
  sample_ids <- hdr[-(1:6)]
  n_ind    <- length(sample_ids)
  if (n_ind < 1L) stop("No sample columns found in dosage_file header.")
  # Variant count = data lines. Counted in R (streamed 64 MB chunks,
  # newline bytes) rather than `wc -l`: portable, and correct for a file
  # without a trailing newline (wc would silently drop the last variant).
  n_lines <- local({
    con <- file(dosage_file, "rb")
    on.exit(close(con))
    n <- 0L
    last <- as.raw(10L)
    repeat {
      chunk <- readBin(con, "raw", n = 64L * 1024L * 1024L)
      if (length(chunk) == 0L) break
      n <- n + sum(chunk == as.raw(10L))
      last <- chunk[length(chunk)]
    }
    if (last != as.raw(10L)) n <- n + 1L   # no trailing newline: count it
    n
  })
  n_var <- as.integer(n_lines) - 1L
  if (n_var < 1L)
    stop("dosage_file has no data lines after the header: ", dosage_file)

  message(sprintf("Ingesting dosages: %d variants x %d individuals (chunk %d)...",
                  n_var, n_ind, chunk_size))

  G <- bigstatsr::FBM.code256(n_ind, n_var, code = code,
                              backingfile = plink_prefix)
  map_chr <- character(n_var); map_id <- character(n_var)
  map_pos <- integer(n_var);   map_a1 <- character(n_var); map_a2 <- character(n_var)

  done <- 0L
  while (done < n_var) {
    take <- min(chunk_size, n_var - done)
    dt <- data.table::fread(dosage_file, skip = 1L + done, nrows = take,
                            header = FALSE, sep = "\t", showProgress = FALSE)
    idx <- (done + 1L):(done + take)
    map_chr[idx] <- as.character(dt[[1]]); map_id[idx] <- as.character(dt[[2]])
    map_pos[idx] <- as.integer(dt[[4]])
    map_a1[idx]  <- as.character(dt[[5]]); map_a2[idx] <- as.character(dt[[6]])
    # dosages: variants x samples (cols 7..) -> encode -> samples x variants
    dose <- as.matrix(dt[, 7:(6 + n_ind)])
    G[, idx] <- matrix(encode(as.vector(t(dose))), n_ind, take)
    done <- done + take
    message(sprintf("  ...%d/%d variants", done, n_var))
  }

  fam <- data.frame(family.ID = sample_ids, sample.ID = sample_ids,
                    paternal.ID = 0L, maternal.ID = 0L, sex = 0L,
                    affection = -9L, stringsAsFactors = FALSE)
  map <- data.frame(chromosome = map_chr, marker.ID = map_id,
                    genetic.dist = 0, physical.pos = map_pos,
                    allele1 = map_a1, allele2 = map_a2,
                    stringsAsFactors = FALSE)
  obj <- structure(list(genotypes = G, fam = fam, map = map),
                   class = "bigSNP")
  saveRDS(obj, rds_file)
  message("Done. Dosage backing files written to: ", plink_prefix, ".{rds,bk}")
  invisible(rds_file)
}
