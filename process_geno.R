## =============================================================================
## process_geno.R
## Utility functions sourced by run_susie_batch_fixed.R.
## Contains genotype extraction and LD matrix validation.
##
## BEFORE RUNNING — one value to set:
##   VCF_PREFIX below must match your actual VCF filenames.
##   Example: if your files are TCD_UU_covid19_chr1.vcf.gz
##            set VCF_PREFIX <- "TCD_UU_covid19_"
##   The chr value from gene_remaining.tsv is appended directly,
##   so check whether your files use "chr1" or "1" format.
## =============================================================================

# ---------------------------------------------------------------------------
# VCF filename prefix — MUST be set before running
# Format: paste0(VCF_PREFIX, chr, ".vcf.gz")
# Example results:
#   VCF_PREFIX = "TCD_UU_covid19_"  + chr = "chr1"  → TCD_UU_covid19_chr1.vcf.gz
#   VCF_PREFIX = "cohort_genotypes_" + chr = "1"     → cohort_genotypes_1.vcf.gz
# ---------------------------------------------------------------------------
VCF_PREFIX <- "TCD_UU_covid19_"   # produces: TCD_UU_covid19_chr1.vcf.gz

# ---------------------------------------------------------------------------
# extract_genotypes
# Reads a VCF for a genomic region and returns a dosage matrix.
# Returns: matrix with rows=variants, cols=samples (dosage 0/1/2)
#          NULL if extraction fails or region is empty
# ---------------------------------------------------------------------------
extract_genotypes <- function(chr, start, end, vcf_dir) {

  # VCF files are NAMED with chr prefix:  TCD_UU_covid19_chr1.vcf.gz
  # VCF internal CHROM column uses plain integers: 1, 2, ... (no chr prefix)
  # → "chr" goes into the FILENAME but NOT into the seqminer region query
  vcf_file <- file.path(vcf_dir, paste0(VCF_PREFIX, "chr", chr, ".vcf.gz"))
  region   <- paste0(chr, ":", start, "-", end)   # matches internal VCF chr naming

  # Validate file exists before calling seqminer — gives clearer error message
  if (!file.exists(vcf_file)) {
    message("VCF file not found: ", vcf_file,
            "\nCheck VCF_PREFIX in process_geno.R and vcf_dir in run_susie_batch_fixed.R")
    return(NULL)
  }

  # Validate tabix index exists — seqminer requires .tbi for range queries
  if (!file.exists(paste0(vcf_file, ".tbi"))) {
    message("Tabix index not found: ", vcf_file, ".tbi",
            "\nRun: tabix -p vcf ", vcf_file)
    return(NULL)
  }

  genotypes <- tryCatch(
    seqminer::readVCFToMatrixByRange(vcf_file, region),
    error = function(e) {
      message("seqminer error | region: ", region,
              " | file: ", basename(vcf_file),
              " | ", e$message)
      NULL
    }
  )

  if (is.null(genotypes) || length(genotypes) == 0 || is.null(genotypes[[1]]))
    return(NULL)

  as.matrix(genotypes[[1]])
}

# ---------------------------------------------------------------------------
# validate_LD_matrix
# Validates and if necessary fixes an LD matrix.
#
# Returns a named list:
#   $fixed_LD_matrix    — the (possibly corrected) matrix
#   $validation_results — named list of all checks
#   $all_valid          — TRUE if no fixes were applied
# ---------------------------------------------------------------------------
validate_LD_matrix <- function(ld_matrix) {
  library(Matrix)

  orig_names <- rownames(ld_matrix)

  is_square         <- nrow(ld_matrix) == ncol(ld_matrix)
  is_symmetric      <- isSymmetric(ld_matrix, tol = 1e-8)
  diag_vals         <- diag(ld_matrix)
  is_diagonal_valid <- all(abs(diag_vals - 1) < 0.01)
  is_within_bounds  <- all(ld_matrix >= -1 & ld_matrix <= 1, na.rm = TRUE)
  eigenvalues       <- eigen(ld_matrix, symmetric = TRUE, only.values = TRUE)$values
  is_psd            <- all(eigenvalues >= -1e-6)

  if (!is_diagonal_valid) {
    diag(ld_matrix) <- 1
    message("  Fixed: diagonal values set to 1")
  }
  if (!is_symmetric) {
    ld_matrix <- (ld_matrix + t(ld_matrix)) / 2
    message("  Fixed: symmetrised LD matrix")
  }

  rownames(ld_matrix) <- colnames(ld_matrix) <- orig_names

  list(
    fixed_LD_matrix    = ld_matrix,
    validation_results = list(
      is_square                 = is_square,
      is_symmetric              = is_symmetric,
      diagonal_close_to_1       = is_diagonal_valid,
      values_within_bounds      = is_within_bounds,
      is_positive_semi_definite = is_psd,
      min_eigenvalue            = min(eigenvalues)
    ),
    all_valid = is_square && is_symmetric && is_diagonal_valid &&
                is_within_bounds && is_psd
  )
}
