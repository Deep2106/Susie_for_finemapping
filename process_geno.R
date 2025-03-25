
# Simplified genotype extraction function
extract_genotypes <- function(chr, start, end, vcf_dir) {
  # Construct file path and genomic region
  vcf_file <- file.path(vcf_dir, paste0("TCD_UU_covid19_", chr, ".vcf.gz"))
  region <- paste0(chr, ":", start, "-", end)

  # Extract genotype data - already in dosage format (0/1/2)
  genotypes <- seqminer::readVCFToMatrixByRange(vcf_file, region)

  # Return the genotype matrix directly (no transformation needed)
  return(as.matrix(genotypes[[1]]))
}


validate_LD_matrix <- function(ld_matrix, save_fixed = TRUE, output_file = "fixed_ld_matrix.txt") {
  # Load required library
  library(Matrix)

  # ✅ 1. Check if matrix is square
  is_square <- nrow(ld_matrix) == ncol(ld_matrix)

  # ✅ 2. Check if matrix is symmetric
  is_symmetric <- all(ld_matrix == t(ld_matrix))

  # ✅ 3. Check if diagonal values are close to 1
  diagonal_values <- diag(ld_matrix)
  is_diagonal_valid <- all(abs(diagonal_values - 1) < 0.01)  # Allow small numerical errors

  # ✅ 4. Check if values are within valid range (-1 to 1)
  is_within_bounds <- all(ld_matrix >= -1 & ld_matrix <= 1)

  # ✅ 5. Check if matrix is positive semi-definite (All eigenvalues ≥ 0)
  eigenvalues <- eigen(ld_matrix)$values
  is_positive_semi_definite <- all(eigenvalues >= -1e-6)

  # 🛠️ Fix diagonal if needed
  if (!is_diagonal_valid) {
    diag(ld_matrix) <- 1
    message("✅ Fixed diagonal values to 1.")
  }

  # 🛠️ Fix symmetry if needed
  if (!is_symmetric) {
    ld_matrix <- (ld_matrix + t(ld_matrix)) / 2
    message("✅ Fixed LD matrix to be symmetric.")
  }

  # ✅ Prepare Validation Results
  validation_results <- list(
    "Is Square" = is_square,
    "Is Symmetric" = is_symmetric,
    "Diagonal Close to 1" = is_diagonal_valid,
    "Values Between -1 and 1" = is_within_bounds,
    "Positive Semi-Definite" = is_positive_semi_definite,
    "Eigenvalues" = eigenvalues
  )

  # ✅ Print only the first 4 validation results
  print(validation_results[1:4])
# output_file = file.path(ld_dir, paste0(gene, "_fixed.ld"))
#if (save_fixed) {
#    tryCatch({
#      write.table(ld_matrix, output_file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
#      message(paste("✅ Fixed LD matrix saved successfully to:", output_file))
#    }, error = function(e) {
#      message("❌ Error saving LD matrix:", e$message)
#    })
#  }

  # Return the fixed matrix and validation results
 # return(list("fixed_LD_matrix" = ld_matrix, "validation_results" = validation_results))
}




