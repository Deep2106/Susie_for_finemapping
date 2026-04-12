## =============================================================================
## run_susie_batch_fixed.R
## STEP 3 of the SuSiE finemapping pipeline (called by submit_susie_batch_fixed.sh)
##
## SAMPLE SIZE NOTE:
##   VCF files contain 433 samples.
##   MatrixQTL eQTL analysis used 430 WGS-matched samples.
##   Z-scores were computed from those 430 samples.
##
##   The LD matrix passed to SuSiE RSS MUST be computed from the same 430
##   samples as the Z-scores — otherwise the Z/LD relationship is broken.
##
##   FIX: after extracting genotypes from VCF (433 samples), rows are
##   immediately subsetted and reordered to match the 430 sample IDs in
##   the covariate matrix (which defines the analysis sample set and order).
##   This is done BEFORE computing M %*% genotype_mat — otherwise the
##   matrix multiplication fails (430×430 × 433×p is non-conformable).
##
## Usage (called by SLURM — not directly):
##   Rscript run_susie_batch_fixed.R <start_index> <end_index>
## =============================================================================

library(dplyr)
library(data.table)
library(susieR)
library(seqminer)
library(MASS)
library(Matrix)
library(filelock)

source("process_geno.R")

# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
start_index <- as.integer(args[1])
end_index   <- as.integer(args[2])
if (is.na(start_index) | is.na(end_index))
  stop("Usage: Rscript run_susie_batch_fixed.R <start_index> <end_index>")

# ---------------------------------------------------------------------------
# Scientific parameters
# ---------------------------------------------------------------------------
SAMPLE_SIZE        <- 430    # WGS-matched samples used in MatrixQTL eQTL
CIS_WINDOW_BP      <- 1e6    # +/- 1 Mb cis window around TSS
SUSIE_L            <- 10     # max causal signals (GTEx standard)
SUSIE_COVERAGE     <- 0.95   # credible set coverage probability
SUSIE_MIN_ABS_CORR <- 0.5    # purity filter: discard CS where max|r| < 0.5
SUSIE_MAX_ITER     <- 1000   # ELBO convergence ceiling
# refine=TRUE: local-optima escape via ELBO refinement
# Required for reliable multi-signal finemapping — do not set FALSE

# ---------------------------------------------------------------------------
# Paths — replace all placeholders before running
# ---------------------------------------------------------------------------

vcf_dir          <- "/path/to/GENOTYPE_DATA/CHROMOSOMES"      # per-chr VCF.gz (433 samples)
cov_file         <- "/path/to/for_PEER60/Covariates_TCD_UU_Peer60_latest.tsv" # 69 rows × 430 cols
z_dir            <- "/path/to/CIS_TRANS_EQTLS_SEPT_2025/COVS_default/Zscore"        # .fminput.z files
ld_dir           <- "/path/to/CIS_TRANS_EQTLS_SEPT_2025/COVS_default/LD"
pip_dir          <- "/path/to/CIS_TRANS_EQTLS_SEPT_2025/COVS_default/pip_results"
credible_set_dir <- "/path/to/CIS_TRANS_EQTLS_SEPT_2025/COVS_default/credible_sets"

gene_info_file   <- "genes.tsv"                 # eGenes only (filter_egenes.R)
error_log        <- "logs/error_log_susie.txt"
summary_file     <- "logs/gene_runtime_summary_susie.txt"

for (d in c("logs", ld_dir, pip_dir, credible_set_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Safe logging
# ---------------------------------------------------------------------------
safe_write_log <- function(msg, log_file) {
  lock_obj <- lock(paste0(log_file, ".lock"), timeout = 2000)
  write(paste(Sys.time(), msg), file = log_file, append = TRUE)
  unlock(lock_obj)
}

# ---------------------------------------------------------------------------
# Load eGene coordinate list
# ---------------------------------------------------------------------------
gene_info           <- read.table(gene_info_file, sep = "\t", header = TRUE,
                                  check.names = FALSE, stringsAsFactors = FALSE)
colnames(gene_info) <- c("gene", "chr", "tss", "end")
gene_info           <- gene_info %>% dplyr::select(chr, tss, gene)

if (end_index > nrow(gene_info)) end_index <- nrow(gene_info)
cat(sprintf("Processing genes %d to %d of %d eGenes\n\n",
            start_index, end_index, nrow(gene_info)))

# ---------------------------------------------------------------------------
# Load covariates and derive the 430-sample ID set
#
# The covariate matrix defines the analysis sample set:
#   rows = covariates (69), cols = samples (430)
#   colnames = sample IDs used in MatrixQTL and Z-score computation
#
# analysis_samples: character vector of 430 sample IDs in analysis order.
# Used to subset and reorder VCF genotypes from 433 → 430.
# ---------------------------------------------------------------------------
cat("Loading covariates...\n")
covariates <- as.matrix(read.table(cov_file, header = TRUE, row.names = 1,
                                   sep = "\t", check.names = FALSE))
covariates      <- t(covariates)          # now: 430 samples × 69 covariates
analysis_samples <- rownames(covariates)  # 430 sample IDs — the ground truth

cat(sprintf("  Covariate matrix: %d samples × %d covariates\n",
            nrow(covariates), ncol(covariates)))
cat(sprintf("  Analysis samples: %d  (these define the 430-sample set)\n\n",
            length(analysis_samples)))

# ---------------------------------------------------------------------------
# Covariate projection matrix M — computed once, reused per gene
#
# M = I - C(C'C)^{-1}C'   (Frisch-Waugh-Lovell)
# Projecting genotypes through M gives covariate-adjusted genotypes.
# LD = cor(M × X) matches the covariate-adjusted space of the Z-scores.
# M is 430×430 — genotype matrix must be 430 rows before multiplication.
# ---------------------------------------------------------------------------
M <- diag(nrow(covariates)) -
     covariates %*% MASS::ginv(t(covariates) %*% covariates) %*% t(covariates)

# ---------------------------------------------------------------------------
# Helper: extract credible sets with SNP IDs and full metadata
#
# susie_res$sets$cs — named list ($L1, $L2, ...) of INTEGER VECTORS.
# Each integer is a POSITIONAL INDEX into common_snps — NOT a SNP ID.
# Resolved here to actual SNP IDs.
# ---------------------------------------------------------------------------
extract_credible_sets <- function(susie_res, snp_ids) {
  cs_list  <- susie_res$sets$cs
  coverage <- susie_res$sets$coverage
  purity   <- susie_res$sets$purity

  if (is.null(cs_list) || length(cs_list) == 0) return(NULL)

  rows <- lapply(names(cs_list), function(cs_name) {
    idx       <- cs_list[[cs_name]]
    snp_names <- snp_ids[idx]
    pip_vals  <- susie_res$pip[idx]

    cs_cov <- if (!is.null(coverage) && cs_name %in% names(coverage))
                coverage[[cs_name]] else NA_real_
    min_r  <- if (!is.null(purity) && cs_name %in% rownames(purity))
                purity[cs_name, "min.abs.corr"]    else NA_real_
    mean_r <- if (!is.null(purity) && cs_name %in% rownames(purity))
                purity[cs_name, "mean.abs.corr"]   else NA_real_
    med_r  <- if (!is.null(purity) && cs_name %in% rownames(purity))
                purity[cs_name, "median.abs.corr"] else NA_real_

    data.table(
      cs_id           = cs_name,
      snp             = snp_names,
      pip             = pip_vals,
      cs_size         = length(idx),
      cs_coverage     = cs_cov,
      min_abs_corr    = min_r,
      mean_abs_corr   = mean_r,
      median_abs_corr = med_r
    )
  })

  rbindlist(rows)
}

# ---------------------------------------------------------------------------
# Helper: CS membership column for PIP table
# ---------------------------------------------------------------------------
annotate_cs_membership <- function(snp_ids, susie_res) {
  cs_col  <- rep(NA_character_, length(snp_ids))
  cs_list <- susie_res$sets$cs
  if (is.null(cs_list)) return(cs_col)
  for (cs_name in names(cs_list)) {
    members              <- snp_ids[cs_list[[cs_name]]]
    cs_col[snp_ids %in% members] <- cs_name
  }
  cs_col
}

# ===========================================================================
# MAIN LOOP
# ===========================================================================
for (i in start_index:end_index) {
  chr             <- gene_info$chr[i]
  tss             <- gene_info$tss[i]
  gene            <- gene_info$gene[i]
  gene_start_time <- Sys.time()

  tryCatch({

    # -----------------------------------------------------------------------
    # 1. Load Z-scores
    #    Columns: snp (variant ID) | zscore
    #    Produced by calc_zscore_fixed.R — all cis-window SNPs present.
    # -----------------------------------------------------------------------
    z_file <- file.path(z_dir, paste0(gene, ".fminput.z"))
    if (!file.exists(z_file)) stop("Missing Z-score file: ", z_file)
    z_data <- fread(z_file)

    if (!all(c("snp", "zscore") %in% colnames(z_data)))
      stop("Z-score file must have columns 'snp' and 'zscore'. Found: ",
           paste(colnames(z_data), collapse = ", "))

    # SNP ID format: MatrixQTL produces IDs in underscore format (e.g. 1_10492_C_T_b38)
    # seqminer returns VCF ID column rownames in the same underscore format
    # DO NOT apply gsub("_",":") — both sides already match, conversion breaks intersection

    # -----------------------------------------------------------------------
    # 2. Extract cis-window genotypes from VCF (returns 433-sample matrix)
    # -----------------------------------------------------------------------
    win_start    <- max(tss - CIS_WINDOW_BP, 0)
    win_end      <- tss + CIS_WINDOW_BP
    genotype_mat <- extract_genotypes(chr, win_start, win_end, vcf_dir)

    if (is.null(genotype_mat) || ncol(genotype_mat) == 0)
      stop("Genotype extraction returned empty matrix")

    # genotype_mat: rows=variants, cols=samples (433 samples from VCF)
    # Transpose to: rows=samples (433), cols=variants
    genotype_mat <- as.data.frame(t(genotype_mat), check.names = FALSE)

    # -----------------------------------------------------------------------
    # 3. Subset genotype matrix from 433 → 430 samples
    #
    # The VCF contains 433 samples; the analysis used 430.
    # We must restrict to the exact 430 sample IDs used in MatrixQTL
    # and in the same order — this is enforced by analysis_samples which
    # is derived from the covariate matrix column names.
    #
    # Without this step:
    #   M (430×430) %*% genotype_mat (433×p) → non-conformable error
    #   LD computed from 433 samples ≠ LD space of Z-scores (430 samples)
    # -----------------------------------------------------------------------
    vcf_samples <- rownames(genotype_mat)

    # Identify which analysis samples are present in the VCF
    samples_in_vcf <- analysis_samples[analysis_samples %in% vcf_samples]

    if (length(samples_in_vcf) == 0)
      stop("No analysis sample IDs found in VCF colnames. ",
           "Check sample ID format between covariate matrix and VCF.")

    if (length(samples_in_vcf) < length(analysis_samples)) {
      n_missing <- length(analysis_samples) - length(samples_in_vcf)
      safe_write_log(
        paste0(gene, " | WARNING: ", n_missing,
               " analysis samples missing from VCF — proceeding with ",
               length(samples_in_vcf), " samples"),
        error_log
      )
    }

    # Subset rows to analysis samples in analysis order
    genotype_mat <- as.matrix(genotype_mat[samples_in_vcf, , drop = FALSE])
    # genotype_mat is now: 430 rows × n_snps cols — matches M dimensions

    # -----------------------------------------------------------------------
    # 4. Match seqminer positions to Z-score SNP IDs
    #
    # seqminer returns CHROM:POS rownames (e.g. "1:10492")
    # Z-score SNP IDs are CHROM_POS_REF_ALT_build format (e.g. "1_10492_C_T_b38")
    # These will never match directly — we extract CHROM:POS from the Z-score
    # IDs as a matching key, align on position, then relabel the genotype matrix
    # columns with the full Z-score SNP IDs so all output files carry allele info.
    # -----------------------------------------------------------------------

    # Extract CHROM:POS key from Z-score SNP IDs
    # "1_10492_C_T_b38" -> "1:10492"  (first two underscore-delimited fields)
    snp_parts      <- strsplit(z_data$snp, "_")
    z_data$snp_pos <- vapply(snp_parts,
                             function(x) paste(x[1], x[2], sep = ":"),
                             character(1))

    # seqminer colnames are CHROM:POS — intersect on position
    common_pos <- intersect(colnames(genotype_mat), z_data$snp_pos)

    if (length(common_pos) == 0)
      stop("No common SNPs between VCF and Z-score file. ",
           "seqminer format: ", paste(head(colnames(genotype_mat), 3), collapse=", "),
           " | Z-score snp_pos format: ", paste(head(z_data$snp_pos, 3), collapse=", "))

    # Subset both to common positions in the same order
    genotype_mat <- genotype_mat[, common_pos, drop = FALSE]
    z_data       <- z_data[match(common_pos, z_data$snp_pos), ]

    # Relabel genotype matrix columns with full Z-score SNP IDs
    # All downstream output (PIPs, credible sets) uses full allele-aware IDs
    colnames(genotype_mat) <- z_data$snp
    common_snps            <- z_data$snp

    # -----------------------------------------------------------------------
    # 6. LD matrix in covariate-adjusted space
    #
    # X_adj = M %*% genotype_mat
    #   M:            430 × 430
    #   genotype_mat: 430 × n_snps  ← correct after step 3
    #   X_adj:        430 × n_snps
    #
    # use="complete.obs": consistent sample set for all pairwise correlations
    # → PSD guaranteed from cor() itself
    # -----------------------------------------------------------------------
    X_adj  <- M %*% genotype_mat
    LD_mat <- cor(X_adj, use = "complete.obs")
    LD_mat <- (LD_mat + t(LD_mat)) / 2    # enforce exact symmetry

    eig            <- eigen(LD_mat, symmetric = TRUE)
    eig$values[eig$values < 0] <- 1e-6
    LD_mat         <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    dimnames(LD_mat) <- list(common_snps, common_snps)

    # nearPD() safety net
    if (any(eigen(LD_mat, only.values = TRUE)$values <= 0)) {
      LD_mat <- as.matrix(nearPD(LD_mat, corr = TRUE)$mat)
      dimnames(LD_mat) <- list(common_snps, common_snps)
      if (any(eigen(LD_mat, only.values = TRUE)$values <= 0))
        stop("LD matrix not PSD after nearPD")
    }

    saveRDS(LD_mat, file = file.path(ld_dir, paste0(gene, ".ld.RDS")))

    # -----------------------------------------------------------------------
    # 7. SuSiE RSS finemapping
    # -----------------------------------------------------------------------
    susie_res <- susie_rss(
      z            = z_data$zscore,
      R            = LD_mat,
      n            = SAMPLE_SIZE,        # 430 — the analysis sample size
      L            = SUSIE_L,
      coverage     = SUSIE_COVERAGE,
      min_abs_corr = SUSIE_MIN_ABS_CORR,
      max_iter     = SUSIE_MAX_ITER,
      refine       = TRUE
    )

    # -----------------------------------------------------------------------
    # 8. Convergence check
    # -----------------------------------------------------------------------
    if (!isTRUE(susie_res$converged)) {
      saveRDS(susie_res,
              file.path(credible_set_dir, paste0(gene, "_susie.RDS")))
      safe_write_log(
        paste0(gene, " | NOT_CONVERGED | ",
               round(difftime(Sys.time(), gene_start_time, units="secs"), 1), " sec"),
        summary_file
      )
      next
    }

    # -----------------------------------------------------------------------
    # 9. Save full SuSiE object (alpha matrix needed for coloc.susie)
    # -----------------------------------------------------------------------
    saveRDS(susie_res,
            file.path(credible_set_dir, paste0(gene, "_susie.RDS")))

    # -----------------------------------------------------------------------
    # 10. PIP output table
    # -----------------------------------------------------------------------
    z_data$PIP           <- susie_res$pip
    z_data$cs_membership <- annotate_cs_membership(common_snps, susie_res)

    # snp_pos is a temporary matching key — exclude from output
    pip_out <- z_data[, c("snp", "zscore", "PIP", "cs_membership")]

    fwrite(pip_out,
           file = file.path(pip_dir, paste0(gene, "_pip.txt")),
           sep  = "\t")

    # -----------------------------------------------------------------------
    # 11. Credible set table — one row per (CS, SNP), SNP IDs resolved
    # -----------------------------------------------------------------------
    cs_table <- extract_credible_sets(susie_res, common_snps)

    if (!is.null(cs_table) && nrow(cs_table) > 0) {
      cs_table[, gene := gene]
      setcolorder(cs_table, c("gene", "cs_id", "snp", "pip",
                               "cs_size", "cs_coverage",
                               "min_abs_corr", "mean_abs_corr", "median_abs_corr"))
      fwrite(cs_table,
             file = file.path(credible_set_dir, paste0(gene, "_credible_sets.txt")),
             sep  = "\t")
    } else {
      safe_write_log(
        paste0(gene, " | INFO | no CS passed purity filter (min_abs_corr >= ",
               SUSIE_MIN_ABS_CORR, ")"),
        summary_file
      )
    }

    # -----------------------------------------------------------------------
    # 12. Runtime log
    # -----------------------------------------------------------------------
    n_cs <- if (!is.null(susie_res$sets$cs)) length(susie_res$sets$cs) else 0
    safe_write_log(
      paste0(gene, " | done | ",
             round(difftime(Sys.time(), gene_start_time, units="secs"), 1), " sec | ",
             "n_SNPs=", length(common_snps), " | n_CS=", n_cs),
      summary_file
    )

  }, error = function(e) {
    safe_write_log(paste0(gene, " | ERROR | ", conditionMessage(e)), error_log)
  })
}

cat("\nBatch complete: genes", start_index, "to", end_index, "\n")
