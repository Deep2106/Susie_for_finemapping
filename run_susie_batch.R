library(dplyr)
library(data.table)
library(susieR)
library(seqminer)
library(caret)
library(MASS)
library(Matrix)
library(filelock)  # For safe file logging

source("process_geno.R")  # Ensure this file is available

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
start_index <- as.integer(args[1])
end_index <- as.integer(args[2])
if (is.na(start_index) | is.na(end_index)) stop("Start or end index not provided.")

sample_size <- 433

# File paths
vcf_dir <- "/home/shared/Deepak/EQTL/SNPDATA_COVID19_LOCAL_UU/GENOTYPE_DATA/CHROMOSOMES"
cov_file <- "/home/shared/Deepak/EQTL/SNPDATA_COVID19_LOCAL_UU/COVARIATES/AUTOSOMES_TMM_LOG2/CovariatesTCDUU.txt"
z_dir <- "/home/shared/Deepak/EQTL/SNPDATA_COVID19_LOCAL_UU/EQTL/TMM_log2/Zscore"
ld_dir <- "/home/shared/Deepak/EQTL/SNPDATA_COVID19_LOCAL_UU/EQTL/TMM_log2/LD"
pip_dir <- "/home/shared/Deepak/EQTL/SNPDATA_COVID19_LOCAL_UU/EQTL/TMM_log2/pip_results"
credible_set_dir <- "/home/shared/Deepak/EQTL/SNPDATA_COVID19_LOCAL_UU/EQTL/TMM_log2/Zscorecredible_sets"
error_log <- "logs/error_log_single.txt"
summary_file <- "logs/gene_runtime_summary_single.txt"

# Create directories if they don't exist
dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create(ld_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pip_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(credible_set_dir, showWarnings = FALSE, recursive = TRUE)

# Define safe logging functions using filelock to avoid simultaneous writes
safe_write_log <- function(msg, log_file) {
  lock_path <- paste0(log_file, ".lock")
  lock_obj <- lock(lock_path, timeout = 2000)
  write(paste(Sys.time(), msg), file = log_file, append = TRUE)
  unlock(lock_obj)
}

# Load gene info
gene_info <- read.table("/home/shared/Deepak/EQTL/SNPDATA_COVID19_LOCAL_UU/GENE_EXPRESSION/gene_coordinates_ordered_17k_autosomes.txt", 
                        sep = "\t", header = TRUE, check.names = FALSE)
colnames(gene_info) <- c("gene", "chr", "tss", "end")
gene_info <- gene_info %>% dplyr::select(chr, tss, gene)

# Load covariates
covariates <- as.matrix(read.table(cov_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
covariates <- t(covariates)
I <- diag(nrow(covariates))
M <- I - covariates %*% MASS::ginv(t(covariates) %*% covariates) %*% t(covariates)

# ========= LOOP ONE GENE AT A TIME ========= #
for (i in start_index:end_index) {
  chr <- gene_info$chr[i]
  tss <- gene_info$tss[i]
  gene <- gene_info$gene[i]
  gene_start <- Sys.time()
  
  tryCatch({
    z_file <- file.path(z_dir, paste0(gene, ".fminput.z"))
    if (!file.exists(z_file)) stop("Missing Z-score file")
    
    start <- max(tss - 1000000, 0)
    end <- tss + 1000000
    genotype_mat <- extract_genotypes(chr, start, end, vcf_dir)
    if (is.null(genotype_mat) || ncol(genotype_mat) == 0) stop("Genotype extraction failed")
    
    genotype_mat <- as.data.frame(t(genotype_mat), check.names = FALSE)
    z_data <- fread(z_file)
    z_data$ids <- gsub("_", ":", z_data$ids)

    # Check if all z_data ids start with "chr"
   if (all(grepl("^chr", z_data$ids))) {
      # If z_data ids have "chr", ensure genotype_mat colnames start with "chr"
      # Identify column names that do not start with "chr"
      idx <- !grepl("^chr", colnames(genotype_mat))
      if (any(idx)) {
         colnames(genotype_mat)[idx] <- paste0("chr", colnames(genotype_mat)[idx])
       }
     } else {
          # If z_data ids do not have "chr", remove "chr" from genotype_mat colnames
          colnames(genotype_mat) <- gsub("^chr", "", colnames(genotype_mat))
   }

    common_snps <- intersect(colnames(genotype_mat), z_data$ids)
    if (length(common_snps) == 0) stop("No common SNPs")
    
    genotype_mat <- as.matrix(genotype_mat[, common_snps, drop = FALSE])
    z_data <- z_data[z_data$ids %in% common_snps, ]
    
    # Apply covariate adjustment
    X_adj <- M %*% genotype_mat
    LD_mat <- cor(X_adj, use = "pairwise.complete.obs")
    LD_mat <- (LD_mat + t(LD_mat)) / 2
    
    # Manual eigenvalue fix
    eig <- eigen(LD_mat)
    eig$values[eig$values < 0] <- 1e-6
    LD_mat <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    colnames(LD_mat) <- rownames(LD_mat) <- common_snps
    
    # Final safety check with nearPD()
    if (any(eigen(LD_mat)$values <= 0)) {
      LD_mat <- as.matrix(nearPD(LD_mat, corr = TRUE)$mat)
      rownames(LD_mat) <- colnames(LD_mat) <- common_snps
      if (any(eigen(LD_mat)$values <= 0)) stop("LD matrix not PSD")
    }
    
    saveRDS(LD_mat, file = file.path(ld_dir, paste0(gene, ".ld.RDS")))
    
    susie_res <- susie_rss(z_data$zscore[match(common_snps, z_data$ids)], R = LD_mat, refine = FALSE, n = sample_size)
    
    z_data$PIP <- susie_res$pip[match(z_data$ids, names(susie_res$pip))]
    fwrite(z_data, file = file.path(pip_dir, paste0(gene, "_pip.txt")), sep = "\t")
    saveRDS(susie_res,file =file.path(credible_set_dir, paste0(gene, "_PIP.RDS")))    
    cs <- susie_res$sets$cs
    if (!is.null(cs)) {
      fwrite(as.data.frame(cs), file = file.path(credible_set_dir, paste0(gene, "_credible_set.txt")), sep = "\t")
    }
    
    safe_write_log(
  paste(
    gene, "done in", 
    round(difftime(Sys.time(), gene_start, units = "secs"), 2), 
    "sec"
  ), 
  summary_file
)
    
  }, error = function(e) {
    safe_write_log(paste(gene, "ERROR:", e$message), error_log)
  })
}

