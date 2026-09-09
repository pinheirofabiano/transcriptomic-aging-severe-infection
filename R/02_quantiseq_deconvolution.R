# Leukocyte deconvolution with quanTIseq via immunedeconv -------------------
# Provenance: archived project scripts used VST expression as input. This is
# recorded explicitly in R/config.R. If the final analysis is rerun with TPM or
# another supported non-log normalized expression scale, update the config and
# rerun all downstream analyses.

source("R/config.R")
source("R/helpers.R")
library(immunedeconv)
library(dplyr)
library(tibble)

message("quanTIseq input scale configured as: ", DECONV_INPUT_SCALE)

normalize_quantiseq_output <- function(x) {
  x <- as.data.frame(x, check.names = FALSE)
  cell_col <- intersect(c("cell_type", "method_cell_type"), colnames(x))
  if (length(cell_col) != 1) stop("Cannot uniquely identify the quanTIseq cell-type column.")

  y <- x %>%
    column_to_rownames(cell_col) %>%
    t() %>%
    as.data.frame(check.names = FALSE)
  y$sample <- rownames(y)

  needed <- c(
    "Neutrophil", "Monocyte", "B cell",
    "T cell CD4+ (non-regulatory)", "T cell regulatory (Tregs)", "T cell CD8+"
  )
  missing <- setdiff(needed, colnames(y))
  if (length(missing)) stop("Expected quanTIseq categories not found: ", paste(missing, collapse = ", "))

  y %>%
    mutate(
      Neutrophil = .data$Neutrophil,
      Monocyte = .data$Monocyte,
      Lymphocytes = .data$`B cell` +
        .data$`T cell CD4+ (non-regulatory)` +
        .data$`T cell regulatory (Tregs)` +
        .data$`T cell CD8+`
    ) %>%
    select(sample, Neutrophil, Monocyte, Lymphocytes)
}

run_quantiseq <- function(expr_path, cohort_name) {
  expr <- read_expression_matrix(expr_path)
  immune <- immunedeconv::deconvolute(
    gene_expression = expr,
    method = "quantiseq"
  )
  out <- normalize_quantiseq_output(immune)

  # No post-hoc renormalization of Neutrophil, Monocyte, or the summed
  # Lymphocytes variable is performed here. Values are taken from quanTIseq;
  # lymphocytes are the arithmetic sum of the listed B/T-cell fractions.
  write.table(
    out,
    file.path(RESULTS_DIR, paste0(cohort_name, "_quantiseq_fractions.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  out
}

thair_immune <- run_quantiseq(THAIR_DECONV_INPUT, "GSE152641")
zhang_immune <- run_quantiseq(ZHANG_DECONV_INPUT, "GSE196399")
