# Transcriptomic age estimation with RNAAgeCalc ----------------------------
# Raw gene counts are supplied directly to RNAAgeCalc. Do not VST-transform
# the matrix before this step; RNAAgeCalc performs the transformation required
# by its pretrained model internally.

source("R/config.R")
source("R/helpers.R")
library(RNAAgeCalc)

run_clock <- function(count_path, meta_path, cohort_name) {
  counts <- read_expression_matrix(count_path)
  meta <- read_metadata(meta_path, SAMPLE_COL)
  aligned <- align_expression_metadata(counts, meta, SAMPLE_COL)
  counts <- aligned$expr
  meta <- aligned$meta

  if (any(counts < 0, na.rm = TRUE)) stop("RNAAgeCalc count input contains negative values.")

  pred <- RNAAgeCalc::predict_age(
    exprdata = counts,
    tissue = RNAAGE_TISSUE,
    exprtype = RNAAGE_EXPRTYPE,
    idtype = RNAAGE_IDTYPE,
    stype = RNAAGE_STYPE,
    signature = RNAAGE_SIGNATURE,
    genelength = NULL,
    chronage = NULL,
    maxp = NULL
  )

  out <- data.frame(
    sample = rownames(pred),
    RNAAge = pred$RNAAge,
    stringsAsFactors = FALSE
  )
  out <- merge(meta, out, by.x = SAMPLE_COL, by.y = "sample", sort = FALSE)
  out$delta <- out$RNAAge - out[[AGE_COL]]

  write.table(
    out,
    file.path(RESULTS_DIR, paste0(cohort_name, "_rnaage_delta.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  out
}

thair_clock <- run_clock(THAIR_COUNTS, THAIR_META, "GSE152641")
zhang_clock <- run_clock(ZHANG_COUNTS, ZHANG_META, "GSE196399")
