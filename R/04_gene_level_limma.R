# Composition-adjusted gene-level age associations -------------------------
# expression ~ age + sex + Neutrophil + Monocyte + Lymphocytes

source("R/config.R")
source("R/helpers.R")
library(limma)

run_gene_model <- function(vst_path, cohort_name) {
  expr <- read_expression_matrix(vst_path)
  clock <- read.delim(file.path(RESULTS_DIR, paste0(cohort_name, "_rnaage_delta.tsv")),
                      check.names = FALSE, stringsAsFactors = FALSE)
  immune <- read.delim(file.path(RESULTS_DIR, paste0(cohort_name, "_quantiseq_fractions.tsv")),
                       check.names = FALSE, stringsAsFactors = FALSE)
  meta <- merge(clock, immune, by.x = SAMPLE_COL, by.y = "sample", sort = FALSE)

  aligned <- align_expression_metadata(expr, meta, SAMPLE_COL)
  expr <- aligned$expr
  meta <- aligned$meta

  keep <- rowMeans(expr, na.rm = TRUE) > MIN_MEAN_VST
  expr <- expr[keep, , drop = FALSE]

  meta[[SEX_COL]] <- factor(meta[[SEX_COL]])
  design <- model.matrix(
    reformulate(c(AGE_COL, SEX_COL, "Neutrophil", "Monocyte", "Lymphocytes")),
    data = meta
  )

  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit)
  res <- limma::topTable(fit, coef = AGE_COL, number = Inf, sort.by = "none")
  res$gene <- rownames(res)

  write.table(res, file.path(RESULTS_DIR, paste0(cohort_name, "_limma_age_adjusted.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  res
}

covid_gene <- run_gene_model(THAIR_VST, "GSE152641")
pneumonia_gene <- run_gene_model(ZHANG_VST, "GSE196399")

merged <- merge(
  covid_gene[, c("gene", "logFC", "P.Value", "adj.P.Val", "t")],
  pneumonia_gene[, c("gene", "logFC", "P.Value", "adj.P.Val", "t")],
  by = "gene", suffixes = c("_covid", "_pneumonia")
)
merged$direction_concordant <- sign(merged$logFC_covid) == sign(merged$logFC_pneumonia)

write.table(merged, file.path(RESULTS_DIR, "gene_level_cross_cohort.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(subset(merged, direction_concordant),
            file.path(RESULTS_DIR, "gene_level_directionally_concordant.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
