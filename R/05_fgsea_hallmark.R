# Hallmark preranked GSEA ---------------------------------------------------
# Genes are ranked by the moderated t statistic for chronological age from
# the composition-adjusted limma model. Leading-edge genes are exported.

source("R/config.R")
library(fgsea)
library(msigdbr)
library(dplyr)

hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

run_fgsea <- function(cohort_name) {
  res <- read.delim(file.path(RESULTS_DIR, paste0(cohort_name, "_limma_age_adjusted.tsv")),
                    check.names = FALSE, stringsAsFactors = FALSE)
  ranks <- res$t
  names(ranks) <- res$gene
  ranks <- ranks[!is.na(ranks)]
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)

  fg <- fgsea::fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = FGSEA_MIN_SIZE,
    maxSize = FGSEA_MAX_SIZE
  )
  fg <- as.data.frame(fg)
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  fg <- fg[order(fg$padj), ]

  write.table(fg, file.path(RESULTS_DIR, paste0(cohort_name, "_hallmark_fgsea.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  fg
}

fg_covid <- run_fgsea("GSE152641")
fg_pneumonia <- run_fgsea("GSE196399")

merged <- merge(
  fg_covid[, c("pathway", "NES", "pval", "padj", "leadingEdge")],
  fg_pneumonia[, c("pathway", "NES", "pval", "padj", "leadingEdge")],
  by = "pathway", suffixes = c("_covid", "_pneumonia")
)
merged$direction_concordant <- sign(merged$NES_covid) == sign(merged$NES_pneumonia)
merged$significant_both <- merged$padj_covid < 0.05 & merged$padj_pneumonia < 0.05

write.table(merged, file.path(RESULTS_DIR, "hallmark_cross_cohort.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
