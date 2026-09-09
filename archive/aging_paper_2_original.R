

GSE152641_dataset <- read.delim("~/Desktop/projetos_em_andamento/aging_signature/healthy_controls/GSE152641_covid_vst.tsv")

thair_metadata_healthy <- read.delim("~/Desktop/projetos_em_andamento/aging_signature/healthy_controls/meta_thair_healthy_final.tsv")

thair_metadata <- read.delim("~/Desktop/projetos_em_andamento/aging_signature/healthy_controls/metadata_GSE152641_covid.tsv")

thair_healthy <- GSE152641_dataset[,c(1:24)]

thair_covid <- GSE152641_dataset[,c(25:86)]

thair_metadata_covid <- thair_metadata[c(25:86),]

model1 <- lm(thair_healthy ~ age + sex)

library(limma)

design <- model.matrix(
  ~ age + sex,
  data = thair_metadata_healthy
)

fit <- lmFit(thair_healthy, design)

fit <- eBayes(fit)

results_thair_healthy <- topTable(
  fit,
  coef = "age",
  number = Inf
)

head(results)

write.table(results_thair_healthy, "~/Desktop/results_thair_healthy.tsv", sep = "\t")

design <- model.matrix(
  ~ age + sex + Neutrophil + Monocyte + Lymphocytes,
  data = thair_metadata_healthy
)

fit <- lmFit(thair_healthy, design)

fit <- eBayes(fit)

results_thair_healthy_with_cell_decomposition <- topTable(
  fit,
  coef = "age",
  number = Inf
)

head(results)
dim(design)

write.table(results_thair_healthy_with_cell_decomposition, "~/Desktop/results_thair_healthy_with_cell_decomposition.tsv", sep = "\t")

candidate_genes <- rownames(
  results[results$P.Value < 0.01, ]
)

covid_expr <- thair_covid[candidate_genes, ]
covid_expr <- covid_expr[,-55]

design <- model.matrix(
  ~ age + sex + Neutrophil + Monocyte + Lymphocytes,
  data = thair_metadata_covid_final
)


fit <- lmFit(covid_expr, design)
fit <- eBayes(fit)

results_thair_covid_with_cell_composition <- topTable(
  fit,
  coef = "age",
  number = Inf
)

write.table(results_thair_covid_with_cell_composition, "~/Desktop/results_thair_covid_with_cell_decomposition.tsv", sep = "\t")
head(results)

# install if needed
remotes::install_github("omnideconv/immunedeconv")

library(immunedeconv)
library(dplyr)
library(tibble)

# choose method
deconvolution_method <- "quanTIseq"

# expression matrix:
# rows = genes
# columns = samples
# values = normalized expression, preferably gene symbols as rownames

zhang_covid <-  GSE196399_dataset[,c(22:77)]

immune_zhang <- deconvolute(
  gene_expression = as.matrix(zhang_covid),
  method = "quantiseq"
)

head(immune_thair)

library(dplyr)
library(tibble)

immune_zhang_t <- immune_zhang

if ("method_cell_type" %in% colnames(immune_zhang_t)) {
  immune_zhang_t <- immune_zhang_t %>%
    column_to_rownames("method_cell_type")
} else if ("cell_type" %in% colnames(immune_zhang_t)) {
  immune_thair_t <- immune_thair_t %>%
    column_to_rownames("cell_type")
} else {
  stop("Cannot find cell-type column. Check colnames(immune_zhang).")
}

immune_zhang_t <- as.data.frame(immune_zhang)

rownames(immune_zhang_t) <- NULL


library(tibble)

if ("method_cell_type" %in% colnames(immune_zhang_t)) {
  
  immune_zhang_t <- immune_zhang_t %>%
    column_to_rownames("method_cell_type")
  
} else if ("cell_type" %in% colnames(immune_zhang_t)) {
  
  immune_zhang_t <- immune_zhang_t %>%
    column_to_rownames("cell_type")
  
} else {
  
  stop("Cannot find cell-type column.")
  
}


immune_zhang_t <- immune_zhang_t %>%
  t() %>%
  as.data.frame()

immune_zhang_t$sample <- rownames(immune_zhang_t)

head(immune_zhang_t)

immune_zhang_t <- immune_zhang_t %>%
  mutate(
    Neutrophil = Neutrophil,
    Monocyte = Monocyte,
    Lymphocytes =
      `B cell` +
      `T cell CD4+ (non-regulatory)` + `T cell regulatory (Tregs)` +
      `T cell CD8+`
  ) %>%
  select(sample, Neutrophil, Monocyte, Lymphocytes)






library(immunedeconv)
library(dplyr)
library(tibble)

# expression matrix:
# rows = genes
# columns = samples

immune_thair <- deconvolute(
  gene_expression = as.matrix(thair_covid),
  method = "quantiseq"
)

immune_thair_t <- immune_thair %>%
  column_to_rownames("cell_type") %>%
  t() %>%
  as.data.frame()

immune_thair_t$sample <- rownames(immune_thair_t)

head(immune_thair_t)

colnames(immune_thair_t)

immune_thair_t <- immune_thair_t %>%
  mutate(
    Neutrophil = Neutrophil,
    Monocyte = Monocyte,
    Lymphocytes =
      `B cell` +
      `T cell CD4+ (non-regulatory)` + `T cell regulatory (Tregs)` +
      `T cell CD8+`
  ) %>%
  select(sample, Neutrophil, Monocyte, Lymphocytes)

metadata_zhang_pneumonia <- metadata_zhang[c(22:77),]

colnames(metadata_zhang_pneumonia)[1] <- "sample"

zhang_metadata_pneumonia_final <- metadata_zhang_pneumonia %>%
  left_join(immune_zhang_t, by = "sample")

colnames(immune_thair_t)[1] <- "sample_ID"

write.table(zhang_metadata_pneumonia_final, "~/Desktop/zhang_metadata_pneumonia_final.tsv", sep = "\t")

GSE196399_dataset <- read.delim("~/Desktop/projetos_em_andamento/aging_signature/datasets/GSE196399_pneumonia_vst.tsv") 

metadata_zhang <- read.delim("~/Desktop/projetos_em_andamento/aging_signature/datasets/metadata_zhang.tsv") 

design <- model.matrix(
  ~ age + sex + Neutrophil + Monocyte + Lymphocytes,
  data = zhang_metadata_pneumonia_final
)

fit <- lmFit(zhang_covid, design)
fit <- eBayes(fit)

zhang_results_pneumonia_with_cell_composition <- topTable(
  fit,
  coef = "age",
  number = Inf,
  sort.by = "none"
)

write.table(zhang_results_pneumonia_with_cell_composition, "~/Desktop/results_zhang_pneumonia_with_cell_decomposition.tsv", sep = "\t")

head(zhang_results_pneumonia_with_cell_composition, 20)

results_thair_covid_with_cell_composition$gene <- rownames(results_thair_covid_with_cell_composition)
zhang_results_pneumonia_with_cell_composition$gene <- rownames(zhang_results_pneumonia_with_cell_composition)

merged <- merge(
  results_thair_covid_with_cell_composition,
  zhang_results_pneumonia_with_cell_composition,
  by = "gene",
  suffixes = c("_covid", "_pneumonia")
)

preserved <- merged[
  sign(merged$logFC_covid) ==
    sign(merged$logFC_pneumonia),
]

preserved$combined_p_fisher <- pchisq(
  -2 * (log(preserved$P.Value_covid) + log(preserved$P.Value_pneumonia)),
  df = 4,
  lower.tail = FALSE
)

preserved$effect_consistency <- 1 - abs(
  preserved$logFC_covid - preserved$logFC_pneumonia
) / (
  abs(preserved$logFC_covid) + abs(preserved$logFC_pneumonia)
)

preserved$preservation_ratio <- pmin(
  abs(preserved$logFC_covid),
  abs(preserved$logFC_pneumonia)
) / pmax(
  abs(preserved$logFC_covid),
  abs(preserved$logFC_pneumonia)
)

preserved$signature_score <-
  -log10(preserved$combined_p_fisher) *
  preserved$effect_consistency *
  preserved$preservation_ratio

preserved_ranked <- preserved[order(-preserved$signature_score), ]

head(preserved_ranked, 30)

write.table(
  preserved_ranked2,
  "~/Desktop/preserved_inflammation_resistant_aging_signature2.tsv",
  sep = "\t"
)

preserved_ranked[, c(
  "gene",
  "logFC_covid",
  "P.Value_covid",
  "logFC_pneumonia",
  "P.Value_pneumonia",
  "combined_p_fisher",
  "effect_consistency",
  "preservation_ratio",
  "signature_score"
)] |> head(30)

preserved$mean_abs_beta <- rowMeans(
  cbind(abs(preserved$logFC_covid), abs(preserved$logFC_pneumonia))
)

preserved$min_p <- pmin(
  preserved$P.Value_covid,
  preserved$P.Value_pneumonia
)

preserved$max_p <- pmax(
  preserved$P.Value_covid,
  preserved$P.Value_pneumonia
)

preserved$balanced_score <-
  -log10(preserved$combined_p_fisher) *
  preserved$preservation_ratio *
  preserved$effect_consistency *
  preserved$mean_abs_beta

preserved_ranked2 <- preserved[order(-preserved$balanced_score), ]

preserved_ranked2[, c(
  "gene",
  "logFC_covid",
  "P.Value_covid",
  "logFC_pneumonia",
  "P.Value_pneumonia",
  "combined_p_fisher",
  "effect_consistency",
  "preservation_ratio",
  "mean_abs_beta",
  "balanced_score"
)] |> head(30)

plot(
  preserved$logFC_covid,
  preserved$logFC_pneumonia,
  pch = 16,
  xlab = "Age coefficient (severe COVID-19 infections)",
  ylab = "Age coefficient (severe non-COVID-19 infections)"
)

abline(0, 1, col = "red", lwd = 2)


library(ggplot2)

# Figure 4B: distribution of preservation ratios
# Requires preserved$preservation_ratio

fig4b <- ggplot(preserved, aes(x = preservation_ratio)) +
  geom_histogram(
    bins = 20,
    color = "black",
    fill = "grey75"
  ) +
  geom_vline(
    xintercept = 0.25,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = 0.50,
    linetype = "dotted",
    linewidth = 0.8
  ) +
  labs(
    x = "Preservation ratio",
    y = "Number of genes",
    title = "Gene-level preservation of aging effects after immune adjustment"
  ) +
  theme_classic(base_size = 14)

fig4b

fig4b_filtered <- ggplot(
  subset(preserved, P.Value_covid < 0.05 | P.Value_pneumonia < 0.05),
  aes(x = preservation_ratio)
) +
  geom_histogram(
    bins = 15,
    color = "black",
    fill = "grey75"
  ) +
  geom_vline(xintercept = 0.25, linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = 0.50, linetype = "dotted", linewidth = 0.8) +
  labs(
    x = "Preservation ratio",
    y = "Number of genes",
    title = "Preserved age-associated genes across severe infections"
  ) +
  theme_classic(base_size = 14)

fig4b_filtered

covid_rank < results_thair_covid_with_cell_composition$logFC
names(covid_rank) <- results_thair_covid_with_cell_composition$gene

covid_rank <- sort(covid_rank, decreasing = TRUE)

covid_rank

pneumonia_rank <- zhang_results_pneumonia_with_cell_composition$logFC
names(pneumonia_rank) <- zhang_results_pneumonia_with_cell_composition$gene

pneumonia_rank <- sort(pneumonia_rank, decreasing = TRUE)

library(msigdbr)
library(fgsea)
library(dplyr)

hallmarks <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

pathways <- hallmarks %>%
  split(x = .$gene_symbol, f = .$gs_name)

fgsea_covid <- fgsea(
  pathways = pathways,
  stats = covid_rank,
  minSize = 15,
  maxSize = 500
)

fgsea_pneumonia <- fgsea(
  pathways = pathways,
  stats = pneumonia_rank,
  minSize = 15,
  maxSize = 500
)

head(fgsea_covid)

head(covid_rank)

length(covid_rank)

sum(is.na(covid_rank))

sum(duplicated(names(covid_rank)))

covid_rank <- covid_results$logFC

names(covid_rank) <- covid_results$gene

covid_rank <- covid_rank[!is.na(covid_rank)]

covid_rank <- covid_rank[!duplicated(names(covid_rank))]

covid_rank <- sort(covid_rank, decreasing = TRUE)

covid_rank

results_thair_covid_with_cell_composition$gene <- rownames(results_thair_covid_with_cell_composition)
zhang_results_pneumonia_with_cell_composition$gene <- rownames(zhang_results_pneumonia_with_cell_composition)

covid_rank <- results_thair_covid_with_cell_composition$t
names(covid_rank) <- results_thair_covid_with_cell_composition$gene

covid_rank <- covid_rank[!is.na(covid_rank)]
covid_rank <- covid_rank[!duplicated(names(covid_rank))]
covid_rank <- sort(covid_rank, decreasing = TRUE)

length(covid_rank)
head(covid_rank)

design_covid <- model.matrix(
  ~ age + sex + Neutrophil + Monocyte + Lymphocytes,
  data = thair_metadata_covid_final
)

thair_covid <- thair_covid[,-55]
fit_covid_all <- lmFit(thair_covid, design_covid)
fit_covid_all <- eBayes(fit_covid_all)

covid_results_all <- topTable(
  fit_covid_all,
  coef = "age",
  number = Inf,
  sort.by = "none"
)

covid_results_all$gene <- rownames(covid_results_all)

covid_rank <- covid_results_all$t
names(covid_rank) <- covid_results_all$gene

covid_rank <- covid_rank[!is.na(covid_rank)]
covid_rank <- covid_rank[!duplicated(names(covid_rank))]
covid_rank <- sort(covid_rank, decreasing = TRUE)

fgsea_covid <- fgsea(
  pathways = pathways,
  stats = covid_rank,
  minSize = 10,
  maxSize = 500
)

head(fgsea_covid[order(fgsea_covid$pval), ])

covid_pathways <- fgsea_covid[, c("pathway", "NES", "padj")]
pneumonia_pathways <- fgsea_pneumonia[, c("pathway", "NES", "padj")]

merged_pathways <- merge(
  covid_pathways,
  pneumonia_pathways,
  by = "pathway",
  suffixes = c("_covid", "_pneumonia")
)

merged_pathways$direction_preserved <-
  sign(merged_pathways$NES_covid) ==
  sign(merged_pathways$NES_pneumonia)

merged_pathways$NES_delta <-
  abs(merged_pathways$NES_covid -
        merged_pathways$NES_pneumonia)

head(
  merged_pathways[
    order(
      -merged_pathways$direction_preserved,
      merged_pathways$NES_delta
    ),
  ]
)

strong_preserved_pathways <- merged_pathways[
  direction_preserved == TRUE &
    padj_covid < 0.05 &
    padj_pneumonia < 0.05
]

strong_preserved_pathways

library(ggplot2)
library(ggrepel)
library(dplyr)

# Clean pathway labels
merged_pathways <- merged_pathways %>%
  mutate(
    pathway_clean = gsub("HALLMARK_", "", pathway),
    pathway_clean = gsub("_", " ", pathway_clean),
    significant_both = padj_covid < 0.05 & padj_pneumonia < 0.05,
    direction = ifelse(direction_preserved, "Same direction", "Opposite direction")
  )

fig5a <- ggplot(
  merged_pathways,
  aes(x = NES_covid, y = NES_pneumonia)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(
    aes(shape = significant_both),
    size = 3,
    alpha = 0.8
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dotted",
    linewidth = 0.7
  ) +
  geom_text_repel(
    data = subset(merged_pathways, significant_both),
    aes(label = pathway_clean),
    size = 3.5,
    max.overlaps = Inf
  ) +
  labs(
    x = "NES for age association in severe COVID-19",
    y = "NES for age association in non-COVID-19 pneumonia",
    title = "Pathway-level preservation of aging associations after immune adjustment",
    shape = "FDR < 0.05 in both cohorts"
  ) +
  theme_classic(base_size = 14)

fig5a

strong_preserved_pathways_plot <- strong_preserved_pathways %>%
  mutate(
    pathway_clean = gsub("HALLMARK_", "", pathway),
    pathway_clean = gsub("_", " ", pathway_clean)
  ) %>%
  select(pathway_clean, NES_covid, NES_pneumonia) %>%
  tidyr::pivot_longer(
    cols = c(NES_covid, NES_pneumonia),
    names_to = "cohort",
    values_to = "NES"
  ) %>%
  mutate(
    cohort = recode(
      cohort,
      NES_covid = "Severe COVID-19",
      NES_pneumonia = "Severe non-COVID-19 pneumonia"
    )
  )

fig5b <- ggplot(
  strong_preserved_pathways_plot,
  aes(x = cohort, y = reorder(pathway_clean, NES), fill = NES)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = round(NES, 2)),
    size = 4
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "NES"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Preserved age-associated pathways after immune adjustment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.ticks = element_blank()
  )

fig5b
