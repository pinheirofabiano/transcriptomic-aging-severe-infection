
zhang_dataset <- read.delim("~/Desktop/projetos_em_andamento/sepsis_clock/healthy_controls/GSE196399_pneumonia_vst.tsv")

thair_dataset <- read.delim("~/Desktop/projetos_em_andamento/sepsis_clock/healthy_controls/GSE152641_covid_vst.tsv")

zhang_metadata <- read.delim("~/Desktop/projetos_em_andamento/sepsis_clock/healthy_controls/metadata_zhang.tsv")

thair_metadata <- read.delim("~/Desktop/projetos_em_andamento/sepsis_clock/healthy_controls/metadata_GSE152641_covid.tsv")

thair_metadata_healthy <- thair_metadata[c(1:24),]

zhang_metadata_healthy <- zhang_metadata[c(1:21),]

thair_metadata_healthy$study <- rep('thair', 24)

zhang_metadata_healthy$study <- rep("zhang", 21)

thair_metadata_healthy <- thair_metadata_healthy[,-c(3,7)]

colnames(zhang_metadata_healthy)
colnames(thair_metadata_healthy)

colnames(thair_metadata_healthy) <- c("GEO_ID", "study_ID", "age", "sex", "class", "study")

thair_metadata_healthy <- thair_metadata_healthy[,-5]

zhang_metadata_healthy <- zhang_metadata_healthy[,-2]

healthy_metadata_merged <- rbind(thair_metadata_healthy, zhang_metadata_healthy)

colnames(thair_metadata_healthy)

zhang_metadata_healthy <- zhang_metadata_healthy[,c(2,1,3,4,5)]

write.table(healthy_metadata_merged, file = "~/Desktop/healthy_metadata_merged_final.tsv", sep = "\t")

healthy_metadata_merged_final <- read.delim("~/Desktop/healthy_metadata_merged_final.tsv")

thair_healthy_dataset <- thair_dataset[,c(1:24)]

zhang_healthy_dataset <- zhang_dataset[,c(1:21)]

zhang_healthy_dataset$symbols <- rownames(zhang_healthy_dataset)

thair_healthy_dataset$symbols <- rownames(thair_healthy_dataset)


healthy_dataset_final <- merge(thair_healthy_dataset, zhang_healthy_dataset, by= "symbols")

colnames(healthy_dataset_final) <- c("symbols", healthy_metadata_merged_final$patient_ID)

write.table(healthy_dataset_final, "~/Desktop/healthy_dataset_final.tsv", sep = "\t")

aggregate(age ~ healthy_dataset_final, healthy_metadata_merged_final, summary)

table(healthy_metadata_merged_final$study)

aggregate(age ~ healthy_dataset_final, healthy_metadata_merged_final, summary)

aggregate(age ~ healthy_metadata_merged_final$study,
          data = healthy_metadata_merged_final,
          FUN = summary)

thair_healthy_dataset <- thair_healthy_dataset[,-25]

library(limma)

design <- model.matrix(~ age + sex, data=thair_metadata_healthy)

fit <- lmFit(thair_healthy_dataset, design)
fit <- eBayes(fit)

results <- topTable(fit,
                    coef="age",
                    number=Inf)
head(results, 30)

keep <- rowMeans(thair_healthy_dataset) > 5
vst_filtered <- thair_healthy_dataset[keep, ]

design <- model.matrix(~ age + sex, data = thair_metadata_healthy)

fit <- lmFit(vst_filtered, design)
fit <- eBayes(fit)

results <- topTable(fit, coef = "age", number = Inf)

ranked_genes <- results$t
names(ranked_genes) <- rownames(results)

ranked_genes <- sort(ranked_genes,
                     decreasing = TRUE)

# Install if needed
# install.packages("BiocManager")
BiocManager::install(c("fgsea", "msigdbr"))

library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)

# 1. Get Hallmark gene sets
hallmark <- msigdbr(species = "Homo sapiens", category = "H")

pathways <- hallmark %>%
  split(x = .$gene_symbol, f = .$gs_name)

# 2. Make sure ranked_genes is clean
ranked_genes <- ranked_genes[!is.na(ranked_genes)]
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Optional: remove duplicated gene names
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]

# 3. Run fgsea
fgsea_res <- fgsea(
  pathways = pathways,
  stats = ranked_genes,
  minSize = 15,
  maxSize = 500,
  nperm = 10000
)

# 4. Order results by adjusted p-value
fgsea_res <- fgsea_res %>%
  arrange(padj)

# 5. View top pathways
head(fgsea_res[, c("pathway", "NES", "pval", "padj", "size")], 20)

top_pathways <- fgsea_res %>%
  filter(padj < 0.25) %>%
  arrange(padj) %>%
  head(20)

ggplot(top_pathways,
       aes(x = reorder(pathway, NES),
           y = NES)) +
  geom_col(fill = "orange") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(x = "",
       y = "Normalized enrichment score",
       title = "Hallmark GSEA: age-associated pathways")

fgsea_res_export <- fgsea_res

fgsea_res_export$leadingEdge <- sapply(
  fgsea_res_export$leadingEdge,
  paste,
  collapse = ";"
)

write.table(fgsea_res_export, "~/Desktop/hallmark_gsea_age_thair_controls.tsv", sep = "\t")


