# Primary disease-association models ---------------------------------------
# Base model:   delta ~ class + age + sex
# Immune model: delta ~ class + age + sex + Neutrophil + Monocyte + Lymphocytes

source("R/config.R")
source("R/helpers.R")

run_models <- function(cohort_name, disease_label) {
  clock_path <- file.path(RESULTS_DIR, paste0(cohort_name, "_rnaage_delta.tsv"))
  immune_path <- file.path(RESULTS_DIR, paste0(cohort_name, "_quantiseq_fractions.tsv"))

  d <- read.delim(clock_path, check.names = FALSE, stringsAsFactors = FALSE)
  immune <- read.delim(immune_path, check.names = FALSE, stringsAsFactors = FALSE)
  d <- merge(d, immune, by.x = SAMPLE_COL, by.y = "sample", all.x = TRUE, sort = FALSE)

  d[[CLASS_COL]] <- factor(d[[CLASS_COL]], levels = c(HEALTHY_LABEL, disease_label))
  d[[SEX_COL]] <- factor(d[[SEX_COL]])

  base_formula <- reformulate(c(CLASS_COL, AGE_COL, SEX_COL), response = "delta")
  immune_formula <- reformulate(
    c(CLASS_COL, AGE_COL, SEX_COL, "Neutrophil", "Monocyte", "Lymphocytes"),
    response = "delta"
  )

  base_fit <- lm(base_formula, data = d)
  immune_fit <- lm(immune_formula, data = d)

  disease_term <- paste0(CLASS_COL, disease_label)
  base_row <- extract_model_row(base_fit, disease_term, cohort_name, "age- and sex-adjusted")
  immune_row <- extract_model_row(
    immune_fit, disease_term, cohort_name,
    "age-, sex-, and quanTIseq-composition-adjusted"
  )

  write.table(
    rbind(base_row, immune_row),
    file.path(RESULTS_DIR, paste0(cohort_name, "_disease_model_summary.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  saveRDS(list(base = base_fit, immune = immune_fit),
          file.path(RESULTS_DIR, paste0(cohort_name, "_models.rds")))

  rbind(base_row, immune_row)
}

# Change labels here only if your curated metadata uses different disease labels.
thair_models <- run_models("GSE152641", "Severe COVID-19")
zhang_models <- run_models("GSE196399", "Severe pneumonia")

all_models <- rbind(thair_models, zhang_models)
write.table(all_models, file.path(RESULTS_DIR, "primary_model_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(all_models)
