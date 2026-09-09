# Install required packages ------------------------------------------------

cran <- c("dplyr", "tibble", "ggplot2", "ggrepel", "remotes")
bioc <- c("DESeq2", "limma", "fgsea", "msigdbr", "RNAAgeCalc")

missing_cran <- setdiff(cran, rownames(installed.packages()))
if (length(missing_cran)) install.packages(missing_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- setdiff(bioc, rownames(installed.packages()))
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

if (!requireNamespace("immunedeconv", quietly = TRUE)) {
  remotes::install_github("omnideconv/immunedeconv")
}

message("Package setup complete. Run R/99_sessionInfo.R to record versions.")
