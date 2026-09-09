# Project configuration ----------------------------------------------------

PROJECT_ROOT <- normalizePath(".", mustWork = FALSE)
RAW_DIR      <- file.path(PROJECT_ROOT, "data", "raw")
DERIVED_DIR  <- file.path(PROJECT_ROOT, "data", "derived")
RESULTS_DIR  <- file.path(PROJECT_ROOT, "results")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Raw-count inputs for RNAAgeCalc.
THAIR_COUNTS <- file.path(RAW_DIR, "GSE152641_raw_counts.tsv")
ZHANG_COUNTS <- file.path(RAW_DIR, "GSE196399_raw_counts.tsv")

# Curated metadata. Required columns are documented below.
THAIR_META <- file.path(DERIVED_DIR, "GSE152641_metadata.tsv")
ZHANG_META <- file.path(DERIVED_DIR, "GSE196399_metadata.tsv")

# Expression matrices used for gene-level limma models.
THAIR_VST <- file.path(DERIVED_DIR, "GSE152641_vst.tsv")
ZHANG_VST <- file.path(DERIVED_DIR, "GSE196399_vst.tsv")

# Expression matrices used for quanTIseq.
# The recovered archived script used the VST matrices as the deconvolution input.
# Keep this explicit for provenance. If a revised analysis uses TPM/non-log normalized
# expression, update these paths, set DECONV_INPUT_SCALE accordingly, and rerun all
# leukocyte-adjusted analyses.
THAIR_DECONV_INPUT <- THAIR_VST
ZHANG_DECONV_INPUT <- ZHANG_VST
DECONV_INPUT_SCALE <- "VST (recovered archived analysis)"

# Metadata column names.
SAMPLE_COL <- "sample"
CLASS_COL  <- "class"
AGE_COL    <- "age"
SEX_COL    <- "sex"

# Healthy reference label used in the regression models.
HEALTHY_LABEL <- "Healthy"

# RNAAgeCalc configuration.
RNAAGE_TISSUE    <- "blood"
RNAAGE_EXPRTYPE  <- "counts"
RNAAGE_IDTYPE    <- "SYMBOL"
RNAAGE_STYPE     <- "all"
RNAAGE_SIGNATURE <- NULL  # NULL uses the package default DESeq2 signature for blood.

# Gene-level filtering on VST expression.
MIN_MEAN_VST <- 5

# fgsea settings.
FGSEA_MIN_SIZE <- 10
FGSEA_MAX_SIZE <- 500
