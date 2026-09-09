read_expression_matrix <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)

  # If the first column is non-numeric, treat it as a gene identifier column.
  if (ncol(x) > 1 && !is.numeric(x[[1]])) {
    genes <- make.unique(as.character(x[[1]]))
    x <- x[, -1, drop = FALSE]
    rownames(x) <- genes
  }

  x <- as.data.frame(lapply(x, as.numeric), check.names = FALSE)
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x
}

read_metadata <- function(path, sample_col = "sample") {
  if (!file.exists(path)) stop("Missing file: ", path)
  m <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!sample_col %in% colnames(m)) {
    stop("Metadata must contain a '", sample_col, "' column: ", path)
  }
  m
}

align_expression_metadata <- function(expr, meta, sample_col = "sample") {
  common <- intersect(colnames(expr), meta[[sample_col]])
  if (length(common) < 2) stop("Fewer than two matched samples between expression and metadata.")
  expr <- expr[, common, drop = FALSE]
  meta <- meta[match(common, meta[[sample_col]]), , drop = FALSE]
  stopifnot(identical(colnames(expr), meta[[sample_col]]))
  list(expr = expr, meta = meta)
}

extract_model_row <- function(fit, term_pattern, cohort, model_name) {
  sm <- summary(fit)$coefficients
  idx <- grep(term_pattern, rownames(sm), fixed = TRUE)
  if (length(idx) != 1) {
    stop("Could not uniquely identify disease coefficient using pattern: ", term_pattern,
         ". Coefficient names: ", paste(rownames(sm), collapse = ", "))
  }
  data.frame(
    cohort = cohort,
    model = model_name,
    term = rownames(sm)[idx],
    beta_years = unname(sm[idx, "Estimate"]),
    SE = unname(sm[idx, "Std. Error"]),
    t = unname(sm[idx, "t value"]),
    p_value = unname(sm[idx, "Pr(>|t|)"]),
    stringsAsFactors = FALSE
  )
}
