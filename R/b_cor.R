#' Correlation with pairwise complete obs and auto numeric selection
#'
#' @description
#' Computes correlation matrix with automatic numeric column selection and
#' flexible missing value handling. Accepts both data frames and matrices.
#'
#' @param x A vector, matrix, data.frame, or tibble.
#' @param y Optional vector/matrix/data.frame/tibble for cross-correlation.
#' @param dim Character string specifying whether to compute correlations between
#'   "rows" or "columns" (default). Accepts any string starting with "r" or "c"
#'   (case-insensitive), so "r", "row", "rows" all work for rows, and "c", "col",
#'   "cols", "columns" all work for columns. This parameter is ignored when y is provided.
#' @param use Missing-data handling; default "pairwise.complete.obs".
#' @param verbose Logical; if TRUE, print messages about data cleaning.
#' @param ... Passed to [stats::cor()] (e.g., method = "spearman").
#' @return Correlation vector/matrix as in [stats::cor()].
#' @details
#' When dim = "rows", the function transposes the input to compute correlations
#' between rows instead of columns. This is equivalent to cor(t(x)) but preserves
#' row names in the output. The dim parameter only affects x and is ignored when
#' y is provided (cross-correlation always uses the data as-is).
#'
#' Non-numeric columns are automatically removed from data frames with a warning.
#'
#' @examples
#' # Create sample data
#' mat <- matrix(rnorm(20), nrow = 5, ncol = 4)
#' colnames(mat) <- paste0("var", 1:4)
#' rownames(mat) <- paste0("obs", 1:5)
#'
#' # Compute correlations between columns (default)
#' bcor(mat)
#' bcor(mat, dim = "columns")  # Same as above
#'
#' # Compute correlations between rows
#' bcor(mat, dim = "rows")
#'
#' # Flexible input formats
#' bcor(mat, dim = "r")     # Rows (shorthand)
#' bcor(mat, dim = "col")   # Columns (shorthand)
#'
#' # With data frame (auto-selects numeric columns)
#' df <- data.frame(
#'   name = c("A", "B", "C", "D", "E"),
#'   x = rnorm(5),
#'   y = rnorm(5),
#'   z = rnorm(5)
#' )
#' bcor(df)  # Correlations between x, y, z columns
#' bcor(df, dim = "rows")  # Correlations between observations
#'
#' @export
bcor <- function(x, y = NULL, dim = "columns", use = "pairwise.complete.obs",
                 verbose = FALSE, ...) {
  # Parse dim argument (only used when y is NULL)
  if (is.null(y)) {
    dim <- tolower(dim)
    if (substr(dim, 1, 1) == "r") {
      compute_rows <- TRUE
    } else if (substr(dim, 1, 1) == "c") {
      compute_rows <- FALSE
    } else {
      stop("dim must be 'rows' or 'columns' (or any string starting with 'r' or 'c')")
    }
  } else {
    # dim is ignored when y is provided (cross-correlation)
    compute_rows <- FALSE
  }

  # Helper function to select numeric columns using standardized input handling
  select_numeric <- function(obj, verbose = FALSE) {
    if (is.data.frame(obj) || is.matrix(obj)) {
      # Use bprepare_data for consistent handling
      obj <- bprepare_data(obj,
                           na_action = "none",  # Let cor() handle NAs
                           numeric_only = TRUE,
                           output_format = "data.frame",
                           verbose = verbose)
      if (ncol(obj) == 0L) stop("No numeric columns in data.", call. = FALSE)
    } else if (is.vector(obj) && !is.numeric(obj)) {
      stop("Vector input must be numeric.", call. = FALSE)
    }
    obj
  }

  # Select numeric columns
  x <- select_numeric(x, verbose)
  if (!is.null(y)) {
    y <- select_numeric(y, verbose)
  }

  # Transpose x if computing correlations between rows (only when y is NULL)
  if (is.null(y) && compute_rows) {
    # Store original row names
    original_names <- rownames(x)

    # Transpose for row-wise correlation
    x <- t(x)

    # The correlation will now be between what were originally rows
    result <- stats::cor(x = x, y = y, use = use, ...)

    # Set names to original row names if available
    if (!is.null(original_names)) {
      rownames(result) <- original_names
      colnames(result) <- original_names
    }

    return(result)
  }

  # Standard correlation (between columns)
  result <- stats::cor(x = x, y = y, use = use, ...)
  return(result)
}
