#' Compute Euclidean distance with missing value imputation
#'
#' @description
#' Computes Euclidean distance (or other distance metrics) between rows or columns
#' of a matrix or data frame. Unlike \code{dist()}, this function handles missing
#' values by imputing them with the mean, automatically ignores non-numeric columns,
#' and defaults to computing distances between columns (similar to \code{cor()}).
#'
#' @param x A matrix or data frame. Non-numeric columns will be automatically removed.
#' @param dim Character string specifying whether to compute distances between
#'   "rows" or "columns" (default). Accepts any string starting with "r" or "c"
#'   (case-insensitive), so "r", "row", "rows" all work for rows, and "c", "col",
#'   "cols", "columns" all work for columns.
#' @param ... Additional arguments passed to \code{\link{dist}}, such as:
#'   \itemize{
#'     \item \code{method}: The distance measure to use. Options include
#'       "euclidean" (default), "maximum", "manhattan", "canberra", "binary",
#'       or "minkowski".
#'     \item \code{diag}: Logical indicating whether the diagonal of the
#'       distance matrix should be printed.
#'     \item \code{upper}: Logical indicating whether the upper triangle of the
#'       distance matrix should be printed.
#'     \item \code{p}: The power of the Minkowski distance.
#'   }
#'
#' @return An object of class \code{"dist"} containing the computed distances.
#'   This is a lower triangular matrix by default (use \code{as.matrix()} to
#'   convert to a full matrix).
#'
#' @details
#' The function automatically removes non-numeric columns with a warning message
#' indicating which columns were excluded.
#'
#' When missing values are present in the numeric data, the function:
#' \enumerate{
#'   \item Issues a warning indicating the number and percentage of missing values
#'   \item Imputes missing values using the mean of the corresponding row or column:
#'     \itemize{
#'       \item If computing distances between columns (\code{dim = "columns"}),
#'         missing values are imputed with column means
#'       \item If computing distances between rows (\code{dim = "rows"}),
#'         missing values are imputed with row means
#'     }
#'   \item If an entire row or column consists of missing values, the global
#'     mean is used for imputation
#' }
#'
#' @examples
#' # Create sample data
#' mat <- matrix(rnorm(20), nrow = 5, ncol = 4)
#' colnames(mat) <- paste0("var", 1:4)
#' rownames(mat) <- paste0("obs", 1:5)
#'
#' # Compute distances between columns (default)
#' beuc(mat)
#' beuc(mat, "columns")  # Same as above
#'
#' # Compute distances between rows
#' beuc(mat, "rows")
#'
#' # Flexible input formats
#' beuc(mat, "r")     # Rows (shorthand)
#' beuc(mat, "col")   # Columns (shorthand)
#'
#' # With missing values
#' mat[2, 3] <- NA
#' mat[4, 1] <- NA
#' beuc(mat)  # Issues warning and imputes
#'
#' # With mixed data types (non-numeric columns are ignored)
#' df <- data.frame(
#'   name = c("A", "B", "C", "D", "E"),
#'   x = rnorm(5),
#'   y = rnorm(5),
#'   category = factor(c("X", "Y", "X", "Y", "X")),
#'   z = rnorm(5)
#' )
#' beuc(df)  # Warns about non-numeric columns, uses only x, y, z
#'
#' # Use different distance methods
#' beuc(mat, "columns", method = "manhattan")
#' beuc(mat, "rows", method = "maximum")
#'
#' # Get full distance matrix
#' as.matrix(beuc(mat))
#'
#' # Return upper triangle and diagonal
#' beuc(mat, "columns", diag = TRUE, upper = TRUE)
#'
#' @seealso
#' \code{\link{dist}} for the base R distance function,
#' \code{\link{cor}} for correlation computation with similar syntax,
#' \code{\link{bcor}} for correlation with similar data handling
#'
#' @export
beuc <- function(x, dim = "columns", ...) {
  # dim: "rows" or "columns" (default), accepts any string starting with "r" or "c"
  # ... passes additional arguments to dist()

  # Check inputs
  if (!is.matrix(x) && !is.data.frame(x)) {
    x <- as.matrix(x)
  }

  # Parse dim argument
  dim <- tolower(dim)
  if (substr(dim, 1, 1) == "r") {
    margin <- 1  # rows
  } else if (substr(dim, 1, 1) == "c") {
    margin <- 2  # columns
  } else {
    stop("dim must be 'rows' or 'columns' (or any string starting with 'r' or 'c')")
  }

  # Work with matrix
  x <- as.matrix(x)

  # Check for non-numeric columns and remove them
  numeric_cols <- apply(x, 2, function(col) {
    is.numeric(col) || all(is.na(col) | !is.na(as.numeric(as.character(col))))
  })

  if (!all(numeric_cols)) {
    removed_cols <- colnames(x)[!numeric_cols]
    if (is.null(removed_cols)) {
      removed_cols <- which(!numeric_cols)
    }
    warning(sprintf(
      "Removed %d non-numeric column(s): %s",
      sum(!numeric_cols),
      paste(removed_cols, collapse = ", ")
    ))
    x <- x[, numeric_cols, drop = FALSE]
  }

  # Convert to numeric matrix (in case some columns were character but convertible)
  x <- apply(x, 2, as.numeric)

  # Check if any columns remain
  if (ncol(x) == 0) {
    stop("No numeric columns found in the input data")
  }

  # Check for missing values
  if (anyNA(x)) {
    n_missing <- sum(is.na(x))
    total_cells <- length(x)
    pct_missing <- round(100 * n_missing / total_cells, 1)

    warning(sprintf(
      "Found %d missing values (%.1f%%). Imputing with mean.",
      n_missing, pct_missing
    ))

    # Impute missing values with mean
    # Calculate mean along the appropriate margin
    if (margin == 2) {
      # Computing distances between columns, so impute by column mean
      for (j in 1:ncol(x)) {
        if (anyNA(x[, j])) {
          x[is.na(x[, j]), j] <- mean(x[, j], na.rm = TRUE)
        }
      }
    } else {
      # Computing distances between rows, so impute by row mean
      for (i in 1:nrow(x)) {
        if (anyNA(x[i, ])) {
          x[i, is.na(x[i, ])] <- mean(x[i, ], na.rm = TRUE)
        }
      }
    }

    # Check if any NAs remain (can happen if entire row/column was NA)
    if (anyNA(x)) {
      # Use global mean for remaining NAs
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      warning("Some rows/columns were entirely NA. Used global mean for imputation.")
    }
  }

  # Transpose if computing distances between columns
  if (margin == 2) {
    x <- t(x)
  }

  # Compute distances (defaults to euclidean if method not specified)
  dist(x, ...)
}
