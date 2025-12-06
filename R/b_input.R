# b_input.R - Standardized Input Data Handling Utilities
#
# This file contains shared utilities for handling input data across all
# borgworld functions, ensuring consistent behavior for:
# - Converting matrices to data frames and vice versa
# - Removing non-numeric columns
# - Handling missing values
# - Validating data inputs

#' Prepare Data for Analysis
#'
#' @description
#' Standardized data preparation utility used by borgworld analysis functions.
#' Converts input to appropriate format, handles missing values, and selects
#' numeric columns.
#'
#' @param data A data frame, matrix, or tibble to prepare
#' @param na_action Character string specifying missing value handling:
#'   \itemize{
#'     \item "complete" (default): Remove rows with any NA (listwise deletion)
#'     \item "pairwise": Keep all rows (for functions that use pairwise complete obs)
#'     \item "none": Keep all data as-is
#'   }
#' @param numeric_only Logical. If TRUE (default), keep only numeric columns.
#'   If FALSE, keep all columns.
#' @param keep_cols Character vector of column names to always keep, even if
#'   not numeric. Useful for preserving ID columns.
#' @param min_rows Minimum number of rows required after cleaning. Default is 2.
#' @param min_cols Minimum number of numeric columns required. Default is 1.
#' @param output_format Character string specifying output format:
#'   \itemize{
#'     \item "data.frame" (default): Return a data.frame
#'     \item "matrix": Return a matrix
#'   }
#' @param verbose Logical. If TRUE, print messages about dropped columns/rows.
#'   Default is FALSE.
#'
#' @return Prepared data in the specified format
#'
#' @details
#' This function provides a consistent data preparation pipeline:
#' \enumerate{
#'   \item Convert matrix/tibble to data.frame if needed
#'   \item Identify and optionally remove non-numeric columns
#'   \item Handle missing values based on na_action
#'   \item Validate that sufficient data remains
#'   \item Convert to requested output format
#' }
#'
#' @examples
#' # Basic usage - clean data, keep only numeric columns
#' clean_data <- bprepare_data(iris)
#'
#' # Keep character/factor columns
#' bprepare_data(iris, numeric_only = FALSE)
#'
#' # Return as matrix for matrix operations
#' bprepare_data(mtcars, output_format = "matrix")
#'
#' # Pairwise deletion (keep rows with some NAs)
#' bprepare_data(data_with_nas, na_action = "pairwise")
#'
#' @export
bprepare_data <- function(data,
                          na_action = c("complete", "pairwise", "none"),
                          numeric_only = TRUE,
                          keep_cols = NULL,
                          min_rows = 2,
                          min_cols = 1,
                          output_format = c("data.frame", "matrix"),
                          verbose = FALSE) {

  # Match arguments

  na_action <- match.arg(na_action)
  output_format <- match.arg(output_format)

  # Input validation
  if (is.null(data)) {
    stop("data cannot be NULL")
  }

  # Store original dimensions for reporting

  original_rows <- nrow(data)
  original_cols <- ncol(data)

  # Convert matrix to data frame if needed
  if (is.matrix(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    if (verbose) {
      message("Converted matrix to data.frame")
    }
  }

  # Convert tibble to data frame if needed
  if (inherits(data, "tbl_df") || inherits(data, "tbl")) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    if (verbose) {
      message("Converted tibble to data.frame")
    }
  }

  # Ensure it's a data.frame at this point
  if (!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  # Store original column names
  original_colnames <- colnames(data)

  # Identify numeric and non-numeric columns
  numeric_cols <- sapply(data, is.numeric)
  non_numeric_cols <- !numeric_cols

  # Handle column selection
  if (numeric_only) {
    # Start with numeric columns
    cols_to_keep <- numeric_cols

    # Add any requested keep_cols back
    if (!is.null(keep_cols)) {
      for (col in keep_cols) {
        if (col %in% names(data)) {
          cols_to_keep[col] <- TRUE
        }
      }
    }

    # Report dropped columns
    dropped_cols <- names(data)[!cols_to_keep]
    if (length(dropped_cols) > 0 && verbose) {
      message("Dropped non-numeric columns: ", paste(dropped_cols, collapse = ", "))
    }

    # Subset the data
    data <- data[, cols_to_keep, drop = FALSE]

  } else {
    # Keep all columns
    dropped_cols <- character(0)
  }

  # Check we have columns left
  n_numeric <- sum(sapply(data, is.numeric))
  if (n_numeric < min_cols) {
    stop("Insufficient numeric columns. Found ", n_numeric,
         ", required at least ", min_cols)
  }

  # Handle missing values
  if (na_action == "complete") {
    # Listwise deletion - only consider numeric columns for NA detection
    numeric_subset <- data[, sapply(data, is.numeric), drop = FALSE]
    complete_rows <- complete.cases(numeric_subset)
    n_dropped <- sum(!complete_rows)

    if (n_dropped > 0) {
      data <- data[complete_rows, , drop = FALSE]
      if (verbose) {
        message("Removed ", n_dropped, " rows with missing values (",
                round(100 * n_dropped / original_rows, 1), "%)")
      }
    }
  }
  # For "pairwise" and "none", we don't drop any rows

  # Validate sufficient data remains
  if (nrow(data) < min_rows) {
    stop("Insufficient data after cleaning. ", nrow(data),
         " rows remaining, required at least ", min_rows)
  }

  # Convert to output format
  if (output_format == "matrix") {
    # Only include numeric columns in matrix output
    numeric_data <- data[, sapply(data, is.numeric), drop = FALSE]
    result <- as.matrix(numeric_data)
  } else {
    result <- data
  }

  # Add attributes with metadata about the cleaning process
  attr(result, "bprepare_info") <- list(
    original_rows = original_rows,
    original_cols = original_cols,
    final_rows = nrow(data),
    final_cols = ncol(data),
    rows_dropped = original_rows - nrow(data),
    cols_dropped = dropped_cols,
    na_action = na_action
  )

  return(result)
}


#' Ensure Data is Numeric Matrix
#'
#' @description
#' Convenience function that prepares data and returns a numeric matrix.
#' Wrapper around bprepare_data() with matrix output.
#'
#' @param data A data frame, matrix, or tibble
#' @param na_action How to handle missing values ("complete", "pairwise", "none")
#' @param verbose Print messages about cleaning
#'
#' @return A numeric matrix
#'
#' @examples
#' X <- bas_matrix(iris)
#'
#' @export
bas_matrix <- function(data, na_action = "complete", verbose = FALSE) {
  bprepare_data(data,
                na_action = na_action,
                numeric_only = TRUE,
                output_format = "matrix",
                verbose = verbose)
}


#' Ensure Data is Clean Data Frame
#'
#' @description
#' Convenience function that prepares data and returns a clean data frame.
#' This is the recommended replacement for internal uses of bclean().
#'
#' @param data A data frame, matrix, or tibble
#' @param na_action How to handle missing values ("complete", "pairwise", "none")
#' @param numeric_only Keep only numeric columns (default TRUE)
#' @param verbose Print messages about cleaning
#'
#' @return A clean data.frame
#'
#' @examples
#' clean_df <- bas_df(iris)
#'
#' @export
bas_df <- function(data, na_action = "complete", numeric_only = TRUE, verbose = FALSE) {
  bprepare_data(data,
                na_action = na_action,
                numeric_only = numeric_only,
                output_format = "data.frame",
                verbose = verbose)
}


#' Check if Input is a Square Matrix
#'
#' @description
#' Validates that input is a square matrix (same number of rows and columns).
#' Used by functions that expect proximity matrices, correlation matrices, etc.
#'
#' @param x A matrix or data frame to check
#' @param name Optional name to use in error messages (default: "Input")
#'
#' @return TRUE if square, otherwise throws an error
#'
#' @examples
#' cor_mat <- cor(mtcars)
#' bcheck_square(cor_mat)
#'
#' @export
bcheck_square <- function(x, name = "Input") {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (!is.matrix(x)) {
    stop(name, " must be a matrix or data frame")
  }

  if (nrow(x) != ncol(x)) {
    stop(name, " must be a square matrix. Got ",
         nrow(x), " rows and ", ncol(x), " columns")
  }

  invisible(TRUE)
}


#' Check if Input Appears to be a Correlation Matrix
#'
#' @description
#' Heuristically checks if a matrix appears to be a correlation matrix
#' (symmetric, diagonal of 1s, values in [-1, 1]).
#'
#' @param x A matrix to check
#' @param tol Tolerance for numeric comparisons (default: 1e-10)
#'
#' @return TRUE if appears to be a correlation matrix, FALSE otherwise
#'
#' @examples
#' bcheck_cormat(cor(mtcars))  # TRUE
#' bcheck_cormat(cov(mtcars))  # FALSE
#'
#' @export
bcheck_cormat <- function(x, tol = 1e-10) {
  # Must be a matrix
  if (!is.matrix(x) && !is.data.frame(x)) {
    return(FALSE)
  }

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  # Must be square
  if (nrow(x) != ncol(x)) {
    return(FALSE)
  }

  # Must be symmetric
  if (!isSymmetric(x, tol = tol)) {
    return(FALSE)
  }

  # Diagonal should be 1s (or very close)
  if (any(abs(diag(x) - 1) > tol)) {
    return(FALSE)
  }

  # Row names should match column names (if present)
  if (!is.null(rownames(x)) && !is.null(colnames(x))) {
    if (!identical(rownames(x), colnames(x))) {
      return(FALSE)
    }
  }

  # Values should be in [-1, 1] range (with tolerance)
  if (any(x < -1 - tol, na.rm = TRUE) || any(x > 1 + tol, na.rm = TRUE)) {
    return(FALSE)
  }

  return(TRUE)
}


#' Extract Labels from Data
#'
#' @description
#' Extracts labels/names from data for use in analysis output.
#' Tries row names, then column names, then generates default labels.
#'
#' @param data A data frame or matrix
#' @param dim Which dimension to extract labels from:
#'   \itemize{
#'     \item "rows" or "r": Extract row names
#'     \item "cols" or "c": Extract column names
#'   }
#' @param default_prefix Prefix for generated labels if none exist
#'
#' @return Character vector of labels
#'
#' @examples
#' bget_labels(mtcars, "rows")
#' bget_labels(mtcars, "cols")
#'
#' @export
bget_labels <- function(data, dim = "rows", default_prefix = "Item") {
  dim <- tolower(substr(dim, 1, 1))

  if (dim == "r") {
    labels <- rownames(data)
    n <- nrow(data)
  } else if (dim == "c") {
    labels <- colnames(data)
    n <- ncol(data)
  } else {
    stop("dim must be 'rows' or 'cols'")
  }

  # If no labels, generate defaults
  if (is.null(labels) || all(labels == "")) {
    labels <- paste0(default_prefix, 1:n)
  }

  return(labels)
}


#' Parse Formula Variables
#'
#' @description
#' Extracts variable names from a formula object or string.
#'
#' @param formula A formula object or character string
#'
#' @return A list with components:
#'   \itemize{
#'     \item response: The response (dependent) variable name
#'     \item predictors: Character vector of predictor variable names
#'     \item all: All variable names
#'   }
#'
#' @examples
#' bparse_formula(mpg ~ wt + hp)
#' bparse_formula("mpg ~ wt * hp")
#'
#' @export
bparse_formula <- function(formula) {
  # Convert string to formula if needed
  if (is.character(formula)) {
    formula <- as.formula(formula)
  }

  if (!inherits(formula, "formula")) {
    stop("Input must be a formula or character string")
  }

  # Get all variable names
  all_vars <- all.vars(formula)

  # Response is the first variable (LHS)
  response <- all_vars[1]

  # Predictors are the rest
  predictors <- all_vars[-1]

  return(list(
    response = response,
    predictors = predictors,
    all = all_vars,
    formula = formula
  ))
}


#' Validate Formula Variables Exist in Data
#'
#' @description
#' Checks that all variables in a formula exist in the data.
#'
#' @param formula A formula object or character string
#' @param data A data frame
#'
#' @return TRUE if all variables exist, otherwise throws an error
#'
#' @examples
#' bvalidate_formula(mpg ~ wt + hp, mtcars)
#'
#' @export
bvalidate_formula <- function(formula, data) {
  parsed <- bparse_formula(formula)

  missing_vars <- parsed$all[!parsed$all %in% names(data)]

  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  invisible(TRUE)
}


#' Detect Binary Variables
#'
#' @description
#' Identifies which variables in a data frame are binary (have exactly 2 unique values).
#'
#' @param data A data frame
#' @param exclude_na If TRUE (default), NA values are not counted as a unique value
#'
#' @return Logical vector indicating which columns are binary
#'
#' @examples
#' # Create sample data
#' df <- data.frame(
#'   binary = c(0, 1, 0, 1),
#'   continuous = rnorm(4),
#'   categorical = c("a", "b", "c", "d")
#' )
#' bdetect_binary(df)
#'
#' @export
bdetect_binary <- function(data, exclude_na = TRUE) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  sapply(data, function(x) {
    if (exclude_na) {
      x <- x[!is.na(x)]
    }
    length(unique(x)) == 2
  })
}


#' Standardize Variable (Z-score)
#'
#' @description
#' Standardizes a numeric vector to have mean 0 and SD 1.
#'
#' @param x Numeric vector
#' @param population If TRUE, use population SD (n divisor). If FALSE (default),
#'   use sample SD (n-1 divisor).
#' @param na.rm Remove NA values before computing mean/SD (default TRUE)
#'
#' @return Standardized numeric vector
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' bzscore(x)
#'
#' @export
bzscore <- function(x, population = FALSE, na.rm = TRUE) {
  m <- mean(x, na.rm = na.rm)

  if (population) {
    # Population SD (divide by n)
    n <- sum(!is.na(x))
    s <- sqrt(sum((x - m)^2, na.rm = na.rm) / n)
  } else {
    # Sample SD (divide by n-1)
    s <- sd(x, na.rm = na.rm)
  }

  if (s == 0) {
    warning("Variable has zero variance, returning zeros")
    return(rep(0, length(x)))
  }

  (x - m) / s
}
