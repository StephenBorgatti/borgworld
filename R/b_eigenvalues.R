#' Extract Eigenvalues from Data or Correlation Matrix
#'
#' Efficiently extracts eigenvalues from a dataset or correlation matrix.
#' For large datasets (e.g., document-by-word matrices), uses SVD on the
#' column-standardized matrix to avoid computing the full correlation matrix.
#'
#' @param x A data frame, matrix, or correlation matrix.
#' @param use Character string specifying how to handle missing values.
#'   Options are "complete" (listwise deletion), "pairwise" (pairwise complete
#'   observations). Default is "complete".
#' @param efficient Logical. If TRUE (default), uses SVD on standardized data
#'   rather than computing the correlation matrix when input is raw data.
#'   This is much faster for wide matrices (many columns).
#' @param tol Tolerance for determining matrix symmetry. Default is 1e-10.
#'
#' @return An object of class "beigenvalues" containing:
#'   \item{eigenvalues}{A data frame with columns: Component, Eigenvalue,
#'     Ratio, Percent, Cumulative}
#'   \item{n}{Number of observations (if computed from raw data)}
#'   \item{p}{Number of variables}
#'   \item{input_type}{Whether input was "correlation" or "data"}
#'   \item{method}{Either "eigen" or "svd"}
#'
#' @examples
#' # From raw data
#' beigenvalues(mtcars)
#'
#' # From correlation matrix
#' beigenvalues(cor(mtcars))
#'
#' @export
beigenvalues <- function(x, use = "complete", efficient = TRUE, tol = 1e-10) {

  # Convert to matrix if data frame
  if (is.data.frame(x)) {
    # Select only numeric columns
    numeric_cols <- sapply(x, is.numeric)
    if (sum(numeric_cols) < ncol(x)) {
      dropped <- names(x)[!numeric_cols]
      message("Dropped non-numeric columns: ", paste(dropped, collapse = ", "))
    }
    x <- as.matrix(x[, numeric_cols, drop = FALSE])
  }

  x <- as.matrix(x)

  # Check if input appears to be a correlation matrix
  is_cormat <- .is_correlation_matrix(x, tol)

  if (is_cormat) {
    # Direct eigendecomposition of correlation matrix
    eig <- eigen(x, symmetric = TRUE, only.values = TRUE)
    eigenvalues <- eig$values
    n <- NA
    p <- nrow(x)
    input_type <- "correlation"
    method <- "eigen"
  } else {
    # Raw data: compute eigenvalues efficiently
    use <- match.arg(use, c("complete", "pairwise"))

    # Handle missing values
    if (use == "complete") {
      complete_rows <- complete.cases(x)
      if (sum(complete_rows) < nrow(x)) {
        n_dropped <- nrow(x) - sum(complete_rows)
        message("Removed ", n_dropped, " rows with missing values")
      }
      x <- x[complete_rows, , drop = FALSE]
    }

    n <- nrow(x)
    p <- ncol(x)

    if (n < 2) {
      stop("Need at least 2 complete observations")
    }

    if (efficient && use == "complete") {
      # SVD approach: faster for wide matrices
      # Standardize columns (z-scores)
      col_means <- colMeans(x, na.rm = TRUE)
      col_sds <- apply(x, 2, sd, na.rm = TRUE)

      # Handle zero-variance columns
      zero_var <- col_sds < tol
      if (any(zero_var)) {
        warning("Removed ", sum(zero_var), " zero-variance columns: ",
                paste(colnames(x)[zero_var], collapse = ", "))
        x <- x[, !zero_var, drop = FALSE]
        col_means <- col_means[!zero_var]
        col_sds <- col_sds[!zero_var]
        p <- ncol(x)
      }

      # Center and scale
      z <- scale(x, center = col_means, scale = col_sds)

      # SVD: eigenvalues of correlation matrix = singular values^2 / (n-1)
      # Use svd with nu=0, nv=0 for efficiency (only need singular values)
      sv <- svd(z, nu = 0, nv = 0)
      eigenvalues <- (sv$d^2) / (n - 1)

      method <- "svd"
    } else {
      # Traditional approach: compute correlation matrix, then eigen
      if (use == "pairwise") {
        R <- cor(x, use = "pairwise.complete.obs")
      } else {
        R <- cor(x)
      }

      # Check for NA in correlation matrix
      if (any(is.na(R))) {
        warning("Correlation matrix contains NA values; some eigenvalues may be unreliable")
        R[is.na(R)] <- 0
      }

      eig <- eigen(R, symmetric = TRUE, only.values = TRUE)
      eigenvalues <- eig$values
      method <- "eigen"
    }
    input_type <- "data"
  }

  # Build results data frame
  k <- length(eigenvalues)

  # Ratio of eigenvalue to next eigenvalue
  ratios <- c(eigenvalues[-k] / eigenvalues[-1], NA)

  # Percent variance
  total_var <- sum(eigenvalues)
  pct_var <- 100 * eigenvalues / total_var

  # Cumulative variance
  cum_var <- cumsum(pct_var)

  results <- data.frame(
    Component  = 1:k,
    Eigenvalue = eigenvalues,
    Ratio      = ratios,
    Percent    = pct_var,
    Cumulative = cum_var
  )

  # Create output object
  out <- list(
    eigenvalues = results,
    n = n,
    p = p,
    input_type = input_type,
    method = method
  )

  class(out) <- "beigenvalues"
  out
}


#' Check if matrix appears to be a correlation matrix
#' @noRd
.is_correlation_matrix <- function(x, tol = 1e-10) {
  # Must be square
  if (nrow(x) != ncol(x)) return(FALSE)

  # Must be symmetric
  if (!isSymmetric(x, tol = tol)) return(FALSE)

  # Diagonal should be 1s (or very close)
  if (any(abs(diag(x) - 1) > tol)) return(FALSE)

  # Row names should match column names (if present)
  if (!is.null(rownames(x)) && !is.null(colnames(x))) {
    if (!identical(rownames(x), colnames(x))) return(FALSE)
  }

  # Values should be in [-1, 1] range (with tolerance)
  if (any(x < -1 - tol) || any(x > 1 + tol)) return(FALSE)

  TRUE
}


#' Print method for beigenvalues
#' @export
print.beigenvalues <- function(x, digits = 4, n = NULL, ...) {

  cat("\nEigenvalue Decomposition\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")

  if (x$input_type == "data") {
    cat("Observations:", x$n, "   Variables:", x$p, "\n")
  } else {
    cat("Variables:", x$p, "   (input was correlation matrix)\n")
  }
  cat("Method:", ifelse(x$method == "svd", "SVD of standardized data", "Eigen decomposition"), "\n\n")

  # Determine how many rows to show
  df <- x$eigenvalues
  k <- nrow(df)

  if (is.null(n)) {
    # Show all if <= 20, otherwise show eigenvalues >= 1 plus a few more
    if (k <= 20) {
      n <- k
    } else {
      n_ge1 <- sum(df$Eigenvalue >= 1)
      n <- min(k, max(n_ge1 + 3, 10))
    }
  }
  n <- min(n, k)

  # Format the table
  show_df <- df[1:n, ]

  # Format columns
  fmt_eigenvalue <- formatC(show_df$Eigenvalue, digits = digits, format = "f")
  fmt_ratio <- ifelse(is.na(show_df$Ratio), "",
                      formatC(show_df$Ratio, digits = 2, format = "f"))
  fmt_percent <- formatC(show_df$Percent, digits = 2, format = "f")
  fmt_cumulative <- formatC(show_df$Cumulative, digits = 2, format = "f")

  # Column widths
  w_comp <- max(nchar("Comp"), nchar(as.character(show_df$Component)))
  w_eig <- max(nchar("Eigenvalue"), nchar(fmt_eigenvalue))
  w_ratio <- max(nchar("Ratio"), nchar(fmt_ratio))
  w_pct <- max(nchar("Percent"), nchar(fmt_percent))
  w_cum <- max(nchar("Cumulative"), nchar(fmt_cumulative))

  # Header
  header <- sprintf("%*s  %*s  %*s  %*s  %*s",
                    w_comp, "Comp",
                    w_eig, "Eigenvalue",
                    w_ratio, "Ratio",
                    w_pct, "Percent",
                    w_cum, "Cumulative")
  cat(header, "\n")
  cat(paste(rep("-", nchar(header)), collapse = ""), "\n")

  # Rows
  for (i in 1:n) {
    row <- sprintf("%*d  %*s  %*s  %*s  %*s",
                   w_comp, show_df$Component[i],
                   w_eig, fmt_eigenvalue[i],
                   w_ratio, fmt_ratio[i],
                   w_pct, fmt_percent[i],
                   w_cum, fmt_cumulative[i])
    cat(row, "\n")
  }

  if (n < k) {
    cat("... (", k - n, " more components)\n", sep = "")
  }

  cat("\n")

  invisible(x)
}
