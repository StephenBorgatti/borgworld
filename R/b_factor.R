#' Stata-style Factor Analysis
#'
#' Performs factor analysis with Stata-style output formatting. Displays complete
#' eigenvalue table and formatted factor loadings with optional suppression of
#' small values in rotated solutions.
#'
#' @param data A data frame or matrix of observed variables
#' @param nfactors Maximum number of factors to retain for display. If NULL,
#'   determined by mineigen criterion. Default is NULL.
#' @param mineigen Minimum eigenvalue threshold for retaining factors. Factors
#'   with eigenvalues below this threshold are not displayed. Default is 1
#'   (Kaiser criterion).
#' @param rotate Rotation method. Default is "varimax". Options include "none",
#'   "varimax", "promax", "oblimin", "quartimin", etc.
#' @param fm Factor extraction method. Default is "pa" (principal axis). Other
#'   options: "ml" (maximum likelihood), "minres" (minimum residual), "wls"
#'   (weighted least squares), "gls" (generalized least squares).
#' @param cut Threshold for suppressing small loadings in rotated solution.
#'   Loadings with absolute value below this are not displayed. Default is 0.3.
#' @param n.obs Number of observations. If NULL, determined from nrow(data).
#' @param ... Additional arguments passed to psych::fa()
#'
#' @return Invisibly returns the psych::fa object for further analysis
#'
#' @examples
#' \dontrun{
#' # Basic usage with default settings (Kaiser criterion)
#' result <- bfactor(mydata)
#'
#' # Extract maximum 3 factors
#' result <- bfactor(mydata, nfactors = 3)
#'
#' # Use different eigenvalue threshold
#' result <- bfactor(mydata, mineigen = 0.7)
#'
#' # With promax rotation and lower cutoff
#' result <- bfactor(mydata, nfactors = 5, rotate = "promax", cut = 0.25)
#'
#' # Maximum likelihood estimation
#' result <- bfactor(mydata, nfactors = 3, fm = "ml")
#' }
#'
#' @export
bfactor <- function(data, nfactors = NULL, mineigen = 1, rotate = "varimax",
                    fm = "pa", cut = 0.5, n.obs = NULL, ...) {
  # Load required package
  if (!require("psych", quietly = TRUE)) {
    stop("Package 'psych' is required. Please install it.")
  }

  # Compute correlation matrix and eigenvalues once
  cor_mat <- cor(data, use = "pairwise.complete.obs")
  eigen_result <- eigen(cor_mat)
  all_eigenvalues <- eigen_result$values

  # Determine maximum factors to extract from fa()
  # Use a conservative threshold to avoid numerical issues
  max_factors <- sum(all_eigenvalues > 0.001)

  # If user specified nfactors, don't extract more than needed
  if (!is.null(nfactors)) {
    max_factors <- min(max_factors, nfactors)
  }

  # Also respect mineigen
  n_eigen <- sum(all_eigenvalues >= mineigen)
  max_factors <- min(max_factors, n_eigen)

  # Ensure at least 1 factor
  max_factors <- max(1, max_factors)

  # Extract factors using fa()
  fa_result <- psych::fa(data, nfactors = max_factors, rotate = rotate,
                         fm = fm, n.obs = n.obs, ...)

  # Determine how many factors to retain for display
  # This is separate from how many we extracted
  n_eigen <- sum(all_eigenvalues >= mineigen)

  if (is.null(nfactors)) {
    retained_factors <- n_eigen
  } else {
    retained_factors <- min(nfactors, n_eigen)
  }

  # Can't display more than we extracted
  retained_factors <- min(retained_factors, max_factors)

  # Ensure at least 1 factor
  if (retained_factors < 1) {
    retained_factors <- 1
    warning(sprintf("No eigenvalues >= %.2f found. Retaining 1 factor.", mineigen))
  }

  # Get number of observations
  if (is.null(n.obs)) {
    n.obs <- nrow(data)
  }

  # Get number of parameters for retained factors
  n_vars <- ncol(data)
  n_params <- retained_factors * n_vars - retained_factors * (retained_factors - 1) / 2

  # Print header
  cat("\nFactor analysis/correlation\n")
  cat(sprintf("    Method: %s\n",
              switch(fm,
                     "pa" = "principal factors",
                     "ml" = "maximum likelihood",
                     "minres" = "minimum residual",
                     "wls" = "weighted least squares",
                     "gls" = "generalized least squares",
                     fm)))

  rotation_name <- if (rotate == "none") "(unrotated)" else paste0("(", rotate, ")")
  cat(sprintf("    Rotation: %s\n", rotation_name))

  cat(sprintf("%50s = %10d\n", "Number of obs", n.obs))
  cat(sprintf("%50s = %10d\n", "Retained factors", retained_factors))
  cat(sprintf("%50s = %10d\n", "Number of params", n_params))

  # Print eigenvalue table for ALL eigenvalues from correlation matrix
  cat("\n")
  cat(strrep("-", 79))
  cat("\n")
  cat(sprintf("    %-10s %10s %12s %15s %12s\n",
              "Factor", "Eigenvalue", "Difference", "Proportion", "Cumulative"))
  cat(strrep("-", 79))
  cat("\n")

  # Use all eigenvalues from correlation matrix
  differences <- c(diff(all_eigenvalues), NA)
  total_var <- sum(all_eigenvalues)
  proportions <- all_eigenvalues / total_var
  cumulative <- cumsum(proportions)

  # Display all eigenvalues, but only up to ncol(data)
  n_display <- length(all_eigenvalues)
  for (i in 1:n_display) {
    diff_str <- if (is.na(differences[i])) "." else sprintf("%10.5f", differences[i])
    cat(sprintf("    Factor%-2d   %10.5f %12s %15.4f %12.4f\n",
                i, all_eigenvalues[i], diff_str, proportions[i], cumulative[i]))
  }
  cat(strrep("-", 79))
  cat("\n")

  # Print LR test if available (for ML method)
  if (fm == "ml" && !is.null(fa_result$STATISTIC)) {
    cat(sprintf("LR test: independent vs. saturated:  chi2(%d) = %7.2f Prob>chi2 = %.4f\n",
                fa_result$dof, fa_result$STATISTIC, fa_result$PVAL))
  }

  # Print unrotated factor loadings (retained factors only)
  cat("\nFactor loadings (pattern matrix) and unique variances\n\n")
  cat(strrep("-", 70))
  cat("\n")

  loading_header <- sprintf("    %-12s", "Variable")
  for (i in 1:retained_factors) {
    loading_header <- paste0(loading_header, sprintf(" %9s", paste0("Factor", i)))
  }
  loading_header <- paste0(loading_header, sprintf(" %12s", "Uniqueness"))
  cat(loading_header, "\n")
  cat(strrep("-", 70))
  cat("\n")

  loadings_mat <- fa_result$loadings[1:n_vars, 1:retained_factors, drop = FALSE]
  uniqueness <- fa_result$uniquenesses

  for (i in 1:n_vars) {
    row_str <- sprintf("    %-12s", rownames(loadings_mat)[i])
    for (j in 1:retained_factors) {
      row_str <- paste0(row_str, sprintf(" %9.4f", loadings_mat[i, j]))
    }
    row_str <- paste0(row_str, sprintf(" %12.4f", uniqueness[i]))
    cat(row_str, "\n")
  }
  cat(strrep("-", 70))
  cat("\n")

  # Print rotated factor loadings (sorted and with cutoff) - retained factors only
  # No uniqueness column for rotated loadings
  if (rotate != "none") {
    cat("\nRotated factor loadings (pattern matrix) and unique variances\n\n")
    cat(strrep("-", 70))
    cat("\n")

    # Header without uniqueness column
    loading_header_rot <- sprintf("    %-12s", "Variable")
    for (i in 1:retained_factors) {
      loading_header_rot <- paste0(loading_header_rot, sprintf(" %9s", paste0("Factor", i)))
    }
    cat(loading_header_rot, "\n")
    cat(strrep("-", 70))
    cat("\n")

    # Sort variables by their maximum absolute loading
    loadings_rot <- fa_result$loadings[1:n_vars, 1:retained_factors, drop = FALSE]
    max_loadings <- apply(abs(loadings_rot), 1, max)
    which_factor <- apply(abs(loadings_rot), 1, which.max)

    # Sort by factor, then by loading within factor
    sort_order <- order(which_factor, -max_loadings)
    loadings_sorted <- loadings_rot[sort_order, , drop = FALSE]

    for (i in 1:n_vars) {
      row_str <- sprintf("    %-12s", rownames(loadings_sorted)[i])
      for (j in 1:retained_factors) {
        val <- loadings_sorted[i, j]
        if (abs(val) < cut) {
          row_str <- paste0(row_str, sprintf(" %9s", ""))
        } else {
          row_str <- paste0(row_str, sprintf(" %9.4f", val))
        }
      }
      cat(row_str, "\n")
    }
    cat(strrep("-", 70))
    cat("\n")
  }

  # Return the fa object invisibly
  invisible(fa_result)
}
