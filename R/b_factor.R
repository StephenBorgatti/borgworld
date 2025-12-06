#' Calculate Factor Reliability Measures
#'
#' Computes reliability measures (alpha, omega, GLB) for factors from a factor
#' analysis solution. Can be called independently or is automatically called by
#' bfactor().
#'
#' @param fa_result A psych::fa object from factor analysis
#' @param data The original data frame used in the factor analysis
#' @param cut Threshold for determining which items load on each factor.
#'   Default is 0.4.
#' @param use_rotated Logical. Use rotated loadings if TRUE (default),
#'   unrotated if FALSE.
#'
#' @return A data frame with reliability measures for each factor
#'
#' @export
breliability <- function(fa_result, data, cut = 0.4, use_rotated = TRUE) {
  if (!requireNamespace("psych", quietly = TRUE)) {
    stop("Package 'psych' is required. Please install it.")
  }

  # Determine which loadings to use
  loadings_mat <- if (use_rotated && !is.null(fa_result$loadings)) {
    as.matrix(fa_result$loadings)
  } else {
    as.matrix(fa_result$loadings)
  }

  n_factors <- ncol(loadings_mat)

  # Initialize results data frame with correct structure
  reliability_df <- data.frame(
    factor = integer(n_factors),
    n_items = integer(n_factors),
    items = character(n_factors),
    alpha = numeric(n_factors),
    std_alpha = numeric(n_factors),
    omega_t = numeric(n_factors),
    stringsAsFactors = FALSE
  )

  # Fill in NAs
  reliability_df$alpha <- NA_real_
  reliability_df$std_alpha <- NA_real_
  reliability_df$omega_t <- NA_real_

  for (i in 1:n_factors) {
    loadings_vec <- loadings_mat[, i]
    # Use absolute values to identify items that load on this factor
    high_loaders <- which(abs(loadings_vec) >= cut)

    reliability_df$factor[i] <- i
    reliability_df$n_items[i] <- length(high_loaders)
    reliability_df$items[i] <- if (length(high_loaders) > 0) {
      paste(names(loadings_vec)[high_loaders], collapse = ", ")
    } else {
      ""
    }

    if (length(high_loaders) >= 2) {
      factor_data <- data[, high_loaders, drop = FALSE]

      # Identify which items need to be reversed (negative loadings)
      items_to_reverse <- loadings_vec[high_loaders] < 0

      # Create reversed data for omega
      factor_data_reversed <- factor_data
      if (any(items_to_reverse)) {
        for (j in which(items_to_reverse)) {
          # Reverse by subtracting from max + min
          factor_data_reversed[, j] <- max(factor_data[, j], na.rm = TRUE) +
            min(factor_data[, j], na.rm = TRUE) -
            factor_data[, j]
        }
      }

      # Calculate alpha with check.keys = TRUE (it handles reversal automatically)
      alpha_result <- tryCatch({
        suppressMessages(suppressWarnings(
          psych::alpha(factor_data, check.keys = TRUE)
        ))
      }, error = function(e) NULL)

      if (!is.null(alpha_result)) {
        reliability_df$alpha[i] <- alpha_result$total$raw_alpha
        reliability_df$std_alpha[i] <- alpha_result$total$std.alpha
      }

      # Calculate omega (if 3+ items) using reversed data
      if (length(high_loaders) >= 3) {
        omega_result <- tryCatch({
          suppressMessages(suppressWarnings({
            capture.output({
              omega_out <- psych::omega(factor_data_reversed, nfactors = 1, plot = FALSE)
            })
            omega_out
          }))
        }, error = function(e) NULL)

        if (!is.null(omega_result)) {
          reliability_df$omega_t[i] <- omega_result$omega.tot
        }
      }
    }
  }

  return(reliability_df)
}


#' Stata-style Factor Analysis
#'
#' Performs factor analysis with Stata-style output formatting. Displays complete
#' eigenvalue table and formatted factor loadings with optional suppression of
#' small values in rotated solutions. Includes reliability measures for each factor.
#'
#' @param data A data frame, matrix, or DocumentTermMatrix of observed variables.
#'   DocumentTermMatrix objects (from tm/tidytext packages) are automatically
#'   converted to matrix format. Non-numeric columns are automatically removed.
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
#'   Loadings with absolute value below this are not displayed. Default is 0.4.
#' @param n.obs Number of observations. If NULL, determined from nrow(data).
#' @param na_action How to handle missing values: "complete" (default) for listwise
#'   deletion, "pairwise" for pairwise complete observations.
#' @param verbose Logical; if TRUE, print messages about data cleaning.
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
#'
#' # Get reliability measures separately
#' fa_result <- bfactor(mydata, nfactors = 3)
#' reliability <- breliability(fa_result, mydata, cut = 0.4)
#' }
#'
#' @export
bfactor <- function(data, nfactors = NULL, mineigen = 1, rotate = "varimax",
                    fm = "pa", cut = 0.4, n.obs = NULL,
                    na_action = "complete", verbose = FALSE, ...) {

  # Load required package
  if (!requireNamespace("psych", quietly = TRUE)) {
    stop("Package 'psych' is required. Please install it.")
  }

  # Handle DocumentTermMatrix and other sparse matrix types
  if (inherits(data, "DocumentTermMatrix") || inherits(data, "simple_triplet_matrix")) {
    data <- as.matrix(data)
    if (verbose) message("DocumentTermMatrix converted to matrix for factor analysis")
  }

  # Use standardized input handling
  data <- bprepare_data(data,
                        na_action = na_action,
                        numeric_only = TRUE,
                        output_format = "data.frame",
                        verbose = verbose)

  n_vars <- ncol(data)

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
  n_params <- retained_factors * n_vars - retained_factors * (retained_factors - 1) / 2

  # --- Print formatted output using shared utilities ---

  bprint_header("FACTOR ANALYSIS")

  # Print method info
  method_name <- switch(fm,
                        "pa" = "principal factors",
                        "ml" = "maximum likelihood",
                        "minres" = "minimum residual",
                        "wls" = "weighted least squares",
                        "gls" = "generalized least squares",
                        fm)
  rotation_name <- if (rotate == "none") "(unrotated)" else paste0("(", rotate, ")")

  bprint_info(
    "Method" = method_name,
    "Rotation" = rotation_name,
    "Number of obs" = n.obs,
    "Retained factors" = retained_factors,
    "Number of params" = n_params
  )

  # Print correlation matrix (only if smaller than 26x26)
  if (n_vars < 26) {
    bprint_section("Correlation Matrix")
    print(round(cor_mat, 3))
  }

  # Print eigenvalue table
  bprint_eigenvalues(all_eigenvalues, threshold = mineigen)

  # Print LR test if available (for ML method)
  if (fm == "ml" && !is.null(fa_result$STATISTIC)) {
    cat(sprintf("\nLR test: independent vs. saturated: chi2(%d) = %.2f, p = %s\n",
                fa_result$dof, fa_result$STATISTIC, bformat_pval(fa_result$PVAL, width = 0)))
  }

  # Print unrotated factor loadings (retained factors only)
  bprint_section("Factor Loadings (Pattern Matrix) and Unique Variances")

  loadings_mat <- fa_result$loadings[1:n_vars, 1:retained_factors, drop = FALSE]
  uniqueness <- fa_result$uniquenesses

  # Build loading table with uniqueness
  var_names <- rownames(loadings_mat)
  max_var_length <- max(nchar(var_names), nchar("Variable"))

  # Header
  header <- sprintf("%-*s", max_var_length, "Variable")
  for (i in 1:retained_factors) {
    header <- paste0(header, sprintf(" %9s", paste0("Factor", i)))
  }
  header <- paste0(header, sprintf(" %12s", "Uniqueness"))
  cat(header, "\n")
  cat(strrep("-", nchar(header)), "\n")

  # Rows
  for (i in 1:n_vars) {
    row_str <- sprintf("%-*s", max_var_length, var_names[i])
    for (j in 1:retained_factors) {
      row_str <- paste0(row_str, sprintf(" %9.4f", loadings_mat[i, j]))
    }
    row_str <- paste0(row_str, sprintf(" %12.4f", uniqueness[i]))
    cat(row_str, "\n")
  }
  cat(strrep("-", nchar(header)), "\n")

  # Print rotated factor loadings (sorted and with cutoff) - retained factors only
  if (rotate != "none") {
    bprint_section("Rotated Factor Loadings (Pattern Matrix)")
    cat("Loadings below", cut, "suppressed\n\n")

    # Sort variables by their maximum absolute loading
    loadings_rot <- fa_result$loadings[1:n_vars, 1:retained_factors, drop = FALSE]
    max_loadings <- apply(abs(loadings_rot), 1, max)
    which_factor <- apply(abs(loadings_rot), 1, which.max)

    # Sort by factor, then by loading within factor
    sort_order <- order(which_factor, -max_loadings)
    loadings_sorted <- loadings_rot[sort_order, , drop = FALSE]

    bprint_loadings(loadings_sorted, cut = cut, max_vars = n_vars)
  }

  # Calculate and print reliability measures
  bprint_header("Factor Reliability Measures")

  reliability_df <- breliability(fa_result, data, cut = cut,
                                 use_rotated = (rotate != "none"))

  # Print formatted table
  cat(sprintf("\n%-8s %8s %10s %12s %10s\n",
              "Factor", "N Items", "Alpha", "Std Alpha", "Omega"))
  cat(strrep("-", 50), "\n")

  for (i in 1:nrow(reliability_df)) {
    cat(sprintf("%-8d %8d %10s %12s %10s\n",
                reliability_df$factor[i],
                reliability_df$n_items[i],
                if (is.na(reliability_df$alpha[i])) "-" else sprintf("%.4f", reliability_df$alpha[i]),
                if (is.na(reliability_df$std_alpha[i])) "-" else sprintf("%.4f", reliability_df$std_alpha[i]),
                if (is.na(reliability_df$omega_t[i])) "-" else sprintf("%.4f", reliability_df$omega_t[i])))
  }
  cat(strrep("-", 50), "\n")

  # Print items for each factor
  cat("\nItems by Factor:\n")
  for (i in 1:nrow(reliability_df)) {
    if (reliability_df$items[i] != "") {
      cat(sprintf("  Factor %d: %s\n", reliability_df$factor[i], reliability_df$items[i]))
    } else {
      cat(sprintf("  Factor %d: No items above cutoff\n", reliability_df$factor[i]))
    }
  }

  cat("\n")
  bprint_sep(char = "=")

  # Return the fa object invisibly
  invisible(fa_result)
}
