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
  if (!require("psych", quietly = TRUE)) {
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
#'
#' # Get reliability measures separately
#' fa_result <- bfactor(mydata, nfactors = 3)
#' reliability <- breliability(fa_result, mydata, cut = 0.4)
#' }
#'
#' @export
bfactor <- function(data, nfactors = NULL, mineigen = 1, rotate = "varimax",
                    fm = "pa", cut = 0.4, n.obs = NULL, ...) {
  # Load required package
  if (!require("psych", quietly = TRUE)) {
    stop("Package 'psych' is required. Please install it.")
  }

  # drop charater variables
  nvar <- ncol(data)
  data <-  data |>
    select(where(is.numeric))
  if (nvar != ncol(data)) {
    print('Character variable removed')
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

  # Print correlation matrix
  cat("\n")
  cat(strrep("=", 79))
  cat("\n")
  cat("Correlation Matrix:\n")
  cat(strrep("=", 79))
  cat("\n\n")
  print(round(cor_mat, 3))

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

  # Calculate maximum variable name length for proper column alignment
  var_names <- rownames(fa_result$loadings)
  max_var_length <- max(nchar(var_names), nchar("Variable"))

  # Print unrotated factor loadings (retained factors only)
  cat("\nFactor loadings (pattern matrix) and unique variances\n\n")

  # Adjust line width based on variable name length
  line_width <- max(70, 4 + max_var_length + retained_factors * 10 + 13)
  cat(strrep("-", line_width))
  cat("\n")

  loading_header <- sprintf(paste0("    %-", max_var_length, "s"), "Variable")
  for (i in 1:retained_factors) {
    loading_header <- paste0(loading_header, sprintf(" %9s", paste0("Factor", i)))
  }
  loading_header <- paste0(loading_header, sprintf(" %12s", "Uniqueness"))
  cat(loading_header, "\n")
  cat(strrep("-", line_width))
  cat("\n")

  loadings_mat <- fa_result$loadings[1:n_vars, 1:retained_factors, drop = FALSE]
  uniqueness <- fa_result$uniquenesses

  for (i in 1:n_vars) {
    row_str <- sprintf(paste0("    %-", max_var_length, "s"), rownames(loadings_mat)[i])
    for (j in 1:retained_factors) {
      row_str <- paste0(row_str, sprintf(" %9.4f", loadings_mat[i, j]))
    }
    row_str <- paste0(row_str, sprintf(" %12.4f", uniqueness[i]))
    cat(row_str, "\n")
  }
  cat(strrep("-", line_width))
  cat("\n")

  # Print rotated factor loadings (sorted and with cutoff) - retained factors only
  # No uniqueness column for rotated loadings
  if (rotate != "none") {
    cat("\nRotated factor loadings (pattern matrix)\n\n")

    # Adjust line width for rotated section (no uniqueness column)
    line_width_rot <- max(70, 4 + max_var_length + retained_factors * 10)
    cat(strrep("-", line_width_rot))
    cat("\n")

    # Header without uniqueness column
    loading_header_rot <- sprintf(paste0("    %-", max_var_length, "s"), "Variable")
    for (i in 1:retained_factors) {
      loading_header_rot <- paste0(loading_header_rot, sprintf(" %9s", paste0("Factor", i)))
    }
    cat(loading_header_rot, "\n")
    cat(strrep("-", line_width_rot))
    cat("\n")

    # Sort variables by their maximum absolute loading
    loadings_rot <- fa_result$loadings[1:n_vars, 1:retained_factors, drop = FALSE]
    max_loadings <- apply(abs(loadings_rot), 1, max)
    which_factor <- apply(abs(loadings_rot), 1, which.max)

    # Sort by factor, then by loading within factor
    sort_order <- order(which_factor, -max_loadings)
    loadings_sorted <- loadings_rot[sort_order, , drop = FALSE]

    for (i in 1:n_vars) {
      row_str <- sprintf(paste0("    %-", max_var_length, "s"), rownames(loadings_sorted)[i])
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
    cat(strrep("-", line_width_rot))
    cat("\n")
  }

  # Calculate and print reliability measures
  cat("\n")
  cat(strrep("=", 79))
  cat("\n")
  cat("Factor Reliability Measures\n")
  cat(strrep("=", 79))
  cat("\n")

  reliability_df <- breliability(fa_result, data, cut = cut,
                                 use_rotated = (rotate != "none"))

  # Print formatted table header
  cat(sprintf("\n    %-8s %8s %10s %12s %10s\n",
              "Factor", "N Items", "Alpha", "Std Alpha", "Omega"))
  cat(strrep("-", 79))
  cat("\n")

  # Print each row
  for (i in 1:nrow(reliability_df)) {
    cat(sprintf("    %-8d %8d %10s %12s %10s\n",
                reliability_df$factor[i],
                reliability_df$n_items[i],
                if (is.na(reliability_df$alpha[i])) "-" else sprintf("%.4f", reliability_df$alpha[i]),
                if (is.na(reliability_df$std_alpha[i])) "-" else sprintf("%.4f", reliability_df$std_alpha[i]),
                if (is.na(reliability_df$omega_t[i])) "-" else sprintf("%.4f", reliability_df$omega_t[i])))
  }
  cat(strrep("-", 79))
  cat("\n")

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
  cat(strrep("=", 79))
  cat("\n")

  # Return the fa object invisibly
  invisible(fa_result)
}
