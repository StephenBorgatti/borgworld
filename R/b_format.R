# b_format.R - Standardized Output Formatting Utilities
#
# This file contains shared utilities for formatting and printing output
# across all borgworld functions, ensuring consistent presentation of:
# - Headers and titles
# - Tables
# - Numbers (coefficients, p-values, etc.)
# - Separators and dividers

# Default line width for output
BFORMAT_DEFAULT_WIDTH <- 79

#' Print Analysis Header
#'
#' @description
#' Prints a formatted header for analysis output.
#'
#' @param title The main title text
#' @param subtitle Optional subtitle text
#' @param width Line width (default: 79)
#' @param char Character to use for border (default: "=")
#'
#' @return Invisibly returns NULL (called for side effects)
#'
#' @examples
#' bprint_header("PCA RESULTS")
#' bprint_header("REGRESSION ANALYSIS", subtitle = "Robust Standard Errors")
#'
#' @export
bprint_header <- function(title, subtitle = NULL, width = BFORMAT_DEFAULT_WIDTH, char = "=") {
  cat("\n")
  cat(strrep(char, width), "\n")
  cat(title, "\n")
  if (!is.null(subtitle)) {
    cat(subtitle, "\n")
  }
  cat(strrep(char, width), "\n")
  invisible(NULL)
}


#' Print Section Header
#'
#' @description
#' Prints a formatted section header (less prominent than main header).
#'
#' @param title The section title text
#' @param width Line width (default: 79)
#' @param char Character to use for border (default: "-")
#'
#' @return Invisibly returns NULL (called for side effects)
#'
#' @examples
#' bprint_section("Eigenvalues and Variance Explained")
#'
#' @export
bprint_section <- function(title, width = BFORMAT_DEFAULT_WIDTH, char = "-") {
  cat("\n")
  cat(title, "\n")
  cat(strrep(char, width), "\n")
  invisible(NULL)
}


#' Print Separator Line
#'
#' @description
#' Prints a separator line.
#'
#' @param width Line width (default: 79)
#' @param char Character to use (default: "-")
#' @param newline_before Add newline before separator (default: FALSE)
#' @param newline_after Add newline after separator (default: FALSE)
#'
#' @return Invisibly returns NULL (called for side effects)
#'
#' @examples
#' bprint_sep()
#' bprint_sep(char = "=", newline_before = TRUE)
#'
#' @export
bprint_sep <- function(width = BFORMAT_DEFAULT_WIDTH, char = "-",
                       newline_before = FALSE, newline_after = FALSE) {
  if (newline_before) cat("\n")
  cat(strrep(char, width), "\n")
  if (newline_after) cat("\n")
  invisible(NULL)
}


#' Print Key-Value Pairs
#'
#' @description
#' Prints formatted key-value information (like "Number of obs: 100").
#'
#' @param ... Named arguments where names are labels and values are the values
#' @param label_width Width for labels (default: auto-calculated)
#' @param value_format Format string for values (default: auto-detect)
#' @param sep Separator between label and value (default: ": ")
#'
#' @return Invisibly returns NULL (called for side effects)
#'
#' @examples
#' bprint_info(
#'   "Number of obs" = 100,
#'   "R-squared" = 0.85,
#'   "Variables" = 5
#' )
#'
#' @export
bprint_info <- function(..., label_width = NULL, value_format = NULL, sep = ": ") {
  args <- list(...)
  labels <- names(args)

  if (is.null(labels) || any(labels == "")) {
    stop("All arguments must be named")
  }

  # Calculate label width if not specified
  if (is.null(label_width)) {
    label_width <- max(nchar(labels))
  }

  for (i in seq_along(args)) {
    label <- labels[i]
    value <- args[[i]]

    # Format value based on type
    if (is.null(value_format)) {
      if (is.numeric(value)) {
        if (abs(value) >= 1000 || (abs(value) < 0.001 && value != 0)) {
          value_str <- formatC(value, digits = 4, format = "g")
        } else if (value == round(value)) {
          value_str <- formatC(value, format = "d", big.mark = ",")
        } else {
          value_str <- sprintf("%.4f", value)
        }
      } else {
        value_str <- as.character(value)
      }
    } else {
      value_str <- sprintf(value_format, value)
    }

    cat(sprintf("%-*s%s%s\n", label_width, label, sep, value_str))
  }

  invisible(NULL)
}


#' Format Number for Display
#'
#' @description
#' Formats a number for display with appropriate precision and width.
#'
#' @param x Numeric value to format
#' @param width Target width (default: 10)
#' @param digits Number of decimal places (default: 4)
#' @param scientific Use scientific notation for very large/small values (default: TRUE)
#'
#' @return Character string with formatted number
#'
#' @examples
#' bformat_num(1234.5678)
#' bformat_num(0.00001)
#' bformat_num(1e10)
#'
#' @export
bformat_num <- function(x, width = 10, digits = 4, scientific = TRUE) {
  if (is.na(x) || is.infinite(x)) {
    return(sprintf("%*s", width, ifelse(is.na(x), "NA", ifelse(x > 0, "Inf", "-Inf"))))
  }

  if (scientific) {
    # Use scientific notation for very small or very large numbers
    if (abs(x) < 1e-4 && x != 0) {
      formatted <- sprintf("%.2e", x)
    } else if (abs(x) >= 1e7) {
      formatted <- sprintf("%.2e", x)
    } else if (abs(x) >= 10000) {
      formatted <- sprintf("%.*f", max(0, digits - 2), x)
    } else if (abs(x) >= 100) {
      formatted <- sprintf("%.*f", max(0, digits - 1), x)
    } else if (abs(x) >= 10) {
      formatted <- sprintf("%.*f", digits, x)
    } else if (abs(x) >= 1) {
      formatted <- sprintf("%.*f", digits, x)
    } else {
      formatted <- sprintf("%.*f", digits, x)
    }

    # If still too wide, force scientific notation
    if (nchar(formatted) > width) {
      formatted <- sprintf("%.2e", x)
    }
  } else {
    formatted <- sprintf("%.*f", digits, x)
  }

  # Right-align in the field
  sprintf("%*s", width, formatted)
}


#' Format P-value for Display
#'
#' @description
#' Formats a p-value with appropriate precision.
#'
#' @param p P-value to format
#' @param width Target width (default: 7)
#' @param digits Decimal places for values >= 0.001 (default: 3)
#' @param min_display Minimum value to display numerically (default: 0.0001)
#'
#' @return Character string with formatted p-value
#'
#' @examples
#' bformat_pval(0.0342)
#' bformat_pval(0.00001)
#' bformat_pval(0.1234)
#'
#' @export
bformat_pval <- function(p, width = 7, digits = 3, min_display = 0.0001) {
  if (is.na(p)) {
    return(sprintf("%*s", width, "NA"))
  }

  if (p < min_display) {
    formatted <- sprintf("<%.*f", digits, min_display)
  } else if (p < 0.001) {
    formatted <- sprintf("%.4f", p)
  } else {
    formatted <- sprintf("%.*f", digits, p)
  }

  sprintf("%*s", width, formatted)
}


#' Get Significance Stars
#'
#' @description
#' Returns significance stars based on p-value.
#'
#' @param p P-value
#' @param levels Named numeric vector of significance levels (default: standard levels)
#'
#' @return Character string with significance stars
#'
#' @examples
#' bsig_stars(0.001)   # "***"
#' bsig_stars(0.03)    # "*"
#' bsig_stars(0.15)    # ""
#'
#' @export
bsig_stars <- function(p, levels = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "." = 0.1)) {
  if (is.na(p)) return("")

  for (stars in names(levels)) {
    if (p < levels[stars]) {
      return(stars)
    }
  }
  return("")
}


#' Print Significance Legend
#'
#' @description
#' Prints the standard significance legend.
#'
#' @param levels Significance levels to display (default: standard levels)
#'
#' @return Invisibly returns NULL (called for side effects)
#'
#' @examples
#' bprint_sig_legend()
#'
#' @export
bprint_sig_legend <- function(levels = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "." = 0.1)) {
  parts <- paste(sprintf("'%s' %g", names(levels), levels), collapse = " ")
  cat("Signif. codes: 0", parts, "\n")
  invisible(NULL)
}


#' Print Formatted Table
#'
#' @description
#' Prints a data frame as a formatted table with aligned columns.
#'
#' @param df Data frame to print
#' @param col_widths Named vector of column widths (optional)
#' @param col_formats Named vector of sprintf format strings (optional)
#' @param row_names Include row names (default: FALSE)
#' @param header Include header row (default: TRUE)
#' @param sep_char Character for header separator (default: "-")
#' @param max_rows Maximum rows to print (default: NULL for all)
#'
#' @return Invisibly returns NULL (called for side effects)
#'
#' @examples
#' df <- data.frame(
#'   Variable = c("wt", "hp"),
#'   Coef = c(-3.878, -0.032),
#'   SE = c(0.633, 0.009)
#' )
#' bprint_table(df)
#'
#' @export
bprint_table <- function(df,
                         col_widths = NULL,
                         col_formats = NULL,
                         row_names = FALSE,
                         header = TRUE,
                         sep_char = "-",
                         max_rows = NULL) {

  if (!is.data.frame(df)) {
    df <- as.data.frame(df)
  }

  # Calculate column widths if not provided
  if (is.null(col_widths)) {
    col_widths <- sapply(names(df), function(col) {
      max(nchar(col), max(nchar(as.character(df[[col]])), na.rm = TRUE))
    })
  }

  # Build header
  if (header) {
    header_str <- ""
    for (i in seq_along(names(df))) {
      col <- names(df)[i]
      w <- col_widths[col]
      header_str <- paste0(header_str, sprintf("%*s", w, col))
      if (i < length(names(df))) header_str <- paste0(header_str, "  ")
    }
    cat(header_str, "\n")
    cat(strrep(sep_char, nchar(header_str)), "\n")
  }

  # Determine rows to print
  if (!is.null(max_rows) && nrow(df) > max_rows) {
    rows_to_print <- 1:max_rows
    truncated <- TRUE
  } else {
    rows_to_print <- 1:nrow(df)
    truncated <- FALSE
  }

  # Print rows
  for (r in rows_to_print) {
    row_str <- ""
    for (i in seq_along(names(df))) {
      col <- names(df)[i]
      w <- col_widths[col]
      val <- df[r, col]

      # Format value
      if (!is.null(col_formats) && col %in% names(col_formats)) {
        formatted <- sprintf(col_formats[col], val)
      } else if (is.numeric(val)) {
        formatted <- bformat_num(val, width = w)
      } else {
        formatted <- sprintf("%*s", w, as.character(val))
      }

      row_str <- paste0(row_str, formatted)
      if (i < length(names(df))) row_str <- paste0(row_str, "  ")
    }
    cat(row_str, "\n")
  }

  if (truncated) {
    cat("... [", nrow(df) - max_rows, " more rows]\n", sep = "")
  }

  invisible(NULL)
}


#' Print Eigenvalue Table
#'
#' @description
#' Prints a formatted eigenvalue table with variance explained.
#'
#' @param eigenvalues Numeric vector of eigenvalues
#' @param max_display Maximum number of eigenvalues to display (default: 20)
#' @param threshold Highlight eigenvalues above this threshold (default: 1)
#'
#' @return Invisibly returns NULL (called for side effects)
#'
#' @examples
#' eigs <- c(3.2, 1.5, 0.8, 0.5)
#' bprint_eigenvalues(eigs)
#'
#' @export
bprint_eigenvalues <- function(eigenvalues, max_display = 20, threshold = 1) {
  k <- length(eigenvalues)
  n_show <- min(k, max_display)

  # Calculate statistics
  total <- sum(eigenvalues)
  pct <- 100 * eigenvalues / total
  cum_pct <- cumsum(pct)
  diff <- c(eigenvalues[-k] - eigenvalues[-1], NA)

  # Build data frame
  df <- data.frame(
    Component = 1:k,
    Eigenvalue = eigenvalues,
    Difference = diff,
    Percent = pct,
    Cumulative = cum_pct
  )

  # Print section header
  bprint_section("Eigenvalues and Variance Explained")

  # Header
  cat(sprintf("    %-10s %12s %12s %12s %12s\n",
              "Component", "Eigenvalue", "Difference", "Percent", "Cumulative"))
  cat(strrep("-", 79), "\n")

  # Rows
  for (i in 1:n_show) {
    diff_str <- if (is.na(df$Difference[i])) "." else sprintf("%12.4f", df$Difference[i])
    marker <- if (df$Eigenvalue[i] >= threshold) "*" else ""

    cat(sprintf("    %-10d %12.4f %12s %12.2f %12.2f%s\n",
                df$Component[i],
                df$Eigenvalue[i],
                diff_str,
                df$Percent[i],
                df$Cumulative[i],
                marker))
  }

  cat(strrep("-", 79), "\n")

  if (k > n_show) {
    cat("... [", k - n_show, " more components]\n", sep = "")
  }

  # Summary
  n_above <- sum(eigenvalues >= threshold)
  cat("\nComponents with eigenvalue >= ", threshold, ": ", n_above, "\n", sep = "")

  invisible(NULL)
}


#' Print Loading Matrix
#'
#' @description
#' Prints a formatted loading matrix (for PCA, factor analysis).
#'
#' @param loadings Matrix of loadings
#' @param cut Cutoff for suppressing small loadings (default: NULL, show all)
#' @param sort Sort variables by primary loading (default: FALSE)
#' @param max_vars Maximum variables to display (default: 20)
#'
#' @return Invisibly returns NULL (called for side effects)
#'
#' @examples
#' # Assuming 'loadings' is a loading matrix from PCA/FA
#' # bprint_loadings(loadings, cut = 0.3)
#'
#' @export
bprint_loadings <- function(loadings, cut = NULL, sort = FALSE, max_vars = 20) {
  if (!is.matrix(loadings)) {
    loadings <- as.matrix(loadings)
  }

  # Get dimensions
  n_vars <- nrow(loadings)
  n_factors <- ncol(loadings)

  # Set column names if missing
  if (is.null(colnames(loadings))) {
    colnames(loadings) <- paste0("PC", 1:n_factors)
  }
  if (is.null(rownames(loadings))) {
    rownames(loadings) <- paste0("Var", 1:n_vars)
  }

  # Sort by primary loading if requested
  if (sort) {
    primary_factor <- apply(abs(loadings), 1, which.max)
    primary_loading <- sapply(1:nrow(loadings), function(i) loadings[i, primary_factor[i]])
    sort_order <- order(primary_factor, -abs(primary_loading))
    loadings <- loadings[sort_order, , drop = FALSE]
  }

  # Calculate column widths
  var_width <- max(nchar(rownames(loadings)), nchar("Variable"))
  factor_width <- 10

  # Header
  header <- sprintf("%-*s", var_width, "Variable")
  for (j in 1:n_factors) {
    header <- paste0(header, sprintf(" %*s", factor_width, colnames(loadings)[j]))
  }
  cat(header, "\n")
  cat(strrep("-", nchar(header)), "\n")

  # Rows (up to max_vars)
  n_show <- min(n_vars, max_vars)
  for (i in 1:n_show) {
    row_str <- sprintf("%-*s", var_width, rownames(loadings)[i])

    for (j in 1:n_factors) {
      val <- loadings[i, j]

      if (!is.null(cut) && abs(val) < cut) {
        val_str <- sprintf("%*s", factor_width, "")
      } else {
        val_str <- sprintf(" %*.4f", factor_width - 1, val)
      }
      row_str <- paste0(row_str, val_str)
    }
    cat(row_str, "\n")
  }

  cat(strrep("-", nchar(header)), "\n")

  if (n_vars > n_show) {
    cat("... [", n_vars - n_show, " more variables]\n", sep = "")
  }

  invisible(NULL)
}


#' Print Coefficient Table (Regression Style)
#'
#' @description
#' Prints a formatted coefficient table for regression models.
#'
#' @param coef_df Data frame with columns: name, estimate, se, statistic, p.value
#' @param ci_lower Optional lower CI bounds
#' @param ci_upper Optional upper CI bounds
#' @param dep_var Name of dependent variable (for header)
#' @param show_sig Show significance stars (default: TRUE)
#'
#' @return Invisibly returns NULL (called for side effects)
#'
#' @export
bprint_coef_table <- function(coef_df, ci_lower = NULL, ci_upper = NULL,
                              dep_var = "Y", show_sig = TRUE) {

  # Determine variable name width
  var_width <- max(nchar(as.character(coef_df[[1]])), nchar(dep_var), 12)

  # Build header
  if (!is.null(ci_lower) && !is.null(ci_upper)) {
    header <- sprintf("%*s | %11s %10s %8s %9s     [95%% CI]",
                      var_width, dep_var,
                      "Coefficient", "Std.Err", "t/z", "P>|t|")
  } else {
    header <- sprintf("%*s | %11s %10s %8s %9s",
                      var_width, dep_var,
                      "Coefficient", "Std.Err", "t/z", "P>|t|")
  }

  total_width <- nchar(header)
  cat(strrep("-", total_width), "\n")
  cat(header, "\n")
  cat(paste0(strrep("-", var_width), "+", strrep("-", total_width - var_width - 1)), "\n")

  # Print each row
  for (i in 1:nrow(coef_df)) {
    name <- as.character(coef_df[i, 1])
    est <- coef_df[i, 2]
    se <- coef_df[i, 3]
    stat <- coef_df[i, 4]
    pval <- coef_df[i, 5]

    # Format values
    est_str <- bformat_num(est, width = 11)
    se_str <- bformat_num(se, width = 10)
    stat_str <- sprintf("%8.2f", stat)
    pval_str <- bformat_pval(pval, width = 9)

    # Significance stars
    stars <- if (show_sig) bsig_stars(pval) else ""

    if (!is.null(ci_lower) && !is.null(ci_upper)) {
      ci_str <- sprintf("[%9.4f, %9.4f]",
                        ci_lower[i], ci_upper[i])
      cat(sprintf("%*s | %s %s %s %s %s %s\n",
                  var_width, name,
                  est_str, se_str, stat_str, pval_str, stars, ci_str))
    } else {
      cat(sprintf("%*s | %s %s %s %s %s\n",
                  var_width, name,
                  est_str, se_str, stat_str, pval_str, stars))
    }
  }

  cat(strrep("-", total_width), "\n")

  if (show_sig) {
    bprint_sig_legend()
  }

  invisible(NULL)
}


#' Print Model Summary Statistics
#'
#' @description
#' Prints summary statistics for a model (R-squared, F-stat, etc.).
#'
#' @param n Number of observations
#' @param r_squared R-squared value
#' @param adj_r_squared Adjusted R-squared
#' @param f_stat F-statistic (optional)
#' @param f_pvalue F-test p-value (optional)
#' @param rmse Root mean squared error (optional)
#' @param aic AIC (optional)
#' @param bic BIC (optional)
#'
#' @return Invisibly returns NULL (called for side effects)
#'
#' @export
bprint_model_stats <- function(n, r_squared = NULL, adj_r_squared = NULL,
                               f_stat = NULL, f_pvalue = NULL,
                               rmse = NULL, aic = NULL, bic = NULL) {

  cat("\n")
  bprint_info("Number of obs" = n)

  if (!is.null(r_squared)) {
    bprint_info("R-squared" = r_squared)
  }
  if (!is.null(adj_r_squared)) {
    bprint_info("Adj R-squared" = adj_r_squared)
  }
  if (!is.null(f_stat) && !is.null(f_pvalue)) {
    cat(sprintf("%-20s = %.4f (p = %s)\n", "F-statistic", f_stat, bformat_pval(f_pvalue, width = 0)))
  }
  if (!is.null(rmse)) {
    bprint_info("Root MSE" = rmse)
  }
  if (!is.null(aic)) {
    bprint_info("AIC" = aic)
  }
  if (!is.null(bic)) {
    bprint_info("BIC" = bic)
  }

  cat("\n")
  invisible(NULL)
}
