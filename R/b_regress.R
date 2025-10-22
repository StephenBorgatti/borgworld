#' Stata-style regression output for linear models
#'
#' Produces regression output that mimics Stata's regress command format
#'
#' @param data A data frame containing the variables
#' @param formula A formula object specifying the model
#' @param robust Logical, whether to use robust standard errors (requires sandwich package)
#'
#' @return Invisibly returns the fitted lm model object
#'
#' @importFrom stats lm coef confint nobs resid pf pt qt
#' @export
#'
#' @examples
#' # Basic regression
#' bregress(mtcars, mpg ~ wt + hp)
#'
#' # With piping
#' library(dplyr)
#' mtcars %>% bregress(mpg ~ wt + hp)
#'
#' # With robust standard errors (if sandwich package is installed)
#' \dontrun{
#' bregress(mtcars, mpg ~ wt + hp, robust = TRUE)
#' }
bregress <- function(data, formula, robust = FALSE) {
  # Convert string to formula if needed
  if (is.character(formula)) {
    formula <- as.formula(formula)
  }

  # Fit the model
  model <- lm(formula, data = data)

  # If robust standard errors requested, check for sandwich package
  if (robust) {
    if (!requireNamespace("sandwich", quietly = TRUE)) {
      stop("Package 'sandwich' needed for robust standard errors. Please install it.",
           call. = FALSE)
    }
    # Get robust standard errors
    robust_se <- sqrt(diag(sandwich::vcovHC(model, type = "HC1")))
    coef_vals <- coef(model)
    t_vals <- coef_vals / robust_se
    p_vals <- 2 * pt(abs(t_vals), df = model$df.residual, lower.tail = FALSE)

    # Override the coefficient summary
    coef_summary <- cbind(
      Estimate = coef_vals,
      `Std. Error` = robust_se,
      `t value` = t_vals,
      `Pr(>|t|)` = p_vals
    )

    # Recalculate confidence intervals with robust SE
    conf_int <- cbind(
      coef_vals - qt(0.975, model$df.residual) * robust_se,
      coef_vals + qt(0.975, model$df.residual) * robust_se
    )
  } else {
    coef_summary <- summary(model)$coefficients
    conf_int <- confint(model, level = 0.95)
  }

  # Extract basic information
  n <- nobs(model)
  k <- length(coef(model)) - 1
  df_model <- k
  df_resid <- model$df.residual
  df_total <- n - 1

  # Calculate sums of squares
  ss_resid <- sum(resid(model)^2)
  ss_total <- sum((model$model[[1]] - mean(model$model[[1]]))^2)
  ss_model <- ss_total - ss_resid

  # Mean squares
  ms_model <- ss_model / df_model
  ms_resid <- ss_resid / df_resid

  # F-statistic and p-value
  f_stat <- ms_model / ms_resid
  f_pvalue <- pf(f_stat, df_model, df_resid, lower.tail = FALSE)

  # R-squared and adjusted R-squared
  r_squared <- ss_model / ss_total
  adj_r_squared <- 1 - (1 - r_squared) * (df_total) / df_resid

  # Root MSE
  root_mse <- sqrt(ms_resid)

  # Get coefficient names
  coef_names <- rownames(coef_summary)
  coef_names[coef_names == "(Intercept)"] <- "_cons"

  # Print header
  if (robust) {
    cat("\nLinear regression with robust standard errors\n")
  }

  # Print ANOVA table with fixed column widths (Stata-like)
  w_src <- 12; w_ss <- 14; w_df <- 9; w_ms <- 14

  header <- sprintf("%*s | %*s %*s %*s   Number of obs   = %9d",
                    w_src, "Source", w_ss, "SS", w_df, "df", w_ms, "MS", n)
  cat("\n", header, "\n", sep = "")
  cat(paste0(rep("-", nchar(header)), collapse = ""), "\n", sep = "")

  row_fmt <- function(src, ss, df, ms, tail = "") {
    ss_s <- formatC(ss, digits = 7, format = "g")
    ms_s <- formatC(ms, digits = 7, format = "g")
    df_s <- formatC(df, format = "d")
    sprintf("%*s | %*s %*s %*s%s",
            w_src, src,
            w_ss, ss_s,
            w_df, df_s,
            w_ms, ms_s,
            tail)
  }

  cat(row_fmt("Model",    ss_model,  df_model, ms_model,
              sprintf("   Prob > F        = %10.4f", f_pvalue)), "\n")
  cat(row_fmt("Residual", ss_resid,  df_resid, ms_resid,
              sprintf("   R-squared       = %10.4f", r_squared)), "\n")

  # Fixed separator line - align "+" with "|" symbols
  sep <- sprintf("%*s+%s", w_src + 1, "",
                 paste0(rep("-", w_ss + 1 + w_df + 1 + w_ms + 1), collapse = ""))
  cat(sep, sprintf("   Adj R-squared   = %10.4f", adj_r_squared), "\n", sep = "")

  cat(row_fmt("Total",    ss_total,  df_total, ss_total/df_total,
              sprintf("   Root MSE        = %9.3f", root_mse)), "\n")

  # Print coefficient table
  cat("\n")

  # Determine dependent variable name
  dep_var <- as.character(formula[[2]])

  # Find the longest variable name to determine column width
  all_var_names <- coef_names
  max_var_length <- max(nchar(all_var_names), nchar(dep_var))
  # Use at least the same width as the "Source" column (12 chars) for consistency
  var_col_width <- max(12, max_var_length)

  # Header for coefficient table with properly aligned columns
  header_line <- sprintf("%*s | Coefficient  Std. err.      t    P>|t|     [95%% conf. interval]",
                         var_col_width, dep_var)

  # Calculate the exact width based on the header line
  total_width <- nchar(header_line)

  # Create separator lines that match the header width exactly
  separator_line <- paste0(rep("-", total_width), collapse = "")
  cat(separator_line, "\n")
  cat(header_line, "\n", sep = "")

  # Separator line between header and data
  sep_line <- paste0(paste0(rep("-", var_col_width), collapse = ""),
                     "+",
                     paste0(rep("-", total_width - var_col_width - 1), collapse = ""))
  cat(sep_line, "\n")

  # Helper function to format numbers that might need scientific notation
  format_coef <- function(x, width = 11) {
    # Determine the best format based on the magnitude
    if (is.na(x) || is.infinite(x)) {
      formatted <- sprintf("%*s", width, "NA")
    } else if (abs(x) < 1e-4 && x != 0) {
      # Very small numbers: use scientific notation
      formatted <- sprintf("%.2e", x)
    } else if (abs(x) >= 1e7) {
      # Very large numbers: use scientific notation
      formatted <- sprintf("%.2e", x)
    } else if (abs(x) >= 10000) {
      # Large numbers: fewer decimal places
      formatted <- sprintf("%.1f", x)
    } else if (abs(x) >= 100) {
      # Medium-large numbers
      formatted <- sprintf("%.3f", x)
    } else if (abs(x) >= 10) {
      # Medium numbers
      formatted <- sprintf("%.4f", x)
    } else if (abs(x) >= 1) {
      # Small-medium numbers
      formatted <- sprintf("%.5f", x)
    } else {
      # Small numbers but not tiny
      formatted <- sprintf("%.6f", x)
    }

    # If it's still too wide, force scientific notation
    if (nchar(formatted) > width) {
      formatted <- sprintf("%.2e", x)
    }

    # Right-align in the field
    sprintf("%*s", width, formatted)
  }

  # Print each coefficient row with better spacing
  for (i in 1:length(coef_names)) {
    var_name <- coef_names[i]

    coef_val <- coef_summary[i, 1]
    se_val <- coef_summary[i, 2]
    t_val <- coef_summary[i, 3]
    p_val <- coef_summary[i, 4]
    ci_low <- conf_int[i, 1]
    ci_high <- conf_int[i, 2]

    # Handle very small p-values
    if (p_val < 0.0005) {
      p_str <- "  0.000"
    } else {
      p_str <- sprintf("%7.3f", p_val)
    }

    # Print with consistent column positions matching the header
    cat(sprintf("%*s", var_col_width, var_name), " | ",
        format_coef(coef_val, 10), "  ",  # Coefficient column
        format_coef(se_val, 9), "  ",      # Std. err. column
        sprintf("%6.2f", t_val), "  ",     # t column
        p_str, "    ",                      # P>|t| column
        format_coef(ci_low, 10), " ",      # Lower CI (reduced spacing)
        format_coef(ci_high, 10),          # Upper CI (ends at the "]")
        "\n", sep = "")
  }
  cat(separator_line, "\n")

  # Calculate VIF and Tolerance for predictors (exclude intercept)
  if (k > 1) {  # Only if there are 2 or more predictors
    # Get the model matrix without intercept
    X <- model.matrix(model)[, -1, drop = FALSE]

    # Calculate VIF for each predictor
    vif_values <- numeric(k)
    tolerance_values <- numeric(k)
    predictor_names <- names(coef(model))[-1]  # Exclude intercept

    for (j in 1:k) {
      # Regress each predictor on all other predictors
      if (k > 1) {
        X_j <- X[, j]
        X_others <- X[, -j, drop = FALSE]
        r_squared_j <- summary(lm(X_j ~ X_others))$r.squared
        vif_values[j] <- 1 / (1 - r_squared_j)
        tolerance_values[j] <- 1 - r_squared_j
      } else {
        # If only one predictor, VIF = 1
        vif_values[j] <- 1
        tolerance_values[j] <- 1
      }
    }

    # Print VIF/Tolerance table
    cat("\n")

    # Create VIF table header with R-squared, Tolerance, VIF, and Sqrt(VIF)
    vif_header <- sprintf("%*s |    R-squared   Tolerance         VIF   Sqrt(VIF)",
                          var_col_width, "Variable")
    vif_table_width <- nchar(vif_header)
    vif_separator <- paste0(rep("-", vif_table_width), collapse = "")

    cat(vif_separator, "\n")
    cat(vif_header, "\n")
    cat(paste0(paste0(rep("-", var_col_width), collapse = ""),
               "+",
               paste0(rep("-", vif_table_width - var_col_width - 1), collapse = "")), "\n")

    # Calculate R-squared values for display
    r_squared_values <- numeric(k)
    for (j in 1:k) {
      if (k > 1) {
        r_squared_values[j] <- 1 - tolerance_values[j]
      } else {
        r_squared_values[j] <- 0
      }
    }

    # Print VIF values for each predictor
    for (j in 1:k) {
      var_name <- predictor_names[j]
      if (nchar(var_name) > var_col_width) {
        var_name <- substr(var_name, 1, var_col_width)
      }

      # Format all values
      rsq_str <- sprintf("%10.4f", r_squared_values[j])
      tol_str <- sprintf("%10.4f", tolerance_values[j])
      vif_str <- sprintf("%10.2f", vif_values[j])
      sqrt_vif_str <- sprintf("%10.4f", sqrt(vif_values[j]))

      cat(sprintf("%*s", var_col_width, var_name), " | ",
          rsq_str, "  ",
          tol_str, "  ",
          vif_str, "  ",
          sqrt_vif_str,
          "\n", sep = "")
    }

    # Add mean VIF at the bottom
    mean_vif <- mean(vif_values)
    mean_rsq <- mean(r_squared_values)
    mean_tol <- mean(tolerance_values)
    mean_sqrt_vif <- sqrt(mean_vif)

    cat(paste0(paste0(rep("-", var_col_width), collapse = ""),
               "+",
               paste0(rep("-", vif_table_width - var_col_width - 1), collapse = "")), "\n")
    cat(sprintf("%*s", var_col_width, "Mean"), " | ",
        sprintf("%10.4f", mean_rsq), "  ",
        sprintf("%10.4f", mean_tol), "  ",
        sprintf("%10.2f", mean_vif), "  ",
        sprintf("%10.4f", mean_sqrt_vif),
        "\n", sep = "")
    cat(vif_separator, "\n")
  }

  # Return the model invisibly
  invisible(model)
}
