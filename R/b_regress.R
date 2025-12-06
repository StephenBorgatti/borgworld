#' Stata-style regression output for linear models
#'
#' Produces regression output that mimics Stata's regress command format.
#' Accepts both data frames and matrices as input.
#'
#' @param data A data frame or matrix containing the variables
#' @param formula A formula object or character string specifying the model
#' @param robust Logical, whether to use robust standard errors (requires sandwich package)
#' @param alpha Significance level for confidence intervals (default: 0.05)
#' @param vif Logical, whether to print VIF table (default: TRUE for multiple predictors)
#' @param verbose Logical; if TRUE, print messages about data cleaning
#'
#' @return Invisibly returns the fitted lm model object
#'
#' @importFrom stats lm coef confint nobs resid pf pt qt fitted effects model.matrix model.frame model.response complete.cases
#' @export
#'
#' @examples
#' # Basic regression
#' bregress(mtcars, mpg ~ wt + hp)
#'
#' # With character formula
#' bregress(mtcars, "mpg ~ wt + hp")
#'
#' # With piping
#' library(dplyr)
#' mtcars %>% bregress(mpg ~ wt + hp)
#'
#' # With robust standard errors (if sandwich package is installed)
#' \dontrun{
#' bregress(mtcars, mpg ~ wt + hp, robust = TRUE)
#' }
bregress <- function(data, formula, robust = FALSE, alpha = 0.05,
                     vif = TRUE, verbose = FALSE) {

  # Handle matrix input - convert to data frame
  if (is.matrix(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    if (verbose) message("Converted matrix to data.frame")
  }

  # Parse and validate formula
  parsed <- bparse_formula(formula)
  bvalidate_formula(formula, data)
  formula <- parsed$formula

  # Fit the model
  model <- lm(formula, data = data, na.action = na.exclude)

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
    t_crit <- qt(1 - alpha/2, model$df.residual)
    conf_int <- cbind(
      coef_vals - t_crit * robust_se,
      coef_vals + t_crit * robust_se
    )
  } else {
    coef_summary <- summary(model)$coefficients
    conf_int <- confint(model, level = 1 - alpha)
  }

  # Extract basic information
  n <- nobs(model)
  k <- length(coef(model)) - 1
  df_model <- k
  df_resid <- model$df.residual
  df_total <- n - 1

  # Calculate sums of squares using model internals
  ss_resid <- sum(model$residuals^2, na.rm = TRUE)
  # Get the response variable from complete cases only
  y_vals <- model.response(model.frame(model))
  ss_total <- sum((y_vals - mean(y_vals))^2)
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

  # --- Print formatted output ---

  # Print header
  if (robust) {
    bprint_header("LINEAR REGRESSION", subtitle = "Robust Standard Errors (HC1)")
  } else {
    bprint_header("LINEAR REGRESSION")
  }

  # Format the SS, DF, and MS values first to determine needed widths
  ss_model_s <- formatC(ss_model, digits = 7, format = "g")
  ss_resid_s <- formatC(ss_resid, digits = 7, format = "g")
  ss_total_s <- formatC(ss_total, digits = 7, format = "g")

  ms_model_s <- formatC(ms_model, digits = 7, format = "g")
  ms_resid_s <- formatC(ms_resid, digits = 7, format = "g")
  ms_total_s <- formatC(ss_total/df_total, digits = 7, format = "g")

  df_model_s <- as.character(df_model)
  df_resid_s <- as.character(df_resid)
  df_total_s <- as.character(df_total)

  # Calculate dynamic column widths based on actual data
  w_src <- 12  # Keep Source column fixed
  w_ss <- max(nchar("SS"), nchar(ss_model_s), nchar(ss_resid_s), nchar(ss_total_s)) + 2
  w_df <- max(nchar("df"), nchar(df_model_s), nchar(df_resid_s), nchar(df_total_s)) + 2
  w_ms <- max(nchar("MS"), nchar(ms_model_s), nchar(ms_resid_s), nchar(ms_total_s)) + 2

  # Print ANOVA table with dynamic column widths
  header <- sprintf("%*s | %*s %*s %*s   Number of obs   = %9d",
                    w_src, "Source", w_ss, "SS", w_df, "df", w_ms, "MS", n)
  cat("\n", header, "\n", sep = "")
  cat(strrep("-", nchar(header)), "\n", sep = "")

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
                 strrep("-", w_ss + 1 + w_df + 1 + w_ms + 1))
  cat(sep, sprintf("   Adj R-squared   = %10.4f", adj_r_squared), "\n", sep = "")

  cat(row_fmt("Total",    ss_total,  df_total, ss_total/df_total,
              sprintf("   Root MSE        = %9.3f", root_mse)), "\n")

  # Print coefficient table
  cat("\n")

  # Determine dependent variable name
  dep_var <- parsed$response

  # Find the longest variable name to determine column width
  max_var_length <- max(nchar(coef_names), nchar(dep_var))
  var_col_width <- max(12, max_var_length)

  # Header for coefficient table with properly aligned columns
  ci_label <- sprintf("[%.0f%% conf. interval]", 100 * (1 - alpha))
  header_line <- sprintf("%*s | Coefficient  Std. err.      t    P>|t|     %s",
                         var_col_width, dep_var, ci_label)

  total_width <- nchar(header_line)
  separator_line <- strrep("-", total_width)
  cat(separator_line, "\n")
  cat(header_line, "\n", sep = "")

  # Separator line between header and data
  sep_line <- paste0(strrep("-", var_col_width), "+", strrep("-", total_width - var_col_width - 1))
  cat(sep_line, "\n")

  # Print each coefficient row
  for (i in 1:length(coef_names)) {
    var_name <- coef_names[i]

    coef_val <- coef_summary[i, 1]
    se_val <- coef_summary[i, 2]
    t_val <- coef_summary[i, 3]
    p_val <- coef_summary[i, 4]
    ci_low <- conf_int[i, 1]
    ci_high <- conf_int[i, 2]

    # Format p-value
    p_str <- bformat_pval(p_val, width = 7)

    # Get significance stars
    stars <- bsig_stars(p_val)

    cat(sprintf("%*s | %s  %s  %6.2f  %s %s  %s %s\n",
                var_col_width, var_name,
                bformat_num(coef_val, width = 10, digits = 5),
                bformat_num(se_val, width = 9, digits = 5),
                t_val,
                p_str,
                stars,
                bformat_num(ci_low, width = 10, digits = 5),
                bformat_num(ci_high, width = 10, digits = 5)))
  }
  cat(separator_line, "\n")
  bprint_sig_legend()

  # Calculate VIF and Tolerance for predictors (exclude intercept)
  if (vif && k > 1) {
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
    bprint_section("Variance Inflation Factors")

    # Create VIF table header with R-squared, Tolerance, VIF, and Sqrt(VIF)
    vif_header <- sprintf("%*s |    R-squared   Tolerance         VIF   Sqrt(VIF)",
                          var_col_width, "Variable")
    vif_table_width <- nchar(vif_header)
    vif_separator <- strrep("-", vif_table_width)

    cat(vif_separator, "\n")
    cat(vif_header, "\n")
    cat(paste0(strrep("-", var_col_width), "+", strrep("-", vif_table_width - var_col_width - 1)), "\n")

    # Calculate R-squared values for display
    r_squared_values <- 1 - tolerance_values

    # Print VIF values for each predictor
    for (j in 1:k) {
      var_name <- predictor_names[j]
      if (nchar(var_name) > var_col_width) {
        var_name <- substr(var_name, 1, var_col_width)
      }

      cat(sprintf("%*s | %10.4f  %10.4f  %10.2f  %10.4f\n",
                  var_col_width, var_name,
                  r_squared_values[j],
                  tolerance_values[j],
                  vif_values[j],
                  sqrt(vif_values[j])))
    }

    # Add mean VIF at the bottom
    mean_vif <- mean(vif_values)
    mean_rsq <- mean(r_squared_values)
    mean_tol <- mean(tolerance_values)
    mean_sqrt_vif <- sqrt(mean_vif)

    cat(paste0(strrep("-", var_col_width), "+", strrep("-", vif_table_width - var_col_width - 1)), "\n")
    cat(sprintf("%*s | %10.4f  %10.4f  %10.2f  %10.4f\n",
                var_col_width, "Mean",
                mean_rsq, mean_tol, mean_vif, mean_sqrt_vif))
    cat(vif_separator, "\n")

    # VIF interpretation guide
    if (max(vif_values) > 10) {
      cat("\nNote: VIF > 10 indicates potentially problematic multicollinearity\n")
    } else if (max(vif_values) > 5) {
      cat("\nNote: VIF > 5 may indicate moderate multicollinearity\n")
    }
  }

  cat("\n")

  # return calculated vars with NAs for cases with NA
  model$fitted.values <- fitted(model)
  model$residuals <- resid(model)
  model$effects <- effects(model)

  # Return the model invisibly
  invisible(model)
}
