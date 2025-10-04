# bregress

#' Stata-style regression output for linear models
#'
#' Produces regression output that mimics Stata's regress command format
#'
#' @param formula A formula object specifying the model
#' @param data A data frame containing the variables
#' @param robust Logical, whether to use robust standard errors (requires sandwich package)
#'
#' @return Invisibly returns the fitted lm model object
#'
#' @importFrom stats lm coef confint nobs resid pf pt qt
#' @export
#'
#' @examples
#' # Basic regression
#' bregress(mpg ~ wt + hp, data = mtcars)
#'
#' # With robust standard errors (if sandwich package is installed)
#' \dontrun{
#' bregress(mpg ~ wt + hp, data = mtcars, robust = TRUE)
#' }
bregress <- function(formula, data, robust = FALSE) {
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
  cat(paste0(rep("-", nchar(header)), collapse = ""), " \n", sep = "")

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

  sep <- sprintf("%*s+%s", w_src, "",
                 paste0(rep("-", w_ss + 1 + w_df + 1 + w_ms + 2), collapse = ""))
  cat(sep, sprintf("   Adj R-squared   = %10.4f", adj_r_squared), "\n", sep = "")

  cat(row_fmt("Total",    ss_total,  df_total, ss_total/df_total,
              sprintf("   Root MSE        = %10.3f", root_mse)), "\n")

  # Print coefficient table
  cat("\n")
  cat("------------------------------------------------------------------------------\n")

  # Determine dependent variable name and truncate if needed
  dep_var <- as.character(formula[[2]])
  if (nchar(dep_var) > 12) {
    dep_var <- substr(dep_var, 1, 12)
  }

  # Header for coefficient table - note the exact spacing
  cat(sprintf("%12s", dep_var), " | Coefficient  Std. err.      t    P>|t|     [95% conf. interval]\n", sep = "")
  cat("-------------+----------------------------------------------------------------\n")

  # Print each coefficient row with exact spacing to match Stata
  for (i in 1:length(coef_names)) {
    # Truncate variable name if too long
    var_name <- coef_names[i]
    if (nchar(var_name) > 12) {
      var_name <- substr(var_name, 1, 12)
    }

    coef_val <- coef_summary[i, 1]
    se_val <- coef_summary[i, 2]
    t_val <- coef_summary[i, 3]
    p_val <- coef_summary[i, 4]
    ci_low <- conf_int[i, 1]
    ci_high <- conf_int[i, 2]

    # Handle very small p-values
    if (p_val < 0.0005) {
      p_str <- "   0.000"
    } else {
      p_str <- sprintf("%8.3f", p_val)
    }

    # Print with exact spacing
    cat(sprintf("%12s", var_name), " | ",
        sprintf("%10.7g", coef_val), "  ",
        sprintf("%10.7g", se_val),
        sprintf("%8.2f", t_val),
        p_str, "    ",
        sprintf("%10.7g", ci_low), "  ",
        sprintf("%10.7g", ci_high),
        "\n", sep = "")
  }
  cat("------------------------------------------------------------------------------\n")

  # Return the model invisibly
  invisible(model)
}
