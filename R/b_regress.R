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

  # [Rest of function remains the same...]
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

  # Print ANOVA table
  cat("\n")
  cat("      Source |       SS           df       MS      Number of obs   =",
      sprintf("%9d", n), "\n")
  cat("-------------+----------------------------------   F(", df_model, ", ",
      df_resid, ")", sprintf("%7s", ""), "=",
      sprintf("%9.2f", f_stat), "\n", sep = "")
  cat("       Model |", sprintf("%11.7g", ss_model),
      sprintf("%10d", df_model),
      sprintf("%11.7g", ms_model),
      "   Prob > F        =",
      sprintf("%10.4f", f_pvalue), "\n")
  cat("    Residual |", sprintf("%11.7g", ss_resid),
      sprintf("%10d", df_resid),
      sprintf("%11.7g", ms_resid),
      "   R-squared       =",
      sprintf("%10.4f", r_squared), "\n")
  cat("-------------+----------------------------------   Adj R-squared   =",
      sprintf("%10.4f", adj_r_squared), "\n")
  cat("       Total |", sprintf("%11.7g", ss_total),
      sprintf("%10d", df_total),
      sprintf("%11.7g", ss_total/df_total),
      "   Root MSE        =",
      sprintf("%10g", root_mse), "\n")

  # Print coefficient table
  cat("\n")
  cat("------------------------------------------------------------------------------\n")

  # Determine dependent variable name
  dep_var <- as.character(formula[[2]])

  # Header for coefficient table
  cat(sprintf("%12s", dep_var), "| Coefficient  Std. err.      t    P>|t|     [95% conf. interval]\n")
  cat("-------------+----------------------------------------------------------------\n")

  # Print each coefficient row
  for (i in 1:length(coef_names)) {
    coef_val <- coef_summary[i, 1]
    se_val <- coef_summary[i, 2]
    t_val <- coef_summary[i, 3]
    p_val <- coef_summary[i, 4]
    ci_low <- conf_int[i, 1]
    ci_high <- conf_int[i, 2]

    cat(sprintf("%12s", substr(coef_names[i], 1, 12)), "|",
        sprintf("%11.7g", coef_val),
        sprintf("%11.7g", se_val),
        sprintf("%8.2f", t_val),
        sprintf("%8.3f", p_val),
        sprintf("%13.7g", ci_low),
        sprintf("%11.7g", ci_high),
        "\n")
  }
  cat("------------------------------------------------------------------------------\n")

  # Return the model invisibly
  invisible(model)
}
