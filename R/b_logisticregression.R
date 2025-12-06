#' Better Logistic Regression
#'
#' @description
#' Performs logistic regression analysis with automatic data cleaning,
#' missing value handling, and Stata-style output formatting.
#' Accepts both data frames and matrices as input.
#'
#' @param data A data frame or matrix containing the variables
#' @param formula A formula or character string specifying the model (e.g., "Y ~ X1 + X2")
#' @param robust Logical; use robust standard errors (default FALSE)
#' @param alpha Significance level for confidence intervals (default 0.05)
#' @param vif Logical; calculate Variance Inflation Factors (default FALSE)
#' @param classification_threshold Numeric; threshold for classification (default 0.5)
#' @param hosmer_lemeshow Logical; perform Hosmer-Lemeshow test (default FALSE)
#' @param roc Logical; calculate ROC curve and AUC (default FALSE)
#' @param print_output Logical; print formatted output (default TRUE)
#' @param verbose Logical; print messages about data cleaning (default FALSE)
#'
#' @return A list containing:
#' \item{model}{The fitted glm model object}
#' \item{coefficients}{Data frame with coefficients, odds ratios, and statistics}
#' \item{pseudo_rsquare}{McFadden's pseudo R-squared}
#' \item{aic}{Akaike Information Criterion}
#' \item{bic}{Bayesian Information Criterion}
#' \item{classification_table}{Classification table (if requested)}
#' \item{vif}{Variance Inflation Factors (if requested)}
#' \item{hosmer_lemeshow}{Hosmer-Lemeshow test results (if requested)}
#' \item{auc}{Area Under the ROC Curve (if requested)}
#' \item{predicted_probs}{Vector of predicted probabilities (aligned with original data)}
#' \item{n_obs}{Number of observations used}
#' \item{n_missing}{Number of observations dropped due to missing values}
#'
#' @importFrom stats glm binomial coef logLik AIC BIC qnorm predict complete.cases as.formula
#' @export
blogisticregression <- function(data, formula,
                                robust = FALSE,
                                alpha = 0.05,
                                vif = FALSE,
                                classification_threshold = 0.5,
                                hosmer_lemeshow = FALSE,
                                roc = FALSE,
                                print_output = TRUE,
                                verbose = FALSE) {

  # Handle matrix input - convert to data frame
  if (is.matrix(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    if (verbose) message("Converted matrix to data.frame")
  }

  # Parse and validate formula
  parsed <- bparse_formula(formula)
  bvalidate_formula(formula, data)
  formula_obj <- parsed$formula
  y_var <- parsed$response
  x_vars <- parsed$predictors

  # Store original data for alignment
  original_data <- data
  n_original <- nrow(data)

  # Select only variables in the model
  data_subset <- data[, parsed$all, drop = FALSE]

  # Count missing before dropping
  n_missing_by_var <- sapply(data_subset, function(x) sum(is.na(x)))

  # Drop missing values
  complete_rows <- complete.cases(data_subset)
  data_clean <- data_subset[complete_rows, , drop = FALSE]

  n_obs <- nrow(data_clean)
  n_missing <- n_original - n_obs

  # Check that Y is binary (or can be made binary)
  y_values <- unique(data_clean[[y_var]])
  if (length(y_values) > 2) {
    stop("Dependent variable must be binary (have exactly 2 unique values)")
  }

  # Convert Y to 0/1 if needed
  if (!all(y_values %in% c(0, 1))) {
    # Convert to 0/1 based on sorted unique values
    y_sorted <- sort(y_values)
    data_clean[[y_var]] <- as.numeric(data_clean[[y_var]] == y_sorted[2])
    message(sprintf("Note: Converted '%s' to 0/1 (0 = %s, 1 = %s)",
                    y_var, y_sorted[1], y_sorted[2]))
  }

  # Fit logistic regression model
  model <- glm(formula_obj, data = data_clean, family = binomial(link = "logit"))

  # Get coefficients and standard errors
  if (robust) {
    if (!requireNamespace("sandwich", quietly = TRUE)) {
      stop("Package 'sandwich' needed for robust standard errors. Please install it.")
    }
    if (!requireNamespace("lmtest", quietly = TRUE)) {
      stop("Package 'lmtest' needed for robust standard errors. Please install it.")
    }
    vcov_matrix <- sandwich::vcovHC(model, type = "HC1")
    coef_test <- lmtest::coeftest(model, vcov. = vcov_matrix)
    se_vector <- coef_test[, "Std. Error"]
    z_vector <- coef_test[, "z value"]
    p_vector <- coef_test[, "Pr(>|z|)"]
  } else {
    model_summary <- summary(model)
    coef_test <- coef(model_summary)
    se_vector <- coef_test[, "Std. Error"]
    z_vector <- coef_test[, "z value"]
    p_vector <- coef_test[, "Pr(>|z|)"]
  }

  # Create coefficient table with odds ratios
  coef_values <- coef(model)
  ci_lower <- coef_values - qnorm(1 - alpha/2) * se_vector
  ci_upper <- coef_values + qnorm(1 - alpha/2) * se_vector

  coef_table <- data.frame(
    Variable = names(coef_values),
    Coef = coef_values,
    Odds_Ratio = exp(coef_values),
    Std_Err = se_vector,
    z = z_vector,
    P_value = p_vector,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    OR_CI_lower = exp(ci_lower),
    OR_CI_upper = exp(ci_upper),
    stringsAsFactors = FALSE
  )
  rownames(coef_table) <- NULL

  # Calculate pseudo R-squared (McFadden)
  null_model <- glm(as.formula(paste0(y_var, " ~ 1")),
                    data = data_clean,
                    family = binomial(link = "logit"))
  pseudo_rsq <- 1 - (logLik(model)[1] / logLik(null_model)[1])

  # Calculate AIC and BIC
  aic_value <- AIC(model)
  bic_value <- BIC(model)

  # Get predicted probabilities (for cleaned data)
  pred_probs_clean <- predict(model, type = "response")

  # Align predictions with original data
  predicted_probs <- rep(NA, n_original)
  predicted_probs[complete_rows] <- pred_probs_clean

  # Create classification table if requested
  classification_table <- NULL
  if (!is.null(classification_threshold)) {
    pred_class <- ifelse(pred_probs_clean > classification_threshold, 1, 0)
    actual_class <- data_clean[[y_var]]

    # Create confusion matrix
    conf_matrix <- table(Actual = actual_class, Predicted = pred_class)

    # Calculate metrics
    if (all(dim(conf_matrix) == c(2, 2))) {
      tn <- conf_matrix[1, 1]
      fp <- conf_matrix[1, 2]
      fn <- conf_matrix[2, 1]
      tp <- conf_matrix[2, 2]

      sensitivity <- tp / (tp + fn)
      specificity <- tn / (tn + fp)
      accuracy <- (tp + tn) / sum(conf_matrix)

      classification_table <- list(
        confusion_matrix = conf_matrix,
        sensitivity = sensitivity,
        specificity = specificity,
        accuracy = accuracy,
        threshold = classification_threshold
      )
    }
  }

  # Calculate VIF if requested
  vif_values <- NULL
  if (vif && length(x_vars) > 1) {
    if (!requireNamespace("car", quietly = TRUE)) {
      warning("Package 'car' needed for VIF calculation. Skipping.")
    } else {
      # Only calculate VIF for numeric predictors
      numeric_vars <- x_vars[sapply(x_vars, function(v) is.numeric(data_clean[[v]]))]
      if (length(numeric_vars) > 1) {
        # Create a temporary linear model for VIF calculation
        temp_lm <- lm(formula_obj, data = data_clean)
        vif_values <- car::vif(temp_lm)
      }
    }
  }

  # Perform Hosmer-Lemeshow test if requested
  hosmer_lemeshow_test <- NULL
  if (hosmer_lemeshow) {
    if (!requireNamespace("ResourceSelection", quietly = TRUE)) {
      warning("Package 'ResourceSelection' needed for Hosmer-Lemeshow test. Skipping.")
    } else {
      hl_test <- ResourceSelection::hoslem.test(data_clean[[y_var]], pred_probs_clean, g = 10)
      hosmer_lemeshow_test <- list(
        statistic = hl_test$statistic,
        p_value = hl_test$p.value,
        df = hl_test$parameter,
        groups = 10
      )
    }
  }

  # Calculate ROC and AUC if requested
  auc_value <- NULL
  if (roc) {
    if (!requireNamespace("pROC", quietly = TRUE)) {
      warning("Package 'pROC' needed for ROC analysis. Skipping.")
    } else {
      roc_obj <- pROC::roc(data_clean[[y_var]], pred_probs_clean, quiet = TRUE)
      auc_value <- pROC::auc(roc_obj)[1]
    }
  }

  # Print formatted output
  if (print_output) {
    if (robust) {
      bprint_header("LOGISTIC REGRESSION", subtitle = "Robust Standard Errors (HC1)")
    } else {
      bprint_header("LOGISTIC REGRESSION")
    }

    bprint_info(
      "Dependent variable" = y_var,
      "Number of obs" = n_obs
    )
    if (n_missing > 0) {
      bprint_info("Obs dropped (missing)" = n_missing)
    }
    bprint_info(
      "Pseudo R-squared" = pseudo_rsq,
      "AIC" = aic_value,
      "BIC" = bic_value,
      "Log Likelihood" = logLik(model)[1]
    )

    # Print coefficient table
    bprint_section("Coefficients and Odds Ratios")

    var_width <- max(nchar(as.character(coef_table$Variable)))

    # Header
    cat(sprintf("%-*s %10s %10s %10s %8s %10s\n",
                var_width, "Variable", "Coef", "OR", "Std.Err", "z", "P>|z|"))
    cat(strrep("-", var_width + 55), "\n")

    # Print each coefficient row
    for (i in 1:nrow(coef_table)) {
      stars <- bsig_stars(coef_table$P_value[i])
      cat(sprintf("%-*s %s %s %s %8.3f %s %s\n",
                  var_width, coef_table$Variable[i],
                  bformat_num(coef_table$Coef[i], width = 10),
                  bformat_num(coef_table$Odds_Ratio[i], width = 10),
                  bformat_num(coef_table$Std_Err[i], width = 10),
                  coef_table$z[i],
                  bformat_pval(coef_table$P_value[i]),
                  stars))
    }

    # Print confidence intervals
    ci_pct <- sprintf("%.0f%%", 100 * (1 - alpha))
    cat(sprintf("\n%s Confidence Intervals:\n", ci_pct))
    cat(sprintf("%-*s %20s %20s\n",
                var_width, "Variable", "CI (Coef)", "CI (OR)"))
    cat(strrep("-", var_width + 45), "\n")

    for (i in 1:nrow(coef_table)) {
      cat(sprintf("%-*s [%8.4f, %8.4f] [%8.4f, %8.4f]\n",
                  var_width, coef_table$Variable[i],
                  coef_table$CI_lower[i], coef_table$CI_upper[i],
                  coef_table$OR_CI_lower[i], coef_table$OR_CI_upper[i]))
    }

    cat(strrep("-", var_width + 45), "\n")
    bprint_sig_legend()

    # Print classification table if available
    if (!is.null(classification_table)) {
      bprint_section(sprintf("Classification Table (threshold = %.2f)",
                             classification_threshold))
      print(classification_table$confusion_matrix)
      cat("\n")
      bprint_info(
        "Sensitivity (Recall)" = classification_table$sensitivity,
        "Specificity" = classification_table$specificity,
        "Accuracy" = classification_table$accuracy
      )
      cat("\nSensitivity: of all actual positives, how many were detected?\n")
      cat("Specificity: of all actual negatives, how many were detected?\n")
    }

    # Print VIF if calculated
    if (!is.null(vif_values)) {
      bprint_section("Variance Inflation Factors")
      vif_df <- data.frame(Variable = names(vif_values), VIF = vif_values)
      bprint_table(vif_df)
    }

    # Print Hosmer-Lemeshow test if calculated
    if (!is.null(hosmer_lemeshow_test)) {
      bprint_section("Hosmer-Lemeshow Goodness of Fit Test")
      bprint_info(
        "Chi-squared" = hosmer_lemeshow_test$statistic,
        "df" = hosmer_lemeshow_test$df,
        "P-value" = hosmer_lemeshow_test$p_value
      )
      if (hosmer_lemeshow_test$p_value > 0.05) {
        cat("(Good fit: p > 0.05 suggests no evidence of poor fit)\n")
      } else {
        cat("(Poor fit: p < 0.05 suggests evidence of poor fit)\n")
      }
    }

    # Print AUC if calculated
    if (!is.null(auc_value)) {
      bprint_section("ROC Analysis")
      bprint_info("Area Under Curve (AUC)" = auc_value)
      if (auc_value >= 0.9) {
        cat("(Excellent discrimination)\n")
      } else if (auc_value >= 0.8) {
        cat("(Good discrimination)\n")
      } else if (auc_value >= 0.7) {
        cat("(Acceptable discrimination)\n")
      } else if (auc_value >= 0.6) {
        cat("(Poor discrimination)\n")
      } else {
        cat("(No discrimination)\n")
      }
    }

    cat("\n")
  }

  # Prepare return object
  results <- list(
    model = model,
    coefficients = coef_table,
    pseudo_rsquare = pseudo_rsq,
    aic = aic_value,
    bic = bic_value,
    log_likelihood = logLik(model)[1],
    classification_table = classification_table,
    vif = vif_values,
    hosmer_lemeshow = hosmer_lemeshow_test,
    auc = auc_value,
    predicted_probs = predicted_probs,
    n_obs = n_obs,
    n_missing = n_missing,
    formula = deparse(formula_obj),
    robust = robust
  )

  invisible(results)
}
