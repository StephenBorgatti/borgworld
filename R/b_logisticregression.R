#' Better Logistic Regression
#'
#' @description
#' Performs logistic regression analysis with automatic data cleaning,
#' missing value handling, and Stata-style output formatting.
#' Follows borgworld package conventions for user-friendly defaults.
#'
#' @param data A data frame containing the variables
#' @param formula A formula or character string specifying the model (e.g., "Y ~ X1 + X2")
#' @param robust Logical; use robust standard errors (default FALSE)
#' @param alpha Significance level for confidence intervals (default 0.05)
#' @param vif Logical; calculate Variance Inflation Factors (default FALSE)
#' @param classification_threshold Numeric; threshold for classification (default 0.5)
#' @param hosmer_lemeshow Logical; perform Hosmer-Lemeshow test (default FALSE)
#' @param roc Logical; calculate ROC curve and AUC (default FALSE)
#' @param print_output Logical; print formatted output (default TRUE)
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
#' @export
blogisticregression <- function(data, formula,
                                robust = FALSE,
                                alpha = 0.05,
                                vif = FALSE,
                                classification_threshold = 0.5,
                                hosmer_lemeshow = FALSE,
                                roc = FALSE,
                                print_output = TRUE) {

  # Load required packages
  require(dplyr)

  # Convert formula to character if needed
  if (inherits(formula, "formula")) {
    formula_str <- deparse(formula)
  } else {
    formula_str <- formula
  }

  # Parse formula to get variable names
  formula_obj <- as.formula(formula_str)
  all_vars <- all.vars(formula_obj)
  y_var <- all_vars[1]
  x_vars <- all_vars[-1]

  # Store original data for alignment
  original_data <- data
  n_original <- nrow(data)

  # Clean data using bclean() approach
  # Select only variables in the model
  data_subset <- data %>%
    select(all_of(all_vars))

  # Count missing before dropping
  n_missing_by_var <- sapply(data_subset, function(x) sum(is.na(x)))

  # Drop missing values
  data_clean <- data_subset %>%
    drop_na()

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
    message(paste0("Note: Converted '", y_var, "' to 0/1 (0 = ", y_sorted[1], ", 1 = ", y_sorted[2], ")"))
  }

  # Fit logistic regression model
  model <- glm(formula_obj, data = data_clean, family = binomial(link = "logit"))

  # Get coefficients and standard errors
  if (robust) {
    require(sandwich)
    require(lmtest)
    vcov_matrix <- vcovHC(model, type = "HC1")
    coef_test <- coeftest(model, vcov. = vcov_matrix)
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

  # Align predictions with original data using balign_predictions
  if (exists("balign_predictions", mode = "function")) {
    predicted_probs <- balign_predictions(original_data, data_clean, pred_probs_clean)
  } else {
    # Fallback alignment method
    predicted_probs <- rep(NA, n_original)
    # Find which rows in original data are in clean data
    if (n_missing > 0) {
      complete_cases <- complete.cases(original_data[, all_vars])
      predicted_probs[complete_cases] <- pred_probs_clean
    } else {
      predicted_probs <- pred_probs_clean
    }
  }

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
    require(car)
    # Only calculate VIF for numeric predictors
    numeric_vars <- x_vars[sapply(x_vars, function(v) is.numeric(data_clean[[v]]))]
    if (length(numeric_vars) > 1) {
      # Create a temporary linear model for VIF calculation
      temp_lm <- lm(formula_obj, data = data_clean)
      vif_values <- vif(temp_lm)
    }
  }

  # Perform Hosmer-Lemeshow test if requested
  hosmer_lemeshow_test <- NULL
  if (hosmer_lemeshow) {
    require(ResourceSelection)
    hl_test <- hoslem.test(data_clean[[y_var]], pred_probs_clean, g = 10)
    hosmer_lemeshow_test <- list(
      statistic = hl_test$statistic,
      p_value = hl_test$p.value,
      df = hl_test$parameter,
      groups = 10
    )
  }

  # Calculate ROC and AUC if requested
  auc_value <- NULL
  if (roc) {
    require(pROC)
    roc_obj <- roc(data_clean[[y_var]], pred_probs_clean, quiet = TRUE)
    auc_value <- auc(roc_obj)[1]
  }

  # Print formatted output
  if (print_output) {
    cat("\n")
    cat("Logistic Regression Results\n")
    cat("=====================================\n")
    cat("Dependent variable:  ", y_var, "\n")
    cat("Number of obs:       ", formatC(n_obs, format = "d", big.mark = ","), "\n")
    if (n_missing > 0) {
      cat("Obs dropped (missing):", formatC(n_missing, format = "d", big.mark = ","), "\n")
    }
    cat("Pseudo R-squared:    ", sprintf("%.4f", pseudo_rsq), "\n")
    cat("AIC:                 ", sprintf("%.2f", aic_value), "\n")
    cat("BIC:                 ", sprintf("%.2f", bic_value), "\n")
    cat("Log Likelihood:      ", sprintf("%.2f", logLik(model)[1]), "\n")
    if (robust) {
      cat("Standard errors:      Robust (HC1)\n")
    }
    cat("\n")

    # Format coefficient table for printing (matching borgworld style)
    cat("Coefficients and Odds Ratios:\n")
    cat("-------------------------------------\n")

    # Add significance stars directly to p-values
    p_stars <- ifelse(coef_table$P_value < 0.001, "***",
                      ifelse(coef_table$P_value < 0.01, "**",
                             ifelse(coef_table$P_value < 0.05, "*",
                                    ifelse(coef_table$P_value < 0.1, ".", ""))))

    # Determine variable name width
    var_width <- max(nchar(as.character(coef_table$Variable)))

    # Create header
    header_line <- sprintf("%-*s %10s %10s %10s %8s %10s",
                           var_width,
                           "Variable",
                           "Coef",
                           "OR",
                           "Std.Err",
                           "z",
                           "P>|z|")
    cat(header_line, "\n")
    cat(paste(rep("-", nchar(header_line)), collapse = ""), "\n")

    # Print each coefficient row
    for (i in 1:nrow(coef_table)) {
      cat(sprintf("%-*s %10.4f %10.4f %10.4f %8.3f %10.4f %s\n",
                  var_width,
                  coef_table$Variable[i],
                  coef_table$Coef[i],
                  coef_table$Odds_Ratio[i],
                  coef_table$Std_Err[i],
                  coef_table$z[i],
                  coef_table$P_value[i],
                  p_stars[i]))
    }

    cat("\n95% Confidence Intervals:\n")

    # CI header
    ci_header <- sprintf("%-*s %20s %20s",
                         var_width,
                         "Variable",
                         "CI (Coef)",
                         "CI (OR)")
    cat(ci_header, "\n")
    cat(paste(rep("-", nchar(ci_header)), collapse = ""), "\n")

    # Print confidence intervals
    for (i in 1:nrow(coef_table)) {
      cat(sprintf("%-*s [%8.4f, %8.4f] [%8.4f, %8.4f]\n",
                  var_width,
                  coef_table$Variable[i],
                  coef_table$CI_lower[i],
                  coef_table$CI_upper[i],
                  coef_table$OR_CI_lower[i],
                  coef_table$OR_CI_upper[i]))
    }

    cat("-------------------------------------\n")
    cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1\n")

    # Print classification table if available
    if (!is.null(classification_table)) {
      cat("\n")
      cat("Classification Table (threshold = ", classification_threshold, "):\n", sep = "")
      cat("-------------------------------------\n")
      print(classification_table$confusion_matrix)
      cat("\n")
      cat(sprintf("Sensitivity (Recall):  %.4f\n", classification_table$sensitivity))
      cat(sprintf("Specificity:           %.4f\n", classification_table$specificity))
      cat(sprintf("Accuracy:              %.4f\n", classification_table$accuracy))
      cat("\n")
      cat("Sensitivity: of all actual positives, how many were detected?\n")
      cat("Specificity: of all actual negatives, how many were detected?")
    }

    # Print VIF if calculated
    if (!is.null(vif_values)) {
      cat("\n")
      cat("Variance Inflation Factors:\n")
      cat("-------------------------------------\n")

      # Determine variable width
      vif_var_width <- max(nchar(names(vif_values)))

      # Header
      vif_header <- sprintf("%-*s %10s", vif_var_width, "Variable", "VIF")
      cat(vif_header, "\n")
      cat(paste(rep("-", nchar(vif_header)), collapse = ""), "\n")

      # Print each VIF
      for (i in 1:length(vif_values)) {
        cat(sprintf("%-*s %10.3f\n",
                    vif_var_width,
                    names(vif_values)[i],
                    vif_values[i]))
      }
    }

    # Print Hosmer-Lemeshow test if calculated
    if (!is.null(hosmer_lemeshow_test)) {
      cat("\n")
      cat("Hosmer-Lemeshow Goodness of Fit Test:\n")
      cat("-------------------------------------\n")
      cat(sprintf("Chi-squared: %.4f\n", hosmer_lemeshow_test$statistic))
      cat(sprintf("df:          %d\n", hosmer_lemeshow_test$df))
      cat(sprintf("P-value:     %.4f\n", hosmer_lemeshow_test$p_value))
      if (hosmer_lemeshow_test$p_value > 0.05) {
        cat("(Good fit: p > 0.05 suggests no evidence of poor fit)\n")
      } else {
        cat("(Poor fit: p < 0.05 suggests evidence of poor fit)\n")
      }
    }

    # Print AUC if calculated
    if (!is.null(auc_value)) {
      cat("\n")
      cat("ROC Analysis:\n")
      cat("-------------------------------------\n")
      cat(sprintf("Area Under Curve (AUC): %.4f\n", auc_value))
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
    formula = formula_str,
    robust = robust
  )

  invisible(results)
}

# Helper function for aligned predictions if not already available
balign_predictions <- function(original_data, cleaned_data, predictions) {
  # Get the row indices from cleaned data
  n_original <- nrow(original_data)
  n_cleaned <- nrow(cleaned_data)

  # Initialize result vector with NAs
  aligned_predictions <- rep(NA, n_original)

  # If no rows were dropped, return predictions as-is
  if (n_original == n_cleaned) {
    return(predictions)
  }

  # Find which rows from original are in cleaned (complete cases)
  # This assumes row order is preserved after dropping NAs
  row_indices <- which(complete.cases(original_data[, names(cleaned_data)]))

  # Assign predictions to appropriate positions
  if (length(row_indices) == length(predictions)) {
    aligned_predictions[row_indices] <- predictions
  }

  return(aligned_predictions)
}
