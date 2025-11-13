#' Align predictions with original data rows
#'
#' NOTE: This helper function is defined but not currently used.
#' The alignment is now done directly using kept_rows indices.
#'
#' @description
#' Internal utility function to align predictions from cleaned data back to original data,
#' inserting NAs where rows were removed during cleaning.
#'
#' @param original_data Original data frame before cleaning
#' @param cleaned_data Cleaned data frame after removing missing values
#' @param predictions Vector of predictions from model fit on cleaned data
#'
#' @return Vector of predictions with same length as original_data, with NAs in removed positions
#' @keywords internal
balign_predictions <- function(original_data, cleaned_data, predictions) {
  # Work with copies to avoid modifying the original data frames
  orig_copy <- original_data
  clean_copy <- cleaned_data

  # Create a key to match rows
  orig_copy$.row_id <- seq_len(nrow(orig_copy))
  clean_copy$.row_id <- as.numeric(rownames(clean_copy))

  # Create full-length result vector
  result <- rep(NA_real_, nrow(original_data))

  # Fill in predictions at correct positions
  result[clean_copy$.row_id] <- predictions

  return(result)
}

#' Logistic Regression with Interaction Analysis
#'
#' @description
#' Fits a logistic regression model with an interaction term and provides comprehensive output
#' including coefficients, odds ratios, Johnson-Neyman style conditional effects analysis,
#' predicted probabilities, and an interaction plot.
#'
#' @param data A data frame
#' @param Y Character string specifying the binary dependent variable name
#' @param X Character string specifying the focal predictor variable name
#' @param M Character string specifying the moderator variable name
#' @param covariates Optional character vector of covariate names to include
#' @param center Character string specifying centering option: "none" (default),
#'   "all" (center all predictors), or "non-binary" (center only non-binary predictors).
#' @param robust Logical; use robust standard errors (default FALSE)
#' @param alpha Significance level for confidence intervals (default 0.05)
#' @param n_points Number of points for Johnson-Neyman analysis (default 50)
#' @param plot Logical; produce interaction plot (default TRUE)
#' @param plot_method Method for selecting moderator values: "sd" (default), "quartile", or "custom"
#' @param plot_custom_values Optional numeric vector of specific moderator values for plotting
#' @param y Alias for Y (for convenience)
#' @param x Alias for X (for convenience)
#' @param m Alias for M (for convenience)
#' @param w Alias for M (for compatibility)
#'
#' @return A list containing:
#' \item{coefficients}{Data frame with coefficients, odds ratios, and statistics}
#' \item{pseudo_rsquare}{McFadden's pseudo R-squared}
#' \item{conditional_effects}{Data frame with Johnson-Neyman style analysis}
#' \item{predicted_probs}{Vector of predicted probabilities (aligned with original data)}
#' \item{model}{The fitted glm model object}
#' \item{plot}{ggplot object (if plot = TRUE)}
#'
#' @export
blogisticinteraction <- function(data, Y = NULL, X = NULL, M = NULL, covariates = NULL,
                                 center = "none",
                                 robust = FALSE, alpha = 0.05, n_points = 50,
                                 plot = TRUE, plot_method = "sd",
                                 plot_custom_values = NULL,
                                 y = NULL, x = NULL, m = NULL, w = NULL) {

  # Handle parameter aliases - lowercase versions take precedence if both are specified
  if (!is.null(y)) Y <- y
  if (!is.null(x)) X <- x
  if (!is.null(m)) M <- m
  if (!is.null(w) && is.null(M) && is.null(m)) M <- w  # w as alias for M

  # Check that required parameters are provided
  if (is.null(Y) || is.null(X) || is.null(M)) {
    stop("Must provide Y (or y), X (or x), and M (or m/w) parameters")
  }

  # Load required packages
  require(dplyr)
  require(ggplot2)

  # Store original data for alignment later
  original_data <- data

  # Select relevant variables - preserve rownames for alignment
  vars_to_keep <- c(Y, X, M, covariates)
  data_subset <- data[, vars_to_keep, drop = FALSE]

  # Drop missing values while preserving rownames
  complete_rows <- complete.cases(data_subset)
  data_subset <- data_subset[complete_rows, , drop = FALSE]

  # Store the row indices for later alignment (before any modifications)
  kept_rows <- as.numeric(rownames(data_subset))

  # Check that Y is binary
  if (length(unique(data_subset[[Y]])) != 2) {
    stop("Dependent variable must be binary")
  }

  # Store original uncentered data for plotting
  original_X <- data_subset[[X]]
  original_M <- data_subset[[M]]
  # Remove NAs for plotting calculations
  original_X_clean <- original_X[!is.na(original_X)]
  original_M_clean <- original_M[!is.na(original_M)]
  X_mean <- 0
  M_mean <- 0
  covariate_means <- list()

  # Apply centering if requested
  if (center == "all") {
    X_mean <- mean(data_subset[[X]])
    M_mean <- mean(data_subset[[M]])
    data_subset[[X]] <- scale(data_subset[[X]], center = TRUE, scale = FALSE)[,1]
    data_subset[[M]] <- scale(data_subset[[M]], center = TRUE, scale = FALSE)[,1]
    if (!is.null(covariates)) {
      for (cov in covariates) {
        if (is.numeric(data_subset[[cov]])) {
          covariate_means[[cov]] <- mean(data_subset[[cov]])
          data_subset[[cov]] <- scale(data_subset[[cov]], center = TRUE, scale = FALSE)[,1]
        } else {
          covariate_means[[cov]] <- 0
        }
      }
    }
  } else if (center == "non-binary") {
    # Center only non-binary variables
    if (length(unique(data_subset[[X]])) > 2) {
      X_mean <- mean(data_subset[[X]])
      data_subset[[X]] <- scale(data_subset[[X]], center = TRUE, scale = FALSE)[,1]
    }
    if (length(unique(data_subset[[M]])) > 2) {
      M_mean <- mean(data_subset[[M]])
      data_subset[[M]] <- scale(data_subset[[M]], center = TRUE, scale = FALSE)[,1]
    }
    if (!is.null(covariates)) {
      for (cov in covariates) {
        if (is.numeric(data_subset[[cov]]) && length(unique(data_subset[[cov]])) > 2) {
          covariate_means[[cov]] <- mean(data_subset[[cov]])
          data_subset[[cov]] <- scale(data_subset[[cov]], center = TRUE, scale = FALSE)[,1]
        } else {
          covariate_means[[cov]] <- 0
        }
      }
    }
  }

  # Create interaction term
  data_subset$XM_interaction <- data_subset[[X]] * data_subset[[M]]

  # Build formula
  if (is.null(covariates)) {
    formula_str <- paste0(Y, " ~ ", X, " + ", M, " + XM_interaction")
  } else {
    formula_str <- paste0(Y, " ~ ", X, " + ", M, " + XM_interaction + ",
                          paste(covariates, collapse = " + "))
  }

  # Fit model
  model <- glm(as.formula(formula_str), data = data_subset, family = binomial(link = "logit"))

  # Get model summary
  if (robust) {
    require(sandwich)
    require(lmtest)
    vcov_matrix <- vcovHC(model, type = "HC1")
    coef_test <- coeftest(model, vcov. = vcov_matrix)
  } else {
    coef_test <- coef(summary(model))
    vcov_matrix <- vcov(model)
  }

  # Create coefficient table with odds ratios
  coef_table <- data.frame(
    Variable = rownames(coef_test),
    Coef = coef_test[, "Estimate"],
    Odds_Ratio = exp(coef_test[, "Estimate"]),
    Std_Err = coef_test[, "Std. Error"],
    z = coef_test[, "z value"],
    P_value = coef_test[, "Pr(>|z|)"],
    CI_lower = coef_test[, "Estimate"] - qnorm(1 - alpha/2) * coef_test[, "Std. Error"],
    CI_upper = coef_test[, "Estimate"] + qnorm(1 - alpha/2) * coef_test[, "Std. Error"]
  )

  # Calculate pseudo R-squared (McFadden)
  null_model <- glm(as.formula(paste0(Y, " ~ 1")), data = data_subset, family = binomial(link = "logit"))
  pseudo_rsq <- 1 - (logLik(model)[1] / logLik(null_model)[1])

  # Johnson-Neyman style conditional effects analysis
  # Effect of X at different values of M: beta_X + beta_XM * M
  beta_X <- coef(model)[X]
  beta_XM <- coef(model)["XM_interaction"]

  M_values <- seq(min(data_subset[[M]]), max(data_subset[[M]]), length.out = n_points)

  conditional_effects <- data.frame(
    M_value = M_values,
    Conditional_Effect = beta_X + beta_XM * M_values
  )

  # Calculate standard errors for conditional effects
  # Var(beta_X + beta_XM * M) = Var(beta_X) + M^2 * Var(beta_XM) + 2 * M * Cov(beta_X, beta_XM)
  var_beta_X <- vcov_matrix[X, X]
  var_beta_XM <- vcov_matrix["XM_interaction", "XM_interaction"]
  cov_X_XM <- vcov_matrix[X, "XM_interaction"]

  conditional_effects$SE <- sqrt(var_beta_X + M_values^2 * var_beta_XM + 2 * M_values * cov_X_XM)
  conditional_effects$z <- conditional_effects$Conditional_Effect / conditional_effects$SE
  conditional_effects$P_value <- 2 * pnorm(-abs(conditional_effects$z))
  conditional_effects$CI_lower <- conditional_effects$Conditional_Effect - qnorm(1 - alpha/2) * conditional_effects$SE
  conditional_effects$CI_upper <- conditional_effects$Conditional_Effect + qnorm(1 - alpha/2) * conditional_effects$SE

  # Get predicted probabilities
  pred_probs_cleaned <- predict(model, type = "response")

  # Align predictions with original data using stored row indices
  predicted_probs <- rep(NA_real_, nrow(original_data))
  predicted_probs[kept_rows] <- pred_probs_cleaned

  # Initialize results object
  results <- list(
    coefficients = coef_table,
    pseudo_rsquare = pseudo_rsq,
    conditional_effects = conditional_effects,
    predicted_probs = predicted_probs,
    model = model
  )

  # Create interaction plot if requested
  plot_obj <- NULL
  if (plot) {
    # Determine M values for plotting (using original scale, cleaned of NAs)
    if (plot_method == "sd") {
      M_plot_mean <- mean(original_M_clean)
      M_plot_sd <- sd(original_M_clean)
      M_plot_values <- c(M_plot_mean - M_plot_sd, M_plot_mean, M_plot_mean + M_plot_sd)
      M_labels <- c("-1 SD", "Mean", "+1 SD")
    } else if (plot_method == "quartile") {
      M_plot_values <- quantile(original_M_clean, probs = c(0.25, 0.5, 0.75))
      M_labels <- c("25th %ile", "Median", "75th %ile")
    } else if (plot_method == "custom" && !is.null(plot_custom_values)) {
      M_plot_values <- plot_custom_values
      M_labels <- as.character(M_plot_values)
    } else {
      stop("Invalid plot_method or missing plot_custom_values")
    }

    # Create grid for plotting (using original scale, cleaned of NAs)
    X_range <- seq(min(original_X_clean), max(original_X_clean), length.out = 100)

    plot_data <- expand.grid(
      X_val = X_range,
      M_val = M_plot_values
    )

    # Set up prediction data frame (need to center for prediction)
    pred_df <- data.frame(
      X_val = plot_data$X_val - X_mean,  # Center for model
      M_val = plot_data$M_val - M_mean,  # Center for model
      XM_interaction = (plot_data$X_val - X_mean) * (plot_data$M_val - M_mean)
    )
    names(pred_df)[1:2] <- c(X, M)

    # Add covariates at their means (centered if needed)
    if (!is.null(covariates)) {
      for (cov in covariates) {
        if (cov %in% names(covariate_means)) {
          pred_df[[cov]] <- 0  # Already centered in model
        } else {
          pred_df[[cov]] <- mean(data_subset[[cov]], na.rm = TRUE)
        }
      }
    }

    # Get predictions
    pred_df$predicted_prob <- predict(model, newdata = pred_df, type = "response")
    pred_df$M_label <- factor(plot_data$M_val, levels = M_plot_values, labels = M_labels)
    pred_df$X_original <- plot_data$X_val  # Store original scale for plotting

    # Create plot (using original scale)
    plot_obj <- ggplot(pred_df, aes(x = X_original, y = predicted_prob, color = M_label)) +
      geom_line(linewidth = 1) +
      labs(x = X, y = "Predicted Probability", color = M,
           title = "Logistic Regression Interaction Plot") +
      theme_minimal() +
      theme(legend.position = "right")
  }

  # Print output using cat/sprintf for better wrapping behavior
  cat("\n")
  cat("Logistic Regression with Interaction\n")
  cat("=====================================\n")
  cat("Dependent variable:", Y, "\n")
  cat("Sample size:", nrow(data_subset), "\n")
  cat("Centering:", center, "\n")
  cat("Pseudo R-squared (McFadden):", round(pseudo_rsq, 4), "\n\n")

  cat("Coefficients:\n")

  # Determine column widths
  var_width <- max(nchar(as.character(coef_table$Variable)))

  # Print header
  header <- sprintf("%-*s       Coef  Odds_Ratio    Std_Err         z    P_value    CI_lower    CI_upper",
                    var_width, "Variable")
  cat(header, "\n")

  # Print separator
  cat(paste(rep("-", nchar(header)), collapse = ""), "\n")

  # Print each row
  for (i in 1:nrow(coef_table)) {
    row <- coef_table[i, ]

    # Format p-value
    if (is.na(row$P_value)) {
      p_str <- "       NA"
    } else if (row$P_value < 0.0001) {
      p_str <- "   <0.0001"
    } else {
      p_str <- sprintf("%10.4f", row$P_value)
    }

    cat(sprintf("%-*s %10.4f %11.4f %10.4f %9.3f %10s %11.4f %11.4f\n",
                var_width,
                row$Variable,
                row$Coef,
                row$Odds_Ratio,
                row$Std_Err,
                row$z,
                p_str,
                row$CI_lower,
                row$CI_upper))
  }

  cat("\n\nConditional Effect of", X, "at values of", M, ":\n")

  # Print conditional effects header
  M_width <- max(8, nchar(M))
  cond_header <- sprintf("%*s  Conditional_Effect        SE         z    P_value    CI_lower    CI_upper",
                         M_width, M)
  cat(cond_header, "\n")
  cat(paste(rep("-", nchar(cond_header)), collapse = ""), "\n")

  # Determine which points to display (subset for readability)
  n_display <- min(20, nrow(conditional_effects))
  display_indices <- round(seq(1, nrow(conditional_effects), length.out = n_display))

  for (i in display_indices) {
    row <- conditional_effects[i, ]

    # Format p-value
    if (is.na(row$P_value)) {
      p_str <- "       NA"
    } else if (row$P_value < 0.0001) {
      p_str <- "   <0.0001"
    } else {
      p_str <- sprintf("%10.4f", row$P_value)
    }

    cat(sprintf("%*.4f %19.4f %9.4f %9.3f %10s %11.4f %11.4f\n",
                M_width,
                row$M_value,
                row$Conditional_Effect,
                row$SE,
                row$z,
                p_str,
                row$CI_lower,
                row$CI_upper))
  }

  if (plot) {
    results$plot <- plot_obj
    print(plot_obj)
  }

  invisible(results)
}
