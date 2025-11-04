#' Moderation analysis with Johnson-Neyman technique
#'
#' Performs moderation analysis including regression output and Johnson-Neyman
#' regions of significance for continuous moderators, or simple effects for
#' binary/factor moderators
#'
#' @param data A data frame containing the variables
#' @param Y Character string naming the outcome variable
#' @param X Character string naming the focal predictor variable
#' @param M Character string naming the moderator variable
#' @param covariates Optional character vector of covariate names to include
#' @param robust Logical, whether to use robust standard errors
#' @param alpha Significance level for confidence intervals (default = 0.05)
#' @param n_points Number of points to evaluate across moderator range (default = 50)
#' @param plot Logical, whether to create an interaction plot (default = TRUE)
#' @param plot_method Method for selecting moderator values in plot (default = "sd")
#' @param plot_custom_values Custom values for moderator if plot_method = "custom"
#'
#' @return Invisibly returns a list containing the model and conditional effects
#'
#' @importFrom stats lm coef vcov qt pt quantile sd
#' @export
#'
#' @examples
#' # Moderation with continuous moderator
#' bmoderation(mtcars, Y = "mpg", X = "wt", M = "hp")
bmoderation <- function(data, Y, X, M, covariates = NULL,
                        robust = FALSE, alpha = 0.05, n_points = 50,
                        plot = TRUE, plot_method = "sd", plot_custom_values = NULL) {

  # Check that required columns exist
  required_cols <- c(Y, X, M, covariates)
  if (!all(required_cols %in% names(data))) {
    missing <- required_cols[!required_cols %in% names(data)]
    stop("Variables not found in data: ", paste(missing, collapse = ", "))
  }

  # Determine if moderator is binary or factor
  mod_values <- data[[M]]
  unique_vals <- unique(mod_values[!is.na(mod_values)])
  n_unique <- length(unique_vals)
  is_binary <- n_unique == 2
  is_factor <- is.factor(mod_values)

  # Check if it's essentially binary (0/1 or similar)
  if (!is_factor && n_unique <= 5) {
    # Might be a categorical variable coded as numeric
    is_categorical <- TRUE
  } else {
    is_categorical <- FALSE
  }

  # Decide on analysis type
  use_jn <- !is_binary && !is_factor && !is_categorical

  # Create formula with interaction
  if (!is.null(covariates)) {
    cov_string <- paste("+", paste(covariates, collapse = " + "))
  } else {
    cov_string <- ""
  }

  formula_string <- paste(Y, "~", X, "*", M, cov_string)
  formula <- as.formula(formula_string)

  # Print header for moderation analysis
  cat("\n==============================================================================\n")
  cat("MODERATION ANALYSIS\n")
  cat("==============================================================================\n")
  cat("Outcome (Y):", Y, "\n")
  cat("Focal predictor (X):", X, "\n")
  cat("Moderator (M):", M, "\n")
  if (!is.null(covariates)) {
    cat("Covariates:", paste(covariates, collapse = ", "), "\n")
  }
  cat("\n")

  # Run regression using bregress for consistent output
  model <- bregress(data, formula, robust = robust)

  # Extract coefficients and variance-covariance matrix
  coefs <- coef(model)

  # Get variance-covariance matrix
  if (robust) {
    if (!requireNamespace("sandwich", quietly = TRUE)) {
      stop("Package 'sandwich' needed for robust standard errors. Please install it.",
           call. = FALSE)
    }
    vcov_matrix <- sandwich::vcovHC(model, type = "HC1")
  } else {
    vcov_matrix <- vcov(model)
  }

  # Find coefficient indices
  focal_index <- which(names(coefs) == X)
  interaction_name <- paste(X, M, sep = ":")
  if (!interaction_name %in% names(coefs)) {
    interaction_name <- paste(M, X, sep = ":")
  }
  interaction_index <- which(names(coefs) == interaction_name)

  if (length(focal_index) == 0 || length(interaction_index) == 0) {
    stop("Could not find focal or interaction coefficients in model")
  }

  # Get coefficients
  b_focal <- coefs[focal_index]
  b_interaction <- coefs[interaction_index]

  # Get variance components
  var_focal <- vcov_matrix[focal_index, focal_index]
  var_interaction <- vcov_matrix[interaction_index, interaction_index]
  cov_focal_interaction <- vcov_matrix[focal_index, interaction_index]

  df_resid <- model$df.residual
  t_crit <- qt(1 - alpha/2, df_resid)

  if (use_jn) {
    # CONTINUOUS MODERATOR: Johnson-Neyman Analysis

    # Calculate Johnson-Neyman points
    a <- b_interaction^2 - t_crit^2 * var_interaction
    b <- 2 * (b_focal * b_interaction - t_crit^2 * cov_focal_interaction)
    c <- b_focal^2 - t_crit^2 * var_focal

    discriminant <- b^2 - 4*a*c

    jn_points <- NULL
    if (discriminant >= 0 && a != 0) {
      jn_point1 <- (-b - sqrt(discriminant)) / (2*a)
      jn_point2 <- (-b + sqrt(discriminant)) / (2*a)
      jn_points <- sort(c(jn_point1, jn_point2))
    }

    # Get range of moderator
    mod_min <- min(mod_values, na.rm = TRUE)
    mod_max <- max(mod_values, na.rm = TRUE)
    mod_range <- mod_max - mod_min

    # Create sequence of moderator values
    eval_points <- seq(mod_min, mod_max, length.out = n_points)

    # Add J-N points if they're within or near the data range
    if (!is.null(jn_points)) {
      extended_min <- mod_min - 0.1 * mod_range
      extended_max <- mod_max + 0.1 * mod_range

      for (jn in jn_points) {
        if (jn >= extended_min && jn <= extended_max) {
          eval_points <- sort(unique(c(eval_points, jn)))
        }
      }
    }

    # Print Johnson-Neyman analysis header
    cat("\n==============================================================================\n")
    cat("JOHNSON-NEYMAN TECHNIQUE\n")
    cat("==============================================================================\n")

    if (!is.null(jn_points)) {
      cat("\nJohnson-Neyman significance transition points:\n")
      for (i in 1:length(jn_points)) {
        if (jn_points[i] >= mod_min && jn_points[i] <= mod_max) {
          cat(sprintf("  Point %d: %s = %.4f (within observed range)\n",
                      i, M, jn_points[i]))
        } else {
          cat(sprintf("  Point %d: %s = %.4f (outside observed range)\n",
                      i, M, jn_points[i]))
        }
      }
    } else {
      cat("\nNo Johnson-Neyman points found within reasonable range.\n")
      cat("The conditional effect may be significant across entire range or never significant.\n")
    }

    cat("\nObserved range of moderator:\n")
    cat(sprintf("  %s: [%.4f, %.4f]\n", M, mod_min, mod_max))

  } else {
    # BINARY/CATEGORICAL MODERATOR: Simple Effects

    # Sort unique values for consistent display
    eval_points <- sort(unique_vals)

    cat("\n==============================================================================\n")
    cat("CONDITIONAL EFFECTS\n")
    cat("==============================================================================\n")

    if (is_binary) {
      cat("\nModerator is binary with values:", paste(eval_points, collapse = ", "), "\n")
    } else if (is_factor) {
      cat("\nModerator is a factor with", n_unique, "levels\n")
    } else {
      cat("\nModerator appears categorical with", n_unique, "unique values:",
          paste(eval_points, collapse = ", "), "\n")
    }
  }

  # Calculate conditional effects at each evaluation point
  results <- data.frame(
    moderator = eval_points,
    effect = NA,
    se = NA,
    t = NA,
    p = NA,
    llci = NA,
    ulci = NA
  )

  for (i in 1:length(eval_points)) {
    W <- eval_points[i]

    # Conditional effect
    effect <- b_focal + b_interaction * W

    # Standard error of conditional effect
    se <- sqrt(var_focal + W^2 * var_interaction + 2 * W * cov_focal_interaction)

    # t-statistic and p-value
    t_stat <- effect / se
    p_value <- 2 * pt(abs(t_stat), df = df_resid, lower.tail = FALSE)

    # Confidence interval
    margin <- t_crit * se
    llci <- effect - margin
    ulci <- effect + margin

    results[i, ] <- c(W, effect, se, t_stat, p_value, llci, ulci)
  }

  # Select points to display
  if (use_jn) {
    # For continuous moderator, show subset of points
    display_indices <- round(seq(1, nrow(results), length.out = min(22, nrow(results))))
    display_results <- results[display_indices, ]
  } else {
    # For binary/categorical, show all levels
    display_results <- results
  }

  # Format and print the conditional effects table
  cat("\nConditional effect of", X, "at values of the moderator:\n")

  # Determine column widths based on variable name length
  mod_name_width <- max(8, nchar(M))

  # Print header
  header <- sprintf("%*s    effect        se         t         p      LLCI      ULCI",
                    mod_name_width, M)
  cat(header, "\n")

  # Print each row
  for (i in 1:nrow(display_results)) {
    row <- display_results[i, ]

    # Format p-value
    if (row$p < 0.0001) {
      p_str <- "   0.0001"
    } else {
      p_str <- sprintf("%9.4f", row$p)
    }

    # Mark J-N points with asterisk if they're in the display (only for continuous)
    marker <- ""
    if (use_jn && !is.null(jn_points)) {
      for (jn in jn_points) {
        if (abs(row$moderator - jn) < 0.0001) {
          marker <- "*"
          break
        }
      }
    }

    # For factor variables, show the level name
    if (is_factor) {
      mod_label <- sprintf("%.4f", row$moderator)
    } else {
      mod_label <- sprintf("%.4f", row$moderator)
    }

    cat(sprintf("%*s %8.4f %9.4f %9.4f %9s %9.4f %9.4f%s\n",
                mod_name_width,
                mod_label,
                row$effect,
                row$se,
                row$t,
                p_str,
                row$llci,
                row$ulci,
                marker))
  }

  if (use_jn && !is.null(jn_points) && any(jn_points >= min(eval_points) & jn_points <= max(eval_points))) {
    cat("\n* Johnson-Neyman significance transition point\n")
  }

  cat("\n")

  # Return results
  result_list <- list(
    model = model,
    coefficients = coefs,
    focal = X,
    moderator = M,
    outcome = Y,
    conditional_effects = results
  )

  if (use_jn) {
    result_list$johnson_neyman_points <- jn_points
  }

  # Create plot if requested
  if (plot) {
    # Determine best method for plotting based on moderator type
    if (!use_jn) {
      # Categorical moderator - will use all levels
      plot_method <- "categorical"
    } else if (plot_method == "sd") {
      # Check if SD method would go out of range
      mod_mean <- mean(data[[M]], na.rm = TRUE)
      mod_sd <- sd(data[[M]], na.rm = TRUE)
      mod_min <- min(data[[M]], na.rm = TRUE)
      mod_max <- max(data[[M]], na.rm = TRUE)

      if ((mod_mean - mod_sd) < mod_min || (mod_mean + mod_sd) > mod_max) {
        cat("\nNote: Defaulting to percentile method for plot as mean +/- 1 SD exceeds data range.\n")
        plot_method <- "percentile"
      }
    }

    # Call the plotting function
    binteraction_plot(result_list, method = plot_method,
                      custom_values = plot_custom_values, print_table = FALSE)
  }

  invisible(result_list)
}
