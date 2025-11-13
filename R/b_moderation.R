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
#' @param center Character string specifying centering option: "none" (default),
#'   "all" (center all predictors), or "non-binary" (center only non-binary predictors).
#'   Note: Centering affects coefficient interpretation but not the location of
#'   Johnson-Neyman points or conditional effects when viewed in original scale.
#' @param display_scale Character string specifying output display: "original" (default)
#'   shows all moderator values in original units, "centered" shows values in
#'   centered scale when centering is applied (matching PROCESS macro behavior)
#' @param robust Logical, whether to use robust standard errors
#' @param alpha Significance level for confidence intervals (default = 0.05)
#' @param n_points Number of points to evaluate across moderator range (default = 50)
#' @param plot Logical, whether to create an interaction plot (default = TRUE)
#' @param plot_method Method for selecting moderator values in plot (default = "sd")
#' @param plot_custom_values Custom values for moderator if plot_method = "custom"
#' @param y Alias for Y (for convenience)
#' @param x Alias for X (for convenience)
#' @param m Alias for M (for convenience)
#' @param w Alias for M (for compatibility)
#'
#' @return Invisibly returns a list containing the model and conditional effects
#'
#' @importFrom stats lm coef vcov qt pt quantile sd
#' @export
#'
#' @examples
#' # Moderation with continuous moderator
#' bmoderation(mtcars, Y = "mpg", X = "wt", M = "hp")
#' # Also works with lowercase:
#' bmoderation(mtcars, y = "mpg", x = "wt", m = "hp")
#'
#' # With centering of non-binary variables (original scale display)
#' bmoderation(mtcars, Y = "mpg", X = "wt", M = "hp", center = "non-binary")
#'
#' # With centering and centered scale display (like PROCESS macro)
#' bmoderation(mtcars, Y = "mpg", X = "wt", M = "hp",
#'             center = "non-binary", display_scale = "centered")
bmoderation <- function(data, Y = NULL, X = NULL, M = NULL, covariates = NULL,
                        center = "none", display_scale = "original",
                        robust = FALSE, alpha = 0.05, n_points = 50,
                        plot = TRUE, plot_method = "sd", plot_custom_values = NULL,
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

  # Check centering argument
  if (!center %in% c("none", "all", "non-binary")) {
    stop("center must be 'none', 'all', or 'non-binary'")
  }

  # Check display_scale argument
  if (!display_scale %in% c("original", "centered")) {
    stop("display_scale must be 'original' or 'centered'")
  }

  # If no centering is applied, display_scale must be original
  if (center == "none" && display_scale == "centered") {
    warning("display_scale='centered' has no effect when center='none'. Using 'original'.")
    display_scale <- "original"
  }

  # Check that required columns exist
  required_cols <- c(Y, X, M, covariates)
  if (!all(required_cols %in% names(data))) {
    missing <- required_cols[!required_cols %in% names(data)]
    stop("Variables not found in data: ", paste(missing, collapse = ", "))
  }

  # Create a working copy of the data to avoid modifying the original
  data_work <- data

  # Helper function to check if a variable is binary
  is_binary_var <- function(x) {
    unique_vals <- unique(x[!is.na(x)])
    return(length(unique_vals) == 2)
  }

  # Store original means for reporting and back-transformation if needed
  original_means <- list()
  centered_vars <- character()

  # Apply centering based on the option chosen
  if (center != "none") {
    # List of variables to potentially center (not including Y)
    vars_to_consider <- c(X, M, covariates)

    for (var in vars_to_consider) {
      should_center <- FALSE

      if (center == "all") {
        should_center <- TRUE
      } else if (center == "non-binary") {
        should_center <- !is_binary_var(data[[var]])
      }

      if (should_center) {
        original_means[[var]] <- mean(data[[var]], na.rm = TRUE)
        data_work[[var]] <- data[[var]] - original_means[[var]]
        centered_vars <- c(centered_vars, var)
      }
    }
  }

  # Determine if moderator is binary or factor (using original data for this check)
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

  # Print centering information
  cat("\nCentering option: ", center, "\n")
  if (length(centered_vars) > 0) {
    cat("Centered variables:", paste(centered_vars, collapse = ", "), "\n")
    cat("Display scale:", display_scale, "\n")
  } else {
    cat("No variables were centered.\n")
  }

  cat("\n")

  # Run regression using lm() (use centered data)
  model <- lm(formula, data = data_work)

  # Print regression results
  cat("REGRESSION RESULTS\n")
  cat("------------------------------------------------------------------------------\n")
  print(summary(model))
  cat("\n")

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

    # Note: If centering was applied, the J-N calculations are performed on the
    # centered scale, but the points will be converted back to the original scale
    # for reporting. The location of significance transitions does not change with
    # centering - only the numerical representation changes.

    # Calculate Johnson-Neyman points (in centered scale if centering was applied)
    a <- b_interaction^2 - t_crit^2 * var_interaction
    b <- 2 * (b_focal * b_interaction - t_crit^2 * cov_focal_interaction)
    c <- b_focal^2 - t_crit^2 * var_focal

    discriminant <- b^2 - 4*a*c

    jn_points <- NULL
    if (!is.na(discriminant) && !is.na(a) && discriminant >= 0 && a != 0) {
      jn_point1 <- (-b - sqrt(discriminant)) / (2*a)
      jn_point2 <- (-b + sqrt(discriminant)) / (2*a)
      jn_points <- sort(c(jn_point1, jn_point2))
    }

    # Get range of moderator (use original data for determining range)
    mod_min <- min(mod_values, na.rm = TRUE)
    mod_max <- max(mod_values, na.rm = TRUE)
    mod_range <- mod_max - mod_min

    # Create sequence of moderator values
    eval_points_original <- seq(mod_min, mod_max, length.out = n_points)

    # Determine display points based on display_scale preference
    if (M %in% centered_vars) {
      eval_points_centered <- eval_points_original - original_means[[M]]

      # Add J-N points if they're within or near the data range
      if (!is.null(jn_points)) {
        extended_min_centered <- (mod_min - original_means[[M]]) - 0.1 * mod_range
        extended_max_centered <- (mod_max - original_means[[M]]) + 0.1 * mod_range

        for (jn in jn_points) {
          if (jn >= extended_min_centered && jn <= extended_max_centered) {
            eval_points_centered <- sort(unique(c(eval_points_centered, jn)))
            eval_points_original <- sort(unique(c(eval_points_original, jn + original_means[[M]])))
          }
        }
      }

      # Set evaluation and display points based on display_scale
      if (display_scale == "centered") {
        eval_points <- eval_points_centered
        display_points <- eval_points_centered
      } else {
        eval_points <- eval_points_centered  # for calculations
        display_points <- eval_points_original  # for display
      }
    } else {
      eval_points <- eval_points_original
      display_points <- eval_points_original

      # Add J-N points if they're within or near the data range
      if (!is.null(jn_points)) {
        extended_min <- mod_min - 0.1 * mod_range
        extended_max <- mod_max + 0.1 * mod_range

        for (jn in jn_points) {
          if (jn >= extended_min && jn <= extended_max) {
            eval_points <- sort(unique(c(eval_points, jn)))
            display_points <- eval_points
          }
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
        # Display J-N points based on display_scale preference
        if (M %in% centered_vars && display_scale == "centered") {
          jn_display <- jn_points[i]
          # Check if within centered range
          centered_min <- mod_min - original_means[[M]]
          centered_max <- mod_max - original_means[[M]]
          if (jn_display >= centered_min && jn_display <= centered_max) {
            cat(sprintf("  Point %d: %s = %.4f (within observed range)\n",
                        i, M, jn_display))
          } else {
            cat(sprintf("  Point %d: %s = %.4f (outside observed range)\n",
                        i, M, jn_display))
          }
        } else {
          # Display in original scale
          if (M %in% centered_vars) {
            jn_display <- jn_points[i] + original_means[[M]]
          } else {
            jn_display <- jn_points[i]
          }

          if (jn_display >= mod_min && jn_display <= mod_max) {
            cat(sprintf("  Point %d: %s = %.4f (within observed range)\n",
                        i, M, jn_display))
          } else {
            cat(sprintf("  Point %d: %s = %.4f (outside observed range)\n",
                        i, M, jn_display))
          }
        }
      }
    } else {
      cat("\nNo Johnson-Neyman points found within reasonable range.\n")
      cat("The conditional effect may be significant across entire range or never significant.\n")
    }

    cat("\nObserved range of moderator:\n")
    if (M %in% centered_vars && display_scale == "centered") {
      cat(sprintf("  %s: [%.4f, %.4f]\n", M,
                  mod_min - original_means[[M]],
                  mod_max - original_means[[M]]))
    } else {
      cat(sprintf("  %s: [%.4f, %.4f]\n", M, mod_min, mod_max))
    }

  } else {
    # BINARY/CATEGORICAL MODERATOR: Simple Effects

    # Sort unique values for consistent display
    eval_points_original <- sort(unique_vals)

    # If moderator was centered, convert to centered scale
    if (M %in% centered_vars) {
      eval_points_centered <- eval_points_original - original_means[[M]]

      # Set evaluation and display points based on display_scale
      if (display_scale == "centered") {
        eval_points <- eval_points_centered
        display_points <- eval_points_centered
      } else {
        eval_points <- eval_points_centered  # for calculations
        display_points <- eval_points_original  # for display
      }
    } else {
      eval_points <- eval_points_original
      display_points <- eval_points_original
    }

    cat("\n==============================================================================\n")
    cat("CONDITIONAL EFFECTS\n")
    cat("==============================================================================\n")

    if (is_binary) {
      if (M %in% centered_vars && display_scale == "centered") {
        cat("\nModerator is binary with values:",
            paste(sprintf("%.4f", display_points), collapse = ", "),
            "(centered)\n")
      } else {
        cat("\nModerator is binary with values:",
            paste(display_points, collapse = ", "), "\n")
      }
    } else if (is_factor) {
      cat("\nModerator is a factor with", n_unique, "levels\n")
    } else {
      if (M %in% centered_vars && display_scale == "centered") {
        cat("\nModerator appears categorical with", n_unique, "unique values:",
            paste(sprintf("%.4f", display_points), collapse = ", "),
            "(centered)\n")
      } else {
        cat("\nModerator appears categorical with", n_unique, "unique values:",
            paste(display_points, collapse = ", "), "\n")
      }
    }
  }

  # Calculate conditional effects at each evaluation point
  results <- data.frame(
    moderator = display_points,  # Store in display scale
    effect = NA,
    se = NA,
    t = NA,
    p = NA,
    llci = NA,
    ulci = NA
  )

  for (i in 1:length(eval_points)) {
    W <- eval_points[i]  # Use centered/calculation value

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

    results[i, ] <- c(display_points[i], effect, se, t_stat, p_value, llci, ulci)
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
    if (is.na(row$p)) {
      p_str <- "       NA"
    } else if (row$p < 0.0001) {
      p_str <- "   0.0001"
    } else {
      p_str <- sprintf("%9.4f", row$p)
    }

    # Mark J-N points with asterisk if they're in the display (only for continuous)
    marker <- ""
    if (use_jn && !is.null(jn_points)) {
      # Check against J-N points in the display scale
      if (M %in% centered_vars && display_scale == "original") {
        jn_display_vals <- jn_points + original_means[[M]]
      } else {
        jn_display_vals <- jn_points
      }

      for (jn_disp in jn_display_vals) {
        if (abs(row$moderator - jn_disp) < 0.0001) {
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
    conditional_effects = results,
    centering = center,
    centered_vars = centered_vars,
    original_means = original_means,
    display_scale = display_scale
  )

  if (use_jn) {
    # Store J-N points in the display scale
    if (M %in% centered_vars && display_scale == "original") {
      result_list$johnson_neyman_points <- jn_points + original_means[[M]]
    } else {
      result_list$johnson_neyman_points <- jn_points
    }
  }

  # Create plot if requested
  if (plot) {
    # Note: The plotting function receives moderator values in the display scale
    # (stored in result_list$conditional_effects$moderator). If display_scale="centered"
    # and centering was applied, plots will show centered values (matching PROCESS).
    # If display_scale="original", plots show original units for interpretability.

    # Determine best method for plotting based on moderator type
    if (!use_jn) {
      # Categorical moderator - will use all levels
      plot_method <- "categorical"
    } else if (plot_method == "sd") {
      # Check if SD method would go out of range (use original data)
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
