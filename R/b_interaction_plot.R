#' Create interaction plots for moderation analysis
#'
#' Plots interaction effects between a focal predictor and moderator on an outcome.
#' Can be used after bmoderation or independently.
#'
#' @param data Data frame containing the variables (or bmoderation result object)
#' @param focal Character string naming the focal predictor (or NULL if data is bmoderation result)
#' @param moderator Character string naming the moderator (or NULL if data is bmoderation result)
#' @param outcome Character string naming the outcome (or NULL if data is bmoderation result)
#' @param model Optional: fitted lm model (extracted from bmoderation result if provided)
#' @param method For continuous moderators: "sd" (+/-1SD), "percentile" (16/50/84),
#'               "minmax", "quartile" (25/50/75), or "custom"
#' @param custom_values If method="custom", numeric vector of moderator values
#' @param plot_points Logical, whether to show data points on plot (default FALSE)
#' @param print_table Logical, whether to print the summary table (default TRUE)
#' @param y_range_method "data" for full data range, "predicted" for predicted values range (default)
#' @param ... Additional arguments passed to plot
#'
#' @return Invisibly returns a data frame with slope information
#'
#' @importFrom stats predict lm quantile sd
#' @importFrom graphics plot lines legend points mtext par
#' @importFrom grDevices rainbow
#' @export
binteraction_plot <- function(data, focal = NULL, moderator = NULL, outcome = NULL,
                              model = NULL, method = "sd", custom_values = NULL,
                              plot_points = FALSE, print_table = TRUE,
                              y_range_method = "predicted", ...) {

  # Check if data is a bmoderation result
  is_bmod_result <- is.list(data) && !is.null(data$model) && !is.null(data$focal)

  if (is_bmod_result) {
    # Extract from bmoderation result
    model <- data$model
    focal <- data$focal
    moderator <- data$moderator
    df <- model$model
    outcome <- names(df)[1]  # First column in model frame is outcome
  } else {
    # Use provided data frame
    df <- data

    # Check required arguments
    if (is.null(focal) || is.null(moderator) || is.null(outcome)) {
      stop("When using raw data, focal, moderator, and outcome must be specified")
    }

    # Fit model if not provided
    if (is.null(model)) {
      formula <- as.formula(paste(outcome, "~", focal, "*", moderator))
      model <- lm(formula, data = df)
    } else {
      # Extract data from provided model
      df <- model$model
    }
  }

  # Get moderator values and determine type
  mod_values <- df[[moderator]]
  unique_vals <- unique(mod_values[!is.na(mod_values)])
  n_unique <- length(unique_vals)

  # Determine if categorical
  is_categorical <- n_unique <= 5 || is.factor(mod_values)

  # Select moderator values for plotting
  if (is_categorical) {
    # Use all unique values for categorical
    plot_mod_values <- sort(unique_vals)
    method_label <- "Categorical levels"
  } else {
    # Continuous moderator - use specified method
    mod_mean <- mean(mod_values, na.rm = TRUE)
    mod_sd <- sd(mod_values, na.rm = TRUE)
    mod_min <- min(mod_values, na.rm = TRUE)
    mod_max <- max(mod_values, na.rm = TRUE)

    if (method == "sd") {
      plot_mod_values <- c(mod_mean - mod_sd, mod_mean, mod_mean + mod_sd)
      method_label <- "Mean +/- 1 SD"

      # Check if values exceed range
      if (plot_mod_values[1] < mod_min || plot_mod_values[3] > mod_max) {
        warning(paste0("Warning: mean +/- 1 SD values exceed observed range of ", moderator,
                       " [", sprintf("%.3f", mod_min), ", ", sprintf("%.3f", mod_max), "]"))
      }

    } else if (method == "percentile") {
      plot_mod_values <- quantile(mod_values, c(0.16, 0.50, 0.84), na.rm = TRUE)
      method_label <- "16th, 50th, 84th percentiles"

    } else if (method == "minmax") {
      plot_mod_values <- c(mod_min, mod_mean, mod_max)
      method_label <- "Min, Mean, Max"

    } else if (method == "quartile") {
      plot_mod_values <- quantile(mod_values, c(0.25, 0.50, 0.75), na.rm = TRUE)
      method_label <- "25th, 50th, 75th percentiles"

    } else if (method == "custom") {
      if (is.null(custom_values)) {
        stop("custom_values must be provided when method='custom'")
      }
      plot_mod_values <- sort(custom_values)
      method_label <- "Custom values"

    } else {
      stop("Invalid method. Choose 'sd', 'percentile', 'minmax', 'quartile', or 'custom'")
    }
  }

  # Get focal predictor range for plotting
  focal_values <- df[[focal]]
  focal_range <- range(focal_values, na.rm = TRUE)
  focal_seq <- seq(focal_range[1], focal_range[2], length.out = 100)

  # Extract coefficients
  coefs <- coef(model)
  b_intercept <- coefs["(Intercept)"]
  b_focal <- coefs[focal]

  # Find interaction term
  interaction_name <- paste(focal, moderator, sep = ":")
  if (!interaction_name %in% names(coefs)) {
    interaction_name <- paste(moderator, focal, sep = ":")
  }
  b_interaction <- coefs[interaction_name]
  b_moderator <- coefs[moderator]

  # Handle other covariates (set to mean for prediction)
  other_vars <- setdiff(names(df), c(outcome, focal, moderator))
  pred_data <- data.frame(matrix(ncol = length(other_vars), nrow = length(focal_seq)))
  names(pred_data) <- other_vars
  for (var in other_vars) {
    if (is.numeric(df[[var]])) {
      pred_data[[var]] <- mean(df[[var]], na.rm = TRUE)
    } else {
      # Use mode for categorical
      pred_data[[var]] <- names(sort(table(df[[var]]), decreasing = TRUE))[1]
    }
  }

  # Store slope information and predicted values
  slope_info <- data.frame(
    moderator_value = plot_mod_values,
    slope = NA,
    intercept = NA,
    z_score = NA,
    percentile = NA,
    label = NA,
    stringsAsFactors = FALSE
  )

  # Store all predicted values to determine y-range
  all_predictions <- list()

  # Calculate slopes and predictions
  for (i in 1:length(plot_mod_values)) {
    mod_val <- plot_mod_values[i]

    # Calculate simple slope
    slope <- b_focal + b_interaction * mod_val
    intercept <- b_intercept + b_moderator * mod_val

    # Create prediction data
    temp_pred_data <- pred_data
    temp_pred_data[[focal]] <- focal_seq
    temp_pred_data[[moderator]] <- mod_val

    # Predict values and store them
    y_pred <- predict(model, newdata = temp_pred_data)
    all_predictions[[i]] <- y_pred

    # Calculate statistics for this moderator value
    z_score <- (mod_val - mean(mod_values, na.rm = TRUE)) / sd(mod_values, na.rm = TRUE)
    percentile <- mean(mod_values <= mod_val, na.rm = TRUE) * 100

    # Create labels
    if (is_categorical) {
      label <- paste0(moderator, " = ", sprintf("%.1f", mod_val))
    } else {
      label <- paste0(moderator, " = ", sprintf("%.2f", mod_val))
    }

    # Store info
    slope_info[i, ] <- list(mod_val, slope, intercept, z_score, percentile, label)
  }

  # Determine y-range based on method
  if (y_range_method == "predicted") {
    # Use range of predicted values with buffer
    all_pred_vals <- unlist(all_predictions)
    y_min <- min(all_pred_vals)
    y_max <- max(all_pred_vals)
    y_buffer <- (y_max - y_min) * 0.1  # 10% buffer
    y_range <- c(y_min - y_buffer, y_max + y_buffer)
  } else {
    # Use full data range
    y_range <- range(df[[outcome]], na.rm = TRUE)
  }

  # Set up color palette
  if (length(plot_mod_values) <= 3) {
    colors <- c("blue", "darkgreen", "red")
  } else {
    colors <- rainbow(length(plot_mod_values))
  }

  # Save current par and adjust bottom margin for legend
  old_par <- par(no.readonly = TRUE)

  # Increase bottom margin to accommodate legend
  par(mar = c(6.5, 4, 4, 2) + 0.1)

  # Create plot
  plot(focal_range, y_range, type = "n",
       xlab = focal, ylab = outcome,
       main = paste("Interaction:", focal, "×", moderator),
       ...)

  # Add method label as subtitle
  mtext(paste("Method:", method_label), side = 3, line = 0.5, cex = 0.8)

  # Add data points if requested
  if (plot_points) {
    # Only show points within the y_range
    point_mask <- df[[outcome]] >= y_range[1] & df[[outcome]] <= y_range[2]
    if (any(point_mask)) {
      points(df[[focal]][point_mask], df[[outcome]][point_mask],
             pch = 16, col = "gray80", cex = 0.5)
    }
  }

  # Plot lines
  for (i in 1:length(plot_mod_values)) {
    lines(focal_seq, all_predictions[[i]], col = colors[i], lwd = 2)
  }

  # Add legend below the plot
  par(xpd = TRUE)  # Allow drawing outside plot region

  # Calculate legend position - below the x-axis label
  usr <- par("usr")
  legend_y <- usr[3] - (usr[4] - usr[3]) * 0.15  # Position below plot

  # Create legend text with line type indicators
  legend_text <- slope_info$label

  # Place legend horizontally centered below plot
  legend(x = mean(usr[1:2]), y = legend_y,
         legend = legend_text,
         col = colors[1:length(plot_mod_values)],
         lwd = 2,
         horiz = TRUE,  # Horizontal layout
         bty = "n",     # No box
         cex = 0.8,
         xjust = 0.5)   # Center horizontally

  # Reset xpd
  par(xpd = FALSE)

  # Restore original par settings
  suppressWarnings(par(old_par))

  # Print summary table
  if (print_table) {
    cat("\n==============================================================================\n")
    cat("SIMPLE SLOPES SUMMARY\n")
    cat("==============================================================================\n")
    cat("Focal predictor:", focal, "\n")
    cat("Moderator:", moderator, "\n")
    cat("Method:", method_label, "\n")
    cat("Y-axis range:", y_range_method, "\n\n")

    # Format table
    cat(sprintf("%-20s %12s %12s %12s\n",
                "Moderator Value", "Simple Slope", "Z-score", "Percentile"))
    cat(paste(rep("-", 68), collapse = ""), "\n")

    for (i in 1:nrow(slope_info)) {
      # Format moderator value based on type
      if (is_categorical) {
        mod_str <- sprintf("%-20.2f", slope_info$moderator_value[i])
      } else {
        mod_str <- sprintf("%-20.4f", slope_info$moderator_value[i])
      }

      # Format percentile
      pct_str <- sprintf("%.1f%%", slope_info$percentile[i])

      cat(sprintf("%s %12.4f %12.3f %12s\n",
                  mod_str,
                  slope_info$slope[i],
                  slope_info$z_score[i],
                  pct_str))
    }

    cat(paste(rep("-", 68), collapse = ""), "\n")

    # Add interpretation note
    cat("\nNote: Z-scores indicate how many standard deviations each moderator value is from the mean.\n")
    cat("      Percentiles show the proportion of observations at or below each moderator value.\n")

    if (!is_categorical && method == "sd") {
      # Check and report if values are within range
      if (any(plot_mod_values < mod_min | plot_mod_values > mod_max)) {
        cat("\nWarning: Some plotted values fall outside the observed data range.\n")
        cat(sprintf("         Observed range of %s: [%.3f, %.3f]\n",
                    moderator, mod_min, mod_max))
      }
    }

    # Report y-axis range used
    if (y_range_method == "predicted") {
      cat(sprintf("\nY-axis scaled to predicted values range: [%.2f, %.2f]\n",
                  y_range[1], y_range[2]))
    }
  }

  invisible(slope_info)
}
