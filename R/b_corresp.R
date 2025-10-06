#' bcorresp: Perform Correspondence Analysis on a Rectangular Table
#'
#' This function performs correspondence analysis on a contingency table and
#' generates various visualization options including biplots. This version
#' has minimal dependencies and implements CA from scratch.
#'
#' @param data A data frame or matrix containing the contingency table
#' @param ncp Number of dimensions to keep (default = min(nrow-1, ncol-1))
#' @param graph Logical, whether to produce graphs (default = TRUE)
#' @param graph.type Type of graph: "biplot", "row", "col", "scree", or "all" (default = "biplot")
#' @param axes Which axes to plot (default = c(1,2))
#' @param title Main title for the plot (optional)
#' @param point.size Size of points in the plot (default = 3)
#' @param text.size Size of text labels (default = 4)
#' @param row.color Color for row points (default = "blue")
#' @param col.color Color for column points (default = "red")
#' @param save.plot Logical, whether to save the plot (default = FALSE)
#' @param file.name Name of file to save plot (default = "bcorresp_plot.pdf")
#' @param width Width of saved plot in inches (default = 10)
#' @param height Height of saved plot in inches (default = 8)
#' @param verbose Logical, whether to print results summary (default = TRUE)
#' @param use.ggplot2 Logical, whether to use ggplot2 if available (default = TRUE)
#'
#' @return A list of class "bcorresp" containing:
#'   \item{row.coords}{Coordinates of rows in factorial space}
#'   \item{col.coords}{Coordinates of columns in factorial space}
#'   \item{eigenvalues}{Eigenvalues}
#'   \item{variance.explained}{Percentage of variance explained by each dimension}
#'   \item{cumulative.variance}{Cumulative percentage of variance}
#'   \item{row.masses}{Row masses (marginal proportions)}
#'   \item{col.masses}{Column masses (marginal proportions)}
#'   \item{row.inertias}{Contribution of each row to total inertia}
#'   \item{col.inertias}{Contribution of each column to total inertia}
#'   \item{row.contrib}{Contribution of rows to dimensions (percentage)}
#'   \item{col.contrib}{Contribution of columns to dimensions (percentage)}
#'   \item{row.cos2}{Quality of representation for rows (cos2)}
#'   \item{col.cos2}{Quality of representation for columns (cos2)}
#'   \item{chi2}{Chi-square statistic and p-value}
#'   \item{total.inertia}{Total inertia}
#'   \item{data}{Original data}
#'   \item{P}{Correspondence matrix}
#'   \item{plots}{List of generated plots (if graph = TRUE)}
#'
#' @importFrom stats pchisq
#' @importFrom graphics abline barplot legend lines plot points text
#' @importFrom grDevices dev.off pdf
#' @importFrom utils write.csv
#'
#' @examples
#' # Example 1: Basic correspondence analysis
#' data <- matrix(c(10, 5, 15, 20,
#'                  25, 30, 10, 5,
#'                  15, 20, 25, 30),
#'                nrow = 3, byrow = TRUE)
#' rownames(data) <- c("Category1", "Category2", "Category3")
#' colnames(data) <- c("TypeA", "TypeB", "TypeC", "TypeD")
#'
#' result <- bcorresp(data)
#'
#' # Example 2: Without ggplot2 (uses base R plotting)
#' result <- bcorresp(data, use.ggplot2 = FALSE)
#'
#' # Example 3: All plot types
#' result <- bcorresp(data, graph.type = "all")
#'
#' @export
bcorresp <- function(data,
                     ncp = NULL,
                     graph = TRUE,
                     graph.type = "biplot",
                     axes = c(1, 2),
                     title = NULL,
                     point.size = 3,
                     text.size = 4,
                     row.color = "blue",
                     col.color = "red",
                     save.plot = FALSE,
                     file.name = "bcorresp_plot.pdf",
                     width = 10,
                     height = 8,
                     verbose = TRUE,
                     use.ggplot2 = TRUE) {

  # Check if ggplot2 is available
  has_ggplot2 <- use.ggplot2 && requireNamespace("ggplot2", quietly = TRUE)
  if (use.ggplot2 && !has_ggplot2) {
    message("ggplot2 not available. Using base R plotting.")
  }

  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Data must be a matrix or data frame")
  }

  # Convert to matrix if necessary
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  if (mode(data) != "numeric") {
    stop("Data matrix is not numeric. Remove non-numeric columns.")
  }

  # Check for negative values
  if (any(data < 0, na.rm = TRUE)) {
    stop("Data cannot contain negative values for correspondence analysis")
  }

  # Remove NA values
  if (any(is.na(data))) {
    warning("NA values found and will be replaced with 0")
    data[is.na(data)] <- 0
  }

  # Check for all-zero rows or columns
  zero_rows <- which(rowSums(data) == 0)
  zero_cols <- which(colSums(data) == 0)

  if (length(zero_rows) > 0 || length(zero_cols) > 0) {
    warning("Removing rows/columns with all zeros")
    if (length(zero_rows) > 0) data <- data[-zero_rows, , drop = FALSE]
    if (length(zero_cols) > 0) data <- data[, -zero_cols, drop = FALSE]
  }

  # Get dimensions
  n_rows <- nrow(data)
  n_cols <- ncol(data)

  # Set default ncp if not specified
  if (is.null(ncp)) {
    ncp <- min(n_rows - 1, n_cols - 1)
  } else {
    ncp <- min(ncp, n_rows - 1, n_cols - 1)
  }

  # Ensure we have row and column names
  if (is.null(rownames(data))) {
    rownames(data) <- paste0("Row", 1:n_rows)
  }
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("Col", 1:n_cols)
  }

  # STEP 1: Calculate the correspondence matrix
  N <- sum(data)  # Grand total
  P <- data / N   # Correspondence matrix

  # STEP 2: Calculate marginal proportions
  row_masses <- rowSums(P)
  col_masses <- colSums(P)

  # STEP 3: Calculate expected proportions under independence
  E <- outer(row_masses, col_masses)

  # STEP 4: Calculate residuals
  R <- (P - E) / sqrt(E)

  # Handle any NaN values (from division by zero)
  R[is.nan(R)] <- 0

  # STEP 5: Weighted SVD
  # Weight matrices
  Dr <- diag(1 / sqrt(row_masses))
  Dc <- diag(1 / sqrt(col_masses))

  # Weighted residual matrix
  S <- Dr %*% (P - E) %*% Dc

  # Perform SVD
  svd_result <- svd(S)

  # STEP 6: Extract results
  # Eigenvalues (squared singular values)
  eigenvalues <- svd_result$d^2

  # Only keep the requested number of dimensions
  n_dim <- min(ncp, length(eigenvalues))
  eigenvalues <- eigenvalues[1:n_dim]

  # Calculate percentage of variance explained
  total_inertia <- sum(eigenvalues)
  variance_explained <- 100 * eigenvalues / total_inertia
  cumulative_variance <- cumsum(variance_explained)

  # STEP 7: Calculate principal coordinates
  # Row coordinates (standard coordinates)
  row_coords <- Dr %*% svd_result$u[, 1:n_dim, drop = FALSE]
  rownames(row_coords) <- rownames(data)
  colnames(row_coords) <- paste0("Dim", 1:n_dim)

  # Column coordinates (standard coordinates)
  col_coords <- Dc %*% svd_result$v[, 1:n_dim, drop = FALSE]
  rownames(col_coords) <- colnames(data)
  colnames(col_coords) <- paste0("Dim", 1:n_dim)

  # Scale coordinates by singular values (principal coordinates)
  for (i in 1:n_dim) {
    row_coords[, i] <- row_coords[, i] * svd_result$d[i]
    col_coords[, i] <- col_coords[, i] * svd_result$d[i]
  }

  # STEP 8: Calculate contributions
  # Row contributions to dimensions
  row_contrib <- matrix(0, n_rows, n_dim)
  for (i in 1:n_dim) {
    row_contrib[, i] <- 100 * (row_masses * row_coords[, i]^2) / eigenvalues[i]
  }
  rownames(row_contrib) <- rownames(data)
  colnames(row_contrib) <- paste0("Dim", 1:n_dim)

  # Column contributions to dimensions
  col_contrib <- matrix(0, n_cols, n_dim)
  for (i in 1:n_dim) {
    col_contrib[, i] <- 100 * (col_masses * col_coords[, i]^2) / eigenvalues[i]
  }
  rownames(col_contrib) <- colnames(data)
  colnames(col_contrib) <- paste0("Dim", 1:n_dim)

  # STEP 9: Calculate quality of representation (cos2)
  # Row cos2
  row_dist2 <- rowSums((P / row_masses - rep(col_masses, each = n_rows))^2 *
                         matrix(rep(1/col_masses, n_rows), n_rows, byrow = TRUE))
  row_cos2 <- row_coords^2 / row_dist2
  row_cos2[is.nan(row_cos2)] <- 0

  # Column cos2
  col_dist2 <- colSums((t(P) / col_masses - rep(row_masses, each = n_cols))^2 *
                         matrix(rep(1/row_masses, n_cols), n_cols, byrow = TRUE))
  col_cos2 <- col_coords^2 / col_dist2
  col_cos2[is.nan(col_cos2)] <- 0

  # Calculate row and column inertias
  row_inertias <- row_masses * rowSums((P / row_masses - rep(col_masses, each = n_rows))^2 *
                                         matrix(rep(1/col_masses, n_rows), n_rows, byrow = TRUE))
  col_inertias <- col_masses * colSums((t(P) / col_masses - rep(row_masses, each = n_cols))^2 *
                                         matrix(rep(1/row_masses, n_cols), n_cols, byrow = TRUE))

  # Chi-square test
  chi2_stat <- N * total_inertia
  df <- (n_rows - 1) * (n_cols - 1)
  p_value <- stats::pchisq(chi2_stat, df, lower.tail = FALSE)

  # Print summary if verbose
  if (verbose) {
    cat("\n=== bcorresp Analysis Results ===\n\n")
    cat("Data dimensions:", n_rows, "rows x", n_cols, "columns\n")
    cat("Total inertia:", round(total_inertia, 4), "\n")
    cat("Chi-square statistic:", round(chi2_stat, 2),
        "(df =", df, ", p-value =", format.pval(p_value), ")\n\n")

    cat("Variance explained by dimensions:\n")
    var_table <- data.frame(
      Dimension = 1:min(5, n_dim),
      Eigenvalue = round(eigenvalues[1:min(5, n_dim)], 4),
      `Variance (%)` = round(variance_explained[1:min(5, n_dim)], 2),
      `Cumulative (%)` = round(cumulative_variance[1:min(5, n_dim)], 2),
      check.names = FALSE
    )
    print(var_table)
    cat("\n")
  }

  # Create result object
  results <- list(
    row.coords = row_coords,
    col.coords = col_coords,
    eigenvalues = eigenvalues,
    variance.explained = variance_explained,
    cumulative.variance = cumulative_variance,
    row.masses = row_masses,
    col.masses = col_masses,
    row.inertias = row_inertias,
    col.inertias = col_inertias,
    row.contrib = row_contrib,
    col.contrib = col_contrib,
    row.cos2 = row_cos2,
    col.cos2 = col_cos2,
    chi2 = list(statistic = chi2_stat, df = df, p.value = p_value),
    total.inertia = total_inertia,
    data = data,
    P = P,
    n.dim = n_dim
  )

  # Create plots if requested
  if (graph) {
    plots <- list()

    # Prepare data for plotting
    row_plot_data <- as.data.frame(row_coords[, axes, drop = FALSE])
    row_plot_data$label <- rownames(row_coords)
    row_plot_data$type <- "Row"

    col_plot_data <- as.data.frame(col_coords[, axes, drop = FALSE])
    col_plot_data$label <- rownames(col_coords)
    col_plot_data$type <- "Column"

    # Axis labels
    x_label <- paste0("Dimension ", axes[1], " (",
                      round(variance_explained[axes[1]], 1), "%)")
    y_label <- paste0("Dimension ", axes[2], " (",
                      round(variance_explained[axes[2]], 1), "%)")

    if (has_ggplot2) {
      # Use ggplot2 for plotting - use :: notation to avoid NOTE

      # Biplot
      if (graph.type %in% c("biplot", "all")) {
        biplot_data <- rbind(row_plot_data, col_plot_data)
        names(biplot_data)[1:2] <- c("Dim1", "Dim2")

        # Avoid global variable NOTEs
        Dim1 <- Dim2 <- type <- label <- NULL

        p_biplot <- ggplot2::ggplot(biplot_data, ggplot2::aes(x = Dim1, y = Dim2)) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
          ggplot2::geom_point(ggplot2::aes(color = type), size = point.size) +
          ggplot2::geom_text(ggplot2::aes(label = label, color = type),
                             size = text.size, vjust = -0.5) +
          ggplot2::scale_color_manual(values = c(Row = row.color, Column = col.color)) +
          ggplot2::labs(x = x_label, y = y_label,
                        title = ifelse(is.null(title),
                                       "Correspondence Analysis Biplot",
                                       title)) +
          ggplot2::theme_minimal() +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
                         legend.title = ggplot2::element_blank(),
                         legend.position = "bottom")

        plots$biplot <- p_biplot
        if (graph.type == "biplot") print(p_biplot)
      }

      # Row plot
      if (graph.type %in% c("row", "all")) {
        names(row_plot_data)[1:2] <- c("Dim1", "Dim2")

        Dim1 <- Dim2 <- label <- NULL

        p_row <- ggplot2::ggplot(row_plot_data, ggplot2::aes(x = Dim1, y = Dim2)) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
          ggplot2::geom_point(color = row.color, size = point.size) +
          ggplot2::geom_text(ggplot2::aes(label = label), color = row.color,
                             size = text.size, vjust = -0.5) +
          ggplot2::labs(x = x_label, y = y_label,
                        title = "Row Points") +
          ggplot2::theme_minimal() +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))

        plots$row <- p_row
        if (graph.type == "row") print(p_row)
      }

      # Column plot
      if (graph.type %in% c("col", "all")) {
        names(col_plot_data)[1:2] <- c("Dim1", "Dim2")

        Dim1 <- Dim2 <- label <- NULL

        p_col <- ggplot2::ggplot(col_plot_data, ggplot2::aes(x = Dim1, y = Dim2)) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
          ggplot2::geom_point(color = col.color, size = point.size) +
          ggplot2::geom_text(ggplot2::aes(label = label), color = col.color,
                             size = text.size, vjust = -0.5) +
          ggplot2::labs(x = x_label, y = y_label,
                        title = "Column Points") +
          ggplot2::theme_minimal() +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))

        plots$col <- p_col
        if (graph.type == "col") print(p_col)
      }

      # Scree plot
      if (graph.type %in% c("scree", "all")) {
        scree_data <- data.frame(
          Dimension = 1:length(eigenvalues),
          Eigenvalue = eigenvalues,
          Variance = variance_explained
        )

        Dimension <- Variance <- NULL

        p_scree <- ggplot2::ggplot(scree_data, ggplot2::aes(x = Dimension, y = Variance)) +
          ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
          ggplot2::geom_line(color = "red", size = 1) +
          ggplot2::geom_point(color = "red", size = 3) +
          ggplot2::labs(x = "Dimension", y = "Percentage of Variance",
                        title = "Scree Plot") +
          ggplot2::theme_minimal() +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))

        plots$scree <- p_scree
        if (graph.type == "scree") print(p_scree)
      }

    } else {
      # Use base R plotting

      # Set up plot parameters
      if (save.plot) {
        grDevices::pdf(file.name, width = width, height = height)
      }

      # Calculate plot limits
      all_x <- c(row_coords[, axes[1]], col_coords[, axes[1]])
      all_y <- c(row_coords[, axes[2]], col_coords[, axes[2]])
      x_range <- range(all_x) * 1.2
      y_range <- range(all_y) * 1.2

      # Biplot
      if (graph.type %in% c("biplot", "all")) {
        graphics::plot(0, 0, type = "n",
                       xlim = x_range, ylim = y_range,
                       xlab = x_label, ylab = y_label,
                       main = ifelse(is.null(title),
                                     "Correspondence Analysis Biplot",
                                     title))
        graphics::abline(h = 0, v = 0, lty = 2, col = "gray")

        # Plot row points
        graphics::points(row_coords[, axes[1]], row_coords[, axes[2]],
                         col = row.color, pch = 19, cex = point.size/2)
        graphics::text(row_coords[, axes[1]], row_coords[, axes[2]],
                       labels = rownames(row_coords),
                       col = row.color, cex = text.size/4, pos = 3)

        # Plot column points
        graphics::points(col_coords[, axes[1]], col_coords[, axes[2]],
                         col = col.color, pch = 17, cex = point.size/2)
        graphics::text(col_coords[, axes[1]], col_coords[, axes[2]],
                       labels = rownames(col_coords),
                       col = col.color, cex = text.size/4, pos = 3)

        graphics::legend("topright",
                         legend = c("Rows", "Columns"),
                         col = c(row.color, col.color),
                         pch = c(19, 17))
      }

      # Row plot
      if (graph.type %in% c("row", "all")) {
        if (graph.type == "all") {
          # Ask user to press enter for next plot
          readline(prompt = "Press [enter] to see next plot")
        }

        graphics::plot(row_coords[, axes[1]], row_coords[, axes[2]],
                       col = row.color, pch = 19, cex = point.size/2,
                       xlim = range(row_coords[, axes[1]]) * 1.2,
                       ylim = range(row_coords[, axes[2]]) * 1.2,
                       xlab = x_label, ylab = y_label,
                       main = "Row Points")
        graphics::abline(h = 0, v = 0, lty = 2, col = "gray")
        graphics::text(row_coords[, axes[1]], row_coords[, axes[2]],
                       labels = rownames(row_coords),
                       col = row.color, cex = text.size/4, pos = 3)
      }

      # Column plot
      if (graph.type %in% c("col", "all")) {
        if (graph.type == "all") {
          readline(prompt = "Press [enter] to see next plot")
        }

        graphics::plot(col_coords[, axes[1]], col_coords[, axes[2]],
                       col = col.color, pch = 17, cex = point.size/2,
                       xlim = range(col_coords[, axes[1]]) * 1.2,
                       ylim = range(col_coords[, axes[2]]) * 1.2,
                       xlab = x_label, ylab = y_label,
                       main = "Column Points")
        graphics::abline(h = 0, v = 0, lty = 2, col = "gray")
        graphics::text(col_coords[, axes[1]], col_coords[, axes[2]],
                       labels = rownames(col_coords),
                       col = col.color, cex = text.size/4, pos = 3)
      }

      # Scree plot
      if (graph.type %in% c("scree", "all")) {
        if (graph.type == "all") {
          readline(prompt = "Press [enter] to see next plot")
        }

        graphics::barplot(variance_explained,
                          names.arg = 1:length(variance_explained),
                          col = "steelblue",
                          xlab = "Dimension",
                          ylab = "Percentage of Variance",
                          main = "Scree Plot",
                          ylim = c(0, max(variance_explained) * 1.1))
        graphics::lines(1:length(variance_explained) - 0.5 + 1.2/2,
                        variance_explained,
                        col = "red", lwd = 2)
        graphics::points(1:length(variance_explained) - 0.5 + 1.2/2,
                         variance_explained,
                         col = "red", pch = 19, cex = 1.5)
      }

      if (save.plot) {
        grDevices::dev.off()
        cat("Plot saved to:", file.name, "\n")
      }
    }

    results$plots <- plots
  }

  # Add class for S3 methods
  class(results) <- c("bcorresp", "list")

  return(results)
}

#' Print method for bcorresp objects
#'
#' @param x A bcorresp object
#' @param ... Additional arguments (not used)
#' @importFrom utils head
#' @export
print.bcorresp <- function(x, ...) {
  cat("\nbcorresp Analysis Results\n")
  cat("===========================\n")
  cat("Data dimensions:", nrow(x$data), "rows x", ncol(x$data), "columns\n")
  cat("Total inertia:", round(x$total.inertia, 4), "\n")
  cat("Chi-square test: X-squared =", round(x$chi2$statistic, 2),
      ", df =", x$chi2$df,
      ", p-value =", format.pval(x$chi2$p.value), "\n\n")
  cat("Dimensions kept:", x$n.dim, "\n")
  cat("Variance explained (first 5 dims):\n")
  n_show <- min(5, x$n.dim)
  var_summary <- data.frame(
    Dim = 1:n_show,
    Eigenvalue = round(x$eigenvalues[1:n_show], 4),
    `Var(%)` = round(x$variance.explained[1:n_show], 2),
    `Cumul(%)` = round(x$cumulative.variance[1:n_show], 2),
    check.names = FALSE
  )
  print(var_summary, row.names = FALSE)
  invisible(x)
}

#' Summary method for bcorresp objects
#'
#' @param object A bcorresp object
#' @param n_display Number of items to display (default = 5)
#' @param ... Additional arguments (not used)
#' @importFrom utils head
#' @export
summary.bcorresp <- function(object, n_display = 5, ...) {
  cat("\n===============================================\n")
  cat("       BCORRESP ANALYSIS SUMMARY\n")
  cat("===============================================\n\n")

  # Basic information
  cat("DATA INFORMATION:\n")
  cat("  Original dimensions:", nrow(object$data), "rows x",
      ncol(object$data), "columns\n")
  cat("  Total observations:", sum(object$data), "\n")
  cat("  Dimensions retained:", object$n.dim, "\n\n")

  # Chi-square test
  cat("INDEPENDENCE TEST:\n")
  cat("  Chi-square statistic:", round(object$chi2$statistic, 2), "\n")
  cat("  Degrees of freedom:", object$chi2$df, "\n")
  cat("  P-value:", format.pval(object$chi2$p.value), "\n")
  cat("  Total inertia:", round(object$total.inertia, 4), "\n\n")

  # Eigenvalues
  cat("DECOMPOSITION:\n")
  n_show <- min(10, object$n.dim)
  eig_df <- data.frame(
    Dimension = 1:n_show,
    Eigenvalue = round(object$eigenvalues[1:n_show], 4),
    `Variance(%)` = round(object$variance.explained[1:n_show], 2),
    `Cumulative(%)` = round(object$cumulative.variance[1:n_show], 2),
    check.names = FALSE
  )
  print(eig_df, row.names = FALSE)
  cat("\n")

  # Top contributing rows
  cat("TOP CONTRIBUTING ROWS (to first 2 dimensions):\n")
  if (object$n.dim >= 2) {
    row_contrib_total <- rowSums(object$row.contrib[, 1:2, drop = FALSE])
  } else {
    row_contrib_total <- object$row.contrib[, 1]
  }
  top_rows <- utils::head(sort(row_contrib_total, decreasing = TRUE), n_display)
  row_display <- data.frame(
    Row = names(top_rows),
    `Contribution(%)` = round(top_rows, 2),
    `Mass(%)` = round(100 * object$row.masses[names(top_rows)], 2),
    check.names = FALSE
  )
  print(row_display, row.names = FALSE)
  cat("\n")

  # Top contributing columns
  cat("TOP CONTRIBUTING COLUMNS (to first 2 dimensions):\n")
  if (object$n.dim >= 2) {
    col_contrib_total <- rowSums(object$col.contrib[, 1:2, drop = FALSE])
  } else {
    col_contrib_total <- object$col.contrib[, 1]
  }
  top_cols <- utils::head(sort(col_contrib_total, decreasing = TRUE), n_display)
  col_display <- data.frame(
    Column = names(top_cols),
    `Contribution(%)` = round(top_cols, 2),
    `Mass(%)` = round(100 * object$col.masses[names(top_cols)], 2),
    check.names = FALSE
  )
  print(col_display, row.names = FALSE)
  cat("\n")

  # Quality of representation
  cat("QUALITY OF REPRESENTATION (cos2) - Best represented:\n")
  if (object$n.dim >= 2) {
    row_quality <- rowSums(object$row.cos2[, 1:2, drop = FALSE])
    col_quality <- rowSums(object$col.cos2[, 1:2, drop = FALSE])
  } else {
    row_quality <- object$row.cos2[, 1]
    col_quality <- object$col.cos2[, 1]
  }

  cat("\nRows (first 2 dims):\n")
  top_row_quality <- utils::head(sort(row_quality, decreasing = TRUE), n_display)
  for (i in 1:length(top_row_quality)) {
    cat("  ", names(top_row_quality)[i], ": ",
        round(top_row_quality[i], 3), "\n", sep = "")
  }

  cat("\nColumns (first 2 dims):\n")
  top_col_quality <- utils::head(sort(col_quality, decreasing = TRUE), n_display)
  for (i in 1:length(top_col_quality)) {
    cat("  ", names(top_col_quality)[i], ": ",
        round(top_col_quality[i], 3), "\n", sep = "")
  }

  invisible(object)
}

#' Extract coordinates from bcorresp results
#'
#' @param ca_result A bcorresp object
#' @param type Type of coordinates to extract: "row", "column", or "both"
#' @param dimensions Which dimensions to extract (default = c(1, 2))
#' @export
bcorresp.coord <- function(ca_result,
                           type = c("row", "column", "both"),
                           dimensions = c(1, 2)) {
  type <- match.arg(type)

  if (!inherits(ca_result, "bcorresp")) {
    stop("Input must be a bcorresp object")
  }

  if (type == "row") {
    return(ca_result$row.coords[, dimensions, drop = FALSE])
  } else if (type == "column") {
    return(ca_result$col.coords[, dimensions, drop = FALSE])
  } else {
    return(list(
      rows = ca_result$row.coords[, dimensions, drop = FALSE],
      columns = ca_result$col.coords[, dimensions, drop = FALSE]
    ))
  }
}

#' Create a contribution plot for bcorresp results
#'
#' @param ca_result A bcorresp object
#' @param dimension Which dimension to plot contributions for (default = 1)
#' @param type Type of contributions to plot: "row" or "column"
#' @param n_display Number of top contributors to display (default = 10)
#' @param main Main title for the plot (optional)
#' @importFrom graphics abline barplot legend par
#' @export
bcorresp.contrib.plot <- function(ca_result,
                                  dimension = 1,
                                  type = c("row", "column"),
                                  n_display = 10,
                                  main = NULL) {
  type <- match.arg(type)

  if (!inherits(ca_result, "bcorresp")) {
    stop("Input must be a bcorresp object")
  }

  if (type == "row") {
    contrib <- ca_result$row.contrib[, dimension]
    if (is.null(main)) {
      main <- paste("Row Contributions to Dimension", dimension)
    }
  } else {
    contrib <- ca_result$col.contrib[, dimension]
    if (is.null(main)) {
      main <- paste("Column Contributions to Dimension", dimension)
    }
  }

  # Sort and get top contributors
  contrib <- sort(contrib, decreasing = TRUE)
  contrib <- utils::head(contrib, n_display)

  # Create bar plot
  graphics::par(mar = c(5, 8, 4, 2))
  graphics::barplot(rev(contrib),
                    horiz = TRUE,
                    las = 1,
                    names.arg = rev(names(contrib)),
                    col = "steelblue",
                    main = main,
                    xlab = "Contribution (%)")

  # Add reference line at average contribution
  avg_contrib <- 100 / length(if(type == "row") ca_result$row.masses else ca_result$col.masses)
  graphics::abline(v = avg_contrib, lty = 2, col = "red")

  # Add legend
  graphics::legend("topright",
                   legend = "Average",
                   lty = 2,
                   col = "red",
                   bty = "n")
}

#' Export bcorresp results to CSV files
#'
#' @param ca_result A bcorresp object
#' @param file_prefix Prefix for exported file names (default = "bcorresp_results")
#' @param include_data Whether to include the original data (default = FALSE)
#' @importFrom utils write.csv
#' @export
bcorresp.export <- function(ca_result,
                            file_prefix = "bcorresp_results",
                            include_data = FALSE) {

  if (!inherits(ca_result, "bcorresp")) {
    stop("Input must be a bcorresp object")
  }

  # Export row coordinates
  utils::write.csv(ca_result$row.coords,
                   file = paste0(file_prefix, "_row_coords.csv"),
                   row.names = TRUE)

  # Export column coordinates
  utils::write.csv(ca_result$col.coords,
                   file = paste0(file_prefix, "_col_coords.csv"),
                   row.names = TRUE)

  # Export eigenvalues and variance
  eig_df <- data.frame(
    Dimension = 1:length(ca_result$eigenvalues),
    Eigenvalue = ca_result$eigenvalues,
    Variance_Explained = ca_result$variance.explained,
    Cumulative_Variance = ca_result$cumulative.variance
  )
  utils::write.csv(eig_df,
                   file = paste0(file_prefix, "_eigenvalues.csv"),
                   row.names = FALSE)

  # Export contributions
  utils::write.csv(ca_result$row.contrib,
                   file = paste0(file_prefix, "_row_contributions.csv"),
                   row.names = TRUE)
  utils::write.csv(ca_result$col.contrib,
                   file = paste0(file_prefix, "_col_contributions.csv"),
                   row.names = TRUE)

  # Export quality (cos2)
  utils::write.csv(ca_result$row.cos2,
                   file = paste0(file_prefix, "_row_cos2.csv"),
                   row.names = TRUE)
  utils::write.csv(ca_result$col.cos2,
                   file = paste0(file_prefix, "_col_cos2.csv"),
                   row.names = TRUE)

  # Export masses and inertias
  masses_df <- data.frame(
    Row = names(ca_result$row.masses),
    Mass = ca_result$row.masses,
    Inertia = ca_result$row.inertias
  )
  utils::write.csv(masses_df,
                   file = paste0(file_prefix, "_row_masses_inertias.csv"),
                   row.names = FALSE)

  masses_df <- data.frame(
    Column = names(ca_result$col.masses),
    Mass = ca_result$col.masses,
    Inertia = ca_result$col.inertias
  )
  utils::write.csv(masses_df,
                   file = paste0(file_prefix, "_col_masses_inertias.csv"),
                   row.names = FALSE)

  # Optionally export original data
  if (include_data) {
    utils::write.csv(ca_result$data,
                     file = paste0(file_prefix, "_data.csv"),
                     row.names = TRUE)
  }

  cat("Results exported with prefix:", file_prefix, "\n")
  cat("Files created:\n")
  cat("  -", paste0(file_prefix, "_row_coords.csv\n"))
  cat("  -", paste0(file_prefix, "_col_coords.csv\n"))
  cat("  -", paste0(file_prefix, "_eigenvalues.csv\n"))
  cat("  -", paste0(file_prefix, "_row_contributions.csv\n"))
  cat("  -", paste0(file_prefix, "_col_contributions.csv\n"))
  cat("  -", paste0(file_prefix, "_row_cos2.csv\n"))
  cat("  -", paste0(file_prefix, "_col_cos2.csv\n"))
  cat("  -", paste0(file_prefix, "_row_masses_inertias.csv\n"))
  cat("  -", paste0(file_prefix, "_col_masses_inertias.csv\n"))
  if (include_data) {
    cat("  -", paste0(file_prefix, "_data.csv\n"))
  }
}

#' Interpret bcorresp results with guidance
#'
#' @param ca_result A bcorresp object
#' @param dimensions Which dimensions to interpret (default = c(1, 2))
#' @param n_points Number of points to highlight (default = 5)
#' @importFrom utils head tail
#' @export
bcorresp.interpret <- function(ca_result,
                               dimensions = c(1, 2),
                               n_points = 5) {

  if (!inherits(ca_result, "bcorresp")) {
    stop("Input must be a bcorresp object")
  }

  cat("\n===============================================\n")
  cat("      INTERPRETATION GUIDE - BCORRESP\n")
  cat("===============================================\n\n")

  # Variance explained
  var_exp <- ca_result$variance.explained
  cum_var <- sum(var_exp[dimensions])

  cat("1. VARIANCE EXPLAINED:\n")
  cat("   Dimension", dimensions[1], ":", round(var_exp[dimensions[1]], 1), "%\n")
  cat("   Dimension", dimensions[2], ":", round(var_exp[dimensions[2]], 1), "%\n")
  cat("   Total (dims", dimensions[1], "&", dimensions[2], "):",
      round(cum_var, 1), "%\n\n")

  if (cum_var > 70) {
    cat("   [checkmark] Excellent representation in these 2 dimensions (>70%)\n")
  } else if (cum_var > 50) {
    cat("   [warning] Moderate representation in these 2 dimensions (50-70%)\n")
    cat("   --> Consider examining additional dimensions\n")
  } else {
    cat("   [x] Poor representation in these 2 dimensions (<50%)\n")
    cat("   --> More dimensions needed for adequate representation\n")
  }

  cat("\n2. CHI-SQUARE TEST:\n")
  if (ca_result$chi2$p.value < 0.001) {
    cat("   [checkmark] Highly significant association (p < 0.001)\n")
  } else if (ca_result$chi2$p.value < 0.05) {
    cat("   [checkmark] Significant association (p < 0.05)\n")
  } else {
    cat("   [x] No significant association (p >= 0.05)\n")
    cat("   --> Variables may be independent\n")
  }

  # Top contributors
  cat("\n3. TOP CONTRIBUTORS:\n")
  for (d in dimensions) {
    cat("\n   Dimension", d, ":\n")

    # Rows
    top_row_contrib <- utils::head(sort(ca_result$row.contrib[, d],
                                        decreasing = TRUE), n_points)
    cat("   Rows:\n")
    for (i in 1:length(top_row_contrib)) {
      cat("     *", names(top_row_contrib)[i], ":",
          round(top_row_contrib[i], 1), "%\n")
    }

    # Columns
    top_col_contrib <- utils::head(sort(ca_result$col.contrib[, d],
                                        decreasing = TRUE), n_points)
    cat("   Columns:\n")
    for (i in 1:length(top_col_contrib)) {
      cat("     *", names(top_col_contrib)[i], ":",
          round(top_col_contrib[i], 1), "%\n")
    }
  }

  # Quality of representation
  cat("\n4. QUALITY OF REPRESENTATION (cos2):\n")

  row_cos2_sum <- rowSums(ca_result$row.cos2[, dimensions, drop = FALSE])
  col_cos2_sum <- rowSums(ca_result$col.cos2[, dimensions, drop = FALSE])

  well_rep_rows <- sum(row_cos2_sum > 0.5)
  well_rep_cols <- sum(col_cos2_sum > 0.5)

  cat("   Well represented (cos2 > 0.5):\n")
  cat("   * Rows:", well_rep_rows, "out of", length(row_cos2_sum), "\n")
  cat("   * Columns:", well_rep_cols, "out of", length(col_cos2_sum), "\n")

  # Extreme coordinates
  cat("\n5. EXTREME POSITIONS (most distinctive):\n")

  for (d in dimensions) {
    cat("\n   Dimension", d, ":\n")

    # Rows
    row_coords_d <- ca_result$row.coords[, d]
    extreme_rows <- c(
      utils::head(sort(row_coords_d, decreasing = TRUE), 2),
      utils::tail(sort(row_coords_d, decreasing = TRUE), 2)
    )
    cat("   Rows:\n")
    cat("     Positive:", paste(names(extreme_rows[1:2]), collapse = ", "), "\n")
    cat("     Negative:", paste(names(extreme_rows[3:4]), collapse = ", "), "\n")

    # Columns
    col_coords_d <- ca_result$col.coords[, d]
    extreme_cols <- c(
      utils::head(sort(col_coords_d, decreasing = TRUE), 2),
      utils::tail(sort(col_coords_d, decreasing = TRUE), 2)
    )
    cat("   Columns:\n")
    cat("     Positive:", paste(names(extreme_cols[1:2]), collapse = ", "), "\n")
    cat("     Negative:", paste(names(extreme_cols[3:4]), collapse = ", "), "\n")
  }

  cat("\n6. INTERPRETATION TIPS:\n")
  cat("   * Points far from origin --> Strong association patterns\n")
  cat("   * Points close together --> Similar profiles\n")
  cat("   * Points in same direction from origin --> Positive association\n")
  cat("   * Points in opposite directions --> Negative association\n")
  cat("   * Check cos2 before interpreting a point's position\n")
  cat("   * Focus on points with high contributions\n")

  cat("\n7. NEXT STEPS:\n")
  cat("   * Create biplot: bcorresp(data, graph.type='biplot')\n")
  cat("   * Check contributions: bcorresp.contrib.plot(result)\n")
  cat("   * Export results: bcorresp.export(result)\n")
  cat("   * Compare axes: try axes=c(1,3) or axes=c(2,3)\n")

  invisible(ca_result)
}

#' Quick example dataset for bcorresp
#' @export
bcorresp.example <- function() {
  # Create example contingency table
  data <- matrix(
    c(10,  5, 15, 20,  8,
      25, 30, 10,  5, 12,
      15, 20, 25, 30, 18,
      8, 12, 20, 15, 22,
      12,  8, 15, 25, 20),
    nrow = 5, byrow = TRUE
  )

  rownames(data) <- c("Region_A", "Region_B", "Region_C", "Region_D", "Region_E")
  colnames(data) <- c("Product_1", "Product_2", "Product_3", "Product_4", "Product_5")

  cat("Example dataset created:\n")
  print(data)
  cat("\nRun: result <- bcorresp(data)\n")

  return(data)
}
