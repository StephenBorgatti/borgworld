#' Johnson's Hierarchical Clustering
#'
#' Implements Johnson's Hierarchical Clustering for an item-by-item proximity matrix
#'
#' @param data A data frame containing a square proximity matrix (similarities or dissimilarities).
#'        Row names and column names are used as labels if available.
#' @param type Character string specifying the type of proximity measure;
#'        "similarities" or "dissimilarities" (can be abbreviated to "sim" or "diss").
#'        This parameter is required with no default to ensure explicit specification.
#' @param method Linkage method: "average" (default), "single", or "complete" (can be abbreviated)
#' @param labels Optional character vector of item labels. If NULL, uses row names from data
#' @param print_dendrogram Logical; whether to print ASCII dendrogram (default: TRUE)
#' @param plot_dendrogram Logical; whether to plot graphical dendrogram (default: TRUE)
#'
#' @return A list of class "bhiclus" containing:
#'   \item{partition_table}{Matrix with items as rows and nested partitions as columns}
#'   \item{hclust_obj}{hclust object for further analysis}
#'   \item{moca}{Data frame of measures of cluster adequacy from bmoca, including:
#'              n_clusters, Corr, Gamma, Modularity, and Silhouette}
#'   \item{cophenetic}{Cophenetic correlation coefficient from bcophenetic}
#'   \item{labels}{Character vector of item labels}
#'   \item{type}{Type of proximity measure used}
#'   \item{method}{Linkage method used}
#'
#' @export
#' @examples
#' # Create example similarity data frame
#' sim_df <- data.frame(
#'   A = c(1.0, 0.8, 0.3, 0.2, 0.1),
#'   B = c(0.8, 1.0, 0.4, 0.3, 0.2),
#'   C = c(0.3, 0.4, 1.0, 0.7, 0.5),
#'   D = c(0.2, 0.3, 0.7, 1.0, 0.6),
#'   E = c(0.1, 0.2, 0.5, 0.6, 1.0)
#' )
#' rownames(sim_df) <- colnames(sim_df)
#'
#' # Note: type parameter is required
#' result <- bhiclus(sim_df, type = "sim")
#'
#' # Or with dissimilarities
#' dissim_df <- as.data.frame(as.matrix(dist(iris[,1:4])))
#' result2 <- bhiclus(dissim_df, type = "diss")
#' @importFrom stats as.dist as.formula cophenetic

bhiclus <- function(data,
                    type,  # No default - force user to specify
                    method = c("average", "single", "complete"),
                    labels = NULL,
                    print_dendrogram = TRUE,
                    plot_dendrogram = TRUE) {

  # Handle data frame input
  if (is.data.frame(data)) {
    # Extract labels from row names if not provided
    if (is.null(labels) && !is.null(rownames(data))) {
      labels <- rownames(data)
    }
    # Convert to matrix
    proximity_matrix <- as.matrix(data)
  } else if (is.matrix(data)) {
    proximity_matrix <- data
    if (is.null(labels) && !is.null(rownames(data))) {
      labels <- rownames(data)
    }
  } else {
    stop("data must be a data frame or matrix")
  }

  # Validate that it's square
  if (nrow(proximity_matrix) != ncol(proximity_matrix)) {
    stop("data must be a square proximity matrix")
  }

  n <- nrow(proximity_matrix)

  # Match arguments with partial matching
  type <- match.arg(type, c("similarities", "dissimilarities"))
  method <- match.arg(method)

  # Set default labels if still not provided
  if (is.null(labels)) {
    if (!is.null(colnames(proximity_matrix))) {
      labels <- colnames(proximity_matrix)
    } else {
      labels <- paste0("Item", 1:n)
    }
  }

  # Check label length
  if (length(labels) != n) {
    stop("Number of labels must match the dimension of the proximity matrix")
  }

  # Store original proximity matrix for MOCA (before conversion)
  # But convert similarities to dissimilarities if needed for MOCA
  if (type == "similarities") {
    # Find maximum value (excluding diagonal for safety)
    temp_mat <- proximity_matrix
    diag(temp_mat) <- NA
    max_val <- max(temp_mat, na.rm = TRUE)
    proximity_matrix <- max_val - proximity_matrix
  }

  # Ensure diagonal is zero
  diag(proximity_matrix) <- 0

  # Now proximity_matrix is in dissimilarity form

  # Convert to dist object and perform hierarchical clustering
  dist_obj <- as.dist(proximity_matrix)
  hclust_obj <- hclust(dist_obj, method = method)

  # Create partition table
  partition_table <- create_partition_table(hclust_obj, n, labels)

  # Calculate measures of cluster adequacy
  # Always pass type = "d" since we converted to dissimilarities above
  moca_results <- bmoca(proximity_matrix, partition_table, type = "d")

  # Calculate cophenetic correlation
  cophenet_corr <- bcophenetic(proximity_matrix, hclust_obj, type = "d")

  # Print ASCII dendrogram if requested
  if (print_dendrogram) {
    print_ascii_dendrogram(hclust_obj, labels)
  }

  # Plot graphical dendrogram if requested
  if (plot_dendrogram) {
    plot_graphical_dendrogram(hclust_obj, labels)
  }

  # Return results
  result <- list(
    partition_table = partition_table,
    hclust_obj = hclust_obj,
    moca = moca_results,
    cophenetic = cophenet_corr,
    labels = labels,
    type = type,
    method = method
  )

  class(result) <- "bhiclus"
  return(result)
}

#' Cophenetic Correlation Coefficient
#'
#' Calculate the cophenetic correlation coefficient for a hierarchical clustering.
#' This measures how faithfully the dendrogram preserves the pairwise distances
#' between the original observations.
#'
#' @param proximity_matrix A square matrix of proximities (either similarities or dissimilarities).
#' @param hclust_obj An hclust object from hierarchical clustering.
#' @param type Character string specifying the type of proximity measure;
#'        "similarities" or "dissimilarities" (can be abbreviated to "sim" or "diss").
#'        This parameter is required with no default to ensure explicit specification.
#'
#' @return A single numeric value: the cophenetic correlation coefficient.
#'        Higher values (closer to 1) indicate that the dendrogram structure
#'        better preserves the original proximity relationships.
#'
#' @details
#' The cophenetic correlation is calculated between:
#' \itemize{
#'   \item The original proximity matrix (lower triangle)
#'   \item The cophenetic distances from the dendrogram (lower triangle)
#' }
#'
#' For dissimilarity data, a positive correlation is expected (small distances
#' cluster together). For similarity data, the correlation is negated so that
#' higher values always indicate better fit.
#'
#' @export
#' @examples
#' # Create a distance matrix
#' dist_mat <- as.matrix(dist(iris[,1:4]))
#'
#' # Perform hierarchical clustering
#' hc <- hclust(as.dist(dist_mat))
#'
#' # Calculate cophenetic correlation
#' bcophenetic(dist_mat, hc, type = "diss")
#'
#' @seealso \code{\link{bhiclus}} which automatically calls bcophenetic
bcophenetic <- function(proximity_matrix, hclust_obj, type) {
  # Input validation
  if (!is.matrix(proximity_matrix) || nrow(proximity_matrix) != ncol(proximity_matrix)) {
    stop("proximity_matrix must be a square matrix")
  }

  if (!inherits(hclust_obj, "hclust")) {
    stop("hclust_obj must be an hclust object")
  }

  # Match type argument
  type <- match.arg(type, c("similarities", "dissimilarities"))

  # Get lower triangles
  prox_lower <- proximity_matrix[lower.tri(proximity_matrix)]
  coph_dist <- cophenetic(hclust_obj)
  coph_lower <- coph_dist[lower.tri(coph_dist)]

  # Check that dimensions match
  if (length(prox_lower) != length(coph_lower)) {
    warning("Dimension mismatch between proximity matrix and cophenetic distances")
    return(NA)
  }

  # Calculate correlation with error handling
  cophenet_corr <- suppressWarnings(cor(prox_lower, coph_lower, use = "complete.obs"))

  # Check for NA
  if (is.na(cophenet_corr)) {
    warning("Could not calculate cophenetic correlation")
    return(NA)
  }

  # Adjust sign based on type
  if (type == "dissimilarities") {
    # For dissimilarities, expect positive correlation
    return(cophenet_corr)
  } else {
    # For similarities, expect negative correlation with cophenetic distances
    # So negate to make higher values better
    return(-1 * cophenet_corr)
  }
}

#' Create partition table from hclust object
#'
#' @noRd
create_partition_table <- function(hclust_obj, n, labels) {
  # Initialize partition table
  partition_table <- matrix(0, nrow = n, ncol = n - 1)
  rownames(partition_table) <- labels
  colnames(partition_table) <- paste0("P", 1:(n-1))

  # Track which cluster each item belongs to at each step
  cluster_membership <- 1:n

  # Process each merge
  for (step in 1:(n-1)) {
    # Get merged clusters
    cluster1 <- hclust_obj$merge[step, 1]
    cluster2 <- hclust_obj$merge[step, 2]

    # Update cluster memberships
    items_to_merge <- c()

    if (cluster1 < 0) {
      items_to_merge <- c(items_to_merge, -cluster1)
    } else {
      # Find all items in this previously merged cluster
      items_to_merge <- c(items_to_merge, which(cluster_membership == n + cluster1))
    }

    if (cluster2 < 0) {
      items_to_merge <- c(items_to_merge, -cluster2)
    } else {
      items_to_merge <- c(items_to_merge, which(cluster_membership == n + cluster2))
    }

    # Update membership for merged items
    cluster_membership[items_to_merge] <- n + step

    # Create partition vector for this step
    unique_clusters <- unique(cluster_membership)
    partition <- numeric(n)
    for (i in seq_along(unique_clusters)) {
      partition[cluster_membership == unique_clusters[i]] <- i
    }

    partition_table[, step] <- partition
  }

  return(partition_table)
}

#' Print ASCII dendrogram
#'
#' @noRd
print_ascii_dendrogram <- function(hclust_obj, labels) {
  n <- length(labels)
  heights <- round(hclust_obj$height, 1)
  # Calculate dendrogram layout
  order <- hclust_obj$order
  # Print header with labels
  cat("\n")
  # Print labels vertically (downwards)
  max_label_len <- max(nchar(labels))
  for (row in 1:max_label_len) {
    cat("        ")  # Space for height column
    for (i in 1:n) {
      label <- labels[order[i]]
      if (nchar(label) >= row) {
        cat(substr(label, row, row))
      } else {
        cat(" ")
      }
      cat(" ")  # Single space between columns
    }
    cat("\n")
  }
  # Print separator
  cat("------  ")
  for (i in 1:n) {
    cat("- ")
  }
  cat("\n")
  # Build cluster hierarchy
  clusters <- list()
  for (i in 1:n) {
    clusters[[i]] <- i
  }
  # Process each merge and display
  for (step in 1:(n-1)) {
    c1 <- hclust_obj$merge[step, 1]
    c2 <- hclust_obj$merge[step, 2]
    # Get items in merging clusters
    if (c1 < 0) {
      items1 <- -c1
    } else {
      items1 <- clusters[[n + c1]]
    }
    if (c2 < 0) {
      items2 <- -c2
    } else {
      items2 <- clusters[[n + c2]]
    }
    # Store new cluster
    clusters[[n + step]] <- c(items1, items2)
    # For display, use cutree to get the actual clustering at this merge level
    num_clusters <- n - step
    current_clustering <- cutree(hclust_obj, k = num_clusters)
    # Create display line with spaces between columns
    # We need 2n-1 positions: n columns and n-1 spaces
    line <- rep(".", 2*n - 1)
    # Mark each cluster that has more than one member
    for (clust_id in unique(current_clustering)) {
      members <- which(current_clustering == clust_id)
      if (length(members) > 1) {
        # Multi-member cluster - mark its range
        positions <- match(members, order)
        min_pos <- min(positions)
        max_pos <- max(positions)
        # Convert to expanded positions (including spaces)
        # Position i in original becomes position 2*i-1 in expanded
        # Space after position i becomes position 2*i
        expanded_min <- 2 * min_pos - 1
        expanded_max <- 2 * max_pos - 1
        # Mark all positions from min to max (including spaces)
        line[expanded_min:expanded_max] <- "X"
      }
    }
    # Print the line
    cat(sprintf("%6.1f  ", heights[step]))
    for (i in 1:(2*n - 1)) {
      cat(line[i])
    }
    cat("\n")
  }
  cat("\n")
}

#' Plot graphical dendrogram
#'
#' @noRd
plot_graphical_dendrogram <- function(hclust_obj, labels) {
  # Set up plot parameters more carefully
  old_par <- par(no.readonly = TRUE)

  # Only try to restore parameters that can be safely restored
  on.exit({
    suppressWarnings(try(par(old_par), silent = TRUE))
  })

  # Set margins
  par(mar = c(5, 4, 4, 2) + 0.1)

  # Create dendrogram plot
  plot(hclust_obj,
       labels = labels,
       main = "Hierarchical Clustering Dendrogram",
       xlab = "Items",
       ylab = "Distance",
       sub = "",
       hang = -1)

  # Add grid for better readability
  grid(ny = NULL, nx = 0, col = "lightgray", lty = "dotted")
}

#' Print method for bhiclus objects
#'
#' @param x A bhiclus object
#' @param ... Additional arguments (unused)
#' @export
print.bhiclus <- function(x, ...) {
  cat("Johnson's Hierarchical Clustering Results\n")
  cat("=========================================\n\n")

  cat("Number of items:", length(x$labels), "\n")
  cat("Linkage method:", x$method, "\n")
  cat("Proximity type:", x$type, "\n")
  cat("Cophenetic correlation:", round(x$cophenetic, 3), "\n\n")

  cat("Partition Table:\n")
  print(x$partition_table)

  cat("\nMeasures of Cluster Adequacy (MOCA):\n")
  # Round the numeric columns for display
  moca_display <- x$moca
  numeric_cols <- c("Corr", "Gamma", "Modularity", "Silhouette")
  for (col in numeric_cols) {
    if (col %in% colnames(moca_display)) {
      moca_display[[col]] <- round(moca_display[[col]], 3)
    }
  }
  print(moca_display)

  invisible(x)
}

#' Summary method for bhiclus objects
#'
#' @param object A bhiclus object
#' @param ... Additional arguments (unused)
#' @export
summary.bhiclus <- function(object, ...) {
  cat("Johnson's Hierarchical Clustering Summary\n")
  cat("=========================================\n\n")

  cat("Number of items:", length(object$labels), "\n")
  cat("Linkage method:", object$method, "\n")
  cat("Proximity type:", object$type, "\n")

  # Get height statistics
  heights <- object$hclust_obj$height
  cat("\nMerge height statistics:\n")
  cat("  Min:", round(min(heights), 2), "\n")
  cat("  1st Qu:", round(quantile(heights, 0.25), 2), "\n")
  cat("  Median:", round(median(heights), 2), "\n")
  cat("  Mean:", round(mean(heights), 2), "\n")
  cat("  3rd Qu:", round(quantile(heights, 0.75), 2), "\n")
  cat("  Max:", round(max(heights), 2), "\n")

  # MOCA summary
  cat("\nCluster Adequacy (Corr) summary:\n")
  corr_values <- object$moca$Corr
  cat("  Min:", round(min(corr_values), 3), "\n")
  cat("  Max:", round(max(corr_values), 3), "\n")
  cat("  Best partition:", rownames(object$moca)[which.max(corr_values)],
      "with", object$moca$n_clusters[which.max(corr_values)], "clusters\n")

  invisible(object)
}
