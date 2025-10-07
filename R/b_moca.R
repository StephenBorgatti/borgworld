#' Measures of Cluster Adequacy (MOCA)
#'
#' Calculate measures of cluster adequacy for evaluating how well a set of partitions
#' fits a proximity matrix. Implements multiple measures including correlation,
#' modularity, Hubert's Gamma, and average silhouette width.
#'
#' @param proximity_matrix A square matrix of proximities (either similarities or dissimilarities).
#' @param partitions A matrix where each column represents a partition (cluster assignment).
#'        Rows correspond to items and values indicate cluster membership.
#'        Can also be a data frame which will be converted to a matrix.
#' @param type Character string specifying the type of proximity measure;
#'        "similarities" or "dissimilarities" (can be abbreviated to "sim" or "diss").
#'        This parameter is required with no default to ensure explicit specification.
#'
#' @return A data frame with rows corresponding to partitions and columns:
#'   \item{n_clusters}{Number of clusters in each partition}
#'   \item{Corr}{Correlation between proximity and same-cluster matrices.
#'         For dissimilarities: negative correlation expected (multiplied by -1 so higher is better).
#'         Higher values indicate better cluster adequacy.
#'         Returns NA when all items are in a single cluster.}
#'   \item{Gamma}{Hubert's Gamma statistic. Normalized sum of products between proximity
#'         and same-cluster matrices. Different from Corr as it uses products rather than correlation.}
#'   \item{Modularity}{Newman's modularity for weighted networks. For dissimilarities,
#'         distances are converted to similarities using 1/(1+dist).
#'         Measures how much the clustering exceeds random expectation.}
#'   \item{Silhouette}{Average silhouette width across all items. Ranges from -1 to 1,
#'         where higher values indicate better-defined clusters.}
#'
#' @details
#' The measures calculated are:
#'
#' \itemize{
#'   \item \strong{Corr}: Correlation between proximity matrix and binary same-cluster matrix
#'   \item \strong{Gamma}: Hubert's Gamma - normalized sum of products (not correlation)
#'   \item \strong{Modularity}: Newman's modularity Q, treating proximity as a weighted network
#'   \item \strong{Silhouette}: Average of individual silhouette values s(i)
#' }
#'
#' For dissimilarity matrices, correlations are multiplied by -1 and distances are
#' converted to similarities for modularity using 1/(1+dist).
#'
#' When all items are in a single cluster, correlation-based measures return NA.
#'
#' @export
#' @examples
#' # Example with dissimilarities - note type is required
#' dissim <- matrix(c(0, 1, 4, 5,
#'                    1, 0, 3, 4,
#'                    4, 3, 0, 1,
#'                    5, 4, 1, 0), nrow = 4)
#'
#' # Create two partitions
#' partitions <- matrix(c(1, 1, 2, 2,  # First partition: two clusters
#'                        1, 2, 3, 4), # Second partition: four clusters
#'                      nrow = 4, ncol = 2)
#' colnames(partitions) <- c("P1", "P2")
#'
#' # Calculate measures of cluster adequacy - type is required
#' bmoca(dissim, partitions, type = "diss")
#'
#' @seealso \code{\link{bhiclus}} which automatically calls bmoca on its partition table
bmoca <- function(proximity_matrix, partitions, type) {  # No default for type
  # Input validation
  if (!is.matrix(proximity_matrix) || nrow(proximity_matrix) != ncol(proximity_matrix)) {
    stop("proximity_matrix must be a square matrix")
  }

  # Convert to matrix if needed
  if (is.data.frame(partitions)) {
    partitions <- as.matrix(partitions)
  }

  if (!is.matrix(partitions)) {
    stop("partitions must be a matrix or data frame")
  }

  if (nrow(proximity_matrix) != nrow(partitions)) {
    stop("Number of items in proximity_matrix must match number of rows in partitions")
  }

  # Match type argument
  type <- match.arg(type, c("similarities", "dissimilarities"))

  n <- nrow(proximity_matrix)
  n_partitions <- ncol(partitions)

  # Get lower triangle of proximity matrix (excluding diagonal)
  prox_lower <- proximity_matrix[lower.tri(proximity_matrix)]

  # Prepare similarity matrix for modularity
  if (type == "dissimilarities") {
    # Convert distances to similarities for modularity
    sim_matrix <- 1 / (1 + proximity_matrix)
    diag(sim_matrix) <- 0  # No self-loops
  } else {
    sim_matrix <- proximity_matrix
    diag(sim_matrix) <- 0
  }

  # Initialize results
  results <- data.frame(
    n_clusters = integer(n_partitions),
    Corr = numeric(n_partitions),
    Gamma = numeric(n_partitions),
    Modularity = numeric(n_partitions),
    Silhouette = numeric(n_partitions)
  )

  # Calculate for each partition
  for (p in 1:n_partitions) {
    partition <- partitions[, p]

    # Number of clusters
    n_clust <- length(unique(partition))
    results$n_clusters[p] <- n_clust

    # Create same-cluster matrix (1 if same cluster, 0 otherwise)
    same_cluster <- outer(partition, partition, FUN = "==")
    same_cluster_lower <- as.numeric(same_cluster[lower.tri(same_cluster)])

    # 1. Correlation (Corr)
    if (n_clust == 1) {
      # All 1s - correlation is undefined
      results$Corr[p] <- NA
      results$Gamma[p] <- NA
    } else {
      # Calculate correlation
      corr_value <- suppressWarnings(cor(prox_lower, same_cluster_lower))

      # For dissimilarities, multiply by -1
      if (type == "dissimilarities") {
        corr_value <- -1 * corr_value
      }
      results$Corr[p] <- corr_value

      # 2. Hubert's Gamma (normalized sum of products, not correlation)
      # This gives a different result than Corr
      gamma_value <- sum(prox_lower * same_cluster_lower) / length(prox_lower)
      # Normalize by the standard deviations
      sd_prox <- sd(prox_lower)
      sd_same <- sd(same_cluster_lower)

      if (sd_prox > 0 && sd_same > 0) {
        gamma_value <- gamma_value / (sd_prox * sd_same)
      } else {
        gamma_value <- NA
      }

      # For dissimilarities, negate so higher is better
      if (type == "dissimilarities") {
        gamma_value <- -gamma_value
      }

      results$Gamma[p] <- gamma_value
    }

    # 3. Newman's Modularity
    # Calculate modularity using similarity matrix
    if (n_clust == 1) {
      results$Modularity[p] <- 0  # No community structure
    } else {
      # Calculate modularity Q
      total_weight <- sum(sim_matrix) / 2  # Each edge counted twice

      if (total_weight > 0) {
        # Degree of each node
        degrees <- rowSums(sim_matrix)

        # Modularity calculation
        Q <- 0
        for (i in 1:(n-1)) {
          for (j in (i+1):n) {
            if (partition[i] == partition[j]) {
              # Same community
              Q <- Q + (sim_matrix[i,j] - (degrees[i] * degrees[j]) / (2 * total_weight))
            }
          }
        }
        results$Modularity[p] <- Q / total_weight
      } else {
        results$Modularity[p] <- NA
      }
    }

    # 4. Average Silhouette Width
    if (n_clust == 1 || n_clust == n) {
      # Silhouette undefined for 1 cluster or n singleton clusters
      results$Silhouette[p] <- NA
    } else {
      # Use the cluster package's silhouette function if available
      # Otherwise use our implementation
      if (requireNamespace("cluster", quietly = TRUE)) {
        # Use robust implementation from cluster package
        if (type == "dissimilarities") {
          sil_obj <- cluster::silhouette(partition, as.dist(proximity_matrix))
        } else {
          # Convert similarities to dissimilarities for silhouette
          dist_for_sil <- as.dist(1 - proximity_matrix)
          sil_obj <- cluster::silhouette(partition, dist_for_sil)
        }
        results$Silhouette[p] <- mean(sil_obj[, "sil_width"])
      } else {
        # Manual implementation with fixes for edge cases
        sil_values <- numeric(n)

        for (i in 1:n) {
          same_clust <- which(partition == partition[i])

          # For singleton clusters, silhouette is 0 by convention
          if (length(same_clust) == 1) {
            sil_values[i] <- 0
            next
          }

          # a(i) = average distance to items in same cluster
          if (type == "dissimilarities") {
            a_i <- mean(proximity_matrix[i, same_clust[same_clust != i]])
          } else {
            # For similarities, use 1-sim as distance
            a_i <- mean(1 - proximity_matrix[i, same_clust[same_clust != i]])
          }

          # b(i) = minimum average distance to items in other clusters
          other_clusts <- unique(partition[partition != partition[i]])
          if (length(other_clusts) > 0) {
            b_values <- numeric(length(other_clusts))
            for (k in seq_along(other_clusts)) {
              other_items <- which(partition == other_clusts[k])
              if (type == "dissimilarities") {
                b_values[k] <- mean(proximity_matrix[i, other_items])
              } else {
                b_values[k] <- mean(1 - proximity_matrix[i, other_items])
              }
            }
            b_i <- min(b_values)

            # Silhouette value
            if (max(a_i, b_i) > 0) {
              sil_values[i] <- (b_i - a_i) / max(a_i, b_i)
            } else {
              sil_values[i] <- 0
            }
          } else {
            # No other clusters
            sil_values[i] <- 0
          }
        }

        results$Silhouette[p] <- mean(sil_values)
      }
    }
  }

  # Add row names if partitions have column names
  if (!is.null(colnames(partitions))) {
    rownames(results) <- colnames(partitions)
  }

  return(results)
}
