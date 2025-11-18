################################################################################
# QAP Correlation Type I Error Test under Network Autocorrelation
################################################################################
#
# Purpose: Verify that QAP correlation maintains proper Type I error control
#          (5% false positive rate at α=0.05) when testing correlation between
#          two independent undirected networks that both exhibit transitivity.
#
# Author: Generated for borgworld package validation
# Date: 2025-11-18
################################################################################

# Set seed for reproducibility
set.seed(42)

# Load required packages
if (!require("sna", quietly = TRUE)) {
  stop("Package 'sna' is required. Install with: install.packages('sna')")
}

################################################################################
# FUNCTIONS
################################################################################

#' QAP Correlation Test
#'
#' Performs Quadratic Assignment Procedure (QAP) correlation test between
#' two matrices using row/column permutations.
#'
#' @param X First matrix (will be permuted)
#' @param Y Second matrix (stays fixed)
#' @param nperm Number of permutations (default 1000)
#' @return List with observed correlation, p-value, and permutation distribution
qap_cor <- function(X, Y, nperm = 1000) {
  # Input validation
  if (!is.matrix(X) || !is.matrix(Y)) {
    stop("X and Y must be matrices")
  }
  if (!identical(dim(X), dim(Y))) {
    stop("X and Y must have the same dimensions")
  }

  n <- nrow(X)

  # Vectorize matrices (typically we exclude diagonal for adjacency matrices)
  # For adjacency matrices, we usually look at off-diagonal elements
  get_off_diag <- function(mat) {
    mat[lower.tri(mat) | upper.tri(mat)]
  }

  # Observed correlation
  x_vec <- get_off_diag(X)
  y_vec <- get_off_diag(Y)
  obs_cor <- cor(x_vec, y_vec)

  # Permutation distribution
  perm_cors <- numeric(nperm)

  for (i in 1:nperm) {
    # Permute rows and columns of X simultaneously
    perm_idx <- sample(1:n)
    X_perm <- X[perm_idx, perm_idx]
    x_perm_vec <- get_off_diag(X_perm)
    perm_cors[i] <- cor(x_perm_vec, y_vec)
  }

  # Two-tailed p-value
  p_value <- mean(abs(perm_cors) >= abs(obs_cor))

  return(list(
    obs_cor = obs_cor,
    p_value = p_value,
    perm_dist = perm_cors,
    nperm = nperm
  ))
}


#' Generate Undirected Network with Transitivity
#'
#' Generates an undirected binary adjacency matrix with specified structural properties
#' using an iterative rewiring approach.
#'
#' @param n Number of nodes
#' @param density Target edge density (proportion of possible edges)
#' @param transitivity Target transitivity/clustering coefficient (0-1)
#' @param max_iter Maximum iterations for optimization (default 1000)
#' @param tol Tolerance for meeting targets (default 0.05)
#' @return n x n symmetric binary adjacency matrix
generate_network <- function(n, density = 0.1,
                             transitivity = 0.3, max_iter = 1000,
                             tol = 0.05) {

  # Start with random undirected graph at target density
  # For undirected graphs, we only count unique edges (upper triangle)
  n_edges <- round(n * (n - 1) / 2 * density)

  # Initialize with random edges
  adj <- matrix(0, n, n)
  diag(adj) <- 0  # no self-loops

  # Randomly place edges in upper triangle only
  possible_edges <- which(upper.tri(adj))
  edge_idx <- sample(possible_edges, n_edges)
  adj[edge_idx] <- 1

  # Make symmetric (undirected)
  adj <- adj + t(adj)

  # Iteratively modify to achieve target transitivity
  for (iter in 1:max_iter) {
    current_trans <- sna::gtrans(adj, mode = "graph")  # Use graph mode for undirected

    # Check if we've met target
    if (abs(current_trans - transitivity) < tol) {
      break
    }

    # Adjust transitivity
    if (current_trans < transitivity - tol) {
      # Add transitive edges (close triads)
      adj <- add_transitive_edge(adj)
    } else if (current_trans > transitivity + tol) {
      # Remove some transitive edges
      adj <- remove_transitive_edge(adj)
    }

    # Maintain approximate density
    current_density <- sum(adj) / (n * (n - 1))  # Count both directions
    if (current_density < density - 0.02) {
      # Add random edge (and its symmetric pair)
      zeros <- which(upper.tri(adj) & adj == 0)
      if (length(zeros) > 0) {
        idx <- sample(zeros, 1)
        adj[idx] <- 1
        # Add symmetric edge
        i <- ((idx - 1) %% n) + 1
        j <- ((idx - 1) %/% n) + 1
        adj[j, i] <- 1
      }
    } else if (current_density > density + 0.02) {
      # Remove random edge (and its symmetric pair)
      ones <- which(upper.tri(adj) & adj == 1)
      if (length(ones) > 0) {
        idx <- sample(ones, 1)
        adj[idx] <- 0
        # Remove symmetric edge
        i <- ((idx - 1) %% n) + 1
        j <- ((idx - 1) %/% n) + 1
        adj[j, i] <- 0
      }
    }
  }

  return(adj)
}


#' Add a Transitive Edge (Close a Triad) for Undirected Graph
#'
#' @param adj Symmetric adjacency matrix
#' @return Modified symmetric adjacency matrix
add_transitive_edge <- function(adj) {
  n <- nrow(adj)

  # Look for open triads: i-j and j-k exist, but not i-k
  for (attempt in 1:10) {  # Try a few times
    i <- sample(1:n, 1)
    j_options <- which(adj[i, ] == 1)

    if (length(j_options) > 0) {
      j <- sample(j_options, 1)
      k_options <- which(adj[j, ] == 1 & adj[i, ] == 0 & (1:n) != i)

      if (length(k_options) > 0) {
        k <- sample(k_options, 1)
        # Add edge i-k (both directions for undirected)
        adj[i, k] <- 1
        adj[k, i] <- 1
        return(adj)
      }
    }
  }

  return(adj)
}


#' Remove a Transitive Edge for Undirected Graph
#'
#' @param adj Symmetric adjacency matrix
#' @return Modified symmetric adjacency matrix
remove_transitive_edge <- function(adj) {
  n <- nrow(adj)

  # Look for transitive triads: i-j, j-k, and i-k all exist
  for (attempt in 1:10) {
    i <- sample(1:n, 1)
    k_options <- which(adj[i, ] == 1)

    if (length(k_options) > 0) {
      k <- sample(k_options, 1)
      # Check if there's a j such that i-j-k triangle exists
      j_options <- which(adj[i, ] == 1 & adj[k, ] == 1 & (1:n) != k)

      if (length(j_options) > 0) {
        # Remove i-k edge (both directions for undirected)
        adj[i, k] <- 0
        adj[k, i] <- 0
        return(adj)
      }
    }
  }

  return(adj)
}


################################################################################
# SIMULATION PARAMETERS
################################################################################

cat("\n")
cat(paste(rep("=", 80), collapse=""), "\n", sep="")
cat("QAP CORRELATION TYPE I ERROR TEST\n")
cat(paste(rep("=", 80), collapse=""), "\n\n", sep="")

# Simulation parameters
n_simulations <- 2000  # Number of independent (X,Y) pairs
n_permutations <- 1000  # QAP permutations per test
network_sizes <- c(30, 50)  # Test different network sizes

# Network generation parameters
target_density <- 0.1
target_transitivity <- 0.3


################################################################################
# RUN SIMULATION
################################################################################

results_summary <- list()

for (n_nodes in network_sizes) {
  cat("\n")
  cat(paste(rep("-", 80), collapse=""), "\n", sep="")
  cat(sprintf("NETWORK SIZE: n = %d nodes\n", n_nodes))
  cat(paste(rep("-", 80), collapse=""), "\n\n", sep="")

  # Storage for results
  p_values <- numeric(n_simulations)
  actual_cors <- numeric(n_simulations)
  transitivity_X <- numeric(n_simulations)
  transitivity_Y <- numeric(n_simulations)
  density_X <- numeric(n_simulations)
  density_Y <- numeric(n_simulations)

  # Track runtime
  start_time <- Sys.time()

  cat(sprintf("Generating %d pairs of networks and running QAP tests...\n",
              n_simulations))
  cat("Progress: ")

  # Run simulations
  for (i in 1:n_simulations) {
    # Progress indicator
    if (i %% 100 == 0) {
      cat(sprintf("%d ", i))
      flush.console()
    }

    # Generate two independent undirected networks with desired properties
    X <- generate_network(n_nodes, target_density, target_transitivity)
    Y <- generate_network(n_nodes, target_density, target_transitivity)

    # Record actual network properties
    transitivity_X[i] <- sna::gtrans(X, mode = "graph")
    transitivity_Y[i] <- sna::gtrans(Y, mode = "graph")
    density_X[i] <- sum(X) / (n_nodes * (n_nodes - 1))
    density_Y[i] <- sum(Y) / (n_nodes * (n_nodes - 1))

    # Run QAP test
    qap_result <- qap_cor(X, Y, nperm = n_permutations)

    # Store results
    p_values[i] <- qap_result$p_value
    actual_cors[i] <- qap_result$obs_cor
  }

  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

  cat("\n\n")

  ############################################################################
  # ANALYZE RESULTS
  ############################################################################

  # Type I error rate (proportion of p < 0.05)
  type1_error_rate <- mean(p_values < 0.05)

  # Store summary
  results_summary[[as.character(n_nodes)]] <- list(
    n = n_nodes,
    type1_error_rate = type1_error_rate,
    p_values = p_values,
    actual_cors = actual_cors,
    transitivity_X = transitivity_X,
    transitivity_Y = transitivity_Y,
    density_X = density_X,
    density_Y = density_Y,
    runtime = runtime
  )

  ############################################################################
  # OUTPUT RESULTS
  ############################################################################

  cat("RESULTS:\n")
  cat("--------\n\n")

  cat(sprintf("Type I Error Rate (p < 0.05): %.4f\n", type1_error_rate))
  cat(sprintf("Expected rate (nominal α):     %.4f\n", 0.05))
  cat(sprintf("Absolute deviation:            %.4f\n\n",
              abs(type1_error_rate - 0.05)))

  # 95% confidence interval for Type I error rate
  # Using normal approximation to binomial
  se <- sqrt(0.05 * 0.95 / n_simulations)
  ci_lower <- 0.05 - 1.96 * se
  ci_upper <- 0.05 + 1.96 * se

  cat(sprintf("Expected 95%% CI for Type I error: [%.4f, %.4f]\n",
              ci_lower, ci_upper))

  if (type1_error_rate >= ci_lower && type1_error_rate <= ci_upper) {
    cat("Type I error rate is within expected range\n\n")
  } else {
    cat("WARNING: Type I error rate is outside expected range!\n\n")
  }

  # P-value distribution statistics
  cat("P-value Distribution:\n")
  cat(sprintf("  Mean:   %.4f (should be ~0.5 if uniform)\n", mean(p_values)))
  cat(sprintf("  Median: %.4f (should be ~0.5 if uniform)\n", median(p_values)))
  cat(sprintf("  Min:    %.4f\n", min(p_values)))
  cat(sprintf("  Max:    %.4f\n\n", max(p_values)))

  # Kolmogorov-Smirnov test for uniformity
  ks_test <- ks.test(p_values, "punif")
  cat("Kolmogorov-Smirnov test for uniformity:\n")
  cat(sprintf("  D = %.4f, p-value = %.4f\n", ks_test$statistic, ks_test$p.value))
  if (ks_test$p.value > 0.05) {
    cat("  P-values are consistent with uniform distribution\n\n")
  } else {
    cat("  P-values deviate from uniform distribution\n\n")
  }

  # Actual correlation between X and Y (should be near 0)
  cat("Actual Correlations between Independent Networks:\n")
  cat(sprintf("  Mean:   %.4f (should be ~0)\n", mean(actual_cors)))
  cat(sprintf("  SD:     %.4f\n", sd(actual_cors)))
  cat(sprintf("  Range:  [%.4f, %.4f]\n\n", min(actual_cors), max(actual_cors)))

  # Network properties achieved
  cat("Achieved Network Properties:\n")
  cat("\n  DENSITY:\n")
  cat(sprintf("    Target: %.2f\n", target_density))
  cat(sprintf("    X: %.3f ± %.3f\n", mean(density_X), sd(density_X)))
  cat(sprintf("    Y: %.3f ± %.3f\n\n", mean(density_Y), sd(density_Y)))

  cat("  TRANSITIVITY:\n")
  cat(sprintf("    Target: %.2f\n", target_transitivity))
  cat(sprintf("    X: %.3f ± %.3f\n", mean(transitivity_X), sd(transitivity_X)))
  cat(sprintf("    Y: %.3f ± %.3f\n\n", mean(transitivity_Y), sd(transitivity_Y)))

  # Runtime
  cat(sprintf("Runtime: %.1f seconds (%.2f minutes)\n", runtime, runtime/60))
  cat(sprintf("Average time per simulation: %.3f seconds\n\n",
              runtime/n_simulations))
}


################################################################################
# FINAL SUMMARY
################################################################################

cat("\n")
cat(paste(rep("=", 80), collapse=""), "\n", sep="")
cat("SUMMARY ACROSS ALL NETWORK SIZES\n")
cat(paste(rep("=", 80), collapse=""), "\n\n", sep="")

summary_df <- data.frame(
  n = network_sizes,
  Type_I_Error = sapply(results_summary, function(x) x$type1_error_rate),
  Mean_Cor = sapply(results_summary, function(x) mean(x$actual_cors)),
  Runtime_sec = sapply(results_summary, function(x) x$runtime)
)

print(summary_df)

cat("\n")
cat("CONCLUSION:\n")
cat("-----------\n")

all_within_ci <- all(abs(summary_df$Type_I_Error - 0.05) < 0.015)

if (all_within_ci) {
  cat("QAP correlation successfully controls Type I error rate at α=0.05\n")
  cat("  even when undirected networks exhibit transitivity.\n\n")
} else {
  cat("WARNING: QAP correlation shows inflated or deflated Type I error\n")
  cat("  when networks have autocorrelation structure.\n\n")
}

cat(sprintf("Total runtime: %.1f minutes\n", sum(summary_df$Runtime_sec)/60))

cat("\n")
cat(paste(rep("=", 80), collapse=""), "\n", sep="")
cat("SIMULATION COMPLETE\n")
cat(paste(rep("=", 80), collapse=""), "\n\n", sep="")

################################################################################
# END
################################################################################
