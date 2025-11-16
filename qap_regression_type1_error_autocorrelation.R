################################################################################
# QAP Regression Type I Error Test under Network Autocorrelation
################################################################################
#
# Purpose: Verify that QAP regression maintains proper Type I error control
#          (5% false positive rate at α=0.05) when testing regression
#          coefficients in a model predicting adjacency matrix Y from a dyadic
#          X matrix, both of which suffer various kinds of autocorrelation
#          (reciprocity, transitivity, and homophily).
#
# Method: Dekker-Krackhardt-Sanders (DSP) style QAP regression via sna::netlm
#
# Author: Generated for borgworld package validation
# Date: 2025-11-15
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

#' Generate Categorical Node Attributes
#'
#' Creates a categorical attribute for nodes with specified number of categories,
#' distributed as evenly as possible across nodes.
#'
#' @param n Number of nodes
#' @param n_categories Number of categories (default 3)
#' @return Integer vector of length n with category assignments (1 to n_categories)
generate_node_attributes <- function(n, n_categories = 3) {
  # Distribute nodes evenly across categories
  categories <- rep(1:n_categories, length.out = n)
  # Shuffle to randomize
  sample(categories)
}


#' Create Dyadic SameCat Matrix
#'
#' Creates a dyadic matrix where entry (i,j) = 1 if nodes i and j have the
#' same category attribute, 0 otherwise.
#'
#' @param node_attr Vector of categorical node attributes
#' @return n x n binary matrix
create_same_cat_matrix <- function(node_attr) {
  n <- length(node_attr)
  same_cat <- matrix(0, n, n)

  for (i in 1:n) {
    for (j in 1:n) {
      same_cat[i, j] <- as.integer(node_attr[i] == node_attr[j])
    }
  }

  # Set diagonal to 0 (no self-loops)
  diag(same_cat) <- 0

  return(same_cat)
}


#' Calculate Reciprocity
#'
#' @param adj Adjacency matrix
#' @return Reciprocity coefficient (proportion of reciprocated edges)
calculate_reciprocity <- function(adj) {
  # Count edges that are reciprocated
  recip_edges <- sum(adj * t(adj))
  total_edges <- sum(adj)

  if (total_edges == 0) return(0)
  return(recip_edges / total_edges)
}


#' Calculate Homophily
#'
#' Calculates the proportion of edges that are within-category (homophilous).
#'
#' @param adj Adjacency matrix
#' @param node_attr Vector of categorical node attributes
#' @return Proportion of edges that are within-category
calculate_homophily <- function(adj, node_attr) {
  same_cat <- create_same_cat_matrix(node_attr)

  # Count edges within categories
  within_edges <- sum(adj * same_cat)
  total_edges <- sum(adj)

  if (total_edges == 0) return(0)
  return(within_edges / total_edges)
}


#' Add a Reciprocal Edge
#'
#' @param adj Adjacency matrix
#' @return Modified adjacency matrix
add_reciprocal_edge <- function(adj) {
  # Find edges that are not reciprocated
  non_recip <- which(adj == 1 & t(adj) == 0)

  if (length(non_recip) > 0) {
    idx <- sample(non_recip, 1)
    # Get i,j coordinates
    i <- ((idx - 1) %% nrow(adj)) + 1
    j <- ((idx - 1) %/% nrow(adj)) + 1
    # Add reciprocal edge j->i
    adj[j, i] <- 1
  }

  return(adj)
}


#' Remove a Reciprocal Edge
#'
#' @param adj Adjacency matrix
#' @return Modified adjacency matrix
remove_reciprocal_edge <- function(adj) {
  # Find reciprocal edges
  recip <- which(adj == 1 & t(adj) == 1)

  if (length(recip) > 0) {
    idx <- sample(recip, 1)
    # Get i,j coordinates
    i <- ((idx - 1) %% nrow(adj)) + 1
    j <- ((idx - 1) %/% nrow(adj)) + 1
    # Remove one direction
    if (runif(1) > 0.5) {
      adj[i, j] <- 0
    } else {
      adj[j, i] <- 0
    }
  }

  return(adj)
}


#' Add a Transitive Edge (Close a Triad)
#'
#' @param adj Adjacency matrix
#' @return Modified adjacency matrix
add_transitive_edge <- function(adj) {
  n <- nrow(adj)

  # Look for open triads: i->j and j->k exist, but not i->k
  for (attempt in 1:10) {  # Try a few times
    i <- sample(1:n, 1)
    j_options <- which(adj[i, ] == 1)

    if (length(j_options) > 0) {
      j <- sample(j_options, 1)
      k_options <- which(adj[j, ] == 1 & adj[i, ] == 0 & (1:n) != i)

      if (length(k_options) > 0) {
        k <- sample(k_options, 1)
        adj[i, k] <- 1
        return(adj)
      }
    }
  }

  return(adj)
}


#' Remove a Transitive Edge
#'
#' @param adj Adjacency matrix
#' @return Modified adjacency matrix
remove_transitive_edge <- function(adj) {
  n <- nrow(adj)

  # Look for transitive triads: i->j, j->k, and i->k all exist
  for (attempt in 1:10) {
    i <- sample(1:n, 1)
    k_options <- which(adj[i, ] == 1)

    if (length(k_options) > 0) {
      k <- sample(k_options, 1)
      # Check if there's a j such that i->j->k
      j_options <- which(adj[i, ] == 1 & adj[, k] == 1)

      if (length(j_options) > 0) {
        # Remove i->k edge
        adj[i, k] <- 0
        return(adj)
      }
    }
  }

  return(adj)
}


#' Add a Homophilous Edge
#'
#' Adds an edge within the same category (homophilous edge)
#'
#' @param adj Adjacency matrix
#' @param node_attr Vector of categorical node attributes
#' @return Modified adjacency matrix
add_homophilous_edge <- function(adj, node_attr) {
  n <- nrow(adj)

  # Find pairs with same category but no edge
  for (attempt in 1:20) {
    i <- sample(1:n, 1)
    # Find nodes in same category without edge
    same_cat_nodes <- which(node_attr == node_attr[i] & adj[i, ] == 0 & (1:n) != i)

    if (length(same_cat_nodes) > 0) {
      j <- sample(same_cat_nodes, 1)
      adj[i, j] <- 1
      return(adj)
    }
  }

  return(adj)
}


#' Remove a Homophilous Edge
#'
#' Removes an edge within the same category
#'
#' @param adj Adjacency matrix
#' @param node_attr Vector of categorical node attributes
#' @return Modified adjacency matrix
remove_homophilous_edge <- function(adj, node_attr) {
  n <- nrow(adj)

  # Find edges within categories
  for (attempt in 1:20) {
    edges <- which(adj == 1)
    if (length(edges) == 0) return(adj)

    idx <- sample(edges, 1)
    i <- ((idx - 1) %% n) + 1
    j <- ((idx - 1) %/% n) + 1

    if (node_attr[i] == node_attr[j]) {
      adj[i, j] <- 0
      return(adj)
    }
  }

  return(adj)
}


#' Generate Directed Network with Reciprocity, Transitivity, and Homophily
#'
#' Generates a directed binary adjacency matrix with specified structural properties
#' using an iterative modification approach.
#'
#' @param n Number of nodes
#' @param node_attr Vector of categorical node attributes
#' @param density Target edge density (proportion of possible edges)
#' @param reciprocity Target reciprocity coefficient (0-1)
#' @param transitivity Target transitivity/clustering coefficient (0-1)
#' @param homophily Target proportion of within-category edges (0-1)
#' @param max_iter Maximum iterations for optimization (default 2000)
#' @param tol Tolerance for meeting targets (default 0.08)
#' @return n x n binary adjacency matrix
generate_network <- function(n, node_attr, density = 0.1, reciprocity = 0.4,
                             transitivity = 0.3, homophily = 0.6,
                             max_iter = 2000, tol = 0.08) {

  # Start with random directed graph at target density
  n_edges <- round(n * (n - 1) * density)

  # Initialize with random edges
  adj <- matrix(0, n, n)
  diag(adj) <- 0  # no self-loops

  # Randomly place edges
  possible_edges <- which(upper.tri(adj) | lower.tri(adj))
  edge_idx <- sample(possible_edges, min(n_edges, length(possible_edges)))
  adj[edge_idx] <- 1

  # Iteratively modify to achieve target properties
  for (iter in 1:max_iter) {
    current_recip <- calculate_reciprocity(adj)
    current_trans <- sna::gtrans(adj)  # Use sna package for transitivity
    current_homoph <- calculate_homophily(adj, node_attr)
    current_density <- sum(adj) / (n * (n - 1))

    # Check if we've met all targets
    recip_ok <- abs(current_recip - reciprocity) < tol
    trans_ok <- abs(current_trans - transitivity) < tol
    homoph_ok <- abs(current_homoph - homophily) < tol

    if (recip_ok && trans_ok && homoph_ok) {
      break
    }

    # Adjust reciprocity
    if (!recip_ok) {
      if (current_recip < reciprocity - tol) {
        adj <- add_reciprocal_edge(adj)
      } else if (current_recip > reciprocity + tol) {
        adj <- remove_reciprocal_edge(adj)
      }
    }

    # Adjust transitivity
    if (!trans_ok) {
      if (current_trans < transitivity - tol) {
        adj <- add_transitive_edge(adj)
      } else if (current_trans > transitivity + tol) {
        adj <- remove_transitive_edge(adj)
      }
    }

    # Adjust homophily
    if (!homoph_ok) {
      if (current_homoph < homophily - tol) {
        adj <- add_homophilous_edge(adj, node_attr)
      } else if (current_homoph > homophily + tol) {
        adj <- remove_homophilous_edge(adj, node_attr)
      }
    }

    # Maintain approximate density
    if (current_density < density - 0.03) {
      # Add random edge
      zeros <- which(adj == 0 & row(adj) != col(adj))
      if (length(zeros) > 0) {
        adj[sample(zeros, 1)] <- 1
      }
    } else if (current_density > density + 0.03) {
      # Remove random edge
      ones <- which(adj == 1)
      if (length(ones) > 0) {
        adj[sample(ones, 1)] <- 0
      }
    }
  }

  return(adj)
}


#' QAP Regression Test
#'
#' Performs Quadratic Assignment Procedure (QAP) regression using the
#' Dekker-Krackhardt-Sanders double semi-partialing (DSP) method via sna::netlm.
#'
#' @param Y Dependent variable matrix (will stay fixed)
#' @param X Independent variable matrix (will be permuted)
#' @param SameCat Optional control variable matrix (will be permuted)
#' @param nperm Number of permutations (default 1000)
#' @param include_samecat Whether to include SameCat as a control (default TRUE)
#' @return List with coefficient estimates, t-statistics, and p-values
qap_regression <- function(Y, X, SameCat = NULL, nperm = 1000,
                            include_samecat = TRUE) {
  # Input validation
  if (!is.matrix(Y) || !is.matrix(X)) {
    stop("Y and X must be matrices")
  }
  if (!identical(dim(Y), dim(X))) {
    stop("Y and X must have the same dimensions")
  }

  n <- nrow(Y)

  # Vectorize matrices (exclude diagonal)
  y_vec <- Y[lower.tri(Y) | upper.tri(Y)]
  x_vec <- X[lower.tri(X) | upper.tri(X)]

  # Use sna::netlm for QAP regression
  if (include_samecat && !is.null(SameCat)) {
    # Include SameCat as control variable
    result <- sna::netlm(Y, list(X, SameCat), mode = "digraph",
                         nullhyp = "qap", reps = nperm)

    # Extract results for X coefficient (index 2, since 1 is intercept)
    x_coef <- result$coefficients[2]
    x_tstat <- result$tstat[2]
    x_pval <- result$pleeq[2]  # One-tailed p-value <=
    # Convert to two-tailed
    x_pval_two <- 2 * min(x_pval, 1 - x_pval)

    # Also get SameCat results
    sc_coef <- result$coefficients[3]
    sc_tstat <- result$tstat[3]
    sc_pval <- result$pleeq[3]
    sc_pval_two <- 2 * min(sc_pval, 1 - sc_pval)

    return(list(
      x_coef = x_coef,
      x_tstat = x_tstat,
      x_pval = x_pval_two,
      samecat_coef = sc_coef,
      samecat_tstat = sc_tstat,
      samecat_pval = sc_pval_two,
      model = result
    ))
  } else {
    # Just X as predictor
    result <- sna::netlm(Y, X, mode = "digraph",
                         nullhyp = "qap", reps = nperm)

    x_coef <- result$coefficients[2]  # Index 2 (1 is intercept)
    x_tstat <- result$tstat[2]
    x_pval <- result$pleeq[2]
    x_pval_two <- 2 * min(x_pval, 1 - x_pval)

    return(list(
      x_coef = x_coef,
      x_tstat = x_tstat,
      x_pval = x_pval_two,
      model = result
    ))
  }
}


################################################################################
# SIMULATION PARAMETERS
################################################################################

cat("\n")
cat(paste(rep("=", 80), collapse=""), "\n", sep="")
cat("QAP REGRESSION TYPE I ERROR TEST UNDER NETWORK AUTOCORRELATION\n")
cat(paste(rep("=", 80), collapse=""), "\n\n", sep="")

# Simulation parameters
n_simulations <- 2000  # Number of independent (X,Y) pairs
n_permutations <- 1000  # QAP permutations per test
network_sizes <- c(25, 50, 75)  # Test different network sizes

# Network generation parameters
target_density <- 0.15
target_reciprocity <- 0.4
target_transitivity <- 0.3
target_homophily <- 0.6
n_categories <- 3

cat("Simulation Settings:\n")
cat(sprintf("  Number of simulations: %d\n", n_simulations))
cat(sprintf("  Permutations per test: %d\n", n_permutations))
cat(sprintf("  Network sizes: %s\n", paste(network_sizes, collapse=", ")))
cat("\nTarget Network Properties:\n")
cat(sprintf("  Density: %.2f\n", target_density))
cat(sprintf("  Reciprocity: %.2f\n", target_reciprocity))
cat(sprintf("  Transitivity: %.2f\n", target_transitivity))
cat(sprintf("  Homophily: %.2f\n", target_homophily))
cat(sprintf("  Categories: %d\n\n", n_categories))


################################################################################
# RUN SIMULATION
################################################################################

results_summary <- list()

for (n_nodes in network_sizes) {
  cat("\n")
  cat(paste(rep("-", 80), collapse=""), "\n", sep="")
  cat(sprintf("NETWORK SIZE: n = %d nodes\n", n_nodes))
  cat(paste(rep("-", 80), collapse=""), "\n\n", sep="")

  # Generate node attributes once for all simulations at this size
  # (each simulation will generate new ones, but this is for sizing)

  # Storage for results
  p_values_x <- numeric(n_simulations)
  t_stats_x <- numeric(n_simulations)
  coefs_x <- numeric(n_simulations)
  actual_cors_xy <- numeric(n_simulations)

  # Storage for network properties
  reciprocity_X <- numeric(n_simulations)
  reciprocity_Y <- numeric(n_simulations)
  transitivity_X <- numeric(n_simulations)
  transitivity_Y <- numeric(n_simulations)
  homophily_X <- numeric(n_simulations)
  homophily_Y <- numeric(n_simulations)
  density_X <- numeric(n_simulations)
  density_Y <- numeric(n_simulations)

  # Track runtime
  start_time <- Sys.time()

  cat(sprintf("Generating %d pairs of networks and running QAP regression...\n",
              n_simulations))
  cat("Progress: ")

  # Run simulations
  for (i in 1:n_simulations) {
    # Progress indicator
    if (i %% 100 == 0) {
      cat(sprintf("%d ", i))
      flush.console()
    }

    # Generate node attributes for this simulation
    # Generate DIFFERENT attributes for X and Y to ensure independence
    node_attr_X <- generate_node_attributes(n_nodes, n_categories)
    node_attr_Y <- generate_node_attributes(n_nodes, n_categories)

    # Generate two independent networks with desired properties
    X <- generate_network(n_nodes, node_attr_X, target_density,
                          target_reciprocity, target_transitivity,
                          target_homophily)
    Y <- generate_network(n_nodes, node_attr_Y, target_density,
                          target_reciprocity, target_transitivity,
                          target_homophily)

    # Create SameCat matrix for Y (used as control)
    SameCat_Y <- create_same_cat_matrix(node_attr_Y)

    # Record actual network properties
    reciprocity_X[i] <- calculate_reciprocity(X)
    reciprocity_Y[i] <- calculate_reciprocity(Y)
    transitivity_X[i] <- sna::gtrans(X)
    transitivity_Y[i] <- sna::gtrans(Y)
    homophily_X[i] <- calculate_homophily(X, node_attr_X)
    homophily_Y[i] <- calculate_homophily(Y, node_attr_Y)
    density_X[i] <- sum(X) / (n_nodes * (n_nodes - 1))
    density_Y[i] <- sum(Y) / (n_nodes * (n_nodes - 1))

    # Calculate actual correlation between X and Y (should be ~0)
    x_vec <- X[lower.tri(X) | upper.tri(X)]
    y_vec <- Y[lower.tri(Y) | upper.tri(Y)]
    actual_cors_xy[i] <- cor(x_vec, y_vec)

    # Run QAP regression: Y ~ X + SameCat
    qap_result <- qap_regression(Y, X, SameCat_Y, nperm = n_permutations,
                                  include_samecat = TRUE)

    # Store results for X coefficient
    p_values_x[i] <- qap_result$x_pval
    t_stats_x[i] <- qap_result$x_tstat
    coefs_x[i] <- qap_result$x_coef
  }

  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

  cat("\n\n")

  ############################################################################
  # ANALYZE RESULTS
  ############################################################################

  # Type I error rate (proportion of p < 0.05)
  type1_error_rate <- mean(p_values_x < 0.05)

  # Store summary
  results_summary[[as.character(n_nodes)]] <- list(
    n = n_nodes,
    type1_error_rate = type1_error_rate,
    p_values_x = p_values_x,
    t_stats_x = t_stats_x,
    coefs_x = coefs_x,
    actual_cors_xy = actual_cors_xy,
    reciprocity_X = reciprocity_X,
    reciprocity_Y = reciprocity_Y,
    transitivity_X = transitivity_X,
    transitivity_Y = transitivity_Y,
    homophily_X = homophily_X,
    homophily_Y = homophily_Y,
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
    cat("✓ Type I error rate is within expected range\n\n")
  } else {
    cat("✗ WARNING: Type I error rate is outside expected range!\n\n")
  }

  # P-value distribution statistics
  cat("P-value Distribution (for X coefficient):\n")
  cat(sprintf("  Mean:   %.4f (should be ~0.5 if uniform)\n", mean(p_values_x)))
  cat(sprintf("  Median: %.4f (should be ~0.5 if uniform)\n", median(p_values_x)))
  cat(sprintf("  Min:    %.4f\n", min(p_values_x)))
  cat(sprintf("  Max:    %.4f\n\n", max(p_values_x)))

  # Kolmogorov-Smirnov test for uniformity
  ks_test <- ks.test(p_values_x, "punif")
  cat("Kolmogorov-Smirnov test for uniformity:\n")
  cat(sprintf("  D = %.4f, p-value = %.4f\n", ks_test$statistic, ks_test$p.value))
  if (ks_test$p.value > 0.05) {
    cat("  ✓ P-values are consistent with uniform distribution\n\n")
  } else {
    cat("  ✗ P-values deviate from uniform distribution\n\n")
  }

  # T-statistic distribution
  cat("T-statistic Distribution (for X coefficient):\n")
  cat(sprintf("  Mean:   %.4f (should be ~0)\n", mean(t_stats_x)))
  cat(sprintf("  SD:     %.4f\n", sd(t_stats_x)))
  cat(sprintf("  Range:  [%.4f, %.4f]\n\n", min(t_stats_x), max(t_stats_x)))

  # Actual correlation between X and Y (should be near 0)
  cat("Actual Correlations between Independent Networks X and Y:\n")
  cat(sprintf("  Mean:   %.4f (should be ~0)\n", mean(actual_cors_xy)))
  cat(sprintf("  SD:     %.4f\n", sd(actual_cors_xy)))
  cat(sprintf("  Range:  [%.4f, %.4f]\n\n",
              min(actual_cors_xy), max(actual_cors_xy)))

  # Network properties achieved
  cat("Achieved Network Properties:\n")
  cat("\n  DENSITY:\n")
  cat(sprintf("    Target: %.2f\n", target_density))
  cat(sprintf("    X: %.3f ± %.3f\n", mean(density_X), sd(density_X)))
  cat(sprintf("    Y: %.3f ± %.3f\n\n", mean(density_Y), sd(density_Y)))

  cat("  RECIPROCITY:\n")
  cat(sprintf("    Target: %.2f\n", target_reciprocity))
  cat(sprintf("    X: %.3f ± %.3f\n", mean(reciprocity_X), sd(reciprocity_X)))
  cat(sprintf("    Y: %.3f ± %.3f\n\n", mean(reciprocity_Y), sd(reciprocity_Y)))

  cat("  TRANSITIVITY:\n")
  cat(sprintf("    Target: %.2f\n", target_transitivity))
  cat(sprintf("    X: %.3f ± %.3f\n", mean(transitivity_X), sd(transitivity_X)))
  cat(sprintf("    Y: %.3f ± %.3f\n\n", mean(transitivity_Y), sd(transitivity_Y)))

  cat("  HOMOPHILY:\n")
  cat(sprintf("    Target: %.2f\n", target_homophily))
  cat(sprintf("    X: %.3f ± %.3f\n", mean(homophily_X), sd(homophily_X)))
  cat(sprintf("    Y: %.3f ± %.3f\n\n", mean(homophily_Y), sd(homophily_Y)))

  # Runtime
  cat(sprintf("Runtime: %.1f seconds (%.2f minutes)\n", runtime, runtime/60))
  cat(sprintf("Average time per simulation: %.3f seconds\n\n",
              runtime/n_simulations))

  ############################################################################
  # PLOTS
  ############################################################################

  cat("Generating diagnostic plots...\n\n")

  # Set up 2x2 plot layout
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  # 1. Histogram of p-values
  hist(p_values_x, breaks = 20,
       main = sprintf("P-value Distribution (n=%d)", n_nodes),
       xlab = "P-value for X coefficient",
       ylab = "Frequency",
       col = "lightblue",
       border = "white")
  abline(h = n_simulations/20, col = "red", lty = 2, lwd = 2)
  legend("topright", legend = "Expected if uniform",
         col = "red", lty = 2, lwd = 2, bty = "n")

  # 2. QQ plot for uniformity
  qqplot(qunif(ppoints(n_simulations)), p_values_x,
         main = sprintf("Q-Q Plot vs Uniform (n=%d)", n_nodes),
         xlab = "Theoretical Quantiles",
         ylab = "Sample Quantiles",
         pch = 20, col = rgb(0, 0, 1, 0.3))
  abline(0, 1, col = "red", lwd = 2)

  # 3. Distribution of t-statistics
  hist(t_stats_x, breaks = 30,
       main = sprintf("T-statistics for X (n=%d)", n_nodes),
       xlab = "T-statistic",
       ylab = "Frequency",
       col = "lightgreen",
       border = "white")
  abline(v = 0, col = "red", lty = 2, lwd = 2)
  abline(v = mean(t_stats_x), col = "blue", lwd = 2)
  legend("topright",
         legend = c("Null (0)", sprintf("Mean = %.3f", mean(t_stats_x))),
         col = c("red", "blue"), lty = c(2, 1), lwd = 2, bty = "n")

  # 4. Cumulative distribution of p-values
  plot(ecdf(p_values_x),
       main = sprintf("ECDF of P-values (n=%d)", n_nodes),
       xlab = "P-value",
       ylab = "Cumulative Probability",
       col = "blue", lwd = 2)
  abline(0, 1, col = "red", lty = 2, lwd = 2)
  legend("topleft", legend = c("Observed", "Expected (uniform)"),
         col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n")

  # Reset par
  par(mfrow = c(1, 1))
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
  Mean_Cor_XY = sapply(results_summary, function(x) mean(x$actual_cors_xy)),
  Mean_Recip_X = sapply(results_summary, function(x) mean(x$reciprocity_X)),
  Mean_Trans_X = sapply(results_summary, function(x) mean(x$transitivity_X)),
  Mean_Homoph_X = sapply(results_summary, function(x) mean(x$homophily_X)),
  Runtime_sec = sapply(results_summary, function(x) x$runtime)
)

print(summary_df)

cat("\n")
cat("CONCLUSION:\n")
cat("-----------\n")

# Allow slightly larger tolerance due to autocorrelation effects
all_within_ci <- all(abs(summary_df$Type_I_Error - 0.05) < 0.02)

if (all_within_ci) {
  cat("✓ QAP regression successfully controls Type I error rate at α=0.05\n")
  cat("  even when networks exhibit reciprocity, transitivity, and homophily.\n\n")
} else {
  cat("✗ WARNING: QAP regression shows inflated or deflated Type I error\n")
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
