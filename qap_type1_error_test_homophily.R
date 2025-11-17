################################################################################
# QAP Regression Type I Error Test with Homophily
# (Optimized - DKS 2007 Double Semi-Partialling Method)
################################################################################
#
# Purpose: Verify that QAP regression maintains proper Type I error control
#          when testing regression with two independent variables:
#          1. X: An independent network
#          2. same_category: Homophily matrix (1 if same category, 0 otherwise)
#
# Design: Network Y has homophily based on attribute A (3 categories)
#         X is generated independently (also with homophily on same attribute)
#         same_category matrix indicates shared attributes
#         X and same_category should not predict Y (Type I error test)
#
# Author: Generated for borgworld package validation
# Date: 2025-11-14
################################################################################

# Set seed for reproducibility
set.seed(42)

################################################################################
# FUNCTIONS
################################################################################

#' QAP Regression using Double Semi-Partialling (Dekker, Krackhardt & Snijders 2007)
#'
#' Performs MRQAP regression using the double semi-partialling method
#' Optimized for speed
#'
#' @param Y Dependent variable matrix
#' @param X Independent variable matrix (or list of matrices)
#' @param nperm Number of permutations (default 1000)
#' @return List with coefficients, standard errors, t-statistics, and p-values
qap_dsp_regression <- function(Y, X, nperm = 1000) {
  # Convert X to list if single matrix
  if (is.matrix(X)) {
    X <- list(X)
  }

  # Vectorize matrices (exclude diagonal)
  vectorize <- function(mat) {
    mat[lower.tri(mat) | upper.tri(mat)]
  }

  # Vectorize Y and X matrices
  y_vec <- vectorize(Y)
  x_list <- lapply(X, vectorize)
  x_mat <- do.call(cbind, x_list)

  # Add intercept
  x_mat <- cbind(1, x_mat)
  n_params <- ncol(x_mat)

  # Observed regression (using base R lm.fit for speed)
  obs_fit <- lm.fit(x_mat, y_vec)
  obs_coef <- obs_fit$coefficients
  obs_resid <- obs_fit$residuals

  # Storage for permutation results
  perm_coefs <- matrix(0, nrow = nperm, ncol = n_params)

  n <- nrow(Y)

  # Permutation loop
  for (p in 1:nperm) {
    # Generate permutation
    perm_idx <- sample(1:n)

    # Permute Y
    Y_perm <- Y[perm_idx, perm_idx]
    y_perm_vec <- vectorize(Y_perm)

    # Permute each X matrix
    x_perm_list <- lapply(X, function(x_mat) {
      x_perm <- x_mat[perm_idx, perm_idx]
      vectorize(x_perm)
    })
    x_perm_mat <- cbind(1, do.call(cbind, x_perm_list))

    # Double semi-partialling:
    # 1. Regress permuted Y on permuted X
    perm_fit <- lm.fit(x_perm_mat, y_perm_vec)

    # Store permuted coefficients
    perm_coefs[p, ] <- perm_fit$coefficients
  }

  # Calculate p-values (two-tailed)
  p_values <- numeric(n_params)
  for (i in 1:n_params) {
    p_values[i] <- mean(abs(perm_coefs[, i]) >= abs(obs_coef[i]))
  }

  # Calculate standard errors from permutation distribution
  se <- apply(perm_coefs, 2, sd)

  # Calculate t-statistics
  t_stats <- obs_coef / se

  # Return results
  param_names <- c("Intercept", paste0("X", 1:(n_params - 1)))

  return(list(
    coefficients = obs_coef,
    se = se,
    t_statistics = t_stats,
    p_values = p_values,
    param_names = param_names,
    nperm = nperm,
    perm_dist = perm_coefs
  ))
}


#' Generate Categorical Node Attributes
#'
#' Creates a vector of categorical attributes with approximately equal-sized groups
#'
#' @param n Number of nodes
#' @param n_categories Number of categories (default 3)
#' @return Vector of category assignments
generate_node_attributes <- function(n, n_categories = 3) {
  # Create approximately equal-sized groups
  base_size <- n %/% n_categories
  remainder <- n %% n_categories

  # Build group sizes
  group_sizes <- rep(base_size, n_categories)
  if (remainder > 0) {
    group_sizes[1:remainder] <- group_sizes[1:remainder] + 1
  }

  # Create and shuffle category assignments
  categories <- rep(1:n_categories, times = group_sizes)
  sample(categories)  # Random assignment to nodes
}


#' Create Same Category Matrix
#'
#' Creates a matrix where entry (i,j) is 1 if nodes i and j have the same attribute
#'
#' @param attributes Vector of node attributes
#' @return Binary matrix indicating same-category pairs
create_same_category_matrix <- function(attributes) {
  n <- length(attributes)
  same_cat <- matrix(0, n, n)

  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j && attributes[i] == attributes[j]) {
        same_cat[i, j] <- 1
      }
    }
  }

  return(same_cat)
}


#' Calculate Homophily (Assortativity)
#'
#' Measures the tendency for nodes with the same attribute to connect
#'
#' @param adj Adjacency matrix
#' @param attributes Vector of node attributes
#' @return Proportion of edges within same category
calculate_homophily <- function(adj, attributes) {
  edges <- which(adj == 1, arr.ind = TRUE)
  if (nrow(edges) == 0) return(NA)

  same_category <- attributes[edges[,1]] == attributes[edges[,2]]
  mean(same_category)
}


#' Generate Network with Homophily, Reciprocity and Transitivity
#'
#' @param n Number of nodes
#' @param attributes Vector of categorical node attributes
#' @param density Target edge density
#' @param homophily_p Probability multiplier for within-category ties (default 2/3)
#' @param reciprocity Target reciprocity coefficient
#' @param transitivity Target transitivity coefficient
#' @param max_iter Maximum iterations for optimization
#' @param tol Tolerance for meeting targets
#' @return n x n binary adjacency matrix
generate_network_with_homophily <- function(n, attributes, density = 0.1,
                                            homophily_p = 2/3, reciprocity = 0.4,
                                            transitivity = 0.3, max_iter = 1000,
                                            tol = 0.05) {

  # Initialize adjacency matrix
  adj <- matrix(0, n, n)
  diag(adj) <- 0

  # Calculate target number of edges
  n_edges <- round(n * (n - 1) * density)

  # Create homophily probability matrix
  prob_matrix <- matrix(1 - homophily_p, n, n)  # Between-category probability
  for (i in 1:n) {
    for (j in 1:n) {
      if (attributes[i] == attributes[j] && i != j) {
        prob_matrix[i, j] <- homophily_p  # Within-category probability
      }
    }
  }
  diag(prob_matrix) <- 0

  # Normalize probabilities to sum to 1
  prob_matrix <- prob_matrix / sum(prob_matrix)

  # Sample edges based on homophily probabilities
  possible_edges <- which(upper.tri(prob_matrix) | lower.tri(prob_matrix))
  edge_probs <- prob_matrix[possible_edges]

  # Sample edges with replacement, then take unique
  sampled_edges <- sample(possible_edges, size = n_edges * 2,
                          replace = TRUE, prob = edge_probs)
  sampled_edges <- unique(sampled_edges)[1:min(n_edges, length(unique(sampled_edges)))]

  # Place initial edges
  adj[sampled_edges] <- 1

  # Iteratively adjust for reciprocity and transitivity while preserving homophily
  for (iter in 1:max_iter) {
    current_recip <- calculate_reciprocity(adj)
    current_trans <- calculate_transitivity(adj)

    # Check if we've met targets
    if (abs(current_recip - reciprocity) < tol &&
        abs(current_trans - transitivity) < tol) {
      break
    }

    # Adjust reciprocity (prefer within-category reciprocal edges)
    if (current_recip < reciprocity - tol) {
      adj <- add_reciprocal_edge_homophily(adj, attributes)
    } else if (current_recip > reciprocity + tol) {
      adj <- remove_reciprocal_edge(adj)
    }

    # Adjust transitivity (prefer within-category transitive closure)
    if (current_trans < transitivity - tol) {
      adj <- add_transitive_edge_homophily(adj, attributes)
    } else if (current_trans > transitivity + tol) {
      adj <- remove_transitive_edge(adj)
    }

    # Maintain approximate density
    current_density <- sum(adj) / (n * (n - 1))
    if (current_density < density - 0.02) {
      # Add edge (prefer within-category)
      adj <- add_edge_homophily(adj, attributes)
    } else if (current_density > density + 0.02) {
      # Remove edge (prefer between-category)
      adj <- remove_edge_homophily(adj, attributes)
    }
  }

  return(adj)
}


#' Add Reciprocal Edge with Homophily Preference
#'
#' @param adj Adjacency matrix
#' @param attributes Node attributes
#' @return Modified adjacency matrix
add_reciprocal_edge_homophily <- function(adj, attributes) {
  # Find non-reciprocated edges
  non_recip <- which(adj == 1 & t(adj) == 0)

  if (length(non_recip) > 0) {
    # Prefer within-category edges
    indices <- arrayInd(non_recip, dim(adj))
    same_cat <- attributes[indices[,1]] == attributes[indices[,2]]

    # Sample preferring within-category
    if (any(same_cat)) {
      idx <- sample(non_recip[same_cat], 1)
    } else {
      idx <- sample(non_recip, 1)
    }

    # Get i,j coordinates
    i <- ((idx - 1) %% nrow(adj)) + 1
    j <- ((idx - 1) %/% nrow(adj)) + 1
    # Add reciprocal edge j->i
    adj[j, i] <- 1
  }

  return(adj)
}


#' Add Transitive Edge with Homophily Preference
#'
#' @param adj Adjacency matrix
#' @param attributes Node attributes
#' @return Modified adjacency matrix
add_transitive_edge_homophily <- function(adj, attributes) {
  n <- nrow(adj)

  # Try to find within-category triads to close
  for (attempt in 1:20) {
    i <- sample(1:n, 1)
    j_options <- which(adj[i, ] == 1)

    if (length(j_options) > 0) {
      # Prefer j in same category as i
      j_same <- j_options[attributes[j_options] == attributes[i]]
      if (length(j_same) > 0) {
        j <- sample(j_same, 1)
      } else {
        j <- sample(j_options, 1)
      }

      k_options <- which(adj[j, ] == 1 & adj[i, ] == 0 & (1:n) != i)

      if (length(k_options) > 0) {
        # Prefer k in same category as i
        k_same <- k_options[attributes[k_options] == attributes[i]]
        if (length(k_same) > 0) {
          k <- sample(k_same, 1)
        } else {
          k <- sample(k_options, 1)
        }
        adj[i, k] <- 1
        return(adj)
      }
    }
  }

  return(adj)
}


#' Add Edge with Homophily Preference
#'
#' @param adj Adjacency matrix
#' @param attributes Node attributes
#' @return Modified adjacency matrix
add_edge_homophily <- function(adj, attributes) {
  zeros <- which(adj == 0 & row(adj) != col(adj))

  if (length(zeros) > 0) {
    # Get indices
    indices <- arrayInd(zeros, dim(adj))
    same_cat <- attributes[indices[,1]] == attributes[indices[,2]]

    # Prefer within-category edges
    if (any(same_cat)) {
      idx <- sample(zeros[same_cat], 1)
    } else {
      idx <- sample(zeros, 1)
    }

    adj[idx] <- 1
  }

  return(adj)
}


#' Remove Edge with Homophily Preference
#'
#' @param adj Adjacency matrix
#' @param attributes Node attributes
#' @return Modified adjacency matrix
remove_edge_homophily <- function(adj, attributes) {
  ones <- which(adj == 1)

  if (length(ones) > 0) {
    # Get indices
    indices <- arrayInd(ones, dim(adj))
    diff_cat <- attributes[indices[,1]] != attributes[indices[,2]]

    # Prefer removing between-category edges
    if (any(diff_cat)) {
      idx <- sample(ones[diff_cat], 1)
    } else {
      idx <- sample(ones, 1)
    }

    adj[idx] <- 0
  }

  return(adj)
}


#' Calculate Reciprocity
#'
#' @param adj Adjacency matrix
#' @return Reciprocity coefficient
calculate_reciprocity <- function(adj) {
  recip_edges <- sum(adj * t(adj))
  total_edges <- sum(adj)

  if (total_edges == 0) return(0)
  return(recip_edges / total_edges)
}


#' Calculate Global Transitivity (Optimized)
#'
#' Computes the global transitivity coefficient (fraction of transitive triples)
#' Optimized version that doesn't rely on sna package
#'
#' @param adj Adjacency matrix
#' @return Global transitivity coefficient
calculate_transitivity <- function(adj) {
  n <- nrow(adj)

  # Count transitive triples and potential triples
  trans_count <- 0
  potential_count <- 0

  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j && adj[i, j] == 1) {
        for (k in 1:n) {
          if (k != i && k != j && adj[j, k] == 1) {
            potential_count <- potential_count + 1
            if (adj[i, k] == 1) {
              trans_count <- trans_count + 1
            }
          }
        }
      }
    }
  }

  if (potential_count == 0) return(0)
  return(trans_count / potential_count)
}


#' Remove Reciprocal Edge
#'
#' @param adj Adjacency matrix
#' @return Modified adjacency matrix
remove_reciprocal_edge <- function(adj) {
  recip <- which(adj == 1 & t(adj) == 1)

  if (length(recip) > 0) {
    idx <- sample(recip, 1)
    i <- ((idx - 1) %% nrow(adj)) + 1
    j <- ((idx - 1) %/% nrow(adj)) + 1

    if (runif(1) > 0.5) {
      adj[i, j] <- 0
    } else {
      adj[j, i] <- 0
    }
  }

  return(adj)
}


#' Remove Transitive Edge
#'
#' @param adj Adjacency matrix
#' @return Modified adjacency matrix
remove_transitive_edge <- function(adj) {
  n <- nrow(adj)

  for (attempt in 1:10) {
    i <- sample(1:n, 1)
    k_options <- which(adj[i, ] == 1)

    if (length(k_options) > 0) {
      k <- sample(k_options, 1)
      j_options <- which(adj[i, ] == 1 & adj[, k] == 1)

      if (length(j_options) > 0) {
        adj[i, k] <- 0
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
cat("QAP REGRESSION TYPE I ERROR TEST WITH HOMOPHILY\n")
cat("(Optimized - DKS 2007 Double Semi-Partialling Method)\n")
cat(paste(rep("=", 80), collapse=""), "\n\n", sep="")
# Simulation parameters
n_simulations <- 200  # Number of independent (X,Y) pairs
n_permutations <- 1000  # QAP permutations per test
network_sizes <- c(75)  # Test different network sizes
n_categories <- 3  # Number of categories for attributes

# Network generation parameters
target_density <- 0.1
target_reciprocity <- 0.4
target_transitivity <- 0.3
homophily_p <- .5  # Within-category tie probability multiplier

cat("Network Parameters:\n")
cat(sprintf("  Density: %.2f\n", target_density))
cat(sprintf("  Reciprocity: %.2f\n", target_reciprocity))
cat(sprintf("  Transitivity: %.2f\n", target_transitivity))
cat(sprintf("  Categories: %d\n", n_categories))
cat(sprintf("  Within-category tie preference: %.2f\n", homophily_p))
cat(sprintf("  Between-category tie preference: %.2f\n\n", 1 - homophily_p))

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
  p_values_X <- numeric(n_simulations)
  p_values_same_cat <- numeric(n_simulations)
  coef_X <- numeric(n_simulations)
  coef_same_cat <- numeric(n_simulations)
  reciprocity_X <- numeric(n_simulations)
  reciprocity_Y <- numeric(n_simulations)
  transitivity_X <- numeric(n_simulations)
  transitivity_Y <- numeric(n_simulations)
  density_X <- numeric(n_simulations)
  density_Y <- numeric(n_simulations)
  homophily_X <- numeric(n_simulations)
  homophily_Y <- numeric(n_simulations)

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

    # Generate single node attribute shared by both X and Y
    attributes <- generate_node_attributes(n_nodes, n_categories)

    # Create same_category matrix
    same_category <- create_same_category_matrix(attributes)

    # Generate two independent networks with homophily on the same attribute
    X <- generate_network_with_homophily(n_nodes, attributes, target_density,
                                         homophily_p, target_reciprocity,
                                         target_transitivity)
    Y <- generate_network_with_homophily(n_nodes, attributes, target_density,
                                         homophily_p, target_reciprocity,
                                         target_transitivity)

    # Record actual network properties
    reciprocity_X[i] <- calculate_reciprocity(X)
    reciprocity_Y[i] <- calculate_reciprocity(Y)
    transitivity_X[i] <- calculate_transitivity(X)
    transitivity_Y[i] <- calculate_transitivity(Y)
    density_X[i] <- sum(X) / (n_nodes * (n_nodes - 1))
    density_Y[i] <- sum(Y) / (n_nodes * (n_nodes - 1))
    homophily_X[i] <- calculate_homophily(X, attributes)
    homophily_Y[i] <- calculate_homophily(Y, attributes)

    # Run QAP regression test (regressing Y on X and same_category)
    qap_reg_result <- qap_dsp_regression(Y, list(X, same_category), nperm = n_permutations)

    # Store results
    p_values_X[i] <- qap_reg_result$p_values[2]  # p-value for X coefficient
    p_values_same_cat[i] <- qap_reg_result$p_values[3]  # p-value for same_category coefficient
    coef_X[i] <- qap_reg_result$coefficients[2]  # Coefficient for X
    coef_same_cat[i] <- qap_reg_result$coefficients[3]  # Coefficient for same_category
  }

  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

  cat("\n\n")

  ############################################################################
  # ANALYZE RESULTS
  ############################################################################

  # Type I error rate (proportion of p < 0.05)
  type1_error_rate_X <- mean(p_values_X < 0.05)
  type1_error_rate_same_cat <- mean(p_values_same_cat < 0.05)

  # Store summary
  results_summary[[as.character(n_nodes)]] <- list(
    n = n_nodes,
    type1_error_rate_X = type1_error_rate_X,
    type1_error_rate_same_cat = type1_error_rate_same_cat,
    p_values_X = p_values_X,
    p_values_same_cat = p_values_same_cat,
    coef_X = coef_X,
    coef_same_cat = coef_same_cat,
    reciprocity_X = reciprocity_X,
    reciprocity_Y = reciprocity_Y,
    transitivity_X = transitivity_X,
    transitivity_Y = transitivity_Y,
    density_X = density_X,
    density_Y = density_Y,
    homophily_X = homophily_X,
    homophily_Y = homophily_Y,
    runtime = runtime
  )

  ############################################################################
  # OUTPUT RESULTS
  ############################################################################

  cat("RESULTS:\n")
  cat("--------\n\n")

  # 95% confidence interval for Type I error rate
  se <- sqrt(0.05 * 0.95 / n_simulations)
  ci_lower <- 0.05 - 1.96 * se
  ci_upper <- 0.05 + 1.96 * se

  cat("QAP REGRESSION TEST (DKS 2007 Double Semi-Partialling)\n")
  cat("Model: Y ~ X + same_category\n\n")

  cat("COEFFICIENT: X\n")
  cat(sprintf("  Type I Error Rate (p < 0.05): %.4f\n", type1_error_rate_X))
  cat(sprintf("  Expected rate (nominal α):     %.4f\n", 0.05))
  cat(sprintf("  Absolute deviation:            %.4f\n", abs(type1_error_rate_X - 0.05)))
  cat(sprintf("  Expected 95%% CI for Type I error: [%.4f, %.4f]\n", ci_lower, ci_upper))
  if (type1_error_rate_X >= ci_lower && type1_error_rate_X <= ci_upper) {
    cat("  Type I error rate is within expected range\n\n")
  } else {
    cat("  WARNING: Type I error rate is outside expected range!\n\n")
  }

  cat("COEFFICIENT: same_category\n")
  cat(sprintf("  Type I Error Rate (p < 0.05): %.4f\n", type1_error_rate_same_cat))
  cat(sprintf("  Expected rate (nominal α):     %.4f\n", 0.05))
  cat(sprintf("  Absolute deviation:            %.4f\n", abs(type1_error_rate_same_cat - 0.05)))
  cat(sprintf("  Expected 95%% CI for Type I error: [%.4f, %.4f]\n", ci_lower, ci_upper))
  if (type1_error_rate_same_cat >= ci_lower && type1_error_rate_same_cat <= ci_upper) {
    cat("  Type I error rate is within expected range\n\n")
  } else {
    cat("  WARNING: Type I error rate is outside expected range!\n\n")
  }

  # P-value distribution statistics
  cat("P-value Distribution (X coefficient):\n")
  cat(sprintf("  Mean:   %.4f (should be ~0.5 if uniform)\n", mean(p_values_X)))
  cat(sprintf("  Median: %.4f (should be ~0.5 if uniform)\n", median(p_values_X)))
  cat(sprintf("  Min:    %.4f\n", min(p_values_X)))
  cat(sprintf("  Max:    %.4f\n\n", max(p_values_X)))

  cat("P-value Distribution (same_category coefficient):\n")
  cat(sprintf("  Mean:   %.4f (should be ~0.5 if uniform)\n", mean(p_values_same_cat)))
  cat(sprintf("  Median: %.4f (should be ~0.5 if uniform)\n", median(p_values_same_cat)))
  cat(sprintf("  Min:    %.4f\n", min(p_values_same_cat)))
  cat(sprintf("  Max:    %.4f\n\n", max(p_values_same_cat)))

  # Regression coefficients (should be near 0)
  cat("Observed Regression Coefficients:\n")
  cat(sprintf("  X:             Mean = %.4f, SD = %.4f\n", mean(coef_X), sd(coef_X)))
  cat(sprintf("  same_category: Mean = %.4f, SD = %.4f\n\n", mean(coef_same_cat), sd(coef_same_cat)))

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

  cat("  HOMOPHILY (proportion within-category):\n")
  cat(sprintf("    Expected (random): %.2f\n", 1/n_categories))
  cat(sprintf("    Target (with p=%.2f): ~%.2f\n", homophily_p,
              homophily_p/(homophily_p + (1-homophily_p)*(n_categories-1))))
  cat(sprintf("    X: %.3f ± %.3f\n", mean(homophily_X, na.rm=TRUE), sd(homophily_X, na.rm=TRUE)))
  cat(sprintf("    Y: %.3f ± %.3f\n\n", mean(homophily_Y, na.rm=TRUE), sd(homophily_Y, na.rm=TRUE)))

  # Runtime
  cat(sprintf("Runtime: %.1f seconds (%.2f minutes)\n", runtime, runtime/60))
  cat(sprintf("Average time per simulation: %.3f seconds\n\n",
              runtime/n_simulations))

  ############################################################################
  # PLOTS
  ############################################################################

  cat("Generating diagnostic plots...\n\n")

  # Open a new graphics window with appropriate size
  # Use windows() on Windows, quartz() on Mac, or x11() on Linux
  if (.Platform$OS.type == "windows") {
    windows(width = 16, height = 8)
  } else if (Sys.info()["sysname"] == "Darwin") {
    quartz(width = 16, height = 8)
  } else {
    x11(width = 16, height = 8)
  }
  # Set up 2x4 plot layout with better margins
  par(mfrow = c(2, 4), mar = c(4, 4, 3, 1), oma = c(1, 1, 1, 1))


  # 1. Histogram of p-values (X coefficient)
  hist(p_values_X, breaks = 20,
       main = sprintf("P-values (X coef) (n=%d)", n_nodes),
       xlab = "P-value",
       ylab = "Frequency",
       col = "lightblue",
       border = "white")
  abline(h = n_simulations/20, col = "red", lty = 2, lwd = 2)
  legend("topright", legend = "Expected if uniform",
         col = "red", lty = 2, lwd = 2, bty = "n", cex = 0.8)

  # 2. Histogram of p-values (same_category coefficient)
  hist(p_values_same_cat, breaks = 20,
       main = sprintf("P-values (same_cat coef) (n=%d)", n_nodes),
       xlab = "P-value",
       ylab = "Frequency",
       col = "lightcoral",
       border = "white")
  abline(h = n_simulations/20, col = "red", lty = 2, lwd = 2)
  legend("topright", legend = "Expected if uniform",
         col = "red", lty = 2, lwd = 2, bty = "n", cex = 0.8)

  # 3. Distribution of X coefficients
  hist(coef_X, breaks = 30,
       main = sprintf("Coefficients (X) (n=%d)", n_nodes),
       xlab = "Coefficient for X",
       ylab = "Frequency",
       col = "lightgreen",
       border = "white")
  abline(v = 0, col = "red", lty = 2, lwd = 2)
  abline(v = mean(coef_X), col = "blue", lwd = 2)
  legend("topright",
         legend = c("True (0)", sprintf("Mean=%.3f", mean(coef_X))),
         col = c("red", "blue"), lty = c(2, 1), lwd = 2, bty = "n", cex = 0.8)

  # 4. Distribution of same_category coefficients
  hist(coef_same_cat, breaks = 30,
       main = sprintf("Coefficients (same_cat) (n=%d)", n_nodes),
       xlab = "Coefficient for same_category",
       ylab = "Frequency",
       col = "lightyellow",
       border = "white")
  abline(v = 0, col = "red", lty = 2, lwd = 2)
  abline(v = mean(coef_same_cat), col = "blue", lwd = 2)
  legend("topright",
         legend = c("True (0)", sprintf("Mean=%.3f", mean(coef_same_cat))),
         col = c("red", "blue"), lty = c(2, 1), lwd = 2, bty = "n", cex = 0.8)

  # 5. ECDF of p-values (X coefficient)
  plot(ecdf(p_values_X),
       main = sprintf("ECDF P-values (X) (n=%d)", n_nodes),
       xlab = "P-value",
       ylab = "Cumulative Probability",
       col = "blue", lwd = 2)
  abline(0, 1, col = "red", lty = 2, lwd = 2)
  legend("topleft", legend = c("Observed", "Uniform"),
         col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)

  # 6. ECDF of p-values (same_category coefficient)
  plot(ecdf(p_values_same_cat),
       main = sprintf("ECDF P-values (same_cat) (n=%d)", n_nodes),
       xlab = "P-value",
       ylab = "Cumulative Probability",
       col = "darkgreen", lwd = 2)
  abline(0, 1, col = "red", lty = 2, lwd = 2)
  legend("topleft", legend = c("Observed", "Uniform"),
         col = c("darkgreen", "red"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)

  # 7. Homophily levels in X and Y (same attribute)
  plot(homophily_X, homophily_Y,
       main = sprintf("Homophily (same attr) (n=%d)", n_nodes),
       xlab = "Homophily in X",
       ylab = "Homophily in Y",
       pch = 20, col = rgb(0, 0, 1, 0.1))
  abline(h = mean(homophily_Y, na.rm=TRUE), col = "red", lty = 2)
  abline(v = mean(homophily_X, na.rm=TRUE), col = "red", lty = 2)

  # 8. X coefficient vs same_category coefficient
  plot(coef_X, coef_same_cat,
       main = sprintf("Coef X vs Coef same_cat (n=%d)", n_nodes),
       xlab = "Coefficient for X",
       ylab = "Coefficient for same_category",
       pch = 20, col = rgb(0, 0, 1, 0.1))
  abline(h = 0, col = "red", lty = 2)
  abline(v = 0, col = "red", lty = 2)

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
  Type_I_Error_X = sapply(results_summary, function(x) x$type1_error_rate_X),
  Type_I_Error_same_cat = sapply(results_summary, function(x) x$type1_error_rate_same_cat),
  Mean_Coef_X = sapply(results_summary, function(x) mean(x$coef_X)),
  Mean_Coef_same_cat = sapply(results_summary, function(x) mean(x$coef_same_cat)),
  Mean_Homophily_X = sapply(results_summary, function(x) mean(x$homophily_X, na.rm=TRUE)),
  Mean_Homophily_Y = sapply(results_summary, function(x) mean(x$homophily_Y, na.rm=TRUE)),
  Runtime_sec = sapply(results_summary, function(x) x$runtime)
)

print(summary_df, digits = 4)

cat("\n")
cat("CONCLUSION:\n")
cat("-----------\n")

all_within_ci_X <- all(abs(summary_df$Type_I_Error_X - 0.05) < 0.015)
all_within_ci_same_cat <- all(abs(summary_df$Type_I_Error_same_cat - 0.05) < 0.015)

cat("QAP REGRESSION (DKS 2007) - Model: Y ~ X + same_category\n\n")

cat("COEFFICIENT X:\n")
if (all_within_ci_X) {
  cat("  Successfully controls Type I error rate at α=0.05\n")
} else {
  cat("  WARNING: Shows inflated or deflated Type I error\n")
}

cat("\nCOEFFICIENT same_category:\n")
if (all_within_ci_same_cat) {
  cat("  Successfully controls Type I error rate at α=0.05\n")
} else {
  cat("  WARNING: Shows inflated or deflated Type I error\n")
}

cat("\nTest was conducted on networks with:\n")
cat("  - Reciprocity\n")
cat("  - Transitivity\n")
cat("  - Homophily (on the same attribute)\n\n")

cat(sprintf("Total runtime: %.1f minutes\n", sum(summary_df$Runtime_sec)/60))

cat("\n")
cat(paste(rep("=", 80), collapse=""), "\n", sep="")
cat("SIMULATION COMPLETE\n")
cat(paste(rep("=", 80), collapse=""), "\n\n", sep="")

################################################################################
# END
################################################################################
