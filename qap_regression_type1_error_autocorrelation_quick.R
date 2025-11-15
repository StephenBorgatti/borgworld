################################################################################
# QUICK TEST VERSION - QAP Regression Type I Error Test
################################################################################
# This is a rapid test version with reduced parameters to verify functionality
# before running the full simulation.

# Set seed for reproducibility
set.seed(42)

# Load required packages
if (!require("sna", quietly = TRUE)) {
  stop("Package 'sna' is required. Install with: install.packages('sna')")
}

# Note: Functions are defined inline below (copied from main script)
# This avoids running the full simulation when sourcing

################################################################################
# HELPER FUNCTIONS
################################################################################

generate_node_attributes <- function(n, n_categories = 3) {
  categories <- rep(1:n_categories, length.out = n)
  sample(categories)
}

create_same_cat_matrix <- function(node_attr) {
  n <- length(node_attr)
  same_cat <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      same_cat[i, j] <- as.integer(node_attr[i] == node_attr[j])
    }
  }
  diag(same_cat) <- 0
  return(same_cat)
}

calculate_reciprocity <- function(adj) {
  recip_edges <- sum(adj * t(adj))
  total_edges <- sum(adj)
  if (total_edges == 0) return(0)
  return(recip_edges / total_edges)
}

calculate_homophily <- function(adj, node_attr) {
  same_cat <- create_same_cat_matrix(node_attr)
  within_edges <- sum(adj * same_cat)
  total_edges <- sum(adj)
  if (total_edges == 0) return(0)
  return(within_edges / total_edges)
}

add_reciprocal_edge <- function(adj) {
  non_recip <- which(adj == 1 & t(adj) == 0)
  if (length(non_recip) > 0) {
    idx <- sample(non_recip, 1)
    i <- ((idx - 1) %% nrow(adj)) + 1
    j <- ((idx - 1) %/% nrow(adj)) + 1
    adj[j, i] <- 1
  }
  return(adj)
}

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

add_transitive_edge <- function(adj) {
  n <- nrow(adj)
  for (attempt in 1:10) {
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

add_homophilous_edge <- function(adj, node_attr) {
  n <- nrow(adj)
  for (attempt in 1:20) {
    i <- sample(1:n, 1)
    same_cat_nodes <- which(node_attr == node_attr[i] & adj[i, ] == 0 & (1:n) != i)
    if (length(same_cat_nodes) > 0) {
      j <- sample(same_cat_nodes, 1)
      adj[i, j] <- 1
      return(adj)
    }
  }
  return(adj)
}

remove_homophilous_edge <- function(adj, node_attr) {
  n <- nrow(adj)
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

generate_network <- function(n, node_attr, density = 0.1, reciprocity = 0.4,
                             transitivity = 0.3, homophily = 0.6,
                             max_iter = 2000, tol = 0.08) {
  n_edges <- round(n * (n - 1) * density)
  adj <- matrix(0, n, n)
  diag(adj) <- 0
  possible_edges <- which(upper.tri(adj) | lower.tri(adj))
  edge_idx <- sample(possible_edges, min(n_edges, length(possible_edges)))
  adj[edge_idx] <- 1

  for (iter in 1:max_iter) {
    current_recip <- calculate_reciprocity(adj)
    current_trans <- sna::gtrans(adj)
    current_homoph <- calculate_homophily(adj, node_attr)
    current_density <- sum(adj) / (n * (n - 1))

    recip_ok <- abs(current_recip - reciprocity) < tol
    trans_ok <- abs(current_trans - transitivity) < tol
    homoph_ok <- abs(current_homoph - homophily) < tol

    if (recip_ok && trans_ok && homoph_ok) {
      break
    }

    if (!recip_ok) {
      if (current_recip < reciprocity - tol) {
        adj <- add_reciprocal_edge(adj)
      } else if (current_recip > reciprocity + tol) {
        adj <- remove_reciprocal_edge(adj)
      }
    }

    if (!trans_ok) {
      if (current_trans < transitivity - tol) {
        adj <- add_transitive_edge(adj)
      } else if (current_trans > transitivity + tol) {
        adj <- remove_transitive_edge(adj)
      }
    }

    if (!homoph_ok) {
      if (current_homoph < homophily - tol) {
        adj <- add_homophilous_edge(adj, node_attr)
      } else if (current_homoph > homophily + tol) {
        adj <- remove_homophilous_edge(adj, node_attr)
      }
    }

    if (current_density < density - 0.03) {
      zeros <- which(adj == 0 & row(adj) != col(adj))
      if (length(zeros) > 0) {
        adj[sample(zeros, 1)] <- 1
      }
    } else if (current_density > density + 0.03) {
      ones <- which(adj == 1)
      if (length(ones) > 0) {
        adj[sample(ones, 1)] <- 0
      }
    }
  }

  return(adj)
}

qap_regression <- function(Y, X, SameCat = NULL, nperm = 1000,
                            include_samecat = TRUE) {
  if (!is.matrix(Y) || !is.matrix(X)) {
    stop("Y and X must be matrices")
  }
  if (!identical(dim(Y), dim(X))) {
    stop("Y and X must have the same dimensions")
  }

  n <- nrow(Y)
  y_vec <- Y[lower.tri(Y) | upper.tri(Y)]
  x_vec <- X[lower.tri(X) | upper.tri(X)]

  if (include_samecat && !is.null(SameCat)) {
    result <- sna::netlm(Y, list(X, SameCat), mode = "digraph",
                         nullhyp = "qap", reps = nperm)
    x_coef <- result$coefficients[2]
    x_tstat <- result$tstat[2]
    x_pval <- result$pleeq[2]
    x_pval_two <- 2 * min(x_pval, 1 - x_pval)

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
    result <- sna::netlm(Y, X, mode = "digraph",
                         nullhyp = "qap", reps = nperm)
    x_coef <- result$coefficients[2]
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
# QUICK TEST
################################################################################

cat("\n")
cat(paste(rep("=", 80), collapse=""), "\n", sep="")
cat("QUICK TEST RUN (Reduced Parameters)\n")
cat(paste(rep("=", 80), collapse=""), "\n\n", sep="")

# Test parameters - much smaller for quick verification
n_simulations <- 50  # Just 50 tests for quick check
n_permutations <- 100  # Fewer permutations
network_sizes <- c(20)  # Just one small size

cat("Testing with:\n")
cat(sprintf("  - %d simulations\n", n_simulations))
cat(sprintf("  - %d permutations per test\n", n_permutations))
cat(sprintf("  - Network size: %d nodes\n\n", network_sizes[1]))

# Network parameters
target_density <- 0.15
target_reciprocity <- 0.4
target_transitivity <- 0.3
target_homophily <- 0.6
n_categories <- 3

# Run quick test
n_nodes <- network_sizes[1]
p_values_x <- numeric(n_simulations)
t_stats_x <- numeric(n_simulations)
actual_cors_xy <- numeric(n_simulations)

start_time <- Sys.time()
cat("Running tests: ")

for (i in 1:n_simulations) {
  if (i %% 10 == 0) cat(sprintf("%d ", i))

  # Generate node attributes
  node_attr_X <- generate_node_attributes(n_nodes, n_categories)
  node_attr_Y <- generate_node_attributes(n_nodes, n_categories)

  # Generate two independent networks
  X <- generate_network(n_nodes, node_attr_X, target_density,
                        target_reciprocity, target_transitivity,
                        target_homophily)
  Y <- generate_network(n_nodes, node_attr_Y, target_density,
                        target_reciprocity, target_transitivity,
                        target_homophily)

  # Create SameCat matrix
  SameCat_Y <- create_same_cat_matrix(node_attr_Y)

  # Calculate actual correlation
  x_vec <- X[lower.tri(X) | upper.tri(X)]
  y_vec <- Y[lower.tri(Y) | upper.tri(Y)]
  actual_cors_xy[i] <- cor(x_vec, y_vec)

  # Run QAP regression
  qap_result <- qap_regression(Y, X, SameCat_Y, nperm = n_permutations,
                                include_samecat = TRUE)
  p_values_x[i] <- qap_result$x_pval
  t_stats_x[i] <- qap_result$x_tstat
}

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat("\n\n")
cat("RESULTS:\n")
cat("--------\n")
cat(sprintf("Type I Error Rate (p < 0.05): %.3f\n", mean(p_values_x < 0.05)))
cat(sprintf("Mean actual correlation X-Y: %.4f (should be ~0)\n",
            mean(actual_cors_xy)))
cat(sprintf("Mean t-statistic: %.4f (should be ~0)\n", mean(t_stats_x)))
cat(sprintf("Runtime: %.1f seconds\n\n", runtime))

# Quick histogram
par(mfrow = c(1, 2))

hist(p_values_x, breaks = 10, main = "P-value Distribution",
     xlab = "P-value", col = "lightblue", border = "white")
abline(h = n_simulations/10, col = "red", lty = 2, lwd = 2)

hist(t_stats_x, breaks = 10, main = "T-statistic Distribution",
     xlab = "T-statistic", col = "lightgreen", border = "white")
abline(v = 0, col = "red", lty = 2, lwd = 2)

par(mfrow = c(1, 1))

cat("Quick test completed successfully!\n")
cat("  Run the full script (qap_regression_type1_error_autocorrelation.R)\n")
cat("  for complete analysis with 2000 simulations.\n\n")
