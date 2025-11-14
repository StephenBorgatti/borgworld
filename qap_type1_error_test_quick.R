################################################################################
# QUICK TEST VERSION - QAP Correlation Type I Error Test
################################################################################
# This is a rapid test version with reduced parameters to verify functionality

# Source the main script functions
#source("qap_type1_error_test.R", echo = FALSE)

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
target_density <- 0.1
target_reciprocity <- 0.4
target_transitivity <- 0.3

# Run quick test
n_nodes <- network_sizes[1]
p_values <- numeric(n_simulations)
actual_cors <- numeric(n_simulations)

start_time <- Sys.time()
cat("Running tests: ")

for (i in 1:n_simulations) {
  if (i %% 10 == 0) cat(sprintf("%d ", i))

  # Generate two independent networks
  X <- generate_network(n_nodes, target_density, target_reciprocity,
                        target_transitivity)
  Y <- generate_network(n_nodes, target_density, target_reciprocity,
                        target_transitivity)

  # Run QAP test
  qap_result <- qap_cor(X, Y, nperm = n_permutations)
  p_values[i] <- qap_result$p_value
  actual_cors[i] <- qap_result$obs_cor
}

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat("\n\n")
cat("RESULTS:\n")
cat("--------\n")
cat(sprintf("Type I Error Rate (p < 0.05): %.3f\n", mean(p_values < 0.05)))
cat(sprintf("Mean actual correlation: %.4f (should be ~0)\n", mean(actual_cors)))
cat(sprintf("Runtime: %.1f seconds\n\n", runtime))

cat("Quick test completed successfully!\n")
cat("  Run the full script (qap_type1_error_test.R) for complete analysis.\n\n")
