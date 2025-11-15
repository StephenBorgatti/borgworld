# QAP Regression Type I Error Test under Network Autocorrelation

## Overview

This script tests whether QAP (Quadratic Assignment Procedure) regression maintains proper Type I error control (5% false positive rate at α=0.05) when testing regression coefficients in the presence of network autocorrelation.

## Background

QAP regression is a permutation-based method for testing hypotheses about relationships between network matrices (adjacency matrices). This test verifies that the method correctly controls Type I error rates when:

1. **X and Y are independent** (true correlation = 0)
2. Both networks exhibit **autocorrelation** through:
   - **Reciprocity**: Mutual ties (if i→j exists, j→i is more likely)
   - **Transitivity**: Triadic closure (if i→j and j→k exist, i→k is more likely)
   - **Homophily**: Within-category ties (nodes with same attribute are more likely to be connected)

## Files

- **`qap_regression_type1_error_autocorrelation.R`**: Full simulation with 2000 iterations
- **`qap_regression_type1_error_autocorrelation_quick.R`**: Quick test with 50 iterations

## Requirements

```r
install.packages("sna")
```

## Usage

### Quick Test (~ 2-3 minutes)

```r
source("qap_regression_type1_error_autocorrelation_quick.R")
```

This runs a reduced version with:
- 50 simulations
- 100 permutations per test
- Network size: 20 nodes

### Full Simulation (~ 30-60 minutes)

```r
source("qap_regression_type1_error_autocorrelation.R")
```

This runs the complete analysis with:
- 2000 simulations
- 1000 permutations per test
- Network sizes: 10, 30, 50 nodes

## Method

### Network Generation

For each simulation, the script generates two **independent** directed binary adjacency matrices (X and Y) with:

1. **Categorical node attributes**: Nodes assigned to one of 3 categories
2. **Target properties**:
   - Density: 0.15 (15% of possible edges)
   - Reciprocity: 0.40 (40% of edges are reciprocated)
   - Transitivity: 0.30 (30% transitivity coefficient)
   - Homophily: 0.60 (60% of edges are within-category)

The generation uses an iterative algorithm that:
- Starts with a random graph at target density
- Iteratively adds/removes reciprocal, transitive, and homophilous edges
- Stops when all targets are met within tolerance (±0.08)

### QAP Regression

For each (X, Y) pair:

1. Create **SameCat** dyadic matrix: `SameCat[i,j] = 1` if nodes i and j share the same category
2. Run QAP regression: `Y ~ X + SameCat`
   - Uses Dekker-Krackhardt-Sanders (DSP) method via `sna::netlm()`
   - 1000 permutations (or 100 for quick test)
3. Extract p-value for **X coefficient's t-statistic**
4. Record whether p < 0.05 (false positive)

### Null Hypothesis

Since X and Y are generated independently:
- True correlation between X and Y is 0
- True regression coefficient for X should be 0
- Any p < 0.05 is a **Type I error** (false positive)

## Expected Results

If QAP regression properly controls Type I error:

1. **Type I error rate** ≈ 0.05 (within 95% CI: [0.040, 0.060] for 2000 simulations)
2. **P-value distribution** should be approximately uniform [0, 1]
3. **T-statistics** should be centered near 0
4. **Actual X-Y correlation** should be near 0 (confirms independence)

## Output

The script provides:

### Summary Statistics

- Type I error rate and 95% confidence interval
- P-value distribution statistics (mean, median, range)
- Kolmogorov-Smirnov test for uniformity
- T-statistic distribution
- Achieved network properties (density, reciprocity, transitivity, homophily)
- Runtime information

### Diagnostic Plots

1. **Histogram of p-values**: Should be approximately flat
2. **Q-Q plot**: Should follow diagonal line if uniform
3. **T-statistic distribution**: Should be centered at 0
4. **ECDF of p-values**: Should follow diagonal if uniform

## Interpretation

### Success Criteria

✓ Type I error rate within [0.04, 0.06]
✓ P-values approximately uniform (KS test p > 0.05)
✓ Mean t-statistic near 0
✓ Mean X-Y correlation near 0

### If Type I Error Rate is Inflated (> 0.06)

This would indicate that QAP regression produces **too many false positives** under these autocorrelation conditions. Possible causes:
- The permutation procedure doesn't adequately account for the autocorrelation structure
- The network properties create dependencies that violate QAP assumptions

### If Type I Error Rate is Deflated (< 0.04)

This would indicate QAP regression is **too conservative** under these conditions.

## Technical Details

### QAP Regression Implementation

The script uses `sna::netlm()` with:
- `mode = "digraph"`: Directed graphs
- `nullhyp = "qap"`: QAP permutation test
- **DSP method**: Double semi-partialing (Dekker-Krackhardt-Sanders)

### Two-Tailed P-values

`sna::netlm()` returns one-tailed p-values (`pleeq`). The script converts to two-tailed:

```r
p_two_tailed = 2 * min(p_one_tailed, 1 - p_one_tailed)
```

### Independence of X and Y

To ensure X and Y are truly independent:
- Different random seeds for each network
- Different node attribute assignments for X and Y
- No shared structure between generation processes

## Customization

To modify network parameters, edit these variables in the script:

```r
# Simulation size
n_simulations <- 2000      # Number of (X,Y) pairs
n_permutations <- 1000     # QAP permutations
network_sizes <- c(10, 30, 50)  # Node counts

# Network properties
target_density <- 0.15
target_reciprocity <- 0.4
target_transitivity <- 0.3
target_homophily <- 0.6
n_categories <- 3
```

## Citation

If using this simulation in research, please cite:

- **Dekker, D., Krackhardt, D., & Snijders, T. A. B.** (2007). Sensitivity of MRQAP tests to collinearity and autocorrelation conditions. *Psychometrika*, 72(4), 563-581.

- **Butts, C. T.** (2008). network: a Package for Managing Relational Data in R. *Journal of Statistical Software*, 24(2).

## Contact

For questions or issues with this simulation, please open an issue in the borgworld package repository.
