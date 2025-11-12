# QAP Correlation Type I Error Test

## Overview

This simulation tests whether the **Quadratic Assignment Procedure (QAP)** correlation test maintains proper Type I error control when testing the correlation between two independent networks that exhibit network autocorrelation (reciprocity and transitivity).

## Research Question

**Does QAP correlation control Type I error at the nominal α = 0.05 level when networks have internal structure?**

Under the null hypothesis (no true correlation between X and Y), a valid statistical test should reject at most 5% of the time when α = 0.05. This simulation verifies that QAP maintains this property even when individual networks exhibit:

- **Reciprocity** (~0.4): Tendency for edges to be bidirectional (if i→j exists, j→i is likely)
- **Transitivity** (~0.3): Clustering/closure effects (if i→j and j→k exist, i→k is likely)
- **Low density** (~0.1): Sparse networks with ~10% of possible edges

## Files

- **`qap_type1_error_test.R`**: Main simulation script (full analysis)
- **`qap_type1_error_test_quick.R`**: Quick test version with reduced parameters
- **`QAP_TYPE1_ERROR_README.md`**: This documentation file

## Requirements

### R Packages
```r
install.packages("sna")
```

The `sna` package is required for calculating network transitivity using the `gtrans()` function.

### System Requirements
- R version 4.1.0 or higher
- At least 2GB RAM recommended
- Runtime: ~30-60 minutes for full simulation (2000 simulations × 2 network sizes)

## Usage

### Quick Test (Recommended First)

Run the quick test version to verify everything works:

```bash
Rscript qap_type1_error_test_quick.R
```

This runs 50 simulations with 100 permutations each on 20-node networks (~1-2 minutes).

### Full Simulation

Run the complete analysis:

```bash
Rscript qap_type1_error_test.R
```

This runs:
- **2000 simulations** per network size
- **1000 QAP permutations** per test
- Two network sizes: **n = 30** and **n = 50** nodes
- Expected runtime: **30-60 minutes**

## What the Script Does

### 1. Network Generation

For each simulation, generates two completely **independent** directed networks (X and Y) with:
- Target edge density: 0.1
- Target reciprocity: 0.4
- Target transitivity: 0.3

The generation uses an iterative rewiring algorithm that:
1. Starts with a random directed graph
2. Iteratively adds/removes reciprocal pairs
3. Adds/removes transitive edges (closes/opens triads)
4. Maintains approximate target density

### 2. QAP Correlation Test

For each (X, Y) pair:
1. Computes observed correlation between off-diagonal elements
2. Permutes rows and columns of X (1000 times)
3. Computes correlation for each permutation
4. Calculates two-tailed p-value

### 3. Analysis

Records:
- **Type I error rate**: Proportion of p-values < 0.05
- **P-value distribution**: Should be approximately uniform [0,1]
- **Actual correlations**: Should be near 0 (networks are independent)
- **Network properties**: Verify achieved reciprocity, transitivity, density

## Output

The script produces:

### Console Output

1. **Type I Error Rate**
   - Observed rate (should be ~0.05)
   - Expected 95% CI
   - Pass/fail assessment

2. **P-value Distribution**
   - Mean, median (should be ~0.5)
   - Kolmogorov-Smirnov test for uniformity

3. **Network Properties**
   - Actual density, reciprocity, transitivity achieved
   - Mean ± SD across all simulations

4. **Runtime Statistics**

### Diagnostic Plots

Four plots per network size:

1. **Histogram of p-values** (should be approximately flat/uniform)
2. **Q-Q plot** vs uniform distribution (should follow diagonal)
3. **Distribution of actual correlations** (should be centered near 0)
4. **Cumulative distribution** of p-values (should match diagonal)

### Summary Table

```
   n Type_I_Error Mean_Cor Runtime_sec
1 30       0.0520   0.0012      1234.5
2 50       0.0485  -0.0023      2345.6
```

## Expected Results

If QAP correctly controls Type I error:

- ✓ Type I error rate ≈ 0.05 (within 95% CI: [0.040, 0.060])
- ✓ P-value distribution is approximately uniform
- ✓ Mean actual correlation ≈ 0
- ✓ K-S test p-value > 0.05 (fails to reject uniformity)

## Interpretation

### Success Criteria

The test **succeeds** if:
1. Type I error rate is within expected 95% CI [0.040, 0.060]
2. P-values are approximately uniformly distributed
3. Mean correlation is close to 0 (true value)

This demonstrates that QAP correlation is a **valid** test that properly controls false positives even when networks have autocorrelation structure.

### Potential Issues

If Type I error rate is significantly different from 0.05:

- **Inflated (> 0.06)**: QAP may be too liberal (too many false positives)
  - Could indicate the permutation test doesn't properly account for network structure
  - May need double-semi-partialling or other adjustments

- **Deflated (< 0.04)**: QAP may be too conservative (too few true positives)
  - Less concerning but reduces statistical power

## Technical Details

### QAP Implementation

The `qap_cor()` function:
- Uses off-diagonal elements only (excludes self-loops)
- Performs **double permutation** (rows and columns simultaneously)
- Computes **two-tailed** p-value: `P(|r_perm| >= |r_obs|)`
- Default: 1000 permutations

### Network Generation Algorithm

The `generate_network()` function uses an iterative approach rather than ERGM because:
- Much faster for moderate networks (n=30-50)
- Simpler implementation with no convergence issues
- Sufficient for this validation study

The algorithm:
1. Initializes with random directed graph at target density
2. Iteratively adjusts edges to meet targets
3. Stops when within tolerance (default: 0.05) or max iterations (1000)

### Statistical Theory

Under H₀ (no correlation), the QAP test should produce:
- **Uniform p-values**: P ~ Uniform(0,1)
- **Type I error**: P(p < α) = α for any α

The simulation tests this empirically by:
1. Generating pairs of truly independent networks (H₀ is true)
2. Running QAP tests
3. Checking if rejection rate matches nominal α

## Reproducibility

- Set seed: `set.seed(42)` for reproducible results
- Save results: Modify script to save `results_summary` object

```r
# Add at end of script:
saveRDS(results_summary, "qap_type1_results.rds")
```

## References

### QAP Method
- Hubert, L., & Schultz, J. (1976). Quadratic assignment as a general data analysis strategy. *British Journal of Mathematical and Statistical Psychology*.
- Krackhardt, D. (1988). Predicting with networks: Nonparametric multiple regression analysis of dyadic data. *Social Networks*.

### Network Autocorrelation
- Dekker, D., Krackhardt, D., & Snijders, T. A. (2007). Sensitivity of MRQAP tests to collinearity and autocorrelation conditions. *Psychometrika*.

## Contact

For questions or issues, please open an issue at:
https://github.com/StephenBorgatti/borgworld/issues

---

**Note**: This simulation is designed for methodological validation. For production QAP analysis, consider using established packages like `sna::qaptest()` or `asnipe::network_permutation()`.
