# b_data.R - Data manipulation utilities

#' Cut a numeric vector into k ordered groups via k-means
#'
#' @param x Numeric vector (coercible).
#' @param k Integer number of groups (default 3).
#' @return Integer vector (1..k) aligned so 1=lowest mean group, k=highest.
#'         Non-finite inputs return NA.
#' @export
#' @examples
#'   # Create example data
#'   x <- rnorm(100)
#'   # Cut into 3 groups
#'   bautocut(x, 3)
bautocut <- function(x, k = 3) {
  x <- as.numeric(x)
  idx <- is.finite(x)
  if (sum(idx) < k || length(unique(x[idx])) < k) {
    stop("Need at least ", k, " distinct finite values.")
  }

  set.seed(1)
  km <- stats::kmeans(x[idx], centers = k)

  centers <- tapply(x[idx], km$cluster, mean)
  ranks <- rank(centers, ties.method = "first")  # 1=lowest center ... k=highest

  out <- rep(NA_integer_, length(x))
  out[idx] <- ranks[km$cluster]
  out
}


#' Clean data for borgworld analysis functions
#'
#' Removes rows with missing values and selects only numeric columns.
#' This is a preprocessing utility used by various borgworld functions.
#' This is a convenience wrapper around bprepare_data().
#'
#' @param data A data frame or matrix
#' @param na_action How to handle missing values: "complete" (default) for listwise
#'   deletion, "pairwise" for keeping all rows, "none" to keep all data
#' @param verbose Print messages about cleaning (default FALSE)
#'
#' @return A data frame with only complete cases and numeric columns
#'
#' @examples
#' # Clean data, removing non-numeric columns and rows with NAs
#' bclean(iris)
#'
#' @export
bclean <- function(data, na_action = "complete", verbose = FALSE) {
  bprepare_data(data,
                na_action = na_action,
                numeric_only = TRUE,
                output_format = "data.frame",
                verbose = verbose)
}


#' Standardize columns of a data frame
#'
#' Standardizes (z-scores) numeric columns after removing NAs and non-numeric columns.
#' By default uses population standard deviation (n), but can optionally use
#' sample standard deviation (n-1).
#'
#' @param data A data frame or matrix
#' @param sample Logical; if FALSE (default), uses population SD (divide by n).
#'   If TRUE, uses sample SD (divide by n-1, R's default).
#' @param center Logical; if TRUE (default), centers the data by subtracting the mean.
#'   If FALSE, only scales by the standard deviation.
#' @param na_action How to handle missing values: "complete" (default) for listwise
#'   deletion, "pairwise" for keeping all rows, "none" to keep all data
#' @param verbose Print messages about cleaning (default FALSE)
#'
#' @return A data frame with standardized numeric columns
#'
#' @details
#' This function first cleans the data using bprepare_data() to remove rows with
#' NAs and select only numeric columns, then standardizes each column to have
#' mean 0 and standard deviation 1.
#'
#' The population standard deviation formula (default) divides by n, while
#' the sample standard deviation divides by n-1. Most statistical software
#' uses the sample SD by default, but population SD is sometimes preferred
#' for descriptive statistics.
#'
#' @examples
#' # Create example data
#' set.seed(123)
#' data <- data.frame(
#'   x = rnorm(100, mean = 10, sd = 2),
#'   y = rnorm(100, mean = 50, sd = 10),
#'   z = letters[1:100]  # Will be removed automatically
#' )
#'
#' # Standardize using population SD (default)
#' data_std <- bstandardize(data)
#'
#' # Verify: should have mean ~0 and SD ~1
#' colMeans(data_std)
#' apply(data_std, 2, sd) * sqrt(99/100)  # Convert back to pop SD
#'
#' # Standardize using sample SD
#' data_std_sample <- bstandardize(data, sample = TRUE)
#'
#' @importFrom stats sd
#' @export
bstandardize <- function(data, sample = FALSE, center = TRUE,
                         na_action = "complete", verbose = FALSE) {
  # Clean the data first using standardized input handling
  data <- bprepare_data(data,
                        na_action = na_action,
                        numeric_only = TRUE,
                        output_format = "data.frame",
                        verbose = verbose)

  # Check if any data remains
  if (nrow(data) == 0 || ncol(data) == 0) {
    stop("No data remaining after cleaning")
  }

  n <- nrow(data)

  # Calculate means for centering
  if (center) {
    means <- colMeans(data)
  } else {
    means <- rep(0, ncol(data))
  }

  # Calculate standard deviations
  if (sample) {
    # Sample SD (n-1), R's default
    sds <- apply(data, 2, sd)
  } else {
    # Population SD (n)
    sds <- apply(data, 2, function(x) {
      sqrt(sum((x - mean(x))^2) / length(x))
    })
  }

  # Check for zero variance columns
  zero_var <- sds == 0
  if (any(zero_var)) {
    warning(paste("Column(s)", paste(names(data)[zero_var], collapse = ", "),
                  "have zero variance and cannot be standardized"))
    sds[zero_var] <- 1  # Don't divide by zero
  }

  # Standardize
  data_std <- as.data.frame(scale(data, center = means, scale = sds))

  return(data_std)
}
