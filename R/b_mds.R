#' Unified interface for MDS methods (Alternative with flexible matching)
#'
#' @description
#' Wrapper function that provides a single interface for both classical and
#' non-metric MDS with flexible method name matching.
#'
#' @param d A distance matrix or dist object
#' @param k Integer specifying the number of dimensions. Default is 2.
#' @param method Character string specifying the MDS method. Accepts:
#'   \itemize{
#'     \item Classical MDS: "classical", "classic", "metric", "cmdscale", or
#'           any abbreviation starting with "c" or "m"
#'     \item Non-metric MDS: "nonmetric", "non-metric", "ordinal", "isoMDS",
#'           or any abbreviation starting with "n", "o", or "i"
#'   }
#'   Default is "classical".
#' @param ... Additional arguments passed to the underlying MDS function.
#'
#' @return An MDS object from the selected method.
#'
#' @export
bmds <- function(d, k = 2, method = "classical", ...) {

  # Convert to lowercase for matching
  method <- tolower(method)

  # Define method patterns
  classical_patterns <- c("^c", "^m", "^classical", "^classic", "^metric", "^cmdscale")
  nonmetric_patterns <- c("^n", "^o", "^i", "^non", "^ordinal", "^iso")

  # Check which pattern matches
  if(any(sapply(classical_patterns, function(p) grepl(p, method)))) {
    result <- bclassicalmds(d, k = k, ...)
    used_method <- "classical"
  } else if(any(sapply(nonmetric_patterns, function(p) grepl(p, method)))) {
    result <- bnonmetricmds(d, k = k, ...)
    used_method <- "nonmetric"
  } else {
    stop("Method '", method, "' not recognized. Use 'classical' or 'nonmetric' (or abbreviations).")
  }

  # Add method attribute
  attr(result, "mds.method") <- used_method

  # Optional: Print method used
  message("Using ", used_method, " MDS")

  return(result)
}
