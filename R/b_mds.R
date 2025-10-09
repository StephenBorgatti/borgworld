#' Unified Interface for MDS Methods
#'
#' @description
#' Wrapper function that provides a single interface for both classical and
#' non-metric MDS with flexible method name matching.
#'
#' @param x A square matrix of similarities or dissimilarities
#' @param type Character string specifying the input type. Must be one of
#'   "similarities", "dissimilarities", or any unambiguous abbreviation
#'   (e.g., "s", "sim", "d", "dis", "diss"). If missing, the user will be
#'   prompted to specify.
#' @param dim Integer specifying the number of dimensions. Default is 2.
#' @param method Character string specifying the MDS method. Accepts:
#'   \itemize{
#'     \item Classical MDS: "classical", "classic", "metric", "cmdscale", or
#'           any abbreviation starting with "c" or "m"
#'     \item Non-metric MDS: "nonmetric", "non-metric", "ordinal", "isoMDS",
#'           or any abbreviation starting with "n", "o", or "i"
#'   }
#'   Default is "classical".
#' @param plot Logical, whether to produce a plot (default = TRUE)
#' @param labels Optional character vector of labels for points
#' @param verbose Logical, whether to print method selection info (default = TRUE)
#' @param ... Additional arguments passed to the underlying MDS function
#'
#' @return An MDS object from the selected method with additional attributes:
#'   \item{mds.method}{The MDS method used ("classical" or "nonmetric")}
#'
#' @examples
#' # Basic usage
#' dist_mat <- as.matrix(dist(USArrests))
#' bmds(dist_mat, "d")  # Will use classical MDS by default
#' bmds(dist_mat, "dissimilarities", method = "n")  # Non-metric MDS
#'
#' # With correlation matrix
#' cor_mat <- cor(mtcars)
#' bmds(cor_mat, "s", dim = 3)
#' bmds(cor_mat, "sim", method = "nonmetric", dim = 2)
#'
#' @export
bmds <- function(x, type, dim = 2, method = "classical",
                 plot = TRUE, labels = NULL, verbose = TRUE, ...) {

  # Store original call
  orig_call <- match.call()

  # Check if type is provided, if not, prompt user
  if (missing(type)) {
    if (interactive()) {
      cat("Please specify the type of your input matrix:\n")
      cat("  1. similarities\n")
      cat("  2. dissimilarities\n")
      response <- readline("Enter your choice (1 or 2, or 's'/'d'): ")

      response <- tolower(trimws(response))
      if (response %in% c("1", "s", "sim", "similarities")) {
        type <- "similarities"
      } else if (response %in% c("2", "d", "dis", "diss", "dissimilarities")) {
        type <- "dissimilarities"
      } else {
        stop("Invalid choice. Please specify 'similarities' or 'dissimilarities'")
      }
      cat("Using type:", type, "\n")
    } else {
      stop("Argument 'type' is missing with no default. Please specify 'similarities' or 'dissimilarities'")
    }
  }

  # Match type argument with abbreviations
  type <- match.arg(type, c("similarities", "dissimilarities"))

  # Handle dist objects
  if (inherits(x, "dist")) {
    x <- as.matrix(x)
    if (verbose) message("Input converted from dist object to matrix")
  }

  # Ensure matrix format
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  # Check for square matrix
  if (nrow(x) != ncol(x)) {
    stop("Input must be a square matrix or dist object")
  }

  # Convert to lowercase for method matching
  method <- tolower(method)

  # Define method patterns (expanded for better matching)
  classical_patterns <- c("^c", "^m", "^p",  # Single letter abbreviations
                          "^class", "^metric", "^cmdscale", "^pcoa", "^pca",
                          "classical", "classic")

  nonmetric_patterns <- c("^n", "^o", "^i", "^k",  # Single letter abbreviations
                          "^non", "^ord", "^iso", "^krus",
                          "nonmetric", "non-metric", "ordinal", "kruskal")

  # Determine which method to use
  if (any(sapply(classical_patterns, function(p) grepl(p, method)))) {
    used_method <- "classical"
    if (verbose) message("Using classical (metric) MDS")

    # Call classical MDS
    result <- bclassicalmds(x = x,
                            type = type,
                            dim = dim,
                            plot = plot,
                            labels = labels,
                            ...)

  } else if (any(sapply(nonmetric_patterns, function(p) grepl(p, method)))) {
    used_method <- "nonmetric"
    if (verbose) message("Using non-metric MDS (Kruskal)")

    # Call non-metric MDS
    result <- bnonmetricmds(x = x,
                            type = type,
                            dim = dim,
                            plot = plot,
                            labels = labels,
                            ...)

  } else {
    # Provide helpful error message with suggestions
    stop("Method '", method, "' not recognized.\n",
         "  Classical MDS: use 'classical', 'metric', 'c', 'm', etc.\n",
         "  Non-metric MDS: use 'nonmetric', 'ordinal', 'n', 'o', etc.")
  }

  # Add attributes to result
  attr(result, "mds.method") <- used_method
  attr(result, "original.call") <- orig_call

  # Add unified class
  class(result) <- c("bmds.unified", class(result))

  return(result)
}


#' Summary Method for Unified bmds Objects
#'
#' @param object A bmds.unified object
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns the object
#'
#' @export
#' @method summary bmds.unified
summary.bmds.unified <- function(object, ...) {
  cat("\nUnified MDS Analysis Summary\n")
  cat("============================\n")

  # Get attributes
  mds_method <- attr(object, "mds.method")
  orig_call <- attr(object, "original.call")

  if (!is.null(orig_call)) {
    cat("Original call:\n")
    print(orig_call)
    cat("\n")
  }

  if (!is.null(mds_method)) {
    cat("MDS Method Used:", toupper(mds_method), "\n\n")
  }

  # Call the appropriate print method based on the underlying class
  if (inherits(object, "bnmds")) {
    print.bnmds(object)
  } else if (inherits(object, "bmds")) {
    print.bmds(object)
  } else {
    print(object)
  }

  invisible(object)
}
