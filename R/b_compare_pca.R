#' Compare borgworld PCA with R's Native Functions
#'
#' @description
#' Compares the results from \code{bpca} with R's built-in
#' PCA functions to demonstrate the differences in terminology and scaling
#' conventions between statistical computing and social science software.
#'
#' @param result An object of class \code{bpca} from \code{\link{bpca}}
#' @param data The original data matrix or data frame used in the PCA
#'
#' @return Invisibly returns a list containing the comparison results.
#'   The function is primarily called for its side effect of printing
#'   a detailed comparison to the console.
#'
#' @details
#' This function demonstrates how to convert between R's PCA conventions
#' and social science conventions:
#'
#' \itemize{
#'   \item R's \code{prcomp} returns eigenvectors as "rotation"
#'   \item Social science "loadings" = R's rotation * sqrt(eigenvalues)
#'   \item R's scores have variance = eigenvalue
#'   \item Standardized scores = R's scores / sqrt(eigenvalue)
#' }
#'
#' @examples
#' # Create example data
#' set.seed(123)
#' data <- as.data.frame(matrix(rnorm(100 * 5), nrow = 100, ncol = 5))
#' colnames(data) <- paste0("Var", 1:5)
#'
#' # Run borgworld PCA
#' result <- bpca(data)  # No n_components parameter
#'
#' # Or if you want to specify which components to plot:
#' result <- bpca(data, choices = c(1, 3))
#'
#' @seealso \code{\link{bpca}}
#' @importFrom stats prcomp
#'
#' @export
bcompare_pca <- function(result, data) {
  # Check input class
  if (!inherits(result, "bpca")) {
    stop("result must be of class 'bpca' from bpca()")
  }

  cat("\n===============================================\n")
  cat("COMPARISON WITH R FUNCTIONS\n")
  cat("===============================================\n\n")

  # Run prcomp
  pca_r <- stats::prcomp(data, scale. = TRUE, center = TRUE)

  cat("What R's prcomp() calls things:\n")
  cat("--------------------------------\n")
  cat("- 'rotation' = eigenvectors (our raw coefficients)\n")
  cat("- 'sdev' = sqrt(eigenvalues)\n")
  cat("- 'x' = raw scores (variance = eigenvalue)\n\n")

  cat("To get social science style output from prcomp:\n")
  cat("-----------------------------------------------\n")
  cat("loadings = pca$rotation %*% diag(pca$sdev)\n")
  cat("standardized_scores = pca$x %*% diag(1/pca$sdev)\n\n")

  # Show the conversion
  n_show <- min(3, ncol(result$loadings))
  r_loadings <- pca_r$rotation[, 1:n_show] %*% diag(pca_r$sdev[1:n_show])

  cat("First 3 variables, PC1 loadings:\n")
  comparison <- data.frame(
    Our_Loadings = round(result$loadings[1:3, 1], 4),
    R_Converted = round(r_loadings[1:3, 1], 4),
    R_Raw = round(pca_r$rotation[1:3, 1], 4)
  )
  rownames(comparison) <- rownames(result$loadings)[1:3]
  print(comparison)

  invisible(list(
    comparison_table = comparison,
    r_pca = pca_r,
    r_loadings_converted = r_loadings
  ))
}
