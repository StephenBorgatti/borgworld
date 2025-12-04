#' Reshape Coder Matrices to Code Matrices
#'
#' Takes multiple coder matrices (documents × codes) and reshapes them into
#' a list of code matrices (documents × coders), suitable for inter-rater
#' reliability analysis with \code{\link{bkappa}}.
#'
#' @param ... Matrices where rows are documents and columns are codes.
#'   Each matrix represents one coder's ratings. All matrices must have
#'   the same dimensions.
#' @param coder_names Optional character vector of coder names. If NULL,
#'   uses the names of the arguments or generates names like "Coder1", "Coder2", etc.
#'
#' @return A named list of matrices, one per code. Each matrix has documents
#'   as rows and coders as columns. The list is given class "breshape_coders"
#'   for pretty printing.
#'
#' @examples
#' # Three coders rating 5 documents on 3 codes
#' alice <- matrix(c(1,1,0,1,0, 0,1,1,0,1, 1,0,0,1,1), nrow=5, ncol=3,
#'                 dimnames=list(paste0("Doc",1:5), c("Theme1","Theme2","Theme3")))
#' bob   <- matrix(c(1,1,0,0,0, 0,1,1,0,1, 1,1,0,1,1), nrow=5, ncol=3,
#'                 dimnames=list(paste0("Doc",1:5), c("Theme1","Theme2","Theme3")))
#' carol <- matrix(c(1,0,0,1,0, 0,1,1,1,1, 1,0,0,1,0), nrow=5, ncol=3,
#'                 dimnames=list(paste0("Doc",1:5), c("Theme1","Theme2","Theme3")))
#'
#' code_matrices <- breshape_coders(alice, bob, carol)
#' bkappa(code_matrices$Theme1)
#'
#' @export
breshape_coders <- function(..., coder_names = NULL) {
  coder_list <- list(...)
  n_coders <- length(coder_list)

  if (n_coders < 2) {
    stop("At least 2 coder matrices are required")
  }

  # Get dimensions from first matrix
  n_docs <- nrow(coder_list[[1]])
  n_codes <- ncol(coder_list[[1]])


  # Validate all matrices have same dimensions

  for (i in seq_along(coder_list)) {
    if (!is.matrix(coder_list[[i]])) {
      coder_list[[i]] <- as.matrix(coder_list[[i]])
    }
    if (nrow(coder_list[[i]]) != n_docs || ncol(coder_list[[i]]) != n_codes) {
      stop(sprintf("Matrix %d has dimensions %d x %d, expected %d x %d",
                   i, nrow(coder_list[[i]]), ncol(coder_list[[i]]), n_docs, n_codes))
    }
  }

  # Get code names from column names of first matrix
  code_names <- colnames(coder_list[[1]])
  if (is.null(code_names)) {
    code_names <- paste0("Code", seq_len(n_codes))
  }

  # Get document names from row names of first matrix
  doc_names <- rownames(coder_list[[1]])
  if (is.null(doc_names)) {
    doc_names <- paste0("Doc", seq_len(n_docs))
  }

  # Get coder names
  if (is.null(coder_names)) {
    coder_names <- as.character(substitute(list(...)))[-1]
    if (length(coder_names) != n_coders || any(coder_names == "")) {
      coder_names <- paste0("Coder", seq_len(n_coders))
    }
  }

  # Reshape: for each code, create a documents × coders matrix
  result <- vector("list", n_codes)
  names(result) <- code_names

  for (j in seq_len(n_codes)) {
    code_matrix <- matrix(NA_real_, nrow = n_docs, ncol = n_coders,
                          dimnames = list(doc_names, coder_names))
    for (i in seq_len(n_coders)) {
      code_matrix[, i] <- coder_list[[i]][, j]
    }
    result[[j]] <- code_matrix
  }

  class(result) <- c("breshape_coders", "list")
  attr(result, "n_docs") <- n_docs
  attr(result, "n_coders") <- n_coders
  attr(result, "n_codes") <- n_codes
  result
}

#' @export
print.breshape_coders <- function(x, ...) {
  cat(sprintf("Reshaped coder data: %d documents, %d coders, %d codes\n\n",
              attr(x, "n_docs"), attr(x, "n_coders"), attr(x, "n_codes")))
  cat("Codes:", paste(names(x), collapse = ", "), "\n")
  cat("Coders:", paste(colnames(x[[1]]), collapse = ", "), "\n\n")
  cat("Access individual code matrices with $CodeName or [[index]]\n")
  invisible(x)
}
