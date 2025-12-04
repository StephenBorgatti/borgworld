#' Compute Inter-Rater Reliability (Kappa Statistics)
#'
#' Computes Fleiss's kappa for multiple raters and optionally Cohen's kappa
#' for all pairwise combinations. Also returns percent agreement and majority
#' codes per document.
#'
#' @param x A matrix where rows are documents (or units) and columns are coders.
#'   Values should be 0/1 for binary coding, or category labels for nominal coding.
#' @param pairwise Logical. If TRUE, compute Cohen's kappa for all pairs of coders.
#'   Default is TRUE.
#'
#' @return A list of class "bkappa" containing:
#'   \item{fleiss}{Fleiss's kappa (overall agreement among all raters)}
#'   \item{pct_agree}{Overall percent agreement (proportion of documents where all coders agree)}
#'   \item{majority}{Vector of majority codes per document (NA if tied)}
#'   \item{cohen}{Matrix of pairwise Cohen's kappa values (if pairwise=TRUE)}
#'   \item{n_docs}{Number of documents}
#'   \item{n_coders}{Number of coders}
#'   \item{categories}{Unique categories found in the data}
#'
#' @details
#' Fleiss's kappa extends Cohen's kappa to multiple raters. Values range from
#' -1 to 1, where 1 indicates perfect agreement, 0 indicates agreement at
#' chance level, and negative values indicate agreement worse than chance.
#'
#' Common interpretation guidelines (Landis & Koch, 1977):
#' \itemize{
#'   \item < 0.00: Poor
#'   \item 0.00-0.20: Slight
#'   \item 0.21-0.40: Fair
#'   \item 0.41-0.60: Moderate
#'   \item 0.61-0.80: Substantial
#'   \item 0.81-1.00: Almost perfect
#' }
#'
#' @examples
#' # 5 documents rated by 3 coders (binary: 0/1)
#' ratings <- matrix(c(
#'   1, 1, 1,   # Doc1: all agree (1)
#'   0, 0, 1,   # Doc2: 2 say 0, 1 says 1
#'   1, 1, 0,   # Doc3: 2 say 1, 1 says 0
#'   0, 0, 0,   # Doc4: all agree (0)
#'   1, 0, 1    # Doc5: 2 say 1, 1 says 0
#' ), nrow = 5, byrow = TRUE,
#'    dimnames = list(paste0("Doc", 1:5), c("Alice", "Bob", "Carol")))
#'
#' result <- bkappa(ratings)
#' result$fleiss
#' result$majority
#'
#' @references
#' Fleiss, J.L. (1971). Measuring nominal scale agreement among many raters.
#' Psychological Bulletin, 76(5), 378-382.
#'
#' Cohen, J. (1960). A coefficient of agreement for nominal scales.
#' Educational and Psychological Measurement, 20(1), 37-46.
#'
#' @export
bkappa <- function(x, pairwise = TRUE) {

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  n_docs <- nrow(x)
  n_coders <- ncol(x)

  if (n_coders < 2) {
    stop("At least 2 coders are required")
  }

  # Get unique categories (excluding NA)
  categories <- sort(unique(as.vector(x[!is.na(x)])))
  n_cats <- length(categories)

  if (n_cats < 2) {
    warning("Only one category found in data; kappa is undefined")
    fleiss_kappa <- NA_real_
  } else {
    # Compute Fleiss's kappa
    fleiss_kappa <- compute_fleiss(x, categories)
  }

  # Percent agreement: proportion of documents where all coders agree
  pct_agree <- mean(apply(x, 1, function(row) {
    row <- row[!is.na(row)]
    if (length(row) < 2) return(NA)
    length(unique(row)) == 1
  }), na.rm = TRUE)

  # Majority code per document
  majority <- apply(x, 1, function(row) {
    row <- row[!is.na(row)]
    if (length(row) == 0) return(NA)
    tab <- table(row)
    max_count <- max(tab)
    winners <- names(tab)[tab == max_count]
    if (length(winners) > 1) {
      return(NA)  # Tie
    }
    # Return as same type as input
    if (is.numeric(x)) {
      return(as.numeric(winners))
    }
    return(winners)
  })

  # Make majority a proper vector
  if (is.list(majority)) {
    majority <- unlist(majority)
  }
  names(majority) <- rownames(x)

  # Pairwise Cohen's kappa
  cohen_matrix <- NULL
  if (pairwise && n_coders >= 2) {
    cohen_matrix <- matrix(NA_real_, nrow = n_coders, ncol = n_coders,
                           dimnames = list(colnames(x), colnames(x)))
    diag(cohen_matrix) <- 1
    for (i in 1:(n_coders - 1)) {
      for (j in (i + 1):n_coders) {
        k <- compute_cohen(x[, i], x[, j])
        cohen_matrix[i, j] <- k
        cohen_matrix[j, i] <- k
      }
    }
  }

  result <- list(
    fleiss = fleiss_kappa,
    pct_agree = pct_agree,
    majority = majority,
    cohen = cohen_matrix,
    n_docs = n_docs,
    n_coders = n_coders,
    categories = categories
  )

  class(result) <- "bkappa"
  result
}

#' Compute Fleiss's Kappa
#' @keywords internal
compute_fleiss <- function(x, categories) {
  n <- nrow(x)
  k <- ncol(x)  # number of raters
  q <- length(categories)

  # For each document, count how many raters assigned each category
  # n_ij = number of raters who assigned category j to document i
  counts <- matrix(0, nrow = n, ncol = q)
  colnames(counts) <- as.character(categories)

  for (i in seq_len(n)) {
    row <- x[i, ]
    row <- row[!is.na(row)]
    n_i <- length(row)  # number of raters for this document
    if (n_i > 0) {
      for (j in seq_along(categories)) {
        counts[i, j] <- sum(row == categories[j])
      }
    }
  }

  # P_i = agreement for document i
  # P_i = (1 / (n_i * (n_i - 1))) * sum_j(n_ij * (n_ij - 1))
  n_per_doc <- rowSums(!is.na(x))
  P_i <- numeric(n)

  for (i in seq_len(n)) {
    n_i <- n_per_doc[i]
    if (n_i < 2) {
      P_i[i] <- NA
    } else {
      P_i[i] <- sum(counts[i, ] * (counts[i, ] - 1)) / (n_i * (n_i - 1))
    }
  }

  P_bar <- mean(P_i, na.rm = TRUE)

  # p_j = proportion of all ratings in category j
  total_ratings <- sum(n_per_doc)
  p_j <- colSums(counts) / total_ratings

  # P_e = expected agreement by chance
  P_e <- sum(p_j^2)

  # Fleiss's kappa
  if (P_e == 1) {
    return(NA_real_)
  }

  kappa <- (P_bar - P_e) / (1 - P_e)
  kappa
}

#' Compute Cohen's Kappa for Two Raters
#' @keywords internal
compute_cohen <- function(r1, r2) {
  # Remove pairs where either is NA
  valid <- !is.na(r1) & !is.na(r2)
  r1 <- r1[valid]
  r2 <- r2[valid]

  if (length(r1) < 1) return(NA_real_)

  # Observed agreement
  p_o <- mean(r1 == r2)

  # Expected agreement
  categories <- unique(c(r1, r2))
  p_e <- 0
  for (cat in categories) {
    p1 <- mean(r1 == cat)
    p2 <- mean(r2 == cat)
    p_e <- p_e + p1 * p2
  }

  if (p_e == 1) return(NA_real_)

  kappa <- (p_o - p_e) / (1 - p_e)
  kappa
}

#' @export
print.bkappa <- function(x, digits = 3, ...) {
  cat("Inter-Rater Reliability Analysis\n")
  cat(strrep("-", 40), "\n")
  cat(sprintf("Documents: %d    Coders: %d    Categories: %s\n\n",
              x$n_docs, x$n_coders, paste(x$categories, collapse = ", ")))

  cat(sprintf("Fleiss's kappa:    %s\n", format(round(x$fleiss, digits), nsmall = digits)))
  cat(sprintf("Percent agreement: %s\n\n", format(round(x$pct_agree * 100, 1), nsmall = 1)))

  if (!is.null(x$cohen)) {
    cat("Pairwise Cohen's Kappa:\n")
    print(round(x$cohen, digits))
    cat("\n")
  }

  cat("Majority codes:\n")
  print(x$majority)

  invisible(x)
}
