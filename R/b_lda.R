#' Latent Dirichlet Allocation (LDA) Topic Modeling
#'
#' A user-friendly wrapper for LDA topic modeling that handles text preprocessing
#' and returns document-topic distributions.
#'
#' @param data A data frame containing the text to analyze
#' @param text_col The name or index of the column containing text documents
#' @param ntopics Number of topics to extract. Default is 10.
#' @param remove_stopwords Logical. Remove common English stopwords? Default is TRUE.
#' @param remove_numbers Logical. Remove numbers from text? Default is TRUE.
#' @param min_word_length Minimum word length to keep. Default is 3.
#' @param seed Random seed for reproducibility. Default is 12345.
#' @param verbose Logical. Print progress messages? Default is TRUE.
#'
#' @return A list of class "blda" containing:
#'   \itemize{
#'     \item \code{theta}: Document-topic probability matrix (documents x topics)
#'     \item \code{phi}: Topic-word probability matrix (topics x words)
#'     \item \code{dtm}: Document-term matrix
#'     \item \code{loglik}: Log-likelihood of the model
#'     \item \code{terms}: Top 10 terms for each topic
#'     \item \code{ntopics}: Number of topics
#'     \item \code{n_terms}: Vocabulary size
#'     \item \code{empty_docs}: Indices of documents with no terms after preprocessing
#'     \item \code{model}: The underlying LDA model object
#'   }
#'
#' @details
#' The function performs the following preprocessing steps:
#' \enumerate{
#'   \item Convert to lowercase
#'   \item Remove punctuation
#'   \item Remove numbers (optional)
#'   \item Remove English stopwords (optional)
#'   \item Remove short words
#'   \item Strip extra whitespace
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with 5 topics
#' result <- blda(mydata, "text_column", ntopics = 5)
#'
#' # Access document-topic matrix
#' result$theta
#'
#' # Access topic-word matrix
#' result$phi
#'
#' # Access document-term matrix
#' result$dtm
#'
#' # Get log-likelihood
#' result$loglik
#'
#' # View top terms
#' result$terms
#'
#' # Access the underlying model
#' result$model
#' }
#'
#' @export
blda <- function(data, text_col, ntopics = 10,
                 remove_stopwords = TRUE,
                 remove_numbers = TRUE,
                 min_word_length = 3,
                 seed = 12345,
                 verbose = TRUE) {

  # Check for required packages
  if (!requireNamespace("tm", quietly = TRUE)) {
    stop("Package 'tm' is required. Install with: install.packages('tm')")
  }
  if (!requireNamespace("topicmodels", quietly = TRUE)) {
    stop("Package 'topicmodels' is required. Install with: install.packages('topicmodels')")
  }

  # Handle text_col specification (name or index)
  if (is.numeric(text_col)) {
    text_col <- names(data)[text_col]
  }

  if (!text_col %in% names(data)) {
    stop(paste("Column", text_col, "not found in data frame"))
  }

  # Extract text
  texts <- as.character(data[[text_col]])
  n_docs <- length(texts)

  # Check for empty texts
  if (all(is.na(texts) | texts == "")) {
    stop("All documents are empty or NA")
  }

  if (verbose) message("Creating corpus from ", n_docs, " documents...")

  # Create corpus (use VCorpus to avoid transformation warnings)
  corpus <- tm::VCorpus(tm::VectorSource(texts))

  # Preprocess
  if (verbose) message("Preprocessing text...")
  corpus <- tm::tm_map(corpus, tm::content_transformer(tolower))
  corpus <- tm::tm_map(corpus, tm::removePunctuation)
  if (remove_numbers) {
    corpus <- tm::tm_map(corpus, tm::removeNumbers)
  }
  if (remove_stopwords) {
    corpus <- tm::tm_map(corpus, tm::removeWords, tm::stopwords("english"))
  }
  corpus <- tm::tm_map(corpus, tm::stripWhitespace)

  # Create document-term matrix
  dtm <- tm::DocumentTermMatrix(corpus)

  # Remove words shorter than min_word_length
  if (min_word_length > 1) {
    terms_to_keep <- nchar(tm::Terms(dtm)) >= min_word_length
    dtm <- dtm[, terms_to_keep]
  }

  # Identify empty documents
  row_totals <- as.numeric(slam::row_sums(dtm))
  empty_docs <- which(row_totals == 0)

  if (length(empty_docs) > 0) {
    warning(length(empty_docs), " document(s) had no terms after preprocessing and will have NA values")
  }

  # For LDA, we need non-empty documents
  dtm_for_lda <- dtm[row_totals > 0, ]

  if (nrow(dtm_for_lda) < 2) {
    stop("Fewer than 2 documents have terms after preprocessing. Cannot fit LDA.")
  }

  if (verbose) message("DTM has ", nrow(dtm_for_lda), " documents and ", ncol(dtm_for_lda), " terms")

  ntopics <- as.integer(ntopics)

  # Validate ntopics
  if (ntopics < 2) {
    stop("ntopics must be at least 2")
  }
  if (ntopics >= nrow(dtm_for_lda)) {
    stop("ntopics must be less than the number of non-empty documents (", nrow(dtm_for_lda), ")")
  }

  # Run LDA
  if (verbose) message("Fitting LDA model with ", ntopics, " topics...")
  set.seed(seed)
  lda_model <- topicmodels::LDA(dtm_for_lda, k = ntopics, method = "Gibbs",
                                control = list(seed = seed))

  # Extract document-topic distributions
  doc_topics <- topicmodels::posterior(lda_model)$topics

  # Extract topic-word distributions (topics as rows, words as columns)
  topic_words <- topicmodels::posterior(lda_model)$terms

  # Get log-likelihood (extract from model slot for Gibbs)
  log_lik <- lda_model@loglikelihood

  # Create document-topic data frame matching original document order
  theta <- as.data.frame(matrix(NA_real_, nrow = n_docs, ncol = ntopics))
  colnames(theta) <- paste0("Topic", 1:ntopics)
  rownames(theta) <- rownames(data)

  # Fill in values for non-empty documents
  non_empty_idx <- which(row_totals > 0)
  theta[non_empty_idx, ] <- doc_topics

  # Convert topic_words to data frame
  phi <- as.data.frame(topic_words)
  rownames(phi) <- paste0("Topic", 1:ntopics)

  # Get top terms
  top_terms <- topicmodels::terms(lda_model, 10)

  # Build result list
  result <- list(
    theta = theta,
    phi = phi,
    dtm = dtm_for_lda,
    loglik = log_lik,
    terms = top_terms,
    ntopics = ntopics,
    n_terms = ncol(dtm_for_lda),
    empty_docs = empty_docs,
    model = lda_model
  )

  class(result) <- c("blda", "list")

  # Print output
  print(result)

  # Return invisibly
  invisible(result)
}


#' Print method for blda objects
#'
#' @param x A blda object
#' @param n_terms Number of top terms to display per topic. Default is 10.
#' @param n_docs Number of document rows to preview. Default is 6.
#' @param ... Additional arguments (ignored)
#'
#' @export
print.blda <- function(x, n_terms = 10, n_docs = 6, ...) {

  n_total <- nrow(x$theta)

  cat("\nLDA Topic Model Results\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")

  cat("Documents:       ", n_total, "\n")
  cat("Topics:          ", x$ntopics, "\n")
  cat("Vocabulary:      ", x$n_terms, " terms\n")
  if (length(x$empty_docs) > 0) {
    cat("Empty docs:      ", length(x$empty_docs), "\n")
  }
  cat("Log-likelihood: ", format(x$loglik, big.mark = ",", scientific = FALSE), "\n")

  cat("\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Top terms per topic:\n")
  cat(paste(rep("-", 50), collapse = ""), "\n\n")

  if (!is.null(x$terms)) {
    n_show <- min(n_terms, nrow(x$terms))
    print(x$terms[1:n_show, , drop = FALSE], quote = FALSE)
  }

  cat("\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Document-topic distributions (first ", min(n_docs, n_total), " rows):\n")
  cat(paste(rep("-", 50), collapse = ""), "\n\n")

  preview <- utils::head(x$theta, n_docs)
  print(preview, digits = 3)

  invisible(x)
}


#' Extract top terms from blda object
#'
#' @param x A blda object
#' @param n Number of top terms per topic. Default is 10.
#'
#' @return A matrix with terms as rows and topics as columns
#'
#' @export
bldaterms <- function(x, n = 10) {
  if (!inherits(x, "blda")) {
    stop("x must be a blda object")
  }

  if (is.null(x$model)) {
    stop("LDA model not found in object")
  }

  topicmodels::terms(x$model, n)
}


#' Get topic assignments for new documents
#'
#' @param object A blda object (containing the fitted model)
#' @param newdata A data frame with new documents
#' @param text_col The name or index of the text column
#' @param ... Additional arguments (ignored)
#'
#' @return A data frame with topic distributions for new documents
#'
#' @export
predict.blda <- function(object, newdata, text_col, ...) {

  if (!requireNamespace("tm", quietly = TRUE)) {
    stop("Package 'tm' is required")
  }
  if (!requireNamespace("topicmodels", quietly = TRUE)) {
    stop("Package 'topicmodels' is required")
  }

  if (is.null(object$model)) {
    stop("LDA model not found in object")
  }

  # Handle text_col
  if (is.numeric(text_col)) {
    text_col <- names(newdata)[text_col]
  }

  texts <- as.character(newdata[[text_col]])

  # Create corpus with same preprocessing
  corpus <- tm::VCorpus(tm::VectorSource(texts))
  corpus <- tm::tm_map(corpus, tm::content_transformer(tolower))
  corpus <- tm::tm_map(corpus, tm::removePunctuation)
  corpus <- tm::tm_map(corpus, tm::removeNumbers)
  corpus <- tm::tm_map(corpus, tm::removeWords, tm::stopwords("english"))
  corpus <- tm::tm_map(corpus, tm::stripWhitespace)

  # Create DTM using same vocabulary as original model
  dtm <- tm::DocumentTermMatrix(corpus,
                                control = list(dictionary = topicmodels::Terms(object$model)))

  # Get posterior
  post <- topicmodels::posterior(object$model, dtm)

  result <- as.data.frame(post$topics)
  colnames(result) <- paste0("Topic", 1:ncol(result))

  return(result)
}
