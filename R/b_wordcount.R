# bwordcount

#' Word frequency analysis with optional grouping
#'
#' Extracts words from a text column, counts their frequencies, and optionally
#' groups results by a categorical variable. Can return raw frequencies,
#' proportions, or percentages.
#'
#' @param df A data frame or tibble containing the text data.
#' @param text_col Unquoted column name containing text to analyze.
#' @param group_col Optional unquoted column name for grouping. When provided,
#'   word frequencies are calculated within each group and results are pivoted
#'   to wide format with a column per group level.
#' @param values Type of values to return: "freq" for raw counts, "prop" for
#'   proportions (0-1), or "pct" for percentages (0-100). Default is "pct".
#' @param remove_stopwords Logical; if TRUE (default), removes common English
#'   stopwords (e.g., "the", "and", "is") using the tidytext stop_words dataset.
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item \code{word}: The extracted word (lowercase).
#'     \item \code{total}: Overall frequency/proportion/percentage across all text.
#'     \item If \code{group_col} is provided: one column per group level containing
#'       the frequency/proportion/percentage within that group.
#'   }
#'   Results are sorted by \code{total} in descending order. Numeric values are
#'   formatted with appropriate decimal places (1 for percentages, 3 for
#'   proportions, 0 for frequencies).
#'
#' @details
#' Words are extracted using a regex pattern that matches sequences of lowercase
#' letters and apostrophes (e.g., "don't" is kept as one word). All text is
#' converted to lowercase before extraction.
#'
#' When \code{group_col} is provided:
#' \itemize{
#'   \item Proportions/percentages in group columns are calculated relative to
#'     the total words within each group.
#'   \item The \code{total} column shows the overall proportion/percentage
#'     relative to all words in the dataset.
#'   \item Groups with zero occurrences of a word show 0.
#' }
#'
#' @examples
#' \dontrun{
#' # Simple word frequency
#' df <- data.frame(
#'   text = c("The quick brown fox", "The lazy dog", "Quick brown dogs")
#' )
#' bwordcount(df, text)
#'
#' # With grouping
#' df <- data.frame(
#'   text = c("I love cats", "Dogs are great", "Cats are cute", "I love dogs"),
#'   person = c("Alice", "Bob", "Alice", "Bob")
#' )
#' bwordcount(df, text, person)
#'
#' # Raw frequencies without stopword removal
#' bwordcount(df, text, values = "freq", remove_stopwords = FALSE)
#' }
#'
#' @importFrom dplyr mutate count arrange desc group_by ungroup summarise
#'   across left_join anti_join
#' @importFrom tidyr unnest pivot_wider
#' @importFrom stringr str_extract_all
#' @importFrom rlang .data
#' @export
bwordcount <- function(df, text_col, group_col = NULL, values = "pct", remove_stopwords = TRUE) {

  result <- df |>
    dplyr::mutate(word = stringr::str_extract_all(tolower({{ text_col }}), "\\b[a-z']+\\b")) |>
    tidyr::unnest(.data$word)

  if (remove_stopwords) {
    result <- result |> dplyr::anti_join(tidytext::stop_words, by = "word")
  }

  total_words <- nrow(result)

  if (missing(group_col) || is.null(substitute(group_col))) {
    # No grouping - overall frequencies
    result <- result |>
      dplyr::count(.data$word, name = "total") |>
      dplyr::arrange(dplyr::desc(.data$total))

    if (values == "prop") {
      result <- result |> dplyr::mutate(total = .data$total / total_words)
    } else if (values == "pct") {
      result <- result |> dplyr::mutate(total = 100 * .data$total / total_words)
    }

  } else {
    # With grouping
    result <- result |>
      dplyr::count(.data$word, {{ group_col }})

    if (values == "prop") {
      result <- result |>
        dplyr::group_by({{ group_col }}) |>
        dplyr::mutate(n = .data$n / sum(.data$n)) |>
        dplyr::ungroup()
    } else if (values == "pct") {
      result <- result |>
        dplyr::group_by({{ group_col }}) |>
        dplyr::mutate(n = 100 * .data$n / sum(.data$n)) |>
        dplyr::ungroup()
    }

    # Get overall total (freq, prop, or pct based on total_words)
    overall <- result |>
      dplyr::group_by(.data$word) |>
      dplyr::summarise(total = sum(.data$n), .groups = "drop")

    if (values != "freq") {
      # Recalculate total as overall proportion/percentage
      word_totals <- df |>
        dplyr::mutate(word = stringr::str_extract_all(tolower({{ text_col }}), "\\b[a-z']+\\b")) |>
        tidyr::unnest(.data$word)

      if (remove_stopwords) {
        word_totals <- word_totals |> dplyr::anti_join(tidytext::stop_words, by = "word")
      }

      overall <- word_totals |>
        dplyr::count(.data$word, name = "total")

      if (values == "prop") {
        overall <- overall |> dplyr::mutate(total = .data$total / total_words)
      } else if (values == "pct") {
        overall <- overall |> dplyr::mutate(total = 100 * .data$total / total_words)
      }
    }

    result <- result |>
      tidyr::pivot_wider(names_from = {{ group_col }}, values_from = .data$n, values_fill = 0) |>
      dplyr::left_join(overall, by = "word") |>
      dplyr::arrange(dplyr::desc(.data$total))
  }

  # Round numeric columns
  decimals <- switch(values,
                     "pct" = 1,
                     "prop" = 3,
                     "freq" = 0)

  result <- result |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., decimals)))

  result <- as.data.frame(result)

 # Add custom class and store formatting info for print method
  class(result) <- c("bwordcount", class(result))
  attr(result, "decimals") <- decimals

  result
}

#' Print method for bwordcount objects
#'
#' Formats numeric columns with consistent decimal places while keeping
#' underlying data numeric.
#'
#' @param x A bwordcount object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns x.
#' @export
print.bwordcount <- function(x, ...) {
  decimals <- attr(x, "decimals") %||% 1

  # Create formatted copy for display
  display <- x
  for (col in names(display)) {
    if (is.numeric(display[[col]])) {
      display[[col]] <- format(display[[col]], nsmall = decimals, trim = TRUE)
    }
  }

  # Remove custom class for printing
  class(display) <- "data.frame"
  print(display, ...)

  invisible(x)
}
