# bstats

#' Column SDs for numeric columns, long format
#' These are "population sds" using the n divisor
#' @param df A data frame or tibble.
#' @return Base data.frame with columns: variable, sd.
#' @export
bsd <- function(df) {
  pop_sd <- function(x) {
    x <- x[!is.na(x)]
    mu <- mean(x)
    sqrt(mean((x - mu)^2))
  }

  df |>
    dplyr::summarise(dplyr::across(where(is.numeric), \(x) pop_sd(x))) |>
    tidyr::pivot_longer(cols = dplyr::everything(),
                        names_to = "variable",
                        values_to = "sd") |>
    as.data.frame()
}

#' Summarize numeric columns
#' @param df data.frame or tibble
#' @return data.frame with columns: variable, min, max, median, mean, sd, n_missing
#' @export
bsummarize <- function(df) {
  df |>
    dplyr::summarise(
      dplyr::across(
        where(is.numeric),
        list(
          min       = \(x) min(x, na.rm = TRUE),
          max       = \(x) max(x, na.rm = TRUE),
          median    = \(x) stats::median(x, na.rm = TRUE),
          mean      = \(x) mean(x, na.rm = TRUE),
          sd        = \(x) stats::sd(x, na.rm = TRUE),
          n_missing = \(x) sum(is.na(x))
        ),
        .names = "{.col}__{.fn}"
      )
    ) |>
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = c("variable", "stat"),
      names_sep = "__",
      values_to = "value"
    ) |>
    tidyr::pivot_wider(
      names_from = "stat",
      values_from = "value"
    ) |>
    as.data.frame()
}

#' Correlation with pairwise complete obs and auto numeric selection
#' @param x A vector, matrix, data.frame, or tibble.
#' @param y Optional vector/matrix/data.frame/tibble for cross-correlation.
#' @param use Missing-data handling; default "pairwise.complete.obs".
#' @param ... Passed to [stats::cor()] (e.g., method = "spearman").
#' @return Correlation vector/matrix as in [stats::cor()].
#' @export
bcor <- function(x, y = NULL, use = "pairwise.complete.obs", ...) {

  select_numeric <- function(obj) {
    if (is.data.frame(obj)) {
      keep <- vapply(obj, is.numeric, logical(1))
      obj  <- obj[ , keep, drop = FALSE]
      if (ncol(obj) == 0L) stop("No numeric columns in data frame.", call. = FALSE)
    } else if (is.vector(obj) && !is.numeric(obj)) {
      stop("Vector input must be numeric.", call. = FALSE)
    }
    obj
  }

  x <- select_numeric(x)
  if (!is.null(y)) y <- select_numeric(y)

  stats::cor(x = x, y = y, use = use, ...)
}

#' Sum with na.rm = TRUE by default
#'
#' Wrapper around base::sum() that always removes NAs unless overridden.
#'
#' @param x A numeric vector.
#' @param ... Additional arguments passed to [base::sum()].
#' @return Numeric sum of x.
#' @export
bsum <- function(x, ...) {
  sum(x, na.rm = TRUE, ...)
}

#' Mean with na.rm = TRUE by default
#'
#' Wrapper around base::mean() that always removes NAs unless overridden.
#'
#' @param x A numeric vector.
#' @param ... Additional arguments passed to [base::mean()].
#' @return Numeric mean of x.
#' @export
bmean <- function(x, ...) {
  mean(x, na.rm = TRUE, ...)
}
