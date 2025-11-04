#' Center Numeric Variables in a Data Frame
#'
#' This function centers all numeric variables in a data frame by subtracting
#' the mean from each numeric column. Non-numeric columns are left unchanged.
#'
#' @param df A data frame
#' @param na.rm Logical. Should missing values be removed when calculating means? Default is TRUE
#' @param exclude Character vector of column names to exclude from centering. Default is NULL
#'
#' @return A data frame with numeric columns centered (mean = 0) and non-numeric columns unchanged
#'
#' @examples
#' # Create sample data
#' df <- data.frame(
#'   id = c("A", "B", "C", "D"),
#'   x = c(1, 2, 3, 4),
#'   y = c(10, 20, 30, 40),
#'   group = factor(c("G1", "G1", "G2", "G2"))
#' )
#'
#' # Center numeric variables
#' centered_df <- bcenter(df)
#'
#' @export
bcenter <- function(df, na.rm = TRUE, exclude = NULL) {
  # Check input
  if (!is.data.frame(df)) {
    stop("Input must be a data frame")
  }

  # Identify numeric columns
  numeric_cols <- sapply(df, is.numeric)
  numeric_names <- names(df)[numeric_cols]

  # Remove excluded columns if specified
  if (!is.null(exclude)) {
    numeric_names <- setdiff(numeric_names, exclude)
  }

  # Center numeric columns
  df_centered <- df
  for (col in numeric_names) {
    df_centered[[col]] <- df[[col]] - mean(df[[col]], na.rm = na.rm)
  }

  # Add attributes to track centering
  attr(df_centered, "centered") <- TRUE
  attr(df_centered, "centered_columns") <- numeric_names
  attr(df_centered, "centering_means") <- sapply(df[numeric_names], mean, na.rm = na.rm)

  return(df_centered)
}
