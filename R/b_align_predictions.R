#' Align predictions with original data rows
#'
#' @description
#' Internal utility function to align predictions from cleaned data back to original data,
#' inserting NAs where rows were removed during cleaning.
#'
#' @param original_data Original data frame before cleaning
#' @param cleaned_data Cleaned data frame after removing missing values
#' @param predictions Vector of predictions from model fit on cleaned data
#'
#' @return Vector of predictions with same length as original_data, with NAs in removed positions
#' @keywords internal
balign_predictions <- function(original_data, cleaned_data, predictions) {
  # Create a key to match rows
  original_data$.row_id <- seq_len(nrow(original_data))
  cleaned_data$.row_id <- as.numeric(rownames(cleaned_data))

  # Create full-length result vector
  result <- rep(NA_real_, nrow(original_data))

  # Fill in predictions at correct positions
  result[cleaned_data$.row_id] <- predictions

  return(result)
}

