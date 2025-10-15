#' Convert a column to row names
#' @param df A data.frame
#' @param col Column to convert to row names. Can be numeric index or name.
#'   If NULL (default), automatically selects based on heuristics.
#' @param remove Logical. Remove the column after conversion. Default TRUE.
#' @return data.frame with modified row names
#' @export
#' @examples
#' \dontrun{
#' # Auto-detect column to use as row names
#' df <- bcoltonames(df)
#'
#' # Specify column by name
#' df <- bcoltonames(df, "ID")
#'
#' # Specify column by index
#' df <- bcoltonames(df, 1)
#'
#' # Keep the column after setting row names
#' df <- bcoltonames(df, "ID", remove = FALSE)
#' }
bcoltonames <- function(df, col = NULL, remove = TRUE) {
  if (!is.data.frame(df)) {
    stop("Input must be a data.frame")
  }

  # If no column specified, auto-detect
  if (is.null(col)) {
    # Check if we have at least 2 columns (need at least one data column remaining)
    if (ncol(df) < 2) {
      message("Cannot auto-convert to row names: data frame has only one column")
      return(df)
    }

    # Common patterns for ID/name columns
    id_patterns <- c("^id$", "^ID$", "^name$", "^Name$",
                     "^row\\.?names?$", "^sample$", "^Sample$",
                     "^label$", "^Label$", "^identifier$", "^key$")

    # Try to find a suitable column
    suitable_col <- NULL

    # First, check for columns matching common ID patterns
    for (pattern in id_patterns) {
      matches <- grep(pattern, names(df), ignore.case = FALSE)
      if (length(matches) > 0) {
        # Check if this column meets our criteria
        col_data <- df[[matches[1]]]
        has_no_duplicates <- !any(duplicated(col_data))
        has_no_nas <- !any(is.na(col_data))

        if (has_no_duplicates && has_no_nas) {
          suitable_col <- matches[1]
          message(sprintf("Auto-selected column %d (%s) for row names (matched ID pattern)",
                          suitable_col, names(df)[suitable_col]))
          break
        }
      }
    }

    # If no ID pattern matched, fall back to original heuristics
    if (is.null(suitable_col)) {
      for (i in 1:ncol(df)) {
        col_data <- df[[i]]
        col_name <- names(df)[i]

        # Check criteria: empty/blank name OR non-numeric content
        is_suitable <- (is.null(col_name) ||
                          trimws(col_name) == "" ||
                          !is.numeric(col_data))

        # Also check for no duplicates and no NAs
        has_no_duplicates <- !any(duplicated(col_data))
        has_no_nas <- !any(is.na(col_data))

        if (is_suitable && has_no_duplicates && has_no_nas) {
          suitable_col <- i
          message(sprintf("Auto-selected column %d (%s) for row names",
                          suitable_col,
                          ifelse(names(df)[i] == "", "[unnamed]", names(df)[i])))
          break  # Use first suitable column found
        }
      }
    }

    if (is.null(suitable_col)) {
      return(df)
    }

    col <- suitable_col
  }

  # Convert col to numeric index if it's a name
  if (is.character(col)) {
    col_idx <- which(names(df) == col)
    if (length(col_idx) == 0) {
      stop(sprintf("Column '%s' not found in data frame", col))
    }
    col <- col_idx[1]
  }

  # Validate column index
  if (!is.numeric(col) || col < 1 || col > ncol(df)) {
    stop(sprintf("Invalid column index: %s", col))
  }

  # Get the column data
  col_data <- df[[col]]

  # Check for duplicates
  if (any(duplicated(col_data))) {
    warning(sprintf("Column %d contains duplicates. Not converting to row names.", col))
    return(df)
  }

  # Check for NAs
  if (any(is.na(col_data))) {
    warning(sprintf("Column %d contains NAs. Not converting to row names.", col))
    return(df)
  }

  # Set row names
  rownames(df) <- as.character(col_data)

  # Remove column if requested
  if (remove) {
    df <- df[-col]
  }

  return(df)
}
