#' Read a tabular block from the clipboard
#' @param header Logical. First row has names. Default TRUE.
#' @param sep Field separator. Default tab.
#' @param auto_rownames Logical. Automatically convert first column to row names
#'   when it's non-numeric or has no header. Default TRUE.
#' @return data.frame
#' @export
#' @examples
#' \dontrun{
#' # Copy a range in Excel, then:
#' bpaste()
#' }
bpaste <- function(header = TRUE, sep = "\t", auto_rownames = TRUE) {
  df <- utils::read.table("clipboard",
                          header = header, sep = sep,
                          stringsAsFactors = FALSE, check.names = FALSE)

  # Auto-detect if first column should be row names
  if (auto_rownames && ncol(df) > 1) {
    first_col <- df[[1]]
    first_col_name <- names(df)[1]

    # Check if first column should become row names:
    # 1. Column name is empty/blank (common when pasting from Excel with row labels)
    # 2. OR column is non-numeric (character/factor)
    should_convert <- (is.null(first_col_name) ||
                         trimws(first_col_name) == "" ||
                         !is.numeric(first_col))

    if (should_convert) {
      # Check for duplicate row names
      if (any(duplicated(first_col))) {
        warning("First column contains duplicates. Not converting to row names.")
      } else {
        # Set row names and remove first column
        rownames(df) <- first_col
        df <- df[-1]
      }
    }
  }

  return(df)
}
