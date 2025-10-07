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

  # Auto-convert first column to row names if requested
  if (auto_rownames) {
    # Call bcoltonames with auto-detection (col = NULL)
    # It will handle all the logic and messaging
    df <- bcoltonames(df, col = NULL, remove = TRUE)
  }

  return(df)
}
