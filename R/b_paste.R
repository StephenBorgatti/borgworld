#' Read tabular data from system clipboard (clipr version)
#'
#' @description
#' Reads tabular data from the system clipboard using the clipr package for
#' cross-platform compatibility. Works seamlessly across Windows, macOS, and Linux.
#'
#' @param header Logical. If TRUE, the first line of the clipboard content is used
#'   as column headers. Default is TRUE.
#' @param sep Character. The field separator character. Default is "\\t" (tab).
#' @param auto_rownames Logical. If TRUE, automatically processes row names using
#'   the bcoltonames function. Default is TRUE.
#'
#' @return A data frame containing the clipboard contents.
#'
#' @details
#' This function uses the clipr package to handle platform-specific clipboard
#' access, providing a cleaner and more reliable solution than manual OS detection.
#' The clipr package automatically handles:
#' \itemize{
#'   \item Windows clipboard access
#'   \item macOS pbpaste command
#'   \item Linux xclip/xsel utilities
#'   \item RStudio Server clipboard functionality
#' }
#'
#' @note
#' Requires the clipr package: \code{install.packages("clipr")}
#'
#' The function assumes that clipboard content is in a tabular format suitable
#' for read.table(). For best results, copy data from spreadsheet applications
#' like Excel or Google Sheets.
#'
#' @examples
#' \dontrun{
#' # First ensure clipr is installed
#' if (!requireNamespace("clipr", quietly = TRUE)) {
#'   install.packages("clipr")
#' }
#'
#' # Copy some data from Excel or a spreadsheet, then:
#' df <- bpaste()
#'
#' # Read without headers
#' df <- bpaste(header = FALSE)
#'
#' # Read CSV-formatted clipboard content
#' df <- bpaste(sep = ",")
#'
#' # Read without automatic row name processing
#' df <- bpaste(auto_rownames = FALSE)
#' }
#'
#' @import clipr
#' @export
#' @seealso
#' \code{\link[clipr]{read_clip}} for raw clipboard reading,
#' \code{\link[utils]{read.table}} for reading tabular data
#'
bpaste <- function(header = TRUE, sep = "\t", auto_rownames = TRUE) {
  if (!requireNamespace("clipr", quietly = TRUE)) {
    stop("Package 'clipr' is required. Install it with: install.packages('clipr')")
  }

  clipboard_text <- clipr::read_clip()
  df <- utils::read.table(text = clipboard_text, header = header, sep = sep,
                          stringsAsFactors = FALSE, check.names = FALSE)
  if (auto_rownames) df <- bcoltonames(df, col = NULL, remove = TRUE)
  df <- bdropmissingrows(df)
  message(paste("Dataset has", nrow(df), "rows and", ncol(df), "columns."))
  df
}
