#' Read Google Spreadsheet with borgworld conventions
#'
#' Loads a Google spreadsheet and optionally converts a column to row names
#' using bcoltonames, similar to how bpaste works.
#'
#' @param sheet_url URL or ID of the Google Sheet
#' @param sheet Name or index of the specific sheet to read (default is NULL for first sheet)
#' @param col Column to convert to row names. Can be numeric index, name, or NULL for auto-detect.
#'   Set to FALSE to skip row name conversion entirely.
#' @param remove_col Logical. Remove the column after setting as row names (default TRUE)
#' @param skip Number of rows to skip before reading data (default is 0)
#' @param na Character vector of strings to interpret as NA (default is c("", "NA"))
#' @param col_types Column specification (default is NULL for auto-detection)
#' @param auth Should googlesheets4 authenticate? Set to FALSE for public sheets (default is TRUE)
#' @param ... Additional arguments passed to googlesheets4::read_sheet
#'
#' @return A data frame with appropriate row names set if requested
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Read a public Google Sheet with auto row name detection
#' df <- breadgoogle("https://docs.google.com/spreadsheets/d/...", auth = FALSE)
#'
#' # Read without row name conversion
#' df <- breadgoogle("sheet_id", col = FALSE)
#'
#' # Read specific sheet and use column "Sample_ID" as row names
#' df <- breadgoogle("sheet_id", sheet = "Data", col = "Sample_ID")
#'
#' # Read and use first column as row names
#' df <- breadgoogle("sheet_id", col = 1)
#' }
breadgoogle <- function(sheet_url,
                        sheet = NULL,
                        col = NULL,
                        remove_col = TRUE,
                        skip = 0,
                        na = c("", "NA"),
                        col_types = NULL,
                        auth = TRUE,
                        ...) {

  # Check if required package is available
  if (!requireNamespace("googlesheets4", quietly = TRUE)) {
    stop("Package 'googlesheets4' is required. Install it with: install.packages('googlesheets4')")
  }

  # Set authentication
  if (!auth) {
    googlesheets4::gs4_deauth()
  } else {
    # Check if already authenticated, if not, it will prompt
    if (!googlesheets4::gs4_has_token()) {
      message("Google Sheets authentication may be required...")
    }
  }

  # Read the Google Sheet
  message(sprintf("Reading Google Sheet: %s",
                  ifelse(nchar(sheet_url) > 50,
                         paste0(substr(sheet_url, 1, 47), "..."),
                         sheet_url)))

  tryCatch({
    df <- googlesheets4::read_sheet(
      ss = sheet_url,
      sheet = sheet,
      skip = skip,
      na = na,
      col_types = col_types,
      ...
    )

    # Convert to regular data frame (not tibble)
    df <- as.data.frame(df, stringsAsFactors = FALSE)

    message(sprintf("Successfully read %d rows x %d columns", nrow(df), ncol(df)))

  }, error = function(e) {
    if (grepl("auth", tolower(e$message))) {
      stop(paste("Authentication error. Try setting auth = FALSE for public sheets, or authenticate with googlesheets4::gs4_auth()\n",
                 "Original error:", e$message))
    } else {
      stop(paste("Error reading Google Sheet:", e$message))
    }
  })

  # Apply column to row names conversion if not explicitly FALSE
  if (!identical(col, FALSE)) {
    df <- bcoltonames(df, col = col, remove = remove_col)
  }

  return(df)
}


#' Short alias for breadgoogle
#'
#' @param ... All arguments passed to breadgoogle
#' @return A data frame
#' @export
bgoogle <- function(...) {
  breadgoogle(...)
}
