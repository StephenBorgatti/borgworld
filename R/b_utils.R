# butils

#' List variables and their data types
#' @param df A data frame or tibble.
#' @return data.frame with columns: variable, type
#' @export
#' @examples
#' bvars(iris)
bvars <- function(df) {
  data.frame(
    variable = names(df),
    type     = vapply(df, function(x) class(x)[1], character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


##' Copy data fram to the clipboard
#' @param df A data frame or tibble.
#' @return data.frame with columns: variable, type
#' @export
#' @examples
#' bclip(iris)
bclip <- function(df) {
  write.table(df,
              file = "clipboard",
              sep = "\t",
              col.names = NA)   # col.names=NA preserves row names

}

#' @rdname bclip
#' @export
bsnip <- bclip


#' Read a tabular block from the clipboard
#' @param header Logical. First row has names. Default TRUE.
#' @param sep Field separator. Default tab.
#' @return data.frame
#' @export
#' @examples
#' \dontrun{
#' # Copy a range in Excel, then:
#' bpaste()
#' }
bpaste <- function(header = TRUE, sep = "\t") {
  utils::read.table("clipboard",
                    header = header, sep = sep,
                    stringsAsFactors = FALSE, check.names = FALSE)
}


#' Normalize Windows paths to forward slashes
#'
#' Accepts one or more paths (e.g., pasted from Explorer) and returns
#' normalized paths with your chosen slash.
#'
#' @param path Character vector of paths.
#' @param slash Slash to use on Windows. Default "/".
#' @return Character vector of normalized paths.
#' @export
#' @examples
#' bpath(r"(C:\Users\sborg2\Documents\UCINET Data)")
#' bpath(c(r"(C:\Temp\a)", r"(C:\Temp\b)"))
bpath <- function(path, slash = "/") {
  path <- as.character(path)
  # Works cross-platform; on Windows converts backslashes to `slash`
  normalizePath(path, winslash = slash, mustWork = FALSE)
}

#' get working directroy
#'
#' Accepts one or more paths (e.g., pasted from Explorer) and returns
#' normalized paths with your chosen slash.
#'
#' @return working direcotry
#' @export
#' @examples
#' wd()
wd <- function() getwd()
