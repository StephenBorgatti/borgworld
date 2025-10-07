#' Display borgworld package help
#'
#' A convenient wrapper to display the borgworld package documentation
#'
#' @param topic Optional character string specifying a specific help topic.
#'        If NULL (default), shows the package help index.
#'
#' @return Invisibly returns NULL. Called for its side effect of displaying help.
#'
#' @export
#' @examples
#' # Show package index
#' bhelp()
#'
#' # Show help for specific function
#' bhelp("bhiclus")
bhelp <- function(topic = NULL) {
  if (is.null(topic)) {
    # Show package help index
    help(package = "borgworld")
  } else {
    # Show help for specific topic within borgworld
    help(topic, package = "borgworld", help_type = "text")
  }
  invisible(NULL)
}
