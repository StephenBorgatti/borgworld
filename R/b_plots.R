# b_plots

#' Scatterplot with optional color and regression line
#'
#' @description Create a scatterplot of two columns of a data frame. Optionally
#' color points by a third variable, and add a global OLS regression line
#' (`geom_smooth(method = "lm")`) for the whole dataset (not per color group).
#'
#' @param df A data frame or tibble.
#' @param y,x Column names for the axes (unquoted).
#' @param color_by Optional column for point colors (unquoted).
#' @param add_lm Logical; if TRUE add a regression line. Default TRUE.
#' @param se Logical; if TRUE include confidence band around regression line. Default TRUE.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' if (requireNamespace("palmerpenguins", quietly = TRUE)) {
#'   data(penguins, package = "palmerpenguins")
#'   # Simple scatter
#'   bscatter(penguins, body_mass_g, flipper_length_mm)
#'
#'   # Colored by species, with regression line
#'   bscatter(penguins, fbody_mass_g, flipper_length_mm, color_by = species)
#'
#'   # Without SE band
#'   bscatter(penguins, fbody_mass_g, flipper_length_mm, species, se = FALSE)
#' }
bscatter <- function(df, y, x, color_by = NULL, add_lm = TRUE, se = TRUE) {
  stopifnot(is.data.frame(df))
  x <- rlang::enquo(x); y <- rlang::enquo(y); col <- rlang::enquo(color_by)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = !!x, y = !!y))
  if (!rlang::quo_is_null(col)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = !!col))
  } else {
    p <- p + ggplot2::geom_point()
  }

  if (add_lm) {
    p <- p + ggplot2::geom_smooth(
      data = df,
      mapping = ggplot2::aes(x = !!x, y = !!y),
      method = "lm", se = se, inherit.aes = FALSE
    )
  }
  p + ggplot2::labs(x = rlang::as_name(x), y = rlang::as_name(y),
                    color = if (!rlang::quo_is_null(col)) rlang::as_name(col) else NULL)
}
