#' Create scatterplot with regression line using ggplot2
#'
#' @param data Data frame containing the variables
#' @param formula Formula specifying the model (e.g., y ~ x + z1 + z2)
#' @param x Character string naming the x-axis variable. If NULL, uses first predictor
#' @param se Logical, whether to show confidence interval (default TRUE)
#' @param ... Additional arguments passed to geom_smooth()
#'
#' @return A ggplot object
#' @export
#' @import ggplot2
#' @examples
#' bscatterplot(mtcars, mpg ~ wt + cyl + hp)
#' bscatterplot(mtcars, mpg ~ poly(wt, 2) + cyl)
bscatterplot <- function(data, formula, x = NULL, se = TRUE, ...) {
  library(ggplot2)

  # Extract variable names from formula
  response_var <- all.vars(formula)[1]
  formula_terms <- attr(terms(formula), "term.labels")

  # Find base predictor variables (before any transformations)
  predictor_vars <- all.vars(formula)[-1]

  # Determine which x variable to use
  if (is.null(x)) {
    x_var <- predictor_vars[1]
  } else {
    x_var <- x
    if (!x_var %in% names(data)) {
      stop(paste("Variable", x_var, "not found in data"))
    }
  }

  # Create the plot using .data pronoun
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[response_var]])) +
    geom_point(color = "gray40") +
    geom_smooth(method = "lm", formula = formula, se = se, ...) +
    theme_minimal()

  return(p)
}
