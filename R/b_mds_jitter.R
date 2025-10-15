bmds.jitter <- function(matrix, type, spread) {
  pos <- position_jitter(width = spread, height = spread, seed = 123)
  bmds(matrix, type, plot = F)$points |>
    as_tibble(rownames = "label") |>
    ggplot(aes(x = Dim1, y = Dim2, label = label)) +
    geom_point(position = pos, color = "gold", alpha = .6, size = 3) +
    geom_text(position = pos, vjust = - 1.2, size = 3) +
    coord_equal() +
    theme_minimal()
}
