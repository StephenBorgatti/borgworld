# ez_data

#' Cut a numeric vector into k ordered groups via k-means
#'
#' @param x Numeric vector (coercible).
#' @param k Integer number of groups (default 3).
#' @return Integer vector (1..k) aligned so 1=lowest mean group, k=highest.
#'         Non-finite inputs return NA.
#' @export
#' @examples
#'   # Simple scatter
#'   bautocut(centrality,3)
bautocut <- function(x, k = 3) {
  x <- as.numeric(x)
  idx <- is.finite(x)
  if (sum(idx) < k || length(unique(x[idx])) < k) {
    stop("Need at least ", k, " distinct finite values.")
  }

  set.seed(1)
  km <- stats::kmeans(x[idx], centers = k)

  centers <- tapply(x[idx], km$cluster, mean)
  ranks <- rank(centers, ties.method = "first")  # 1=lowest center ... k=highest

  out <- rep(NA_integer_, length(x))
  out[idx] <- ranks[km$cluster]
  out
}

#' Standardize columns of df
#'
#' @param df data frame
#' @return data frame of z-scores
#'         Non-finite inputs return NA.
#' @export
#' @examples
#'   # Simple scatter
#'   bstandardize(hdi)
bstandardize <- function(df, k = 3) {
  df |> mutate(across(all_of(vars), z_pop, .names = "z_{.col}"))

}
