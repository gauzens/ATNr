#' @title Calculate trophic level of species
#'
#' @param fw numeric matrix, the matrix of the food web.
#'
#' @export
#'
#' @return A numeric vector of species' trophic level.

TroLev <- function(fw) {
  fw <- t(fw)
  nn <- rowSums(fw); nn[nn == 0] <- 1
  ww <- diag(1 / nn)
  L1 <- ww %*% fw
  L2 <- L1 - diag(rep(1, length(nn)))
  b <- -1 * rep(1, length(nn))
  Tro.lev <- solve(L2) %*% b
  return(Tro.lev)
}
