#' @title Sort custom input
#'
#' @param BM numeric vector, body mass of species.
#' @param fw adjacency matrix of the food web.
#' 
#' @export
#'
#' @details Body masses and food web matrix should be arranged with the first
#'   elements/columns being for basal species. This is a requirement for the Cpp
#'   class and must be enforced before initializing the Rcpp_Schneider and
#'   Rcpp_Delmas objects.
#'
#' @return A list with sorted body masses (\emph{body.mass}) and food web
#'   matrix (\emph{food.web}).
#'
#' @examples
#' bm <- runif(10, 10, 50)
#' fw <- matrix(as.numeric(runif(100) > .9), 10, 10)
#' sort_input(bm, fw)
sort_input <- function(BM, fw) {
  stopifnot(nrow(fw) == ncol(fw))
  stopifnot(length(BM) == nrow(fw))
  basals <- which(colSums(fw) == 0)
  consumers <- which(colSums(fw) > 0)
  adj <- fw[c(basals, consumers), c(basals, consumers)]
  bm <- BM[c(basals, consumers)]
  return(list("body.mass" = bm, "food.web" = adj))
}
