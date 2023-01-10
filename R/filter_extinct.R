#' Filter Extinct Species
#'
#' @param df deSolve matrix as returned from lsoda_wrapper().
#' @param model ATNr model object, from which extinction threshold is extracted.
#'
#' @details Set to zero species biomass that are below the extinction threshold.
#'
#' @return df with values below th set to zero.
.filter_extinct <- function(df, model) {
  ext <- which(df[, -1] < model$ext, arr.ind = TRUE)
  for (j in unique(ext[, 2])) {
    df[min(ext[ext[, 2] == j, 1]):nrow(df), j + 1] <- 0
  }
  return(df)
}
