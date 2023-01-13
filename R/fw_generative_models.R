#' @title Create a food web based on the niche model
#'
#' @description Function to generate a food web based on the niche model
#'   (Williams and Martinez, 2000) based on the number of species and
#'   connectance. Corrections from Allesina et al. (2008) are used.
#' @details If at least one species has not resource or consumer (i.e. it is an
#'   isolated species), another food web is generated, until a maximum of 100
#'   iterations.
#' @param S integer, number of species.
#' @param C numeric, connectance i.e. the number of realized links over the all
#'   possible links.
#' @export
#' @return A (square) matrix with zeros (no interaction) and ones (species j
#'   consume species i).
#' @references Williams, R. J., & Martinez, N. D. (2000). Simple rules yield
#'   complex food webs. Nature, 404(6774), 180-183.
#'
#'   Allesina, S., Alonso, D., & Pascual, M. (2008). A general model for food
#'   web structure. science, 320(5876), 658-661.
#' @examples
#' set.seed(123)
#' web_niche <- create_niche_model(30, .1)
#' image(t(web_niche))
create_niche_model <- function(S, C) {
  stopifnot(S > 0 && C > 0)
  niche_model <- function(S, C) {
    # niches of species
    niche <- sort(stats::runif(S))
    # feeding ranges, using the correction from Allesina et al. (2008)
    if (((S - 1) / (2 * S * C)) - 1 < 0) {
      stop("Beta distribution parameter < 0. Try to decrease C.")
    }
    diet <- stats::rbeta(S, 1, ((S - 1) / (2 * S * C)) - 1) * niche
    # feeding center, using the correction from Allesina et al. (2008)
    center <- sapply(seq_len(S),
                     function(i) {
                       n <- niche[i]
                       r <- diet[i]
                       ifelse(n + r / 2 <= 1,
                              stats::runif(1, r / 2, n),
                              stats::runif(1, r / 2, 1 - r / 2)
                       )
                     })
    species <- seq_len(S)
    # create food web adjacency matrix
    fw <- matrix(rep(0, S ^ 2), S, S)
    for (sp in species) {
      preys <- (center[sp] - diet[sp] / 2 <= niche) &
        (niche <= center[sp] + diet[sp] / 2)
      fw[preys, sp] <- 1
    }
    return(fw)
  }
  fw <- niche_model(S, C)
  # check for isolated species
  isolated <- ifelse(any(colSums(fw) + rowSums(fw) == 0), TRUE, FALSE)
  # check if trophic levels can be calculated
  tro_lev <- tryCatch(ATNr::TroLev(fw), error = function(e) NULL)
  # check is fw is connected
  connected <- is_connected(fw)
  i <- 0
  while((isolated | is.null(tro_lev) | !connected) & i < 100) {
    fw <- niche_model(S, C)
    # check first if isolated then TL calculation and then detection of connected components
    # no need to make the 3 of them each time as one is enough to reject
    isolated <- ifelse (any(colSums(fw) + rowSums(fw) == 0), TRUE, FALSE)
    if (!isolated) {
      tro_lev <- tryCatch(ATNr::TroLev(fw), error = function(e) NULL)
    }
    if (!isolated & !is.null(tro_lev)){
      connected <- is_connected(fw)
    }
    i <- i + 1
  }
  if (isolated) warning("Presence of an isolated species after 100 iterations.")
  if (is.null(tro_lev)) warning("Trophic levels cannot be calcualted after 100 iterations.")
  if (!is_connected(fw)) warning("Several connected components detected")

  # reorder matrix to put basal species first
  basals <- which(colSums(fw) == 0)
  consumers <- which(colSums(fw) > 0)
  fw <- fw[c(basals, consumers), c(basals, consumers)]
  return(fw)
}

#' @title Make L matrix
#'
#' @param BM float vector, body mass of species.
#' @param nb_b integer, number of basal species.
#' @param Ropt numeric, consumer/resource optimal body mass ratio.
#' @param gamma numeric, code for the width of the Ricker function.
#' @param th float, the threshold below which attack rates are considered = 0.
#' @export
#' @details The L matrix contains the probability for an attack event to be
#'   successful based on allometric rules and a Ricker function defined by
#'   \emph{Ropt} and \emph{gamma}. If at least one species has not resource or
#'   consumer (i.e. it is an isolated species), another food web is generated,
#'   until a maximum of 100 iterations.
#'
#' @return A numeric matrix with the probability for an attack event between two
#'   species to be successful.
#'
#' @examples
#' set.seed(123)
#' mass <- sort(10 ^ runif(30, 2, 6))
#' L <- create_Lmatrix(mass, nb_b = 10, Ropt = 100)
#' image(t(L))
create_Lmatrix <- function(
  BM,
  nb_b,
  Ropt = 100,
  gamma = 2,
  th = 0.01
) {
  stopifnot(all(BM > 0) && nb_b >= 0 && Ropt > 0 && gamma > 0 && th >= 0)
  Lmatrix <- function(BM, nb_b, Ropt, gamma, th) {
    s <- length(BM)
    L <- matrix(rep(BM, s), s, s, byrow = TRUE) /
      (matrix(rep(BM, s), s, s) * Ropt)
    L <- (L * exp(1 - L)) ^ gamma
    L[L < th] <- 0
    L[, 1:nb_b] <- 0
    return(L)
  }
  s <- length(BM)
  L <- Lmatrix(BM, nb_b, Ropt, gamma, th)
  # check for isolated species
  isolated <- ifelse(any(colSums(L) + rowSums(L) == 0), TRUE, FALSE)
  # check for isolated consumers
  cons_no_prey <- ifelse(any(colSums(L[, (nb_b + 1) : s]) == 0), TRUE, FALSE)
  # check if trophic levels can be calculated
  tro_lev <- tryCatch(ATNr::TroLev(L), error = function(e) NULL)
  # check for different connected components
  connected <- ATNr::is_connected(L) #BUG - add namespace to avoid conflicts with igraph::is_connected()
  if (is.null(tro_lev)) warning("Cannot compute trophic levels.")
  if (isolated) warning("Presence of an isolated species.")
  if (cons_no_prey) warning("Presence of consumer without prey.")
  if (!connected) warning("Several conected component detected")
  return(L)
}
