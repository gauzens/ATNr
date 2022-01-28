#' @title Initialize an ATN model, following Schneider et al. 2016, Nature Communication
#'
#' @param nb_s integer, number of total species.
#' @param nb_b integer, number of basal species.
#' @param nb_n integer, number of nutrients.
#' @param BM float vector, body mass of species.
#' @param fw binary adjacency matrix of the food web.
#' 
#' @export
#'
#' @details A model is defined by the total number of species
#' (\emph{nb_s}), the number of basal species (\emph{nb_b}),
#' the number of nutrients (\emph{nb_n}), the body masses
#'  (\emph{BM}) of species, and the adjacency matrix (\emph{fw})
#'  representing species interactions.
#' Nutrients are not counted as species.
#'
#' @return An object of class \emph{ATN (Rcpp_parameters_prefs)}.
#'
#' @examples
#' library(ATNr)
#' set.seed(123)
#' n_species <- 50
#' n_basal <- 20
#' n_nutrients <- 2
#' masses <- sort(10^runif(n_species, 2, 6)) #body mass of species
#' L <- create_Lmatrix(masses, n_basal)
#' fw <- L
#' fw[fw > 0] <- 1
#' mod <- create_model_Unscaled_nuts(n_species, n_basal, n_nutrients, masses, fw)
create_model_Unscaled_nuts <- function(
  nb_s,
  nb_b,
  nb_n = 2,
  BM,
  fw
) {
  # check input is correct
  if (any(c(nb_s, nb_b, nb_n) %% 1 != 0)) {
    stop("nb_s, nb_b and nb_n must all be integers")
  }
  if (length(BM) != nb_s) {
    stop("BM should have length equal to nb_s (", nb_s, ")")
  }
  if (dim(fw)[1] != dim(fw)[2]) {
    stop("Food web matrix is not a square matrix")
  }
  if (nb_s != ncol(fw)) {
    stop("Number of species and food web matrix do not match")
  }

  model <- methods::new(Unscaled_nuts, nb_s, nb_b, nb_n)

  # THIS WE CAN EVEN PUT IN THE CONSTRUCTOR, PERHAPS?
  model[["BM"]] <- BM
  model[["log_BM"]] <- log10(BM)
  model[["fw"]] <- fw
  return(model)
}

#' @title Initialize an ATN model, following Delmas et al. 2017, Methods in Ecology and Evolution
#'
#' @param nb_s integer, number of total species.
#' @param nb_b integer, number of basal species.
#' @param BM float vector, body mass of species.
#' @param fw binary adjacency matrix of the food web.
#' 
#' @export
#'
#' @details A model is defined by the total number of species
#' (\emph{nb_s}), the number of basal species (\emph{nb_b}),
#' the number of nutrients (\emph{nb_n}), the body masses
#' (\emph{BM}) of species, and the adjacency matrix (\emph{fw})
#' representing species interactions.
#'
#' @return An object of class \emph{ATN (Rcpp_parameters_prefs)}.
#'
#' @references
#' @references Delmas, E., Brose, U., Gravel, D., Stouffer, D.B. and Poisot, T.
#'   (2017), Simulations of biomass dynamics in community food webs. Methods
#'   Ecol Evol, 8: 881-886. https://doi.org/10.1111/2041-210X.12713
#'
#' @examples
#' library(ATNr)
#' set.seed(123)
#' n_species <- 50
#' n_basal <- 20
#' masses <- sort(10^runif(n_species, 2, 6)) #body mass of species
#' L <- create_Lmatrix(masses, n_basal)
#' fw <- L
#' fw[fw > 0] <- 1
#' mod <- create_model_Scaled(n_species, n_basal, masses, fw)
create_model_Scaled <- function(
  nb_s,
  nb_b,
  BM,
  fw
) {
  # check input is correct
  if (any(c(nb_s, nb_b) %% 1 != 0)) {
    stop("nb_s, nb_b and nb_n must all be integers")
  }
  if (length(BM) != nb_s) {
    stop("BM should have length equal to nb_s (", nb_s, ")")
  }
  if (dim(fw)[1] != dim(fw)[2]) {
    stop("Food web matrix is not a square matrix")
  }
  if (nb_s != ncol(fw)) {
    stop("Number of species and food web matrix do not match")
  }

   model <- methods::new(Scaled, nb_s, nb_b)

  # THIS WE CAN EVEN PUT IN THE CONSTRUCTOR, PERHAPS?
  model[["BM"]] <- BM
  model[["log_BM"]] <- log10(BM)
  model[["fw"]] <- fw
  return(model)
}


#' @title Initialize an ATN model, following Binzer et al. 201, Global Change Biology
#'
#' @param nb_s integer, number of total species.
#' @param nb_b integer, number of basal species.
#' @param BM float vector, body mass of species.
#' @param fw binary adjacency matrix of the food web.
#' 
#' @export
#'
#' @details A model is defined by the total number of species
#' (\emph{nb_s}), the number of basal species (\emph{nb_b}),
#' the number of nutrients (\emph{nb_n}), the body masses
#'  (\emph{BM}) of species, and the adjacency matrix (\emph{fw})
#'  representing species interactions.
#'
#' @return An object of class \emph{ATN (Rcpp_parameters_prefs)}.
#'
#' @references
#' @references Binzer, A., Guill, C., Rall, B.C. and Brose, U. (2016),
#' Interactive effects of warming, eutrophication and size structure: impacts on biodiversity and food-web structure.
#' Glob Change Biol, 22: 220-227. https://doi.org/10.1111/gcb.13086
#' Gauzens, B., Rall, B.C., Mendonca, V. et al.
#' Biodiversity of intertidal food webs in response to warming across latitudes.
#' Nat. Clim. Chang. 10, 264-269 (2020). https://doi.org/10.1038/s41558-020-0698-z
#'
#' @examples
#' library(ATNr)
#' set.seed(123)
#' n_species <- 50
#' n_basal <- 20
#' masses <- sort(10^runif(n_species, 1, 6)) #body mass of species
#' L <- create_Lmatrix(masses, n_basal)
#' fw <- L
#' fw[fw > 0] <- 1
#' mod <- create_model_Unscaled(n_species, n_basal, masses, fw)
create_model_Unscaled <- function(
  nb_s,
  nb_b,
  BM,
  fw
) {
  # check input is correct
  if (any(c(nb_s, nb_b) %% 1 != 0)) {
    stop("nb_s, nb_b and nb_n must all be integers")
  }
  if (length(BM) != nb_s) {
    stop("BM should have length equal to nb_s (", nb_s, ")")
  }
  if (dim(fw)[1] != dim(fw)[2]) {
    stop("Food web matrix is not a square matrix")
  }
  if (nb_s != ncol(fw)) {
    stop("Number of species and food web matrix do not match")
  }

  model <- methods::new(Unscaled, nb_s, nb_b)

  # THIS WE CAN EVEN PUT IN THE CONSTRUCTOR, PERHAPS?
  model[["BM"]] <- BM
  model[["log_BM"]] <- log10(BM)
  model[["fw"]] <- fw
  return(model)
}
