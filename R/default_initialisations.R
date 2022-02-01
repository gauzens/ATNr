#' @title Make parameter matrix
#'
#' @param BM float vector, body mass of species.
#' @param b0 const
#' @param bprey const
#' @param bpred const
#' @param E const
#' @param T.K, Celsius to Kelvin conversion
#' @param T0, Default temperature in Kelvin
#' @param k, Boltzmann constant
#'
#' @return A matrix filled with estimated values 
#' for a model parameter that depends on prey and predator body masses (see details)
#'  
#' @export
#'
#' @details Make a parameter matrix that depends on both predators
#' and prey and that is used to define attack rates and handling
#' times based on the general allometric equation:
#' \deqn{p_{i,j} = b_0 * BM_i^{bprey} * BM_j^{bpred} * exp(-E * (T0-T.K) / (k * T.K * T0))}
create_matrix_parameter <- function(
  BM,
  b0,
  bprey,
  bpred,
  E,
  T.K,
  T0,
  k
) {
  # create a matrix of link specific parameters (i.e that depend on both prey and predators).
  # Classically, this is needed to define attack rates or handling times.
  # based on the general allometric equation:
  # p_{i,j} = b_0*BM_i^{bprey}*BM_j^{bpred}*exp(-E*(T0-T.K)/(k*T.K*T0))
  nb_s <- length(BM)
  M <- matrix(1, nrow = nb_s, ncol = nb_s)
  return(b0 * ((M * BM[, 1])^bprey * t(M * BM[, 1])^bpred) *
               exp(-E * (T0 - T.K) / (k * T.K * T0)))
}
#' @title Default model parameters as in Schneider et al. 2016
#'
#' @description Initialise the default parametrisation for the model for
#'   Schneider et al. (2016).
#' @param model an object of class \emph{ATN (Rcpp_Unscaled_nuts}.
#' @param L.mat numeric matrix, probability of a consumer to attack and capture an encountered resource. See \code{\link{create_Lmatrix}}.
#' @param temperature numeric, ambient temperature of the ecosystem in Celsius.
#'
#' @export
#' 
#' @references Schneider, F. D., Brose, U., Rall, B. C., & Guill, C. (2016).
#'   Animal diversity and ecosystem functioning in dynamic food webs. Nature
#'   Communications, 7(1), 1-8.
#'
#' @return An object of class \emph{ATN (Rcpp_Unscaled_nuts)} with default
#'   parameters as in Schneider et al. (2016).
initialise_default_Unscaled_nuts <- function(
  model,
  L.mat,
  temperature = 20
) {
  utils::data("schneider", envir = environment())
  schneider[["nb_s"]] <- model$nb_s
  schneider[["nb_b"]] <- model$nb_b
  schneider[["nb_n"]] <- model$nb_n
  schneider[["BM"]] <- model$BM
  schneider[["T.K"]] <- temperature + 273.15
  # w parameter: how a predator splits its foraging time on the different
  # species. by default a predator split its foraging time equally between all
  # the prey sum of w values should be equal to 1 for a given predator.
  model$ext <- 1e-6
  model$w <- sweep(x = model$fw, MARGIN = 2, FUN = "/", colSums(model$fw))
  model$w <- model$w[, (model$nb_b + 1):model$nb_s]

  # Plant nutrient uptake efficiency
  model$K <- with(schneider,
                  matrix(stats::runif(nb_b * nb_n, nut_up_min, nut_up_max),
                         nrow = nb_n, ncol = nb_b))

  # turnover rate of the nutrients
  model$D <- schneider$D

  # maximal nutrient level
  model$S <- with(schneider, stats::rnorm(nb_n, mu_nut, sd_nut))
  # growth rate of the basal species
  model$r <- with(schneider, BM[1 : nb_b] ^ -0.25 * exp(-0.22 * (T0 - T.K) / (k * T.K * T0)))
  # per gram metabolic rate
  model$X <- with(schneider,
                  c(rep(x_P, nb_b), rep(x_A, nb_s - nb_b)) *
                    BM^-0.25 * exp(-0.69 * (T0 - T.K) / (k * T.K * T0)))

  # species efficiencies
  model$e <- with(schneider, c(rep(e_P, nb_b), rep(e_A, nb_s - nb_b)))
  # species specific capture rate (encounter rate * predaion success)
  model$b <- with(schneider, create_matrix_parameter(BM, b0, bprey, bpred, E.b, T.K, T0, k)) * L.mat
  model$b <- model$b[, (model$nb_b + 1):model$nb_s]
  # specific values for plants: BM^beta = 20
  model$b[1:model$nb_b, ] <-  with(schneider,
    t(replicate(model$nb_b, 20 * model$BM[(model$nb_b + 1):nrow(model$BM), 1]^bpred)) *
      L.mat[1:model$nb_b, (model$nb_b + 1):model$nb_s]
  )

  # interference competition
  model$c <- with(schneider, stats::rnorm(nb_s - nb_b, mu_c, sd_c) * exp(-0.65 * (T0 - T.K) / (k * T.K * T0)))
  # handling time
  model$h <- with(schneider, create_matrix_parameter(BM, h0, hprey, hpred, E.h, T.K, T0, k))
  model$h <- model$h[, (model$nb_b + 1):model$nb_s]
  # Hill exponent
  model$q <- stats::rnorm(1, 1.5, 0.2)
  # plant stoichiometry (relative content in  the nutrients) !!!!!!!!!! to update. here assume 2 nutrients only !!!!!!!
  model$V <- with(schneider,
                  matrix(stats::runif(nb_b * nb_n, 1, 2), nrow = nb_n, ncol = nb_b))
  model$V <- sweep(x = model$V, MARGIN = 2, FUN = "/", colSums(model$V))
  # growth rate of plant species !!!!!!!!!!! temperature independant right now !!!!!!!!!!!!!!!
  model$r <- with(schneider, BM[1:nb_b]^-0.25)
  # initialisation of the matrix of feeding rates.
  # all values are 0 for now. Updated at each call of the ODEs estimations.
  model$F <- with(schneider, matrix(0.0, nrow = nb_s, ncol = nb_s - nb_b))

  return(model)
}

#' @title Default parameters for the scaled version of ATN as in Delmas et al.
#'   2016
#'
#' @description Initialise the default parametrisation for the scaled version of
#'   the ATN model as in Delmas et al. (2016).
#'
#' @param model an object of class \emph{Rcpp_Scaled}.
#' 
#' @export
#'
#' @references Delmas, E., Brose, U., Gravel, D., Stouffer, D.B. and Poisot, T.
#'   (2017), Simulations of biomass dynamics in community food webs. Methods
#'   Ecol Evol, 8: 881-886. https://doi.org/10.1111/2041-210X.12713
#'
#' @return An object of class \emph{Rcpp_Scaled} with default
#'   parameters as in Delmas et al. (2017).
#'
initialise_default_Scaled <- function(model) {

  utils::data("schneider", envir = environment())
  schneider[["nb_s"]] <- model$nb_s
  schneider[["nb_b"]] <- model$nb_b
  schneider[["BM"]] <- model$BM

  # allometric constant for growth rate:
  ar <- 1
  # Carrying capacity for basal species i.e. nutrient pool accessible to ALL species simultaneously
  K <- 10

  # w parameter: how a predator splits its foraging time on the different
  # species. by default a predator split its foraging time equally between all
  # the prey sum of w values should be equal to 1 for a given predator.
  model$w <- sweep(x = model$fw, MARGIN = 2, FUN = "/", colSums(model$fw))
  model$w <- model$w[, (model$nb_b + 1):model$nb_s]

  # per gram metabolic rate
  model$X <- with(schneider, 0.314 * BM^-0.25 / BM[1]^-0.25)
  model$X[1:schneider$nb_b] <- 0.0
  # species efficiencies
  model$e <- with(schneider, c(rep(e_P, nb_b), rep(e_A, nb_s - nb_b)))
  # interference competition
  model$c <- rep(0.8, model$nb_s - model$nb_b)
  # max feeding rate
  model$max_feed <- rep(8, model$nb_s - model$nb_b)
  # half sturation density:
  model$B0 <- rep(0.5, model$nb_s - model$nb_b)
  # Hill exponent
  model$q <- stats::rnorm(1, 1.2, 0.2)
  # max growth rate of plant species
  model$r <- with(schneider, (ar * BM[1:nb_b]^-0.25) / (ar * BM[1]^-0.25))
  # max carrying capacity of all plant species
  model$K <- 10
  # initialisation of the matrix of feeding rates.
  # all values are 0 for now. Updated at each call of the ODEs estimations.
  model$F <- with(schneider, matrix(0.0, nrow = model$nb_s, ncol = model$nb_s - model$nb_b))
  
  # plant resource competition, matrix should be symetric by default
  # model$alpha = matrix(runif(model$nb_b*model$nb_b, 0.5, 1), nrow = model$nb_b, ncol = model$nb_b)
  # model$alpha[lower.tri(model$alpha)] = t(model$alpha)[lower.tri(model$alpha)]

  model$alpha <- matrix(0, nrow = model$nb_b, ncol = model$nb_b)
  diag(model$alpha) = 1
  
  return(model)
}



#' @title Default parameters for the scaled version of ATN as in Binzer et al.
#'   2016, with updates from Gauzens et al. 2020
#'
#' @description Initialise the default parametrisation for the scaled version of
#'   the ATN model as in Binzer et al. (2016), with updates from Gauzens et al. 2020
#'
#' @param model an object of class \emph{ATN (Rcpp_Unscaled)}.
#' @param temperature numeric, ambient temperature of the ecosystem in Celsius.
#' 
#' @export
#' 
#' @references Binzer, A., Guill, C., Rall, B. C. & Brose, U.
#' Interactive effects of warming, eutrophication and size structure: impacts on biodiversity and food-web structure.
#' Glob. Change Biol. 22, 220-227 (2016).
#' Gauzens, B., Rall, B.C., Mendonca, V. et al.
#' Biodiversity of intertidal food webs in response to warming across latitudes.
#' Nat. Clim. Chang. 10, 264-269 (2020). https://doi.org/10.1038/s41558-020-0698-z
#'
#' @return An object of class \emph{ATN (Rcpp_Unscaled)} with default
#'   parameters as in Delmas et al. (2017).

initialise_default_Unscaled <- function(model, temperature = 20){
  k <- 8.6173324e-5
  T0 <- 293.15
  T.K <- temperature + 273.15
  model$X <- exp(-16.54) * model$BM^-0.31 * exp(-0.69 * (T0 - T.K) / (k * T.K * T0))
  e_P <- 0.545
  e_A <- 0.906
  model$e <- c(rep(e_P, model$nb_b), rep(e_A, model$nb_s - model$nb_b))
  model$c <- rep(0.8, model$nb_s - model$nb_b) * exp(-0.65 * (T0 - T.K) / (k * T.K * T0))
  model$r <- exp(-15.68) * model$BM[1:model$nb_b]^-0.25 * exp(-0.84 * (T0 - T.K) / (k * T.K * T0))
  model$K <- 40 * model$BM[1:model$nb_b]^0.28 * exp(0.71 * (T0 - T.K) / (k * T.K * T0))
  # here attack rate decrease with consumer BM, and handling time increase with consumer BM
  # Is it an error in Binzer et al. ?
  model$a <- create_matrix_parameter(model$BM, exp(-13.1), 0.25, -0.8, -0.38, T.K, T0, k)
  model$a <- model$a * model$fw
  model$a <- model$a[, (model$nb_b + 1):model$nb_s]

  model$h <- create_matrix_parameter(model$BM, exp(9.66), -0.45, 0.47, -0.26, T.K, T0, k)
  model$h <- model$h * model$fw
  model$h <- model$h[, (model$nb_b + 1):model$nb_s]
  model$q <- 1.2
  
  # plant resource competition, matrix should be symetric by default
  # model$alpha = matrix(runif(model$nb_b*model$nb_b, 0.5, 1), nrow = model$nb_b, ncol = model$nb_b)
  # model$alpha[lower.tri(model$alpha)] = t(model$alpha)[lower.tri(model$alpha)]
  
  model$alpha <- matrix(0, nrow = model$nb_b, ncol = model$nb_b)
  diag(model$alpha) = 1

  return(model)
}
